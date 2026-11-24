#!/usr/bin/env python3
"""
hybrid_tree_bench.py -- run AC17 (FAME), TKN20 and CoverCrypt of the Hybrid-PQ-CP-ABE-Library on the
policies of a tree file (the same tree_med32.txt used by tree_bench for the lattice scheme).

Usage:   python3 hybrid_tree_bench.py <tree_file> <path to Hybrid-PQ-CP-ABE-Library> [reps=15]
Output:  hybrid_results_<time>.csv   (policy, scheme, keygen_ms, enc_ms, dec_ms, ct_bytes, decrypt_ok, wrong_user_rejected)

What is measured (one policy at a time, one ciphertext for the WHOLE policy, OR/AND handled natively):
  keygen = key for a user who holds every attribute of the policy (median of 3, file based: includes writing the key)
  enc    = encrypt a 16-byte secret (buffer API, in memory), median of <reps>
  dec    = decrypt (buffer API, in memory), median of <reps>
  ct     = ciphertext size in bytes (includes the library's AES-GCM wrapper and headers)
"""
import base64, json, os, re, statistics, sys, tempfile, time

if len(sys.argv) < 3:
    print(__doc__); sys.exit(1)
TREE, LIB = sys.argv[1], sys.argv[2]
REPS = int(sys.argv[3]) if len(sys.argv) > 3 else 15
sys.path.insert(0, os.path.join(LIB, "bindings", "python"))
import hybrid_cpabe as abe          # the library's own Python binding

# ------------------------------------------------------------------ tree file
attrs, policies, ell = [], [], None
for line in open(TREE, encoding="utf-8"):
    line = line.split("#")[0].strip()
    if not line: continue
    key, _, rest = line.partition(" ")
    if key == "ATTRS": attrs = rest.split()
    elif key == "POLICY":
        name, _, expr = rest.partition(":")
        policies.append((name.strip(), expr.strip()))

def short(a):                       # library-safe attribute name: lower case letters and digits only
    return re.sub(r"[^a-z0-9]", "", a.lower())

# policy parser: AND, OR, parentheses, k-of-n  ->  nested tuples
def tokenize(s):
    return re.findall(r"[A-Za-z0-9_]+\(|[A-Za-z0-9_]+|\(|\)|,", s)
def parse(tokens):
    pos = 0
    def peek(): return tokens[pos] if pos < len(tokens) else None
    def eat():
        nonlocal pos; pos += 1; return tokens[pos - 1]
    def p_or():
        n = [p_and()]
        while peek() == "OR": eat(); n.append(p_and())
        return n[0] if len(n) == 1 else ("or", n)
    def p_and():
        n = [p_atom()]
        while peek() == "AND": eat(); n.append(p_atom())
        return n[0] if len(n) == 1 else ("and", n)
    def p_atom():
        t = eat()
        if t == "(":
            n = p_or(); assert eat() == ")"; return n
        if re.fullmatch(r"\d+of\(", t):
            k = int(t[:-3]); parts = [p_or()]
            while peek() == ",": eat(); parts.append(p_or())
            assert eat() == ")"
            return ("kof", k, parts)
        return ("leaf", t)
    return p_or()
def leaves(n):
    if n[0] == "leaf": return {n[1]}
    return set().union(*[leaves(c) for c in (n[1] if n[0] != "kof" else n[2])])
def render(n, op_and, op_or):
    """policy text; every AND/OR is binary and fully parenthesised (the form the library's own tests use)"""
    from itertools import combinations
    def fold(op, items):
        out = items[0]
        for x in items[1:]: out = "(" + out + " " + op + " " + x + ")"
        return out
    if n[0] == "leaf": return short(n[1])
    if n[0] == "kof":
        k, parts = n[1], n[2]
        alts = [fold(op_and, [render(p, op_and, op_or) for p in c]) for c in combinations(parts, k)]
        return fold(op_or, alts)
    op = op_and if n[0] == "and" else op_or
    return fold(op, [render(c, op_and, op_or) for c in n[1]])

parsed = []
used_all = set()
for name, expr in policies:
    tree = parse(tokenize(expr)); lv = leaves(tree)
    parsed.append((name, tree, sorted(lv))); used_all |= lv

# ------------------------------------------------------------------ helpers
def load_b64(path):
    data = ""
    for l in open(path):
        l = l.strip()
        if l and not l.startswith("-----"): data += l
    return base64.b64decode(data)
def ms(f):
    t = time.perf_counter(); r = f(); return (time.perf_counter() - t) * 1000.0, r
med = statistics.median

work = tempfile.mkdtemp(prefix="hybrid_bench_")
# CoverCrypt needs every attribute declared as a dimension (present/absent)
cc_conf = os.path.join(work, "covercrypt.json")
json.dump({"dimensions": [{"name": short(a), "attributes": ["present", "absent"]} for a in sorted(used_all)]}, open(cc_conf, "w"))
os.environ["HCPABE_COVERCRYPT_CONFIG"] = cc_conf
os.environ.setdefault("HCPABE_COVERCRYPT_MODE", "hybridized")   # or: post_quantum

SCHEMES = [("AC17_FAME", abe.CPABE_SCHEME_AC17), ("TKN20", abe.CPABE_SCHEME_TKN20),
           ("CoverCrypt_" + os.environ["HCPABE_COVERCRYPT_MODE"], abe.CPABE_SCHEME_COVERCRYPT)]
out = open("hybrid_results_%d.csv" % int(time.time()), "w")
head = "policy,scheme,setup_ms,keygen_ms,enc_ms,dec_ms,ct_bytes,decrypt_ok,wrong_user_rejected"
out.write(head + "\n"); print(head)

secret = os.urandom(16)
for sname, scheme in SCHEMES:
    d = os.path.join(work, sname); os.makedirs(d)
    t_setup, r = ms(lambda: abe.call_setup_with_scheme(d, scheme))
    if r != abe.HCPABE_SUCCESS:
        print("setup failed for", sname); continue
    msk, pkf = os.path.join(d, "cpabe_msk.key"), os.path.join(d, "cpabe_pk.key")
    pk = load_b64(pkf)
    for pname, tree, lv in parsed:
        # AC17 / TKN20 take 'and' / 'or'; CoverCrypt takes '&&' / '||'
        policy = render(tree, "&&", "||") if scheme == abe.CPABE_SCHEME_COVERCRYPT else render(tree, "and", "or")
        skf = os.path.join(d, "user_%s.key" % pname)
        kg = []
        for i in range(3):
            t, r = ms(lambda: abe.call_generate_secret_key_with_scheme(msk, " ".join(short(a) for a in lv), skf, scheme))
            kg.append(t)
        ok_key = (r == abe.HCPABE_SUCCESS)
        sk = load_b64(skf) if ok_key else None
        enc, dec, ct = [], [], None
        ok = False
        try:
            for i in range(-1, REPS):
                t, ct = ms(lambda: abe.call_hybrid_cpabe_encryptBuffer_with_scheme(pk, secret, policy, scheme))
                if i >= 0: enc.append(t)
            for i in range(-1, REPS):
                t, rec = ms(lambda: abe.call_hybrid_cpabe_decryptBuffer(sk, ct))
                if i >= 0: dec.append(t)
            ok = (rec == secret)
        except Exception as e:
            print("  %s / %s failed: %s" % (sname, pname, e))
        # a user who holds only an attribute that is NOT in the policy must not decrypt
        rejected = ""
        other = [a for a in sorted(used_all) if a not in lv]
        if other and ct is not None:
            bad = os.path.join(d, "bad_%s.key" % pname)
            if abe.call_generate_secret_key_with_scheme(msk, short(other[0]), bad, scheme) == abe.HCPABE_SUCCESS:
                try: rejected = "yes" if abe.call_hybrid_cpabe_decryptBuffer(load_b64(bad), ct) != secret else "NO"
                except Exception: rejected = "yes"
        row = "%s,%s,%.1f,%.1f,%s,%s,%s,%s,%s" % (
            pname, sname, t_setup, med(kg),
            "%.2f" % med(enc) if enc else "", "%.2f" % med(dec) if dec else "",
            len(ct) if ct else "", "yes" if ok else "NO", rejected)
        out.write(row + "\n"); out.flush(); print(row)
print("\nSaved:", out.name)
