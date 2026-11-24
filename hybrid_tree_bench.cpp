// hybrid_tree_bench.cpp -- AC17 (FAME 2017), TKN20 (2020) and CoverCrypt (ML-KEM) of the Hybrid-PQ-CP-ABE-Library,
// run on the policies of a tree file, in C++ against the library's in-memory (buffer) C API.
// Same tree file, same machine, same language as tree_bench (the lattice benchmark): no Python, no file I/O in the timings.
//
// Build (from the folder of this file, library already built in ~/Hybrid-PQ-CP-ABE-Library):
//   g++ -O2 -std=c++17 -I$LIB/include hybrid_tree_bench.cpp -L$LIB/lib -lhybrid-pq-cp-abe -Wl,-rpath,$LIB/lib -o hybrid_tree_bench
// Run:   ./hybrid_tree_bench tree_med32.txt 15          (tree file, repetitions)
// Output: hybrid_results_<time>.csv
//   keygen = key for a user who holds every attribute of the policy (median of 3)
//   enc / dec = median of <reps>, encrypt/decrypt of a 16-byte secret, ONE ciphertext for the whole policy
//   ct_bytes = ciphertext size (includes the library's AES-GCM wrapper and headers)
#include "hybrid_pq_cp_abe/hybrid-pq-cp-abe.h"
#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <iostream>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>
using namespace std;
using Clock = chrono::high_resolution_clock;
static double Ms(Clock::time_point a, Clock::time_point b) { return chrono::duration<double, milli>(b - a).count(); }
static double Median(vector<double> v) { sort(v.begin(), v.end()); return v[v.size() / 2]; }

// ---------------------------------------------------------------- policy tree
struct Node { string kind; string name; int k = 0; vector<Node> kids; };   // kind: leaf, and, or, kof
struct Parser {
  vector<string> t; size_t p = 0;
  explicit Parser(const string& s) {
    size_t i = 0;
    while (i < s.size()) {
      if (isspace((unsigned char)s[i])) { i++; continue; }
      if (s[i] == '(' || s[i] == ')' || s[i] == ',') { t.push_back(string(1, s[i])); i++; continue; }
      size_t b = i;
      while (i < s.size() && (isalnum((unsigned char)s[i]) || s[i] == '_')) i++;
      if (i < s.size() && s[i] == '(' && i > b + 2 && s.compare(i - 2, 2, "of") == 0 && isdigit((unsigned char)s[b])) i++;  // "2of("
      if (i == b) throw runtime_error("bad character in policy: " + s);
      t.push_back(s.substr(b, i - b));
    }
  }
  const string* peek() { return p < t.size() ? &t[p] : nullptr; }
  string eat() { return t.at(p++); }
  Node orExpr() {
    vector<Node> n{andExpr()};
    while (peek() && *peek() == "OR") { eat(); n.push_back(andExpr()); }
    if (n.size() == 1) return n[0];
    Node r; r.kind = "or"; r.kids = n; return r;
  }
  Node andExpr() {
    vector<Node> n{atom()};
    while (peek() && *peek() == "AND") { eat(); n.push_back(atom()); }
    if (n.size() == 1) return n[0];
    Node r; r.kind = "and"; r.kids = n; return r;
  }
  Node atom() {
    string x = eat();
    if (x == "(") { Node n = orExpr(); if (eat() != ")") throw runtime_error("missing )"); return n; }
    if (x.size() > 3 && x.back() == '(') {
      Node r; r.kind = "kof"; r.k = atoi(x.c_str());
      r.kids.push_back(orExpr());
      while (peek() && *peek() == ",") { eat(); r.kids.push_back(orExpr()); }
      if (eat() != ")") throw runtime_error("missing ) in k-of-n");
      return r;
    }
    Node r; r.kind = "leaf"; r.name = x; return r;
  }
};
static string Short(const string& a) {  // library-safe name: lower case letters and digits
  string o; for (char c : a) if (isalnum((unsigned char)c)) o += (char)tolower((unsigned char)c); return o;
}
static void Leaves(const Node& n, set<string>& out) {
  if (n.kind == "leaf") out.insert(n.name); else for (auto& k : n.kids) Leaves(k, out);
}
static string Fold(const string& op, const vector<string>& items) {  // binary, fully parenthesised
  string out = items[0];
  for (size_t i = 1; i < items.size(); i++) out = "(" + out + " " + op + " " + items[i] + ")";
  return out;
}
static string Render(const Node& n, const string& opAnd, const string& opOr) {
  if (n.kind == "leaf") return Short(n.name);
  vector<string> kids;
  for (auto& k : n.kids) kids.push_back(Render(k, opAnd, opOr));
  if (n.kind == "kof") {  // all k-combinations of the children, OR-ed
    vector<string> alts; int m = kids.size();
    vector<int> idx(n.k); 
    function<void(int,int)> rec;
    rec = [&](int start, int depth) {
      if (depth == n.k) { vector<string> c; for (int i : idx) c.push_back(kids[i]); alts.push_back(Fold(opAnd, c)); return; }
      for (int i = start; i < m; i++) { idx[depth] = i; rec(i + 1, depth + 1); }
    };
    rec(0, 0);
    return Fold(opOr, alts);
  }
  return Fold(n.kind == "and" ? opAnd : opOr, kids);
}

// ---------------------------------------------------------------- library helpers
struct Buf {
  unsigned char* p = nullptr; size_t n = 0;
  Buf() {} Buf(const Buf&) = delete; Buf& operator=(const Buf&) = delete;
  ~Buf() { if (p) freeBuffer(p); }
};
static void Check(int rc, const char* what) {
  if (rc != HCPABE_SUCCESS) { cerr << what << " failed: " << getErrorMessage(rc) << " (" << rc << ")\n"; exit(2); }
}

int main(int argc, char** argv) {
  if (argc < 2) { cerr << "usage: " << argv[0] << " <tree file> [reps=15]\n"; return 1; }
  string treeFile = argv[1];
  int reps = argc > 2 ? atoi(argv[2]) : 15;
  int ccMode = getenv("HCPABE_COVERCRYPT_MODE") && string(getenv("HCPABE_COVERCRYPT_MODE")) == "post_quantum" ? 1
             : getenv("HCPABE_COVERCRYPT_MODE") && string(getenv("HCPABE_COVERCRYPT_MODE")) == "classic" ? 0 : 2;
  const char* ccName = ccMode == 1 ? "CoverCrypt_post_quantum" : ccMode == 0 ? "CoverCrypt_classic" : "CoverCrypt_hybridized";

  ifstream in(treeFile);
  if (!in) { cerr << "cannot open " << treeFile << "\n"; return 1; }
  vector<pair<string, Node>> pols; set<string> usedAll;
  string line;
  while (getline(in, line)) {
    size_t h = line.find('#'); if (h != string::npos) line = line.substr(0, h);
    istringstream is(line); string key; is >> key;
    if (key == "POLICY") {
      string name, colon; is >> name >> colon; string rest; getline(is, rest);
      Parser P(rest); Node n = P.orExpr(); Leaves(n, usedAll); pols.push_back({name, n});
    }
  }
  // CoverCrypt: every attribute is a dimension with the values present / absent
  string conf = "{\"dimensions\":[";
  bool first = true;
  for (auto& a : usedAll) { conf += string(first ? "" : ",") + "{\"name\":\"" + Short(a) + "\",\"attributes\":[\"present\",\"absent\"]}"; first = false; }
  conf += "]}";

  string outName = "hybrid_results_" + to_string((long long)time(nullptr)) + ".csv";
  ofstream csv(outName);
  string head = "policy,scheme,setup_ms,keygen_ms,enc_ms,dec_ms,ct_bytes,decrypt_ok,wrong_user_rejected";
  csv << head << "\n"; cout << head << endl;

  mt19937 rng(12345);
  vector<unsigned char> secret(16); for (auto& b : secret) b = rng() & 255;

  struct S { const char* name; CPABEScheme scheme; };
  vector<S> schemes{{"AC17_FAME", CPABE_SCHEME_AC17}, {"TKN20", CPABE_SCHEME_TKN20}, {ccName, CPABE_SCHEME_COVERCRYPT}};
  for (auto& sc : schemes) {
    Buf pk, msk;
    auto t0 = Clock::now();
    if (sc.scheme == CPABE_SCHEME_COVERCRYPT)
      Check(hybrid_cpabe_covercrypt_setupBuffer_from_json(conf.c_str(), ccMode, &pk.p, &pk.n, &msk.p, &msk.n), "covercrypt setup");
    else
      Check(hybrid_cpabe_setupBuffer_with_scheme(&pk.p, &pk.n, &msk.p, &msk.n, sc.scheme), "setup");
    double setupMs = Ms(t0, Clock::now());
    bool cc = sc.scheme == CPABE_SCHEME_COVERCRYPT;

    for (auto& pr : pols) {
      set<string> lv; Leaves(pr.second, lv);
      string policy = cc ? Render(pr.second, "&&", "||") : Render(pr.second, "and", "or");
      string attrs; for (auto& a : lv) attrs += (attrs.empty() ? "" : " ") + Short(a);

      vector<double> kg; Buf sk;
      for (int i = 0; i < 3; i++) {
        Buf tmp; auto a = Clock::now();
        Check(hybrid_cpabe_genkeyBuffer_with_scheme(msk.p, msk.n, attrs.c_str(), &tmp.p, &tmp.n, sc.scheme), "keygen");
        kg.push_back(Ms(a, Clock::now()));
        if (i == 2) { sk.p = tmp.p; sk.n = tmp.n; tmp.p = nullptr; }
      }
      vector<double> te, td; Buf ct; bool ok = false;
      for (int r = -1; r < reps; r++) {
        Buf c; auto a = Clock::now();
        int rc = hybrid_cpabe_encryptBuffer_with_scheme(pk.p, pk.n, secret.data(), secret.size(), policy.c_str(), &c.p, &c.n, sc.scheme);
        double d = Ms(a, Clock::now());
        if (rc != HCPABE_SUCCESS) { cerr << "  " << sc.name << " / " << pr.first << " encrypt failed: " << getErrorMessage(rc) << "\n"; break; }
        if (r >= 0) te.push_back(d);
        if (ct.p) freeBuffer(ct.p);               // keep only the last ciphertext
        ct.p = c.p; ct.n = c.n; c.p = nullptr;
      }
      if (ct.p) {
        for (int r = -1; r < reps; r++) {
          Buf pt; auto a = Clock::now();
          int rc = hybrid_cpabe_decryptBuffer(sk.p, sk.n, ct.p, ct.n, &pt.p, &pt.n);
          double d = Ms(a, Clock::now());
          if (rc != HCPABE_SUCCESS) { cerr << "  " << sc.name << " / " << pr.first << " decrypt failed: " << getErrorMessage(rc) << "\n"; break; }
          if (r >= 0) td.push_back(d);
          ok = (pt.n == secret.size() && memcmp(pt.p, secret.data(), secret.size()) == 0);
        }
      }
      // a user who holds only an attribute that is NOT in the policy must not decrypt
      string rejected;
      for (auto& a : usedAll) if (!lv.count(a) && ct.p) {
        Buf bad; Buf pt;
        if (hybrid_cpabe_genkeyBuffer_with_scheme(msk.p, msk.n, Short(a).c_str(), &bad.p, &bad.n, sc.scheme) == HCPABE_SUCCESS) {
          int rc = hybrid_cpabe_decryptBuffer(bad.p, bad.n, ct.p, ct.n, &pt.p, &pt.n);
          rejected = (rc != HCPABE_SUCCESS || pt.n != secret.size() || memcmp(pt.p, secret.data(), secret.size()) != 0) ? "yes" : "NO";
        }
        break;
      }
      char row[512], encS[32] = "", decS[32] = "";
      if (!te.empty()) snprintf(encS, sizeof encS, "%.3f", Median(te));
      if (!td.empty()) snprintf(decS, sizeof decS, "%.3f", Median(td));
      snprintf(row, sizeof row, "%s,%s,%.1f,%.1f,%s,%s,%zu,%s,%s", pr.first.c_str(), sc.name, setupMs, Median(kg),
               encS, decS, ct.n, ok ? "yes" : "NO", rejected.c_str());
      csv << row << "\n"; cout << row << endl;
    }
  }
  cout << "\nSaved: " << outName << endl;
  return 0;
}
