// tree_bench.cpp -- run the BASELINE lattice CP-ABE and the OPTIMIZED (group-wise) one on the policies of tree.txt.
// BASELINE  = the paper's algorithm = the same code with G=1 (one relation over all attributes, no default relation).
// OPTIMIZED = group-wise additive sharing (G groups). Scheme code is copied unchanged from lattice_daemon3.cpp
//             (derived from 02_improved_cpabe.cpp); only the default relation is skipped when G=1.
// PROTOTYPE. No security proof is claimed.
//
// Build (one line, same PALISADE paths as before):
//   g++ -O2 -std=c++17 -o tree_bench -I/usr/local/include/palisade-abe -I/usr/local/include/palisade -I/usr/local/include/palisade/core -L/usr/local/lib tree_bench.cpp -lPALISADEabe -lPALISADEbinfhe -lPALISADEpke -lPALISADEcore
// Run:   ./tree_bench tree.txt 5        (tree file, repetitions; defaults: tree.txt, 5)
#include "abecontext.h"
#include <algorithm>
#include <cctype>
#include <set>
#include <fstream>
#include <functional>
#include <ctime>
#include <cstdio>
#include <stdexcept>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <vector>
#include <sstream>
#include <map>
#include <string>

using namespace lbcrypto;
using namespace std;
using namespace chrono;

typedef NativePoly Elem;
typedef shared_ptr<Matrix<Elem> > MatPtr;

static const usint kBits = 128;

static long long Us(const high_resolution_clock::time_point& a,
                    const high_resolution_clock::time_point& b) {
  return duration_cast<microseconds>(b - a).count();
}

// ---------------------------------------------------------- parameters -----
// Copy of ABEContext::ParamsGenCPABE. CHANGED: `terms` is a parameter. The
// paper uses ell+1. The shared decryption has up to ell+G noise terms, so we use
// terms = ell + G (in our runs log2 q stayed 36 either way).
static shared_ptr<CPABEParams<Elem> > MakeParams(usint ringsize, usint ell,
                                                 usint base, usint terms) {
  double sigma = SIGMA;
  auto qCorrectness = [&](uint32_t n, uint32_t m) -> double {
    return 256 * sigma * SPECTRAL_BOUND(n, m - 2, base) * sqrt(m * n * terms);
  };

  double qPrev = 1e6;
  double q = 0;
  usint k = 0;
  usint m = 0;

  k = floor(log2(qPrev - 1.0) + 1.0);
  m = ceil(k / log2(base)) + 2;
  q = qCorrectness(ringsize, m);

  while (std::abs(q - qPrev) > 0.001 * q) {
    qPrev = q;
    k = floor(log2(qPrev - 1.0) + 1.0);
    m = ceil(k / log2(base)) + 2;
    q = qCorrectness(ringsize, m);
  }
  usint sm = ringsize * 2;
  typename Elem::DggType dgg(sigma);
  typename Elem::DugType dug;
  typename Elem::Integer smodulus;
  typename Elem::Integer srootOfUnity;

  smodulus = FirstPrime<typename Elem::Integer>(floor(log2(q - 1.0)) + 1.0,
                                                2 * ringsize);
  srootOfUnity = RootOfUnity(sm, smodulus);
  dug.SetModulus(smodulus);
  ILParamsImpl<typename Elem::Integer> ilParams(sm, smodulus, srootOfUnity);

  ChineseRemainderTransformFTT<typename Elem::Vector>::PreCompute(
      srootOfUnity, sm, smodulus);
  DiscreteFourierTransform::PreComputeTable(sm);
  EncodingParams eparams(std::make_shared<EncodingParamsImpl>(2));
  auto silparams =
      std::make_shared<ILParamsImpl<typename Elem::Integer> >(ilParams);
  RLWETrapdoorParams<Elem> tparams(silparams, dgg, sigma, base);
  return std::make_shared<CPABEParams<Elem> >(
      std::make_shared<RLWETrapdoorParams<Elem> >(tparams), ell, dug, eparams);
}

// ------------------------------------------------------ data containers ----
// CHANGED: one set of key pieces PER GROUP (paper: one set in total).
struct ImprovedKey {
  vector<MatPtr> skA;     // per group: m x 1        (group relation)
  vector<MatPtr> skB;     // per group: m x per
  vector<MatPtr> skAdef;  // per group: m x 1        (NEW: default relation)
  vector<MatPtr> skBdef;  // per group: m x 1        (NEW)
};

struct ImprovedCt {
  MatPtr ctW, cPos, cNeg;
  MatPtr cDef;  // NEW: one row per group the policy does not mention
  Elem c1;
  vector<bool> touched;  // NEW: which groups the policy mentions
};

// ------------------------------------------------------- the scheme --------
struct ImprovedScheme {
  shared_ptr<CPABEParams<Elem> > params;
  MatPtr A, Bpos, Bneg;
  MatPtr Bdef;  // NEW: G x m public vectors
  Elem D;
  shared_ptr<RLWETrapdoorPair<Elem> > TA;
  usint ell, G, per, m, k, N, radix;

  ImprovedScheme(usint ringsize, usint radix_, usint ell_, usint G_)
      : ell(ell_), G(G_), per(ell_ / G_), radix(radix_) {
    params = MakeParams(ringsize, ell, radix, ell + G);
    N = params->GetTrapdoorParams()->GetN();
    k = params->GetTrapdoorParams()->GetK();
    m = k + 2;
  }

  // ---- 1. Setup: the paper's Setup, plus the Bdef vectors (NEW) ----
  void Setup() {
    typename Elem::DugType& dug = params->GetDUG();
    auto ep = params->GetTrapdoorParams()->GetElemParams();
    auto zero_alloc = Elem::Allocator(ep, Format::COEFFICIENT);

    Matrix<Elem> pubElemBPos(zero_alloc, ell, m);
    Matrix<Elem> pubElemBNeg(zero_alloc, ell, m);
    Elem pubElemD(pubElemBPos(0, 0));

    if (pubElemD.GetFormat() != Format::COEFFICIENT) pubElemD.SwitchFormat();
    pubElemD.SetValues(dug.GenerateVector(N), Format::COEFFICIENT);
    pubElemD.SwitchFormat();

    for (usint i = 0; i < pubElemBPos.GetRows(); i++)
      for (usint j = 0; j < pubElemBPos.GetCols(); j++) {
        if (pubElemBPos(i, j).GetFormat() != Format::COEFFICIENT)
          pubElemBPos(i, j).SwitchFormat();
        pubElemBPos(i, j).SetValues(dug.GenerateVector(N), Format::COEFFICIENT);
        pubElemBPos(i, j).SwitchFormat();
      }
    for (usint i = 0; i < pubElemBNeg.GetRows(); i++)
      for (usint j = 0; j < pubElemBNeg.GetCols(); j++) {
        if (pubElemBNeg(i, j).GetFormat() != Format::COEFFICIENT)
          pubElemBNeg(i, j).SwitchFormat();
        pubElemBNeg(i, j).SetValues(dug.GenerateVector(N), Format::COEFFICIENT);
        pubElemBNeg(i, j).SwitchFormat();
      }
    Bpos = MatPtr(new Matrix<Elem>(pubElemBPos));
    Bneg = MatPtr(new Matrix<Elem>(pubElemBNeg));
    D = pubElemD;

    std::pair<Matrix<Elem>, RLWETrapdoorPair<Elem> > keypair =
        RLWETrapdoorUtility<Elem>::TrapdoorGen(ep, SIGMA, radix, false);
    A = MatPtr(new Matrix<Elem>(keypair.first));
    TA = shared_ptr<RLWETrapdoorPair<Elem> >(
        new RLWETrapdoorPair<Elem>(keypair.second));

    // NEW: G extra uniform public vectors, one per group
    Matrix<Elem> bdef(zero_alloc, G, m);
    for (usint i = 0; i < G; i++)
      for (usint j = 0; j < m; j++) {
        Elem& e = bdef(i, j);
        if (e.GetFormat() != Format::COEFFICIENT) e.SwitchFormat();
        e.SetValues(dug.GenerateVector(N), Format::COEFFICIENT);
        e.SwitchFormat();
      }
    Bdef = MatPtr(new Matrix<Elem>(bdef));
  }

  // ---- 2. KeyGen: CHANGED (random shares + one relation per group) ----
  ImprovedKey KeyGen(const vector<usint>& s) {
    auto ep = params->GetTrapdoorParams()->GetElemParams();
    double sb = SPECTRAL_BOUND(N, k, radix);
    typename Elem::DugType& dug = params->GetDUG();

    // NEW: fresh random shares r_1..r_G of D, different for every user
    vector<Elem> r(G);
    Elem rest = D;
    for (usint g = 0; g + 1 < G; g++) {
      Elem t(dug, ep, Format::COEFFICIENT);
      t.SwitchFormat();
      r[g] = t;
      rest = rest - t;
    }
    r[G - 1] = rest;

    ImprovedKey key;
    for (usint g = 0; g < G; g++) {
      // group relation (the paper's relation, restricted to the group, target r_g)
      Matrix<Elem> skB(Elem::MakeDiscreteGaussianCoefficientAllocator(
                           ep, Format::COEFFICIENT, sb),
                       m, per);
      for (usint j = 0; j < per; j++)
        for (usint i = 0; i < m; i++) skB(i, j).SwitchFormat();

      Elem y(ep, Format::EVALUATION, true);
      for (usint a = 0; a < per; a++) {
        usint idx = g * per + a;
        const Matrix<Elem>& B = (s[idx] == 1) ? *Bpos : *Bneg;
        Elem z = B(idx, 0) * skB(0, a);
        for (usint j = 1; j < m; j++) z += B(idx, j) * skB(j, a);
        y += z;
      }
      y = r[g] - y;  // paper: y = D - ...
      Matrix<Elem> skA = RLWETrapdoorUtility<Elem>::GaussSamp(
          N, k, *A, *TA, y, params->GetTrapdoorParams()->GetDGG(),
          params->GetTrapdoorParams()->GetDGGLargeSigma(), radix);

      key.skA.push_back(MatPtr(new Matrix<Elem>(skA)));
      key.skB.push_back(MatPtr(new Matrix<Elem>(skB)));
      if (G > 1) {  // baseline (G=1) = the paper: no default relation
      // NEW: default relation to the SAME share r_g
      Matrix<Elem> skBd(Elem::MakeDiscreteGaussianCoefficientAllocator(
                            ep, Format::COEFFICIENT, sb),
                        m, 1);
      for (usint i = 0; i < m; i++) skBd(i, 0).SwitchFormat();
      Elem zd = (*Bdef)(g, 0) * skBd(0, 0);
      for (usint j = 1; j < m; j++) zd += (*Bdef)(g, j) * skBd(j, 0);
      Elem yd = r[g] - zd;
      Matrix<Elem> skAd = RLWETrapdoorUtility<Elem>::GaussSamp(
          N, k, *A, *TA, yd, params->GetTrapdoorParams()->GetDGG(),
          params->GetTrapdoorParams()->GetDGGLargeSigma(), radix);

      key.skAdef.push_back(MatPtr(new Matrix<Elem>(skAd)));
      key.skBdef.push_back(MatPtr(new Matrix<Elem>(skBd)));
      }
    }
    return key;
  }

  Elem MakePtxt(const vector<int64_t>& bits) {
    auto ep = params->GetTrapdoorParams()->GetElemParams();
    Plaintext pt = PlaintextFactory::MakePlaintext(
        CoefPacked, ep, params->GetEncodingParams(), bits);
    Elem ptxt = pt->GetElement<Elem>();
    return ptxt;
  }

  // ---- 3. Encrypt: CHANGED for groups the policy does not mention ----
  ImprovedCt Encrypt(const vector<int>& w, Elem ptxt) {
    auto ep = params->GetTrapdoorParams()->GetElemParams();
    const Matrix<Elem>& pubTA = *A;
    const Matrix<Elem>& BPos = *Bpos;
    const Matrix<Elem>& BNeg = *Bneg;

    ImprovedCt ct;
    ct.touched.assign(G, false);
    usint lenW = 0, wild = 0, untouched = 0;
    for (usint g = 0; g < G; g++) {
      usint nz = 0;
      for (usint a = 0; a < per; a++)
        if (w[g * per + a] != 0) nz++;
      if (nz > 0) {
        ct.touched[g] = true;
        lenW += nz;
        wild += per - nz;
      } else {
        untouched++;
      }
    }
    usint ncols = 1 + lenW + 2 * wild + untouched;

    Matrix<Elem> err(Elem::MakeDiscreteGaussianCoefficientAllocator(
                         ep, Format::COEFFICIENT, SIGMA),
                     m, ncols);
    for (usint i = 0; i < m; i++)
      for (usint j = 0; j < ncols; j++) err(i, j).SwitchFormat();

    typename Elem::DugType& dug = params->GetDUG();
    Elem s(dug, ep, Format::COEFFICIENT);
    s.SwitchFormat();

    Matrix<Elem> ctW(Elem::Allocator(ep, Format::EVALUATION), lenW + 1, m);
    Matrix<Elem> cPos(Elem::Allocator(ep, Format::EVALUATION), wild, m);
    Matrix<Elem> cNeg(Elem::Allocator(ep, Format::EVALUATION), wild, m);
    Matrix<Elem> cDef(Elem::Allocator(ep, Format::EVALUATION), untouched, m);

    usint iNoise = 0;
    for (usint j = 0; j < m; j++) ctW(0, j) = pubTA(0, j) * s + err(j, iNoise);
    iNoise++;

    usint iW = 0, iAW = 0, iD = 0;
    for (usint g = 0; g < G; g++) {
      if (ct.touched[g]) {
        // exactly the paper's loop, restricted to the group's positions
        for (usint a = 0; a < per; a++) {
          usint idx = g * per + a;
          if (w[idx] == 1) {
            for (usint j = 0; j < m; j++)
              ctW(iW + 1, j) = BPos(idx, j) * s + err(j, iNoise);
            iNoise++;
            iW++;
          } else if (w[idx] == -1) {
            for (usint j = 0; j < m; j++)
              ctW(iW + 1, j) = BNeg(idx, j) * s + err(j, iNoise);
            iNoise++;
            iW++;
          } else {
            for (usint j = 0; j < m; j++) {
              cPos(iAW, j) = BPos(idx, j) * s + err(j, iNoise);
              cNeg(iAW, j) = BNeg(idx, j) * s + err(j, iNoise + 1);
            }
            iNoise += 2;
            iAW++;
          }
        }
      } else {
        // NEW: ONE vector for the whole group instead of 2*per vectors
        for (usint j = 0; j < m; j++)
          cDef(iD, j) = (*Bdef)(g, j) * s + err(j, iNoise);
        iNoise++;
        iD++;
      }
    }

    // c1: unchanged
    Elem qHalf(ep, Format::COEFFICIENT, true);
    typename Elem::Integer m_q = ep->GetModulus();
    qHalf += (m_q >> 1);
    qHalf.SwitchFormat();
    qHalf.AddILElementOne();

    Elem pubElemD = D;
    Elem err1(ep, Format::COEFFICIENT, true);
    err1.SetValues(params->GetTrapdoorParams()->GetDGG().GenerateVector(N, m_q),
                   Format::COEFFICIENT);
    err1.SwitchFormat();
    if (ptxt.GetFormat() != Format::EVALUATION) ptxt.SwitchFormat();
    if (pubElemD.GetFormat() != Format::EVALUATION) pubElemD.SwitchFormat();
    Elem ctC1(params->GetDUG(), ep, Format::EVALUATION);
    ctC1 = s * pubElemD + err1 + ptxt * qHalf;

    ct.ctW = MatPtr(new Matrix<Elem>(ctW));
    ct.cPos = MatPtr(new Matrix<Elem>(cPos));
    ct.cNeg = MatPtr(new Matrix<Elem>(cNeg));
    ct.cDef = MatPtr(new Matrix<Elem>(cDef));
    ct.c1 = ctC1;
    return ct;
  }

  // ---- 4. Decrypt: CHANGED (sum of per-group results) ----
  // keys[g] = the key used for group g. A normal user passes the same key for
  // all groups; the tests pass keys of different users to check mixing.
  vector<int64_t> Decrypt(const vector<int>& w, const vector<usint>& s,
                          const vector<const ImprovedKey*>& keys,
                          const ImprovedCt& ct) {
    auto ep = params->GetTrapdoorParams()->GetElemParams();
    const Matrix<Elem>& ctW = *ct.ctW;
    const Matrix<Elem>& cPos = *ct.cPos;
    const Matrix<Elem>& cNeg = *ct.cNeg;
    const Matrix<Elem>& cDef = *ct.cDef;

    Elem d(ep, Format::EVALUATION, true);
    usint iW = 0, iAW = 0, iD = 0;
    for (usint g = 0; g < G; g++) {
      const ImprovedKey& K = *keys[g];
      if (ct.touched[g]) {
        const Matrix<Elem>& skA = *K.skA[g];
        const Matrix<Elem>& skB = *K.skB[g];
        for (usint j = 0; j < m; j++) d += ctW(0, j) * skA(j, 0);
        for (usint a = 0; a < per; a++) {
          usint idx = g * per + a;
          if (w[idx] == 1 || w[idx] == -1) {
            for (usint j = 0; j < m; j++) d += ctW(iW + 1, j) * skB(j, a);
            iW++;
          } else {
            if (s[idx] == 1) {
              for (usint j = 0; j < m; j++) d += cPos(iAW, j) * skB(j, a);
            } else {
              for (usint j = 0; j < m; j++) d += cNeg(iAW, j) * skB(j, a);
            }
            iAW++;
          }
        }
      } else {
        // NEW: default relation, 2 products
        const Matrix<Elem>& skAd = *K.skAdef[g];
        const Matrix<Elem>& skBd = *K.skBdef[g];
        for (usint j = 0; j < m; j++) d += ctW(0, j) * skAd(j, 0);
        for (usint j = 0; j < m; j++) d += cDef(iD, j) * skBd(j, 0);
        iD++;
      }
    }

    d = ct.c1 - d;
    d.SwitchFormat();

    typename Elem::Integer m_q = ep->GetModulus();
    typename Elem::Integer dec, threshold = m_q >> 2, qHalf = m_q >> 1;
    vector<int64_t> out(kBits, 0);
    for (usint i = 0; i < N; i++) {
      dec = d.at(i);
      if (dec > qHalf) dec = m_q - dec;
      if (i < kBits) out[i] = (dec > threshold) ? 1 : 0;
    }
    return out;
  }
};

// ------------------------------------------------------------ tree parser ---
// Policy syntax (one policy per POLICY line in tree.txt):
//   expr   := term ( OR term )*
//   term   := atom ( AND atom )*
//   atom   := NAME | ( expr ) | <k>of( expr , expr , ... )
// A policy is expanded into clauses (DNF): every clause is an AND of attributes
// and becomes ONE lattice ciphertext. OR / k-of-n therefore cost several clauses.
typedef std::set<int> TBClause;
typedef std::vector<TBClause> TBDnf;

struct TreeParser {
  std::string src; size_t p = 0;
  std::map<std::string,int>* attrs;
  void ws() { while (p < src.size() && isspace((unsigned char)src[p])) p++; }
  bool isId(char c) { return isalnum((unsigned char)c) || c == '_'; }
  std::string ident() {
    ws(); size_t b = p;
    while (p < src.size() && isId(src[p])) p++;
    return src.substr(b, p - b);
  }
  static TBDnf Or(const TBDnf& a, const TBDnf& b) { TBDnf r = a; r.insert(r.end(), b.begin(), b.end()); return r; }
  static TBDnf And(const TBDnf& a, const TBDnf& b) {
    TBDnf r;
    for (auto& x : a) for (auto& y : b) { TBClause c = x; c.insert(y.begin(), y.end()); r.push_back(c); }
    return r;
  }
  static TBDnf Simplify(TBDnf d) {  // remove duplicates and clauses that contain another clause
    std::sort(d.begin(), d.end(), [](const TBClause& a, const TBClause& b){
      return a.size() != b.size() ? a.size() < b.size() : a < b; });
    d.erase(std::unique(d.begin(), d.end()), d.end());
    TBDnf out;
    for (auto& c : d) {
      bool covered = false;
      for (auto& o : out) if (std::includes(c.begin(), c.end(), o.begin(), o.end())) { covered = true; break; }
      if (!covered) out.push_back(c);
    }
    return out;
  }
  TBDnf expr() {
    TBDnf r = term();
    for (;;) { ws(); size_t save = p; std::string w = ident();
      if (w == "OR") r = Simplify(Or(r, term())); else { p = save; break; } }
    return r;
  }
  TBDnf term() {
    TBDnf r = atom();
    for (;;) { ws(); size_t save = p; std::string w = ident();
      if (w == "AND") r = Simplify(And(r, atom())); else { p = save; break; } }
    return r;
  }
  TBDnf atom() {
    ws();
    if (p < src.size() && src[p] == '(') { p++; TBDnf r = expr(); ws(); if (src[p] != ')') throw std::runtime_error("missing )"); p++; return r; }
    std::string id = ident();
    if (id.empty()) throw std::runtime_error("syntax error near: " + src.substr(p));
    ws();
    if (p < src.size() && src[p] == '(' && id.size() > 2 && id.substr(id.size() - 2) == "of" && isdigit((unsigned char)id[0])) {
      int k = atoi(id.c_str()); p++;
      std::vector<TBDnf> parts;
      for (;;) { parts.push_back(expr()); ws(); if (src[p] == ',') { p++; continue; } if (src[p] == ')') { p++; break; } throw std::runtime_error("bad threshold"); }
      int n = parts.size(); TBDnf all;
      if (k < 1 || k > n) throw std::runtime_error("bad k in " + id);
      std::vector<int> idx(k);
      std::function<void(int,int,TBDnf)> rec = [&](int start, int depth, TBDnf cur) {
        if (depth == k) { all = Or(all, cur); return; }
        for (int i = start; i < n; i++) rec(i + 1, depth + 1, depth == 0 ? parts[i] : And(cur, parts[i]));
      };
      rec(0, 0, TBDnf());
      return Simplify(all);
    }
    auto it = attrs->find(id);
    if (it == attrs->end()) throw std::runtime_error("unknown attribute: " + id);
    TBDnf r; r.push_back(TBClause{it->second}); return r;
  }
};
static vector<int64_t> RandomBits(mt19937& rng) { vector<int64_t> m(kBits); for (usint i = 0; i < kBits; i++) m[i] = rng() & 1; return m; }

// ------------------------------------------------------------ benchmark ----
static double MsMedian(vector<double> v) { sort(v.begin(), v.end()); return v[v.size() / 2]; }

struct TBResult {
  double setup_ms = 0, keygen_ms = 0, enc_ms = 0, dec_ms = 0;
  size_t vectors = 0;  // rows of the ciphertext (each row = m ring elements)
  double ct_kb = 0;
  bool ok = false, neg_ok = false;
  usint touched = 0;
};

// Encrypts every clause of the policy (that is the real encryption cost of the policy),
// decrypts the first clause with a user who holds all attributes of the policy.
static TBResult RunPolicy(ImprovedScheme& S, const TBDnf& dnf, const set<int>& used,
                        int reps, mt19937& rng) {
  usint ell = S.ell;
  TBResult R;
  vector<usint> s(ell, 0);
  for (int a : used) s[a] = 1;
  vector<usint> none(ell, 0);

  vector<double> tk;
  ImprovedKey key;
  for (int i = 0; i < 3; i++) {
    auto t0 = high_resolution_clock::now();
    key = S.KeyGen(s);
    auto t1 = high_resolution_clock::now();
    tk.push_back(Us(t0, t1) / 1000.0);
  }
  R.keygen_ms = MsMedian(tk);
  ImprovedKey keyNone = S.KeyGen(none);

  vector<vector<int>> ws;
  for (auto& c : dnf) {
    vector<int> w(ell, 0);
    for (int a : c) w[a] = 1;
    ws.push_back(w);
  }
  vector<int64_t> msg = RandomBits(rng);
  Elem ptxt = S.MakePtxt(msg);

  vector<double> te, td;
  vector<ImprovedCt> cts;
  for (int r = -1; r < reps; r++) {  // r = -1 is a warm-up
    vector<ImprovedCt> cur;
    auto t0 = high_resolution_clock::now();
    for (auto& w : ws) cur.push_back(S.Encrypt(w, ptxt));
    auto t1 = high_resolution_clock::now();
    if (r >= 0) te.push_back(Us(t0, t1) / 1000.0);
    cts = cur;
  }
  R.enc_ms = MsMedian(te);

  vector<const ImprovedKey*> use(S.G, &key), useNone(S.G, &keyNone);
  vector<int64_t> out;
  for (int r = -1; r < reps; r++) {
    auto t0 = high_resolution_clock::now();
    out = S.Decrypt(ws[0], s, use, cts[0]);
    auto t1 = high_resolution_clock::now();
    if (r >= 0) td.push_back(Us(t0, t1) / 1000.0);
  }
  R.dec_ms = MsMedian(td);
  R.ok = (out == msg);
  vector<int64_t> bad = S.Decrypt(ws[0], none, useNone, cts[0]);
  R.neg_ok = (bad != msg);

  for (auto& c : cts) {
    R.vectors += c.ctW->GetRows() + c.cPos->GetRows() + c.cNeg->GetRows() + c.cDef->GetRows();
  }
  R.ct_kb = ((double)R.vectors * S.m + dnf.size()) * S.N * 8.0 / 1024.0;
  for (bool t : cts[0].touched) if (t) R.touched++;
  return R;
}

int main(int argc, char* argv[]) {
  string treeFile = argc > 1 ? argv[1] : "tree.txt";
  int reps = argc > 2 ? atoi(argv[2]) : 5;
  ifstream in(treeFile);
  if (!in) { cerr << "cannot open " << treeFile << "\n"; return 1; }
  usint ell = 0, G = 0;
  vector<string> names;
  map<string, int> attrs;
  vector<pair<string, string> > policies;
  string line;
  while (getline(in, line)) {
    size_t h = line.find('#'); if (h != string::npos) line = line.substr(0, h);
    istringstream is(line); string key; is >> key;
    if (key == "ELL") is >> ell;
    else if (key == "GROUPS") is >> G;
    else if (key == "ATTRS") { string a; while (is >> a) { attrs[a] = names.size(); names.push_back(a); } }
    else if (key == "POLICY") {
      string name, colon; is >> name >> colon;
      string rest; getline(is, rest);
      policies.push_back(make_pair(name, rest));
    }
  }
  if (names.size() != ell || ell % G != 0 || G < 3 ||
      CPABEMinRingSizeMap.count(make_pair(HEStd_128_classic, ell)) == 0) {
    cerr << "tree.txt: need ELL names in ATTRS (got " << names.size() << "), ELL in 6,8,16,20,32, GROUPS>=3 dividing ELL\n";
    return 1;
  }
  pair<usint, usint> rb = CPABEMinRingSizeMap.at(make_pair(HEStd_128_classic, ell));
  mt19937 rng(12345);

  cout << "Building schemes (ell=" << ell << ", G=" << G << ", reps=" << reps << ") ...\n";
  ImprovedScheme base(rb.first, rb.second, ell, 1);  // G=1 : the paper's algorithm
  ImprovedScheme opt(rb.first, rb.second, ell, G);   // group-wise sharing
  auto a0 = high_resolution_clock::now(); base.Setup(); auto a1 = high_resolution_clock::now();
  opt.Setup(); auto a2 = high_resolution_clock::now();
  cout << "setup_ms baseline=" << Us(a0, a1) / 1000.0 << " optimized=" << Us(a1, a2) / 1000.0 << "\n";

  time_t stamp = time(nullptr);
  string outName = "tree_results_" + to_string((long long)stamp) + ".csv";
  ofstream csv(outName);
  string head = "policy,clauses,attrs_used,groups_in_clause0,scheme,keygen_ms,enc_ms_all_clauses,dec_ms,ct_vectors,ct_KB,decrypt_ok,wrong_user_rejected";
  csv << head << "\n";
  cout << head << "\n";
  double sums[2][4] = {{0}};
  int counted = 0;
  for (auto& pr : policies) {
    TreeParser P; P.src = pr.second; P.attrs = &attrs;
    TBDnf dnf;
    try { dnf = P.expr(); } catch (const exception& e) { cerr << pr.first << ": " << e.what() << "\n"; continue; }
    if (dnf.size() > 64) { cout << pr.first << ",FALLBACK(" << dnf.size() << " clauses > 64 -> pairing CP-ABE)\n"; continue; }
    set<int> used;
    for (auto& c : dnf) used.insert(c.begin(), c.end());
    for (int sch = 0; sch < 2; sch++) {
      ImprovedScheme& S = sch == 0 ? base : opt;
      TBResult R = RunPolicy(S, dnf, used, reps, rng);
      usint g0 = 0;
      { set<int> gs; for (int a : dnf[0]) gs.insert(a / (ell / S.G)); g0 = gs.size(); }
      char buf[512];
      snprintf(buf, sizeof buf, "%s,%zu,%zu,%u,%s,%.1f,%.1f,%.2f,%zu,%.0f,%s,%s",
               pr.first.c_str(), dnf.size(), used.size(), g0, sch == 0 ? "BASELINE" : "OPTIMIZED",
               R.keygen_ms, R.enc_ms, R.dec_ms, R.vectors, R.ct_kb, R.ok ? "yes" : "NO", R.neg_ok ? "yes" : "NO");
      cout << buf << endl;
      csv << buf << "\n";
      sums[sch][0] += R.keygen_ms; sums[sch][1] += R.enc_ms; sums[sch][2] += R.dec_ms; sums[sch][3] += R.ct_kb;
    }
    counted++;
  }
  if (counted) {
    cout << "\nMEAN over " << counted << " policies (baseline -> optimized, speed-up = baseline/optimized)\n";
    const char* nm[4] = {"keygen_ms", "enc_ms", "dec_ms", "ct_KB"};
    for (int i = 0; i < 4; i++)
      printf("  %-10s %10.1f -> %10.1f   x%.2f\n", nm[i], sums[0][i] / counted, sums[1][i] / counted,
             sums[0][i] / sums[1][i]);
  }
  cout << "\nSaved: " << outName << "\n";
  return 0;
}
