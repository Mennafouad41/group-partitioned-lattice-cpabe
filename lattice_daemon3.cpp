// lattice_daemon3.cpp
//
// Long-running helper process for the Java framework. It holds the improved
// lattice CP-ABE (group-wise additive sharing, PALISADE) in memory and answers
// one text command per line on stdin. Replies go to stdout, one line each:
// "OK ..." or "ERR ...". Any other line (library noise) must be ignored by the
// client. The scheme code below is copied UNCHANGED from 02_improved_cpabe.cpp.
//
// STATUS: PROTOTYPE. No security proof is claimed.
// v2: setup (public parameters + master trapdoor), user keys and ciphertexts can be
// saved to / loaded from binary files, so they survive a restart of this process.
// The files are written with plain PALISADE-independent code: for every ring element
// the format flag and the N coefficients as 64-bit integers.
//
// Commands:
//   INIT <ell> <G>                  -> OK <setup_us>
//   KEYGEN <keyId> <s>              s = ell chars of '0'/'1' (user has attr i or not)
//                                   -> OK <us>
//   ENC <ctId> <w> <hex32>          w = ell chars of '1','-','0' (needs / must not
//                                   have / don't care); hex32 = 16 bytes to protect
//                                   -> OK <us>
//   DEC <ctId> <keyId>              -> OK <hex32> <us>   (garbage if no access)
//   SAVESETUP <path>               -> OK <bytes>   (public params + master trapdoor: keep private!)
//   LOADSETUP <path>               -> OK <ell> <G>  (replaces INIT; same keys as when saved)
//   SAVEKEY <keyId> <path>         -> OK <bytes>
//   LOADKEY <keyId> <path>         -> OK <bytes>
//   SAVECT <ctId> <path>           -> OK <bytes>   (<bytes> = real ciphertext size)
//   LOADCT <ctId> <path>           -> OK <bytes>
//   QUIT
// Paths must not contain spaces.
//
// Build (one line):
//   g++ -O2 -o lattice_daemon3 -I/usr/local/include/palisade-abe -I/usr/local/include/palisade -I/usr/local/include/palisade/core -L/usr/local/lib lattice_daemon3.cpp -lPALISADEabe -lPALISADEbinfhe -lPALISADEpke -lPALISADEcore

#include "abecontext.h"

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

      key.skA.push_back(MatPtr(new Matrix<Elem>(skA)));
      key.skB.push_back(MatPtr(new Matrix<Elem>(skB)));
      key.skAdef.push_back(MatPtr(new Matrix<Elem>(skAd)));
      key.skBdef.push_back(MatPtr(new Matrix<Elem>(skBd)));
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

// ------------------------------------------------------------ protocol -----
static string ToHex(const vector<int64_t>& bits) {
  static const char* d = "0123456789abcdef";
  string out;
  for (usint b = 0; b < kBits / 8; b++) {
    int v = 0;
    for (int i = 0; i < 8; i++) v = (v << 1) | (bits[b * 8 + i] & 1);
    out += d[v >> 4];
    out += d[v & 15];
  }
  return out;
}

static bool FromHex(const string& h, vector<int64_t>& bits) {
  if (h.size() != kBits / 4) return false;
  bits.assign(kBits, 0);
  for (usint b = 0; b < kBits / 8; b++) {
    int v = 0;
    for (int c = 0; c < 2; c++) {
      char ch = h[2 * b + c];
      int x;
      if (ch >= '0' && ch <= '9') x = ch - '0';
      else if (ch >= 'a' && ch <= 'f') x = 10 + ch - 'a';
      else if (ch >= 'A' && ch <= 'F') x = 10 + ch - 'A';
      else return false;
      v = (v << 4) | x;
    }
    for (int i = 0; i < 8; i++) bits[b * 8 + i] = (v >> (7 - i)) & 1;
  }
  return true;
}


// ------------------------------------------------------- file format -------
#include <fstream>
static void W32(ostream& o, uint32_t v) { o.write(reinterpret_cast<const char*>(&v), 4); }
static void W64(ostream& o, uint64_t v) { o.write(reinterpret_cast<const char*>(&v), 8); }
static uint32_t R32(istream& i) {
  uint32_t v = 0;
  i.read(reinterpret_cast<char*>(&v), 4);
  if (!i) throw std::runtime_error("file too short");
  return v;
}
static uint64_t R64(istream& i) {
  uint64_t v = 0;
  i.read(reinterpret_cast<char*>(&v), 8);
  if (!i) throw std::runtime_error("file too short");
  return v;
}

static void PutElem(ostream& o, Elem e, usint N) {
  uint8_t f = (e.GetFormat() == Format::EVALUATION) ? 1 : 0;
  o.write(reinterpret_cast<const char*>(&f), 1);
  for (usint i = 0; i < N; i++) W64(o, e.at(i).ConvertToInt());
}

template <class EP>
static Elem GetElem(istream& in, EP ep, usint N) {
  uint8_t f = 0;
  in.read(reinterpret_cast<char*>(&f), 1);
  if (!in) throw std::runtime_error("file too short");
  Format fmt = f ? Format::EVALUATION : Format::COEFFICIENT;
  typename Elem::Vector v(N, ep->GetModulus());
  for (usint i = 0; i < N; i++) v[i] = typename Elem::Integer(R64(in));
  Elem e(ep, fmt, true);
  e.SetValues(v, fmt);
  return e;
}

static void PutMat(ostream& o, const Matrix<Elem>& M, usint N) {
  W32(o, M.GetRows());
  W32(o, M.GetCols());
  for (usint i = 0; i < M.GetRows(); i++)
    for (usint j = 0; j < M.GetCols(); j++) PutElem(o, M(i, j), N);
}

template <class EP>
static MatPtr GetMat(istream& in, EP ep, usint N) {
  usint r = R32(in), c = R32(in);
  Matrix<Elem> M(Elem::Allocator(ep, Format::EVALUATION), r, c);
  for (usint i = 0; i < r; i++)
    for (usint j = 0; j < c; j++) M(i, j) = GetElem(in, ep, N);
  return MatPtr(new Matrix<Elem>(M));
}

static const uint32_t kMagicSetup = 0x4C534554, kMagicKey = 0x4C4B4559, kMagicCt = 0x4C435431;

static long long FileSize(const string& p) {
  ifstream f(p, ios::binary | ios::ate);
  return f ? (long long)f.tellg() : -1;
}

static void SaveSetup(ImprovedScheme& S, const string& path) {
  ofstream o(path, ios::binary);
  if (!o) throw std::runtime_error("cannot write " + path);
  W32(o, kMagicSetup);
  W32(o, S.ell);
  W32(o, S.G);
  PutMat(o, *S.A, S.N);
  PutMat(o, *S.Bpos, S.N);
  PutMat(o, *S.Bneg, S.N);
  PutMat(o, *S.Bdef, S.N);
  PutElem(o, S.D, S.N);
  PutMat(o, S.TA->m_r, S.N);
  PutMat(o, S.TA->m_e, S.N);
}

static void SaveKey(ImprovedScheme& S, const ImprovedKey& K, const vector<usint>& s,
                    const string& path) {
  ofstream o(path, ios::binary);
  if (!o) throw std::runtime_error("cannot write " + path);
  W32(o, kMagicKey);
  W32(o, S.G);
  W32(o, s.size());
  for (usint v : s) W32(o, v);
  for (usint g = 0; g < S.G; g++) {
    PutMat(o, *K.skA[g], S.N);
    PutMat(o, *K.skB[g], S.N);
    PutMat(o, *K.skAdef[g], S.N);
    PutMat(o, *K.skBdef[g], S.N);
  }
}

static void SaveCt(ImprovedScheme& S, const ImprovedCt& ct, const vector<int>& w,
                   const string& path) {
  ofstream o(path, ios::binary);
  if (!o) throw std::runtime_error("cannot write " + path);
  W32(o, kMagicCt);
  W32(o, w.size());
  for (int v : w) W32(o, (uint32_t)(v + 1));  // -1,0,1 -> 0,1,2
  W32(o, ct.touched.size());
  for (bool b : ct.touched) W32(o, b ? 1 : 0);
  PutMat(o, *ct.ctW, S.N);
  PutMat(o, *ct.cPos, S.N);
  PutMat(o, *ct.cNeg, S.N);
  PutMat(o, *ct.cDef, S.N);
  PutElem(o, ct.c1, S.N);
}

struct StoredKey { ImprovedKey key; vector<usint> s; };
struct StoredCt  { ImprovedCt ct; vector<int> w; };

int main() {
  unique_ptr<ImprovedScheme> S;
  map<string, StoredKey> keys;
  map<string, StoredCt> cts;
  usint ell = 0, G = 0;

  cout << "READY" << endl;
  string line;
  while (getline(cin, line)) {
    istringstream is(line);
    string cmd;
    is >> cmd;
    if (cmd.empty()) continue;
    try {
      if (cmd == "QUIT") break;

      if (cmd == "INIT") {
        is >> ell >> G;
        if (G < 3 || ell == 0 || ell % G != 0 ||
            CPABEMinRingSizeMap.count(make_pair(HEStd_128_classic, ell)) == 0) {
          cout << "ERR bad ell/G (ell in 6,8,16,20,32; G>=3; G divides ell)" << endl;
          continue;
        }
        pair<usint, usint> rb =
            CPABEMinRingSizeMap.at(make_pair(HEStd_128_classic, ell));
        S.reset(new ImprovedScheme(rb.first, rb.second, ell, G));
        keys.clear();
        cts.clear();
        auto t0 = high_resolution_clock::now();
        S->Setup();
        auto t1 = high_resolution_clock::now();
        cout << "OK " << Us(t0, t1) << endl;
        continue;
      }


      if (cmd == "LOADSETUP") {
        string path;
        is >> path;
        ifstream in(path, ios::binary);
        if (!in) { cout << "ERR cannot open " << path << endl; continue; }
        if (R32(in) != kMagicSetup) { cout << "ERR not a setup file" << endl; continue; }
        ell = R32(in);
        G = R32(in);
        if (G < 3 || ell == 0 || ell % G != 0 ||
            CPABEMinRingSizeMap.count(make_pair(HEStd_128_classic, ell)) == 0) {
          cout << "ERR bad ell/G in file" << endl;
          continue;
        }
        pair<usint, usint> rb = CPABEMinRingSizeMap.at(make_pair(HEStd_128_classic, ell));
        S.reset(new ImprovedScheme(rb.first, rb.second, ell, G));
        auto ep = S->params->GetTrapdoorParams()->GetElemParams();
        S->A = GetMat(in, ep, S->N);
        S->Bpos = GetMat(in, ep, S->N);
        S->Bneg = GetMat(in, ep, S->N);
        S->Bdef = GetMat(in, ep, S->N);
        S->D = GetElem(in, ep, S->N);
        Matrix<Elem> r = *GetMat(in, ep, S->N);
        Matrix<Elem> e = *GetMat(in, ep, S->N);
        S->TA = shared_ptr<RLWETrapdoorPair<Elem> >(new RLWETrapdoorPair<Elem>(r, e));
        keys.clear();
        cts.clear();
        cout << "OK " << ell << " " << G << endl;
        continue;
      }

      if (!S) { cout << "ERR not initialised" << endl; continue; }

      if (cmd == "SAVESETUP") {
        string path;
        is >> path;
        SaveSetup(*S, path);
        cout << "OK " << FileSize(path) << endl;
        continue;
      }
      if (cmd == "SAVEKEY") {
        string id, path;
        is >> id >> path;
        if (!keys.count(id)) { cout << "ERR no such key" << endl; continue; }
        SaveKey(*S, keys.at(id).key, keys.at(id).s, path);
        cout << "OK " << FileSize(path) << endl;
        continue;
      }
      if (cmd == "SAVECT") {
        string id, path;
        is >> id >> path;
        if (!cts.count(id)) { cout << "ERR no such ciphertext" << endl; continue; }
        SaveCt(*S, cts.at(id).ct, cts.at(id).w, path);
        cout << "OK " << FileSize(path) << endl;
        continue;
      }
      if (cmd == "FREECT") {
        string id;
        is >> id;
        cts.erase(id);
        cout << "OK" << endl;
        continue;
      }
      if (cmd == "LOADKEY") {
        string id, path;
        is >> id >> path;
        ifstream in(path, ios::binary);
        if (!in) { cout << "ERR cannot open " << path << endl; continue; }
        if (R32(in) != kMagicKey) { cout << "ERR not a key file" << endl; continue; }
        if (R32(in) != S->G) { cout << "ERR key file has different G" << endl; continue; }
        auto ep = S->params->GetTrapdoorParams()->GetElemParams();
        StoredKey sk;
        usint n = R32(in);
        sk.s.resize(n);
        for (usint i = 0; i < n; i++) sk.s[i] = R32(in);
        for (usint g = 0; g < S->G; g++) {
          sk.key.skA.push_back(GetMat(in, ep, S->N));
          sk.key.skB.push_back(GetMat(in, ep, S->N));
          sk.key.skAdef.push_back(GetMat(in, ep, S->N));
          sk.key.skBdef.push_back(GetMat(in, ep, S->N));
        }
        keys.erase(id);
        keys.emplace(id, sk);
        cout << "OK " << FileSize(path) << endl;
        continue;
      }
      if (cmd == "LOADCT") {
        string id, path;
        is >> id >> path;
        ifstream in(path, ios::binary);
        if (!in) { cout << "ERR cannot open " << path << endl; continue; }
        if (R32(in) != kMagicCt) { cout << "ERR not a ciphertext file" << endl; continue; }
        auto ep = S->params->GetTrapdoorParams()->GetElemParams();
        StoredCt sc;
        usint nw = R32(in);
        sc.w.resize(nw);
        for (usint i = 0; i < nw; i++) sc.w[i] = (int)R32(in) - 1;
        usint nt = R32(in);
        sc.ct.touched.resize(nt);
        for (usint i = 0; i < nt; i++) sc.ct.touched[i] = (R32(in) != 0);
        sc.ct.ctW = GetMat(in, ep, S->N);
        sc.ct.cPos = GetMat(in, ep, S->N);
        sc.ct.cNeg = GetMat(in, ep, S->N);
        sc.ct.cDef = GetMat(in, ep, S->N);
        sc.ct.c1 = GetElem(in, ep, S->N);
        cts.erase(id);
        cts.emplace(id, sc);
        cout << "OK " << FileSize(path) << endl;
        continue;
      }

      if (cmd == "KEYGEN") {
        string id, sbits;
        is >> id >> sbits;
        if (sbits.size() != ell) { cout << "ERR s length" << endl; continue; }
        vector<usint> s(ell);
        for (usint i = 0; i < ell; i++) s[i] = (sbits[i] == '1') ? 1 : 0;
        auto t0 = high_resolution_clock::now();
        ImprovedKey k = S->KeyGen(s);
        auto t1 = high_resolution_clock::now();
        StoredKey sk;
        sk.key = k;
        sk.s = s;
        keys.erase(id);
        keys.emplace(id, sk);
        cout << "OK " << Us(t0, t1) << endl;
        continue;
      }

      if (cmd == "ENC") {
        string id, wstr, hex;
        is >> id >> wstr >> hex;
        if (wstr.size() != ell) { cout << "ERR w length" << endl; continue; }
        vector<int64_t> bits;
        if (!FromHex(hex, bits)) { cout << "ERR hex" << endl; continue; }
        vector<int> w(ell, 0);
        for (usint i = 0; i < ell; i++)
          w[i] = (wstr[i] == '1') ? 1 : (wstr[i] == '-' ? -1 : 0);
        Elem ptxt = S->MakePtxt(bits);
        auto t0 = high_resolution_clock::now();
        ImprovedCt ct = S->Encrypt(w, ptxt);
        auto t1 = high_resolution_clock::now();
        StoredCt sc;
        sc.ct = ct;
        sc.w = w;
        cts.erase(id);
        cts.emplace(id, sc);
        cout << "OK " << Us(t0, t1) << endl;
        continue;
      }

      if (cmd == "DEC") {
        string cid, kid;
        is >> cid >> kid;
        if (!cts.count(cid)) { cout << "ERR no such ciphertext" << endl; continue; }
        if (!keys.count(kid)) { cout << "ERR no such key" << endl; continue; }
        const StoredCt& sc = cts.at(cid);
        const StoredKey& sk = keys.at(kid);
        vector<const ImprovedKey*> use(G, &sk.key);
        auto t0 = high_resolution_clock::now();
        vector<int64_t> out = S->Decrypt(sc.w, sk.s, use, sc.ct);
        auto t1 = high_resolution_clock::now();
        cout << "OK " << ToHex(out) << " " << Us(t0, t1) << endl;
        continue;
      }

      cout << "ERR unknown command" << endl;
    } catch (const std::exception& e) {
      cout << "ERR exception: " << e.what() << endl;
    }
  }
  return 0;
}
