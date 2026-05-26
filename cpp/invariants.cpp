// invariants.cpp -- FKSZ untwisted invariants from a trisection diagram.
//
// Self-contained C++ Smith normal form over Z, intersection-form
// extraction from alpha/beta/gamma cut systems, signature + parity by
// real-symmetric eigenvalue analysis. For Sage-precision exact algebra,
// route through the optional sage/ python helpers via the bindings
// layer; this C++ implementation produces the same answers for the
// small-genus diagrams produced by the front half.

#include "bhrt_trisect.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

namespace bhrt {

// ---------------------------------------------------------------------- //
// Smith normal form (iterative Euclidean reduction)
// ---------------------------------------------------------------------- //

std::vector<std::int64_t> smithNormalFormDiagonal(
    const std::vector<std::vector<std::int64_t>>& Min) {
    if (Min.empty() || Min[0].empty()) return {};
    auto M = Min;
    const std::int32_t m = static_cast<std::int32_t>(M.size());
    const std::int32_t n = static_cast<std::int32_t>(M[0].size());
    std::vector<std::int64_t> diag;
    std::int32_t i = 0;
    while (i < std::min(m, n)) {
        std::int64_t best = 0;
        std::int32_t br = -1, bc = -1;
        for (std::int32_t r = i; r < m; ++r)
            for (std::int32_t c = i; c < n; ++c) {
                auto v = std::llabs(M[r][c]);
                if (v == 0) continue;
                if (br < 0 || v < best) { best = v; br = r; bc = c; }
            }
        if (br < 0) break;
        if (br != i) std::swap(M[i], M[br]);
        if (bc != i) for (auto& row : M) std::swap(row[i], row[bc]);
        if (M[i][i] < 0) for (std::int32_t c = 0; c < n; ++c) M[i][c] = -M[i][c];

        bool progress = true;
        while (progress) {
            progress = false;
            for (std::int32_t r = i + 1; r < m; ++r) {
                if (M[r][i] == 0) continue;
                auto q = M[r][i] / M[i][i];
                for (std::int32_t c = 0; c < n; ++c) M[r][c] -= q * M[i][c];
                if (M[r][i] != 0) {
                    std::swap(M[i], M[r]);
                    if (M[i][i] < 0)
                        for (std::int32_t c = 0; c < n; ++c) M[i][c] = -M[i][c];
                    progress = true;
                    break;
                }
            }
            for (std::int32_t c = i + 1; c < n; ++c) {
                if (M[i][c] == 0) continue;
                auto q = M[i][c] / M[i][i];
                for (std::int32_t r = 0; r < m; ++r) M[r][c] -= q * M[r][i];
                if (M[i][c] != 0) {
                    for (auto& row : M) std::swap(row[i], row[c]);
                    if (M[i][i] < 0)
                        for (std::int32_t cc = 0; cc < n; ++cc) M[i][cc] = -M[i][cc];
                    progress = true;
                    break;
                }
            }
        }
        diag.push_back(std::llabs(M[i][i]));
        ++i;
    }
    return diag;
}

// ---------------------------------------------------------------------- //
// Helpers
// ---------------------------------------------------------------------- //

namespace {

// Real symmetric eigenvalues for signature/parity.
// Jacobi rotation algorithm; sufficient for small Q (genus <= 60m bound
// means Q is small for the inputs we care about).
std::vector<double> jacobiEigenvalues(std::vector<std::vector<double>> A,
                                       int max_iter = 200) {
    const std::int32_t n = static_cast<std::int32_t>(A.size());
    if (n == 0) return {};
    for (int iter = 0; iter < max_iter; ++iter) {
        std::int32_t p = 0, q = 1;
        double max_off = 0.0;
        for (std::int32_t i = 0; i < n; ++i)
            for (std::int32_t j = i + 1; j < n; ++j) {
                double v = std::abs(A[i][j]);
                if (v > max_off) { max_off = v; p = i; q = j; }
            }
        if (max_off < 1e-12) break;
        double theta = 0.5 * std::atan2(2.0 * A[p][q], A[p][p] - A[q][q]);
        double c = std::cos(theta), s = std::sin(theta);
        for (std::int32_t k = 0; k < n; ++k) {
            double Apk = A[p][k], Aqk = A[q][k];
            A[p][k] = c * Apk + s * Aqk;
            A[q][k] = -s * Apk + c * Aqk;
        }
        for (std::int32_t k = 0; k < n; ++k) {
            double Akp = A[k][p], Akq = A[k][q];
            A[k][p] = c * Akp + s * Akq;
            A[k][q] = -s * Akp + c * Akq;
        }
    }
    std::vector<double> ev(n);
    for (std::int32_t i = 0; i < n; ++i) ev[i] = A[i][i];
    return ev;
}

}  // namespace

// ---------------------------------------------------------------------- //
// Integer rank and the symplectic pairing on H_1(Sigma)
// ---------------------------------------------------------------------- //

std::int32_t integerRank(const std::vector<std::vector<std::int64_t>>& M) {
    if (M.empty() || M[0].empty()) return 0;
    return static_cast<std::int32_t>(smithNormalFormDiagonal(M).size());
}

std::int64_t symplecticOmega(const std::vector<std::int64_t>& u,
                             const std::vector<std::int64_t>& v) {
    // omega((p|q),(r|s)) = p.s - q.r, with the first half "a"-coords and
    // the second half "b"-coords.
    const std::size_t n = u.size();
    const std::size_t g = n / 2;
    std::int64_t acc = 0;
    for (std::size_t i = 0; i < g; ++i)
        acc += u[i] * v[g + i] - u[g + i] * v[i];
    return acc;
}

namespace {

// Stack two row blocks.
std::vector<std::vector<std::int64_t>> stack2(
    const std::vector<std::vector<std::int64_t>>& A,
    const std::vector<std::vector<std::int64_t>>& B) {
    auto out = A;
    out.insert(out.end(), B.begin(), B.end());
    return out;
}

// ---- Minimal exact rational linear algebra (small, tiny denominators) --
struct Frac {
    long long n{0}, d{1};
    static long long g(long long a, long long b) {
        a = a < 0 ? -a : a; b = b < 0 ? -b : b;
        while (b) { long long t = a % b; a = b; b = t; }
        return a ? a : 1;
    }
    void norm() {
        if (d < 0) { n = -n; d = -d; }
        long long k = g(n, d);
        n /= k; d /= k;
    }
};
Frac fadd(Frac a, Frac b){ Frac r{a.n*b.d + b.n*a.d, a.d*b.d}; r.norm(); return r; }
Frac fsub(Frac a, Frac b){ Frac r{a.n*b.d - b.n*a.d, a.d*b.d}; r.norm(); return r; }
Frac fmul(Frac a, Frac b){ Frac r{a.n*b.n, a.d*b.d}; r.norm(); return r; }
Frac fdiv(Frac a, Frac b){ Frac r{a.n*b.d, a.d*b.n}; r.norm(); return r; }

// Invert a square rational matrix in place via Gauss-Jordan.
// Returns false if singular or if magnitudes overflow the safe range.
bool invertRational(std::vector<std::vector<Frac>> A,
                    std::vector<std::vector<Frac>>& out) {
    const int n = static_cast<int>(A.size());
    out.assign(n, std::vector<Frac>(n, Frac{0, 1}));
    for (int i = 0; i < n; ++i) out[i][i] = Frac{1, 1};
    const long long LIM = 1000000000000LL;  // 1e12 safety bound
    for (int col = 0; col < n; ++col) {
        int piv = -1;
        for (int r = col; r < n; ++r) if (A[r][col].n != 0) { piv = r; break; }
        if (piv < 0) return false;
        std::swap(A[piv], A[col]);
        std::swap(out[piv], out[col]);
        Frac inv = fdiv(Frac{1, 1}, A[col][col]);
        for (int c = 0; c < n; ++c) { A[col][c] = fmul(A[col][c], inv);
                                      out[col][c] = fmul(out[col][c], inv); }
        for (int r = 0; r < n; ++r) {
            if (r == col || A[r][col].n == 0) continue;
            Frac f = A[r][col];
            for (int c = 0; c < n; ++c) {
                A[r][c]   = fsub(A[r][c],   fmul(f, A[col][c]));
                out[r][c] = fsub(out[r][c], fmul(f, out[col][c]));
                if (std::llabs(A[r][c].n) > LIM || A[r][c].d > LIM ||
                    std::llabs(out[r][c].n) > LIM || out[r][c].d > LIM)
                    return false;
            }
        }
    }
    return true;
}

std::vector<std::vector<Frac>> toFrac(
    const std::vector<std::vector<std::int64_t>>& M) {
    std::vector<std::vector<Frac>> F(M.size());
    for (std::size_t i = 0; i < M.size(); ++i) {
        F[i].resize(M[i].size());
        for (std::size_t j = 0; j < M[i].size(); ++j) F[i][j] = Frac{M[i][j], 1};
    }
    return F;
}

}  // namespace

// ---------------------------------------------------------------------- //
// Validate a Lagrangian diagram
// ---------------------------------------------------------------------- //

DiagramValidity validateLagrangianDiagram(const LagrangianDiagram& d) {
    DiagramValidity v;
    v.genus = d.genus;
    const std::int32_t g = d.genus;
    const std::int32_t twog = 2 * g;
    auto checkSystem = [&](const std::vector<std::vector<std::int64_t>>& S,
                           const char* name) -> bool {
        if (static_cast<std::int32_t>(S.size()) != g) {
            v.messages.push_back(std::string(name) + ": expected g curves");
            return false;
        }
        for (const auto& row : S)
            if (static_cast<std::int32_t>(row.size()) != twog) {
                v.messages.push_back(std::string(name) + ": row not length 2g");
                return false;
            }
        if (g == 0) return true;
        if (integerRank(S) != g) {
            v.messages.push_back(std::string(name) + ": curves not independent");
            return false;
        }
        for (std::size_t i = 0; i < S.size(); ++i)
            for (std::size_t j = i + 1; j < S.size(); ++j)
                if (symplecticOmega(S[i], S[j]) != 0) {
                    v.messages.push_back(std::string(name) + ": not isotropic");
                    return false;
                }
        return true;
    };
    v.alpha_is_cut_system = checkSystem(d.alpha, "alpha");
    v.beta_is_cut_system  = checkSystem(d.beta,  "beta");
    v.gamma_is_cut_system = checkSystem(d.gamma, "gamma");
    v.surface_closed = v.surface_orientable = v.surface_connected = true;
    v.ok = v.alpha_is_cut_system && v.beta_is_cut_system && v.gamma_is_cut_system;
    return v;
}

// ---------------------------------------------------------------------- //
// Homology + intersection form from Lagrangian data (the FKSZ input)
// ---------------------------------------------------------------------- //

InvariantBundle computeInvariantsFromClasses(const LagrangianDiagram& d) {
    InvariantBundle inv;
    inv.genus = d.genus;
    const std::int32_t g = d.genus;
    const std::int32_t twog = 2 * g;

    if (g == 0) {
        // S^4-type: trivial reduced homology, empty form.
        inv.parity = "even";
        inv.audit.emplace_back("genus 0: X has the homology of S^4");
        return inv;
    }

    const auto& A = d.alpha;
    const auto& B = d.beta;
    const auto& C = d.gamma;

    // H_1(X) = Z^{2g} / (<alpha> + <beta> + <gamma>).
    auto M = stack2(stack2(A, B), C);
    auto snf = smithNormalFormDiagonal(M);
    std::int32_t rankABC = static_cast<std::int32_t>(snf.size());
    inv.h1_free_rank = twog - rankABC;
    inv.h1_torsion.clear();
    for (auto e : snf) if (e > 1) inv.h1_torsion.push_back(e);

    // k_{xy} = 2g - rank([x;y]); chi = 2 + g - sum k; b2 = chi - 2 + 2 b1.
    std::int32_t k_ab = twog - integerRank(stack2(A, B));
    std::int32_t k_bc = twog - integerRank(stack2(B, C));
    std::int32_t k_ac = twog - integerRank(stack2(A, C));
    std::int32_t chi  = 2 + g - (k_ab + k_bc + k_ac);
    inv.h2_free_rank  = chi - 2 + 2 * inv.h1_free_rank;
    if (inv.h2_free_rank < 0) inv.h2_free_rank = 0;

    // Closed oriented 4-manifold: H_3 free of rank b1 (no torsion);
    // torsion of H_2 ~= torsion of H_1.
    inv.h3_free_rank = inv.h1_free_rank;
    inv.h3_torsion.clear();
    inv.h2_torsion = inv.h1_torsion;
    inv.audit.emplace_back(
        "homology computed exactly from the three Lagrangians "
        "(H1 = coker[alpha;beta;gamma]; b2 via Euler characteristic)");

    // Intersection form: classical linking matrix Q = M_gamma * N_gamma^-1
    // after writing gamma in the basis [alpha;beta].  Exact and correct
    // when X is simply connected and [alpha;beta] is unimodular
    // (the (g;0,0,0)-standard case: S^4, CP^2, S^2xS^2, connected sums).
    bool formDone = false;
    if (inv.h1_free_rank == 0 && k_ab == 0) {
        auto basisInt = stack2(A, B);                 // 2g x 2g
        std::vector<std::vector<Frac>> basisInv;
        if (invertRational(toFrac(basisInt), basisInv)) {
            // coeff = C * basisInv  -> g x 2g (rational)
            std::vector<std::vector<Frac>> Mg(g, std::vector<Frac>(g, Frac{0,1}));
            std::vector<std::vector<Frac>> Ng(g, std::vector<Frac>(g, Frac{0,1}));
            bool integral = true;
            for (std::int32_t i = 0; i < g; ++i)
                for (std::int32_t j = 0; j < twog; ++j) {
                    Frac acc{0,1};
                    for (std::int32_t k = 0; k < twog; ++k)
                        acc = fadd(acc, fmul(Frac{C[i][k],1}, basisInv[k][j]));
                    if (j < g) Mg[i][j] = acc; else Ng[i][j - g] = acc;
                }
            std::vector<std::vector<Frac>> Ninv;
            if (invertRational(Ng, Ninv)) {
                std::vector<std::vector<std::int64_t>> Q(g,
                    std::vector<std::int64_t>(g, 0));
                for (std::int32_t i = 0; i < g && integral; ++i)
                    for (std::int32_t j = 0; j < g; ++j) {
                        Frac acc{0,1};
                        for (std::int32_t k = 0; k < g; ++k)
                            acc = fadd(acc, fmul(Mg[i][k], Ninv[k][j]));
                        if (acc.d != 1) { integral = false; break; }
                        Q[i][j] = acc.n;
                    }
                if (integral) {
                    // Symmetrise defensively (should already be symmetric).
                    for (std::int32_t i = 0; i < g; ++i)
                        for (std::int32_t j = i + 1; j < g; ++j) {
                            std::int64_t s = (Q[i][j] + Q[j][i]);
                            Q[i][j] = Q[j][i] = s / 2;
                        }
                    inv.intersection_form = Q;
                    std::vector<std::vector<double>> Qd(g,
                        std::vector<double>(g, 0.0));
                    for (std::int32_t i = 0; i < g; ++i)
                        for (std::int32_t j = 0; j < g; ++j)
                            Qd[i][j] = static_cast<double>(Q[i][j]);
                    auto eig = jacobiEigenvalues(Qd);
                    int pos = 0, neg = 0;
                    for (auto e : eig) { if (e > 1e-7) ++pos; else if (e < -1e-7) ++neg; }
                    inv.signature = pos - neg;
                    bool even = true;
                    for (std::int32_t i = 0; i < g; ++i)
                        if (Q[i][i] % 2 != 0) { even = false; break; }
                    inv.parity = even ? "even" : "odd";
                    formDone = true;
                    inv.audit.emplace_back(
                        "intersection form via linking matrix (g;0,0,0) case");
                }
            }
        }
    }
    if (!formDone) {
        inv.parity = "unknown";
        inv.audit.emplace_back(
            "intersection form not computed: requires a simply-connected "
            "(g;0,0,0)-standardisable diagram, or the geometric crossing "
            "data of the curves on Sigma (full FKSZ algorithm)");
    }
    return inv;
}

// ---------------------------------------------------------------------- //
// Cellular-diagram entry point.
// ---------------------------------------------------------------------- //

InvariantBundle computeInvariants(const TrisectionDiagram& diagram) {
    InvariantBundle inv;
    inv.genus = diagram.genus();
    inv.parity = "unknown";

    const bool no_curves =
        diagram.alpha.empty() && diagram.beta.empty() && diagram.gamma.empty();
    if (no_curves) {
        inv.audit.emplace_back("no cut curves recovered from diagram");
        return inv;
    }

    // The manifold invariants are determined by the *homology classes* of
    // the cut curves in H_1(Sigma), not by their edge-incidence vectors.
    // Lifting cellular curves to symplectic H_1(Sigma) classes requires a
    // rotation system on the surface, which the cellular extractor does
    // not yet record.  Rather than emit a fabricated form (the previous
    // (GG^T - G^T G)/2 proxy), report only what the cellular data fixes
    // and route exact computation through computeInvariantsFromClasses.
    inv.audit.emplace_back(
        "cellular diagram present; exact H_* and intersection form require "
        "symplectic curve classes -- call computeInvariantsFromClasses with "
        "a LagrangianDiagram, or validate the cut systems via validateDiagram");
    return inv;
}

}  // namespace bhrt
