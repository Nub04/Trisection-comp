// search_cpu.cpp -- beam search baseline (C++ core).
//
// Enumerates Pachner/edge-collapse candidates, scores them with the
// CPU reference scorer, commits the top-K via the safe executor
// (Regina-backed when BHRT_HAS_REGINA, conservative built-in otherwise),
// and maintains a Pareto front of (pentachora, trisection_genus). The
// optional CUDA path replaces scoreCandidateCPU with the kernel in
// cuda/score_moves.cu; both code paths produce identical scores.

#include "bhrt_trisect.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <utility>

#if BHRT_HAS_CUDA
// Host launcher defined in cuda/score_moves.cu. Returns 0 on success; non-zero
// if no GPU is available (caller then falls back to the CPU scorer).
extern "C" int bhrt_score_moves_host(const float* records, int n_records,
                                     int n_cols, const float* weights,
                                     float* scores);
#endif

// When built with AddressSanitizer, bake in LeakSanitizer suppressions for the
// NVIDIA CUDA driver/runtime. Those libraries hold internal allocations until
// process teardown that LSan reports as leaks; they are not bugs in this
// project. This keeps leak detection ON for bhrt's own code while ignoring the
// driver, with no need to set LSAN_OPTIONS.
#if defined(__SANITIZE_ADDRESS__)
extern "C" const char* __lsan_default_suppressions() {
    return "leak:libcuda\n"
           "leak:libcudart\n"
           "leak:libcuda.so\n";
}
#endif

namespace bhrt {

// ---------------------------------------------------------------------- //
// Move enumeration
// ---------------------------------------------------------------------- //

std::vector<MoveCandidate> enumerateMoves(const Triangulation& T) {
    std::vector<MoveCandidate> out;
    // 2-4 candidates: every tetrahedron equivalence class (the face the
    // 4D move acts on).
    std::set<std::int32_t> seen_tet;
    for (std::int32_t p = 0; p < (std::int32_t)T.size(); ++p)
        for (std::uint8_t f = 0; f < 5; ++f) {
            auto ttid = T.tetraOf(p, f);
            if (!seen_tet.insert(ttid).second) continue;
            out.push_back({MoveType::Pachner_2_3, {ttid}, +2});
        }
    // 3-3 candidates: every triangle equivalence class.
    std::set<std::int32_t> seen_tri;
    for (std::int32_t p = 0; p < (std::int32_t)T.size(); ++p)
        for (std::uint8_t a = 0; a < 5; ++a)
            for (std::uint8_t b = a + 1; b < 5; ++b)
                for (std::uint8_t c = b + 1; c < 5; ++c) {
                    auto tid = T.triangleOf(p, a, b, c);
                    if (!seen_tri.insert(tid).second) continue;
                    out.push_back({MoveType::Pachner_3_3, {tid}, 0});
                }
    // 4-4 / 4-2 / edge-collapse candidates: every edge equivalence class.
    std::set<std::int32_t> seen_edge;
    for (std::int32_t p = 0; p < (std::int32_t)T.size(); ++p)
        for (std::uint8_t a = 0; a < 5; ++a)
            for (std::uint8_t b = a + 1; b < 5; ++b) {
                auto eid = T.edgeOf(p, a, b);
                if (!seen_edge.insert(eid).second) continue;
                out.push_back({MoveType::Pachner_4_4, {eid}, 0});
                out.push_back({MoveType::Pachner_4_2, {eid}, -2});
                out.push_back({MoveType::EdgeCollapse, {eid}, -1});
            }
    // 1-5 expansion at every pentachoron.
    for (std::int32_t p = 0; p < (std::int32_t)T.size(); ++p)
        out.push_back({MoveType::Pachner_1_5, {p}, +4});
    // 5-1 collapse candidate at every degree-5 vertex (locator = vertex id).
    // Pentachoron-degree of each global vertex:
    std::map<std::int32_t, int> vdeg;
    for (std::int32_t p = 0; p < (std::int32_t)T.size(); ++p) {
        std::set<std::int32_t> here;
        for (std::uint8_t v = 0; v < 5; ++v) here.insert(T.vertexOf(p, v));
        for (auto gv : here) ++vdeg[gv];
    }
    for (auto& kv : vdeg)
        if (kv.second == 5)
            out.push_back({MoveType::Pachner_5_1, {kv.first}, -4});
    return out;
}

// ---------------------------------------------------------------------- //
// Built-in combinatorial move executor.
//
// Implements the two Pachner moves that are well-defined *purely
// combinatorially* on a closed 4-manifold triangulation:
//
//   * 1->5  (oneFiveMove):  always legal.  Cone a fresh central vertex
//            over the five tetrahedral facets of one pentachoron.  The
//            external gluings of the old pentachoron are inherited by the
//            five new "cone" pentachora; internal gluings follow the
//            transposition pattern.  This is the standard stellar
//            subdivision of a pentachoron and never changes the PL type.
//
//   * 5->1  (fiveOneMove):  the exact inverse.  Detected at a degree-5
//            internal vertex and *self-verified*: the candidate result U
//            is accepted only if it is valid and oneFiveMove(U, .)
//            reproduces the original isomorphism signature.  A wrong
//            reconstruction therefore degrades to "skip", never to a
//            corrupted triangulation.
//
// The lateral 4D bistellar moves (2-4, 3-3, 4-4, 4-2) and edge collapses
// require committing through a topology-aware executor; when
// BHRT_HAS_REGINA is defined the Regina-backed path (applyMoveRegina)
// supplies all of them.  Without Regina those candidates return nullopt
// and are skipped by the beam search.
// ---------------------------------------------------------------------- //

namespace {

// Cone a new central vertex over pentachoron p (1->5 stellar subdivision).
// Always legal; result has size()+4 pentachora.
Triangulation oneFiveMove(const Triangulation& T, std::int32_t p) {
    const std::int32_t n = static_cast<std::int32_t>(T.size());
    // Compact the kept pentachora (everything except p) to the front.
    std::vector<std::int32_t> old2new(n, -1);
    std::int32_t idx = 0;
    for (std::int32_t q = 0; q < n; ++q)
        if (q != p) old2new[q] = idx++;
    const std::int32_t kept = idx;            // = n - 1
    auto cone = [&](int i) { return kept + i; };  // cone pentachoron P_i

    Triangulation U(T.label());
    U.addPentachora(static_cast<std::size_t>(kept) + 5);

    const auto& pents = T.pentachora();

    // (1) Kept<->kept gluings (everything not touching p), realised once.
    for (std::int32_t q = 0; q < n; ++q) {
        if (q == p) continue;
        for (std::uint8_t f = 0; f < 5; ++f) {
            const auto& g = pents[q].gluings[f];
            if (!g.has_value()) continue;
            const std::int32_t dst = g->dst_pent;
            if (dst == p) continue;            // handled from the cone side
            if (dst == q) { if (f < g->dst_facet) U.glue(old2new[q], f, old2new[q], g->perm); }
            else if (q < dst) U.glue(old2new[q], f, old2new[dst], g->perm);
        }
    }
    // (2) Internal cone<->cone gluings: P_i facet j <-> P_j facet i,
    //     permutation = transposition (i j).
    for (int i = 0; i < 5; ++i)
        for (int j = i + 1; j < 5; ++j) {
            Perm5 tr = perm_identity;
            tr[i] = static_cast<std::uint8_t>(j);
            tr[j] = static_cast<std::uint8_t>(i);
            U.glue(cone(i), static_cast<std::uint8_t>(j), cone(j), tr);
        }
    // (3) Cone external facet i inherits p's original facet-i gluing.
    for (std::uint8_t i = 0; i < 5; ++i) {
        const auto& g = pents[p].gluings[i];
        if (!g.has_value()) continue;          // boundary facet stays boundary
        const std::int32_t dst = g->dst_pent;
        const std::uint8_t df = g->dst_facet;
        if (dst != p) {
            U.glue(cone(i), i, old2new[dst], g->perm);
        } else {                               // self-gluing of p
            if (i < df)       U.glue(cone(i), i, cone(df), g->perm);
            else if (i == df) U.glue(cone(i), i, cone(i),  g->perm);
        }
    }
    return U;
}

// Pentachora incident to a global vertex.
std::vector<std::int32_t> starOfVertex(const Triangulation& T, std::int32_t v) {
    std::vector<std::int32_t> star;
    for (std::int32_t q = 0; q < (std::int32_t)T.size(); ++q)
        for (std::uint8_t loc = 0; loc < 5; ++loc)
            if (T.vertexOf(q, loc) == v) { star.push_back(q); break; }
    return star;
}

// Attempt the inverse 5->1 collapse at a degree-5 internal vertex c.
// Reconstructs the merged pentachoron, then *verifies* by re-expanding;
// returns nullopt unless the round-trip reproduces T exactly.
std::optional<Triangulation> fiveOneMove(const Triangulation& T, std::int32_t c) {
    const std::int32_t n = static_cast<std::int32_t>(T.size());
    auto star = starOfVertex(T, c);
    if (star.size() != 5) return std::nullopt;
    std::set<std::int32_t> starSet(star.begin(), star.end());

    // Local index of c, and the global vertices opposite it, per star pent.
    std::map<std::int32_t, std::uint8_t> cLocal;
    std::set<std::int32_t> Wset;
    for (auto q : star) {
        std::uint8_t cl = 255;
        for (std::uint8_t loc = 0; loc < 5; ++loc)
            if (T.vertexOf(q, loc) == c) { cl = loc; break; }
        if (cl == 255) return std::nullopt;
        cLocal[q] = cl;
        for (std::uint8_t loc = 0; loc < 5; ++loc)
            if (loc != cl) Wset.insert(T.vertexOf(q, loc));
    }
    if (Wset.size() != 5) return std::nullopt;          // not a 1-5 pattern
    std::vector<std::int32_t> W(Wset.begin(), Wset.end());
    auto Windex = [&](std::int32_t gv) -> int {
        for (int k = 0; k < 5; ++k) if (W[k] == gv) return k;
        return -1;
    };

    // Compact kept (non-star) pentachora; the merged pentachoron is last.
    std::vector<std::int32_t> old2new(n, -1);
    std::int32_t idx = 0;
    for (std::int32_t q = 0; q < n; ++q)
        if (!starSet.count(q)) old2new[q] = idx++;
    const std::int32_t kept = idx;
    const std::int32_t Pnew = kept;                     // merged pentachoron

    Triangulation U(T.label());
    U.addPentachora(static_cast<std::size_t>(kept) + 1);
    const auto& pents = T.pentachora();

    // Kept<->kept gluings.
    for (std::int32_t q = 0; q < n; ++q) {
        if (starSet.count(q)) continue;
        for (std::uint8_t f = 0; f < 5; ++f) {
            const auto& g = pents[q].gluings[f];
            if (!g.has_value()) continue;
            const std::int32_t dst = g->dst_pent;
            if (starSet.count(dst)) continue;           // handled from Pnew side
            if (dst == q) { if (f < g->dst_facet) U.glue(old2new[q], f, old2new[q], g->perm); }
            else if (q < dst) U.glue(old2new[q], f, old2new[dst], g->perm);
        }
    }
    // Pnew external facets = the star pentachora's facets opposite c.
    for (auto q : star) {
        const std::uint8_t cl = cLocal[q];
        const auto& g = pents[q].gluings[cl];
        if (!g.has_value()) return std::nullopt;        // c not internal -> bail
        const std::int32_t dst = g->dst_pent;
        if (starSet.count(dst)) return std::nullopt;    // outer facet must be external
        // Relabel q-locals -> Pnew-locals (= Windex), apex c -> omitted slot m.
        Perm5 Pperm{};
        bool used[5] = {false,false,false,false,false};
        for (std::uint8_t v = 0; v < 5; ++v) {
            if (v == cl) continue;
            int k = Windex(T.vertexOf(q, v));
            if (k < 0) return std::nullopt;
            Pperm[k] = g->perm[v];
            used[k] = true;
        }
        int m = -1;
        for (int k = 0; k < 5; ++k) if (!used[k]) { m = k; break; }
        if (m < 0) return std::nullopt;
        Pperm[m] = g->perm[cl];                         // apex -> dst facet
        U.glue(Pnew, static_cast<std::uint8_t>(m), old2new[dst], Pperm);
    }

    // Self-verify: U must re-expand (1->5 at Pnew) back to T.
    if (!U.isValid()) return std::nullopt;
    if (U.isClosed() != T.isClosed()) return std::nullopt;
    if (oneFiveMove(U, Pnew).isoSig() != T.isoSig()) return std::nullopt;
    return U;
}

}  // namespace

#if BHRT_HAS_REGINA
// Defined in regina_bridge.cpp: lateral 2-4 / 3-3 / 4-4 / 4-2 moves and
// edge collapses via Regina.
std::optional<Triangulation> applyMoveRegina(const Triangulation& T,
                                             const MoveCandidate& c);
#endif

std::optional<Triangulation> applyMove(const Triangulation& T,
                                        const MoveCandidate& c) {
    switch (c.move) {
        case MoveType::Pachner_1_5:
            if (c.locator.empty()) return std::nullopt;
            if (c.locator[0] < 0 || c.locator[0] >= (std::int32_t)T.size())
                return std::nullopt;
            return oneFiveMove(T, c.locator[0]);
        case MoveType::Pachner_5_1:
            if (c.locator.empty()) return std::nullopt;
            return fiveOneMove(T, c.locator[0]);
        default:
#if BHRT_HAS_REGINA
            return applyMoveRegina(T, c);
#else
            // Lateral moves need the topology-aware executor; skip safely.
            return std::nullopt;
#endif
    }
}

// ---------------------------------------------------------------------- //
// Flatten / score
// ---------------------------------------------------------------------- //

std::array<float, 5> flattenCandidate(const SearchState& parent,
                                       const MoveCandidate& c) {
    return {
        static_cast<float>(static_cast<int>(c.move)),
        c.locator.empty() ? 0.0f : static_cast<float>(c.locator[0]),
        static_cast<float>(c.delta_pent_estimate),
        static_cast<float>(parent.pentachora),
        parent.ts_feasible ? 1.0f : 0.0f,
    };
}

float scoreCandidateCPU(const std::array<float, 5>& r,
                         const std::array<float, 5>& w) {
    float delta = -r[2];
    float size_proxy = -r[3] / 100.0f;
    float excess = r[2] > 0.0f ? r[2] : 0.0f;
    return w[0] * delta + w[1] * r[4] + w[2] * size_proxy + w[3] * excess;
}

// ---------------------------------------------------------------------- //
// Pareto front
// ---------------------------------------------------------------------- //

void ParetoFront::consider(const SearchState& s) {
    auto key = [](const SearchState& x) -> std::pair<int, int> {
        return {x.pentachora, x.best_genus.value_or(1000000000)};
    };
    auto s_k = key(s);
    std::vector<SearchState> survivors;
    bool dominated = false;
    for (auto& x : items) {
        auto x_k = key(x);
        if (x_k.first <= s_k.first && x_k.second <= s_k.second && x_k != s_k)
            dominated = true;
        if (!(s_k.first <= x_k.first && s_k.second <= x_k.second && s_k != x_k))
            survivors.push_back(x);
    }
    if (!dominated) survivors.push_back(s);
    items = std::move(survivors);
}

// ---------------------------------------------------------------------- //
// Beam search
// ---------------------------------------------------------------------- //

static SearchState makeSeed(const Triangulation& T) {
    SearchState s;
    s.triangulation = T;
    s.pentachora = static_cast<int>(T.size());
    s.isosig = T.edgeDegreeIsoSig();
    auto ts = enumerateTsTricolourings(T, 10000, true);
    s.ts_feasible = !ts.empty();
    if (s.ts_feasible) {
        auto d = extractDiagram(T, ts.front());
        s.best_genus = d.genus();
    }
    return s;
}

ParetoFront beamSearch(const Triangulation& start, const SearchConfig& cfg) {
    // GPU scorer is opt-in: SearchConfig flag, or the BHRT_USE_GPU env var.
    bool use_gpu = cfg.use_gpu_scorer;
#if BHRT_HAS_CUDA
    if (const char* e = std::getenv("BHRT_USE_GPU"))
        if (e[0] && e[0] != '0') use_gpu = true;
    {
        static bool once = false;
        if (!once) {
            once = true;
            std::fprintf(stderr, "[bhrt] CUDA-enabled build; GPU scorer %s\n",
                         use_gpu ? "ON" : "off (set BHRT_USE_GPU=1 to enable)");
        }
    }
#else
    (void)use_gpu;
#endif
    ParetoFront pareto;
    auto seed = makeSeed(start);
    pareto.consider(seed);

    std::vector<SearchState> beam{seed};
    std::set<std::string> visited{seed.isosig};

    std::mt19937_64 rng(cfg.seed);
    auto t0 = std::chrono::steady_clock::now();

    const std::array<float, 5> weights{1.0f, 0.5f, 0.3f, -0.2f, -1.0f};

    for (int iter = 0; iter < cfg.max_iterations; ++iter) {
        double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();
        if (elapsed > cfg.time_limit_seconds) break;

        std::vector<std::pair<SearchState*, MoveCandidate>> pairs;
        for (auto& s : beam) {
            for (auto& c : enumerateMoves(s.triangulation))
                pairs.emplace_back(&s, c);
        }
        if (pairs.empty()) break;

        std::vector<float> scores(pairs.size());
        bool scored_on_gpu = false;
#if BHRT_HAS_CUDA
        if (use_gpu && !pairs.empty()) {
            std::vector<float> recs(pairs.size() * 5);
            for (std::size_t i = 0; i < pairs.size(); ++i) {
                auto r = flattenCandidate(*pairs[i].first, pairs[i].second);
                for (int k = 0; k < 5; ++k) recs[i * 5 + k] = r[k];
            }
            // Same scoring function as the CPU path; only the arithmetic is
            // batched (GPU may fuse a*b+c, so allow the last ulp to differ).
            if (bhrt_score_moves_host(recs.data(), static_cast<int>(pairs.size()),
                                      5, weights.data(), scores.data()) == 0) {
                scored_on_gpu = true;
                static bool announced = false;
                if (!announced) {
                    announced = true;
                    std::fprintf(stderr, "[bhrt] move scoring running on GPU (CUDA)\n");
                }
                // Self-verify the FIRST successful GPU batch against the CPU
                // reference. Any disagreement beyond FP-contraction noise
                // means a broken GPU pipeline (bad arch/JIT, stale objects,
                // ABI drift) -> permanently fall back to the CPU scorer for
                // this run rather than search on garbage scores.
                static bool verified = false;
                if (!verified) {
                    verified = true;
                    bool agree = true;
                    for (std::size_t i = 0; i < pairs.size() && agree; ++i) {
                        auto r = flattenCandidate(*pairs[i].first,
                                                  pairs[i].second);
                        const float ref = scoreCandidateCPU(r, weights);
                        const float tol = 1e-5f *
                            (std::fabs(ref) > 1.0f ? std::fabs(ref) : 1.0f);
                        if (std::fabs(scores[i] - ref) > tol) agree = false;
                    }
                    if (!agree) {
                        std::fprintf(stderr,
                            "[bhrt] GPU scores disagree with CPU reference; "
                            "disabling GPU scorer for this run\n");
                        use_gpu = false;
                        scored_on_gpu = false;   // rescore this batch on CPU
                    } else {
                        std::fprintf(stderr,
                            "[bhrt] GPU scorer verified against CPU reference "
                            "(%zu candidates)\n", pairs.size());
                    }
                }
            }
        }
#endif
        if (!scored_on_gpu) {
            for (std::size_t i = 0; i < pairs.size(); ++i) {
                auto r = flattenCandidate(*pairs[i].first, pairs[i].second);
                scores[i] = scoreCandidateCPU(r, weights);
            }
        }
        std::vector<std::size_t> order(pairs.size());
        for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
        std::sort(order.begin(), order.end(),
                   [&](auto a, auto b) { return scores[a] > scores[b]; });

        std::vector<SearchState> new_beam;
        int considered = 0;
        for (auto idx : order) {
            if (considered >= cfg.beam_width) break;
            auto& parent = *pairs[idx].first;
            auto& cand = pairs[idx].second;
            auto T2 = applyMove(parent.triangulation, cand);
            if (!T2) continue;
            int new_excess = std::max(0,
                parent.excess_height +
                ((int)T2->size() - parent.pentachora));
            if (new_excess > cfg.excess_height) continue;
            auto sig = T2->edgeDegreeIsoSig();
            if (visited.count(sig)) continue;
            visited.insert(sig);
            SearchState child;
            child.triangulation = std::move(*T2);
            child.pentachora = (int)child.triangulation.size();
            child.isosig = sig;
            // The GPU enumerator is a pure filter (survivors re-verified on
            // CPU), so feasibility is identical on either path.
#if BHRT_HAS_CUDA
            auto ts = use_gpu
                ? enumerateTsTricolouringsGPU(child.triangulation, 10000, true)
                : enumerateTsTricolourings(child.triangulation, 10000, true);
#else
            auto ts = enumerateTsTricolourings(child.triangulation, 10000, true);
#endif
            child.ts_feasible = !ts.empty();
            if (child.ts_feasible) {
                auto d = extractDiagram(child.triangulation, ts.front());
                child.best_genus = d.genus();
            }
            child.excess_height = new_excess;
            child.move_history = parent.move_history;
            child.move_history.push_back(std::to_string(static_cast<int>(cand.move)));
            child.score = scores[idx];
            new_beam.push_back(std::move(child));
            pareto.consider(new_beam.back());
            ++considered;
        }
        if (new_beam.empty()) break;
        beam = std::move(new_beam);
    }
    return pareto;
}

}  // namespace bhrt
