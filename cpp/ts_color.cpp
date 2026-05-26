// ts_color.cpp -- direct ts-tricolouring algorithm (Spreer-Tillmann).
//
// Full, self-contained C++ implementation. Mirrors the Python reference
// in python/bhrt_trisect/ts_color.py; both implementations produce the
// same TSColouring results on the same input.

#include "bhrt_trisect.hpp"

#include <algorithm>
#include <cstdio>
#include <deque>
#include <map>
#include <queue>
#include <set>

namespace bhrt {

// ---------------------------------------------------------------------- //
// TwoComplex bit-packed helpers
// ---------------------------------------------------------------------- //

std::int64_t TwoComplex::packEdge(std::int32_t a, std::int32_t b) noexcept {
    if (a > b) std::swap(a, b);
    return (static_cast<std::int64_t>(a) << 21) | static_cast<std::int64_t>(b);
}

std::int64_t TwoComplex::packTri(std::int32_t a, std::int32_t b,
                                  std::int32_t c) noexcept {
    std::int32_t v[3]{a, b, c};
    std::sort(v, v + 3);
    return (static_cast<std::int64_t>(v[0]) << 42)
         | (static_cast<std::int64_t>(v[1]) << 21)
         |  static_cast<std::int64_t>(v[2]);
}

void TwoComplex::addEdge(std::int32_t a, std::int32_t b) {
    if (a == b) return;
    auto e = packEdge(a, b);
    edges.insert(e);
    edge_to_faces.try_emplace(e);
    vertices.insert(a);
    vertices.insert(b);
}

void TwoComplex::addFace(std::int32_t a, std::int32_t b, std::int32_t c) {
    if (a == b || b == c || a == c) return;
    auto key = packTri(a, b, c);
    if (!faces.insert(key).second) return;
    vertices.insert(a);
    vertices.insert(b);
    vertices.insert(c);
    std::int32_t v[3]{a, b, c};
    for (int i = 0; i < 3; ++i) for (int j = i + 1; j < 3; ++j) {
        auto e = packEdge(v[i], v[j]);
        edges.insert(e);
        edge_to_faces[e].insert(key);
    }
}

bool TwoComplex::isConnected() const {
    if (vertices.empty()) return true;
    std::map<std::int32_t, std::vector<std::int32_t>> adj;
    for (auto e : edges) {
        std::int32_t a = static_cast<std::int32_t>(e >> 21);
        std::int32_t b = static_cast<std::int32_t>(e & ((1LL << 21) - 1));
        adj[a].push_back(b);
        adj[b].push_back(a);
    }
    std::set<std::int32_t> seen;
    std::deque<std::int32_t> q;
    auto start = *vertices.begin();
    seen.insert(start);
    q.push_back(start);
    while (!q.empty()) {
        auto x = q.front(); q.pop_front();
        for (auto y : adj[x])
            if (seen.insert(y).second) q.push_back(y);
    }
    return seen.size() == vertices.size();
}

bool TwoComplex::greedyCollapseTo1(int budget) {
    int steps = 0;
    while (!faces.empty()) {
        if (steps++ >= budget) return false;
        std::int64_t free_edge = -1;
        std::int64_t free_face = -1;
        for (auto& kv : edge_to_faces) {
            if (kv.second.size() == 1) {
                free_edge = kv.first;
                free_face = *kv.second.begin();
                break;
            }
        }
        if (free_edge == -1) return false;
        // Remove the free face from all incidences.
        std::int32_t a = static_cast<std::int32_t>(free_face >> 42);
        std::int32_t b = static_cast<std::int32_t>(
            (free_face >> 21) & ((1LL << 21) - 1));
        std::int32_t c = static_cast<std::int32_t>(
            free_face & ((1LL << 21) - 1));
        std::int32_t v[3]{a, b, c};
        for (int i = 0; i < 3; ++i) for (int j = i + 1; j < 3; ++j) {
            auto e = packEdge(v[i], v[j]);
            auto it = edge_to_faces.find(e);
            if (it != edge_to_faces.end()) it->second.erase(free_face);
        }
        faces.erase(free_face);
        edges.erase(free_edge);
        edge_to_faces.erase(free_edge);
    }
    return true;
}

// ---------------------------------------------------------------------- //
// Helpers
// ---------------------------------------------------------------------- //

namespace {

bool everyTriangleHasTwoDistinctVertices(const Triangulation& T) {
    // For each triangle equivalence class, look up its three global
    // vertex ids from one representative.
    std::map<std::int32_t, bool> seen;
    for (std::int32_t p = 0; p < (std::int32_t)T.size(); ++p) {
        for (std::uint8_t a = 0; a < 5; ++a)
            for (std::uint8_t b = a + 1; b < 5; ++b)
                for (std::uint8_t c = b + 1; c < 5; ++c) {
                    auto tid = T.triangleOf(p, a, b, c);
                    if (seen.count(tid)) continue;
                    seen[tid] = true;
                    std::set<std::int32_t> gvs{
                        T.vertexOf(p, a),
                        T.vertexOf(p, b),
                        T.vertexOf(p, c)};
                    if (gvs.size() < 2) return false;
                }
    }
    return true;
}

bool isType221(const Triangulation& T,
               const std::vector<std::int32_t>& colour) {
    for (std::int32_t p = 0; p < (std::int32_t)T.size(); ++p) {
        int counts[3]{0, 0, 0};
        for (std::uint8_t v = 0; v < 5; ++v)
            ++counts[colour[T.vertexOf(p, v)]];
        int sorted[3]{counts[0], counts[1], counts[2]};
        std::sort(sorted, sorted + 3);
        if (!(sorted[0] == 1 && sorted[1] == 2 && sorted[2] == 2))
            return false;
    }
    return true;
}

bool hasMonochromaticTriangle(const Triangulation& T,
                               const std::vector<std::int32_t>& colour) {
    std::set<std::int32_t> seen;
    for (std::int32_t p = 0; p < (std::int32_t)T.size(); ++p) {
        for (std::uint8_t a = 0; a < 5; ++a)
            for (std::uint8_t b = a + 1; b < 5; ++b)
                for (std::uint8_t c = b + 1; c < 5; ++c) {
                    auto tid = T.triangleOf(p, a, b, c);
                    if (!seen.insert(tid).second) continue;
                    int ca = colour[T.vertexOf(p, a)];
                    int cb = colour[T.vertexOf(p, b)];
                    int cc = colour[T.vertexOf(p, c)];
                    if (ca == cb && cb == cc) return true;
                }
    }
    return false;
}

std::vector<std::set<std::int32_t>> monochromaticComponents(
    const Triangulation& T,
    const std::vector<std::int32_t>& colour,
    int k) {
    std::set<std::int32_t> P_k;
    for (std::int32_t v = 0; v < (std::int32_t)colour.size(); ++v)
        if (colour[v] == k) P_k.insert(v);
    std::map<std::int32_t, std::set<std::int32_t>> adj;
    for (auto v : P_k) adj[v];
    std::set<std::int32_t> seen_edges;
    for (std::int32_t p = 0; p < (std::int32_t)T.size(); ++p) {
        for (std::uint8_t a = 0; a < 5; ++a)
            for (std::uint8_t b = a + 1; b < 5; ++b) {
                auto eid = T.edgeOf(p, a, b);
                if (!seen_edges.insert(eid).second) continue;
                int va = T.vertexOf(p, a), vb = T.vertexOf(p, b);
                if (colour[va] == k && colour[vb] == k) {
                    adj[va].insert(vb);
                    adj[vb].insert(va);
                }
            }
    }
    std::set<std::int32_t> seen;
    std::vector<std::set<std::int32_t>> comps;
    for (auto v : P_k) {
        if (seen.count(v)) continue;
        std::set<std::int32_t> comp;
        std::deque<std::int32_t> q{v};
        seen.insert(v);
        while (!q.empty()) {
            auto x = q.front(); q.pop_front();
            comp.insert(x);
            for (auto y : adj[x])
                if (seen.insert(y).second) q.push_back(y);
        }
        comps.push_back(std::move(comp));
    }
    return comps;
}

TwoComplex buildGamma(const Triangulation& T,
                       const std::vector<std::int32_t>& colour,
                       int i, int j) {
    TwoComplex gamma;
    std::set<std::int32_t> seen_tri, seen_edge;
    for (std::int32_t p = 0; p < (std::int32_t)T.size(); ++p) {
        for (std::uint8_t a = 0; a < 5; ++a)
            for (std::uint8_t b = a + 1; b < 5; ++b)
                for (std::uint8_t c = b + 1; c < 5; ++c) {
                    auto tid = T.triangleOf(p, a, b, c);
                    if (!seen_tri.insert(tid).second) continue;
                    int va = T.vertexOf(p, a);
                    int vb = T.vertexOf(p, b);
                    int vc = T.vertexOf(p, c);
                    std::set<int> cs{colour[va], colour[vb], colour[vc]};
                    bool keep = false;
                    if (cs.size() == 2 && cs.count(i) && cs.count(j)) keep = true;
                    if (cs.size() == 1
                        && (*cs.begin() == i || *cs.begin() == j)) keep = true;
                    if (keep) gamma.addFace(va, vb, vc);
                }
        for (std::uint8_t a = 0; a < 5; ++a)
            for (std::uint8_t b = a + 1; b < 5; ++b) {
                auto eid = T.edgeOf(p, a, b);
                if (!seen_edge.insert(eid).second) continue;
                int va = T.vertexOf(p, a), vb = T.vertexOf(p, b);
                if ((colour[va] == i || colour[va] == j)
                    && (colour[vb] == i || colour[vb] == j))
                    gamma.addEdge(va, vb);
            }
    }
    return gamma;
}

}  // namespace

TSColouring isTsTricolouring(const Triangulation& T,
                              const std::vector<std::int32_t>& colour,
                              int collapse_budget) {
    TSColouring tc;
    tc.colour = colour;
    if (T.nVertices() < 3) {
        tc.audit.push_back("precheck: triangulation has fewer than 3 vertices");
        return tc;
    }
    if (!everyTriangleHasTwoDistinctVertices(T)) {
        tc.audit.push_back("precheck: at least one triangle has < 2 distinct vertices");
        return tc;
    }
    if (!isType221(T, colour)) {
        tc.audit.push_back("step3: a pentachoron is not type (2,2,1)");
        return tc;
    }
    if (hasMonochromaticTriangle(T, colour)) {
        tc.audit.push_back("step3: monochromatic triangle present");
        return tc;
    }
    for (int k = 0; k < 3; ++k) {
        auto comps = monochromaticComponents(T, colour, k);
        tc.monochromatic_components[k] = comps;
        if (comps.size() != 1) {
            tc.audit.push_back("step4: Gamma_k disconnected");
            return tc;
        }
    }
    tc.is_c = true;
    int pairs[3][2]{{0, 1}, {0, 2}, {1, 2}};
    for (auto& p : pairs) {
        int i = p[0], j = p[1];
        TwoComplex gamma = buildGamma(T, colour, i, j);
        if (!gamma.isConnected()) {
            tc.audit.push_back("step5: gamma disconnected");
            tc.gamma_2complexes[{i, j}] = gamma;
            return tc;
        }
        TwoComplex scratch = gamma;
        bool ok = scratch.greedyCollapseTo1(collapse_budget);
        tc.gamma_2complexes[{i, j}] = gamma;
        if (!ok) {
            tc.audit.push_back("step5: gamma did not collapse (inconclusive)");
            return tc;
        }
    }
    tc.is_ts = true;
    return tc;
}

// ---------------------------------------------------------------------- //
// Enumerate ts-tricolourings
// ---------------------------------------------------------------------- //

namespace {

// Yield ordered partitions of {0..n-1} into 3 non-empty parts with
// vertex 0 pinned to part 0.
template <class F>
void enumeratePartitions(int n, F&& yield) {
    if (n < 3) return;
    std::vector<int> assignment(n - 1, 0);
    while (true) {
        std::vector<int> p0{0}, p1, p2;
        for (int i = 0; i < n - 1; ++i) {
            int a = assignment[i];
            if (a == 0) p0.push_back(i + 1);
            else if (a == 1) p1.push_back(i + 1);
            else p2.push_back(i + 1);
        }
        if (!p1.empty() && !p2.empty()) {
            yield(p0, p1, p2);
        }
        // Increment assignment in base 3.
        int i = 0;
        while (i < n - 1) {
            if (assignment[i] < 2) { assignment[i]++; break; }
            assignment[i] = 0;
            ++i;
        }
        if (i == n - 1) break;
    }
}

}  // namespace

std::vector<TSColouring> enumerateTsTricolourings(
    const Triangulation& T,
    int  collapse_budget,
    bool stop_at_first) {
    std::vector<TSColouring> out;
    if (T.nVertices() < 3) return out;
    if (!everyTriangleHasTwoDistinctVertices(T)) return out;
    enumeratePartitions(T.nVertices(),
        [&](const std::vector<int>& p0,
            const std::vector<int>& p1,
            const std::vector<int>& p2) {
            if (stop_at_first && !out.empty()) return;
            std::vector<std::int32_t> colour(T.nVertices(), 0);
            for (auto v : p0) colour[v] = 0;
            for (auto v : p1) colour[v] = 1;
            for (auto v : p2) colour[v] = 2;
            if (!isType221(T, colour)) return;
            if (hasMonochromaticTriangle(T, colour)) return;
            auto tc = isTsTricolouring(T, colour, collapse_budget);
            if (tc.is_ts) out.push_back(std::move(tc));
        });
    return out;
}

// ---------------------------------------------------------------------- //
// GPU-accelerated enumeration (optional)
// ---------------------------------------------------------------------- //

#if BHRT_HAS_CUDA
// Host launcher defined in cuda/ts_scan.cu. The GPU exhaustively filters
// all 3^(n-1) candidate colourings down to those passing the cheap
// prechecks; returns non-zero when the GPU path is unavailable.
extern "C" int bhrt_ts_scan_host(const unsigned char* pent_verts, int n_pents,
                                 const unsigned char* tri_verts,  int n_tris,
                                 const unsigned char* edge_verts, int n_edges,
                                 int n_vertices,
                                 unsigned long long total,
                                 unsigned long long* survivors,
                                 unsigned long long cap,
                                 unsigned long long* n_survivors);

std::vector<TSColouring> enumerateTsTricolouringsGPU(
    const Triangulation& T,
    int  collapse_budget,
    bool stop_at_first) {
    // Fall back to the CPU enumerator whenever the GPU filter cannot be
    // used; results are identical either way (the GPU only prunes, and
    // every survivor is re-verified below by the CPU reference check).
    auto cpu = [&] {
        return enumerateTsTricolourings(T, collapse_budget, stop_at_first);
    };
    const int nv = T.nVertices();
    if (nv < 3) return {};
    if (nv > 32) return cpu();
    if (!everyTriangleHasTwoDistinctVertices(T)) return {};

    // total = 3^(nv-1), with overflow guard (nv <= 32 keeps this < 2^50).
    unsigned long long total = 1;
    for (int i = 0; i < nv - 1; ++i) total *= 3ULL;

    // Pack the skeleton into flat device-friendly arrays.
    std::vector<unsigned char> pent(T.size() * 5);
    for (std::int32_t p = 0; p < (std::int32_t)T.size(); ++p)
        for (std::uint8_t v = 0; v < 5; ++v)
            pent[p * 5 + v] = static_cast<unsigned char>(T.vertexOf(p, v));

    std::vector<unsigned char> tris;
    {
        std::set<std::int32_t> seen;
        for (std::int32_t p = 0; p < (std::int32_t)T.size(); ++p)
            for (std::uint8_t a = 0; a < 5; ++a)
                for (std::uint8_t b = a + 1; b < 5; ++b)
                    for (std::uint8_t c = b + 1; c < 5; ++c) {
                        if (!seen.insert(T.triangleOf(p, a, b, c)).second)
                            continue;
                        tris.push_back((unsigned char)T.vertexOf(p, a));
                        tris.push_back((unsigned char)T.vertexOf(p, b));
                        tris.push_back((unsigned char)T.vertexOf(p, c));
                    }
    }
    std::vector<unsigned char> edges;
    {
        std::set<std::int32_t> seen;
        for (std::int32_t p = 0; p < (std::int32_t)T.size(); ++p)
            for (std::uint8_t a = 0; a < 5; ++a)
                for (std::uint8_t b = a + 1; b < 5; ++b) {
                    if (!seen.insert(T.edgeOf(p, a, b)).second) continue;
                    edges.push_back((unsigned char)T.vertexOf(p, a));
                    edges.push_back((unsigned char)T.vertexOf(p, b));
                }
    }

    const unsigned long long cap =
        total < (1ULL << 20) ? total : (1ULL << 20);
    std::vector<unsigned long long> surv(static_cast<std::size_t>(cap));
    unsigned long long n_surv = 0;
    int rc = bhrt_ts_scan_host(
        pent.data(), (int)T.size(),
        tris.data(), (int)(tris.size() / 3),
        edges.data(), (int)(edges.size() / 2),
        nv, total, surv.data(), cap, &n_surv);
    if (rc != 0) return cpu();   // no GPU, error, overflow or unsupported

    static bool announced = false;
    if (!announced) {
        announced = true;
        std::fprintf(stderr,
            "[bhrt] ts-colouring precheck running on GPU (CUDA)\n");
    }

    // CPU-verify every survivor with the full certified check (including
    // the gamma 2-complex collapsibility step the GPU does not attempt).
    std::vector<TSColouring> out;
    for (unsigned long long s = 0; s < n_surv; ++s) {
        std::vector<std::int32_t> colour(nv, 0);
        unsigned long long x = surv[s];
        for (int v = 1; v < nv; ++v) {
            colour[v] = static_cast<std::int32_t>(x % 3ULL);
            x /= 3ULL;
        }
        auto tc = isTsTricolouring(T, colour, collapse_budget);
        if (tc.is_ts) {
            out.push_back(std::move(tc));
            if (stop_at_first) break;
        }
    }
    return out;
}
#endif  // BHRT_HAS_CUDA

ScanCounters countAdmitTrisection(const std::vector<Triangulation>& tris,
                                   int collapse_budget) {
    ScanCounters c;
    for (const auto& T : tris) {
        ++c.n_total;
        if (T.nVertices() < 3) continue;
        if (!everyTriangleHasTwoDistinctVertices(T)) continue;
        ++c.n_pass_precheck;
        bool has_tri = false, has_c = false;
        std::vector<TSColouring> ts_for_T;
        enumeratePartitions(T.nVertices(),
            [&](const std::vector<int>& p0,
                const std::vector<int>& p1,
                const std::vector<int>& p2) {
                std::vector<std::int32_t> colour(T.nVertices(), 0);
                for (auto v : p0) colour[v] = 0;
                for (auto v : p1) colour[v] = 1;
                for (auto v : p2) colour[v] = 2;
                if (!isType221(T, colour)) return;
                if (hasMonochromaticTriangle(T, colour)) return;
                has_tri = true;
                auto tc = isTsTricolouring(T, colour, collapse_budget);
                if (tc.is_c) has_c = true;
                if (tc.is_ts) ts_for_T.push_back(std::move(tc));
            });
        if (has_tri) ++c.n_with_tricolouring;
        if (has_c) ++c.n_with_c_tricolouring;
        if (!ts_for_T.empty()) {
            ++c.n_with_ts_tricolouring;
            c.n_ts_tricolourings_total += static_cast<std::int64_t>(ts_for_T.size());
        }
    }
    return c;
}

}  // namespace bhrt
