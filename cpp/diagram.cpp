// diagram.cpp -- central surface + cut-system extractor (C++ core).
//
// Implements the triangulation-supported trisection diagram. Each
// pentachoron of a ts-tricoloured input contributes one square Sigma-
// face; the four Sigma-vertices sit on the bi-coloured edges of the
// pent that cross between the two pair-colours. Edges of Sigma are
// reconstructed from face cycles; cut curves are read off as meridians
// of non-tree edges of the spine graph of each gamma_{i,j}.

#include "bhrt_trisect.hpp"

#include <algorithm>
#include <deque>
#include <set>
#include <map>
#include <stdexcept>

namespace bhrt {

std::int32_t CellularSurface::euler() const noexcept {
    return num_vertices
         - static_cast<std::int32_t>(edges.size())
         + static_cast<std::int32_t>(faces.size());
}

std::int32_t CellularSurface::genus() const noexcept {
    const auto chi = euler();
    return std::max(0, (2 - chi) / 2);
}

namespace {

// Spine graph of a TwoComplex: collapse to 1-complex, then return (V, E).
struct SpineGraph {
    std::vector<std::int32_t> vertices;
    std::vector<std::pair<std::int32_t, std::int32_t>> edges;
};

SpineGraph spineOfGamma(const TwoComplex& gamma) {
    TwoComplex scratch = gamma;
    scratch.greedyCollapseTo1();
    SpineGraph g;
    std::set<std::int32_t> vs(scratch.vertices.begin(), scratch.vertices.end());
    g.vertices.assign(vs.begin(), vs.end());
    for (auto e : scratch.edges) {
        std::int32_t a = static_cast<std::int32_t>(e >> 21);
        std::int32_t b = static_cast<std::int32_t>(e & ((1LL << 21) - 1));
        g.edges.emplace_back(a, b);
    }
    return g;
}

std::set<std::size_t> maximalSpanningTree(const SpineGraph& g) {
    std::map<std::int32_t, std::int32_t> parent;
    for (auto v : g.vertices) parent[v] = v;
    auto find = [&](std::int32_t x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    std::set<std::size_t> tree;
    for (std::size_t i = 0; i < g.edges.size(); ++i) {
        auto [a, b] = g.edges[i];
        auto ra = find(a), rb = find(b);
        if (ra != rb) { parent[ra] = rb; tree.insert(i); }
    }
    return tree;
}

}  // namespace

TrisectionDiagram extractDiagram(const Triangulation& T,
                                  const TSColouring&  colouring) {
    TrisectionDiagram diag;
    if (!colouring.is_ts) {
        diag.audit.push_back("input colouring is not certified ts; "
                              "diagram extraction will fall back to best-effort");
    }

    // Build Sigma-vertices: one per (pent, bi-coloured edge) pair.
    std::map<std::int32_t, std::int32_t> sigma_vid_of_edge;
    std::int32_t next_v = 0;

    for (std::int32_t p = 0; p < (std::int32_t)T.size(); ++p) {
        int loc_c[5];
        int counts[3]{0, 0, 0};
        for (std::uint8_t v = 0; v < 5; ++v) {
            loc_c[v] = colouring.colour[T.vertexOf(p, v)];
            ++counts[loc_c[v]];
        }
        int singleton = -1;
        for (int c = 0; c < 3; ++c) if (counts[c] == 1) { singleton = c; break; }
        if (singleton < 0) continue;
        int pair_a = (singleton + 1) % 3;
        int pair_b = (singleton + 2) % 3;
        std::vector<std::uint8_t> a_verts, b_verts;
        for (std::uint8_t v = 0; v < 5; ++v) {
            if (loc_c[v] == pair_a) a_verts.push_back(v);
            if (loc_c[v] == pair_b) b_verts.push_back(v);
        }
        if (a_verts.size() != 2 || b_verts.size() != 2) continue;
        std::array<std::int32_t, 4> face4{};
        int slot = 0;
        for (auto a : a_verts) for (auto b : b_verts) {
            auto eid = T.edgeOf(p, a, b);
            auto it = sigma_vid_of_edge.find(eid);
            if (it == sigma_vid_of_edge.end()) {
                sigma_vid_of_edge[eid] = next_v;
                face4[slot++] = next_v++;
            } else {
                face4[slot++] = it->second;
            }
        }
        // 4-cycle ordering: (a0,b0) -> (a0,b1) -> (a1,b1) -> (a1,b0)
        std::vector<std::int32_t> cycle{face4[0], face4[1], face4[3], face4[2]};
        diag.surface.faces.push_back(std::move(cycle));
    }
    diag.surface.num_vertices = next_v;

    // Reconstruct edges + edge_index from face cycles.
    for (const auto& f : diag.surface.faces) {
        for (std::size_t i = 0; i < f.size(); ++i) {
            std::int32_t a = f[i], b = f[(i + 1) % f.size()];
            if (a > b) std::swap(a, b);
            std::pair<std::int32_t, std::int32_t> key{a, b};
            auto [it, inserted] = diag.surface.edge_index.try_emplace(
                key, static_cast<std::int32_t>(diag.surface.edges.size()));
            if (inserted) diag.surface.edges.emplace_back(a, b);
        }
    }

    // Extract cut curves from each gamma_{i,j} spine.
    static const std::vector<std::tuple<int, int, std::string>> roles{
        {0, 1, "alpha"}, {0, 2, "beta"}, {1, 2, "gamma"}};
    for (auto& [i, j, role] : roles) {
        auto it = colouring.gamma_2complexes.find({i, j});
        if (it == colouring.gamma_2complexes.end()) continue;
        auto spine = spineOfGamma(it->second);
        auto tree = maximalSpanningTree(spine);
        for (std::size_t e_idx = 0; e_idx < spine.edges.size(); ++e_idx) {
            if (tree.count(e_idx)) continue;
            auto [u, v] = spine.edges[e_idx];
            // Walk Sigma-faces containing both u and v to define the
            // meridian.
            for (const auto& face : diag.surface.faces) {
                bool has_u = false, has_v = false;
                for (auto x : face) {
                    if (x == u) has_u = true;
                    if (x == v) has_v = true;
                }
                if (!(has_u && has_v)) continue;
                CutCurve c;
                c.role = role;
                for (std::size_t k = 0; k < face.size(); ++k) {
                    std::int32_t a = face[k], b = face[(k + 1) % face.size()];
                    if (a > b) std::swap(a, b);
                    c.edge_sequence.push_back(
                        diag.surface.edge_index.at({a, b}));
                }
                if (role == "alpha") diag.alpha.push_back(std::move(c));
                else if (role == "beta") diag.beta.push_back(std::move(c));
                else diag.gamma.push_back(std::move(c));
                break;
            }
        }
    }
    return diag;
}

// ---------------------------------------------------------------------- //
// Validation: is this cellular diagram a genuine trisection diagram?
// ---------------------------------------------------------------------- //

namespace {

// Union-find with parity (relative to the component root).
struct ParityUF {
    std::vector<int> parent, rankv, par;  // par = parity to parent
    void init(int n) {
        parent.resize(n); rankv.assign(n, 0); par.assign(n, 0);
        for (int i = 0; i < n; ++i) parent[i] = i;
    }
    // returns {root, parity-to-root}
    std::pair<int,int> find(int x) {
        if (parent[x] == x) return {x, 0};
        auto [r, p] = find(parent[x]);
        parent[x] = r;
        par[x] ^= p;
        return {r, par[x]};
    }
    // constrain: parity(a) xor parity(b) == want.  Returns false on conflict.
    bool unite(int a, int b, int want) {
        auto [ra, pa] = find(a);
        auto [rb, pb] = find(b);
        if (ra == rb) return ((pa ^ pb) == want);
        if (rankv[ra] < rankv[rb]) { std::swap(ra, rb); std::swap(pa, pb); }
        parent[rb] = ra;
        par[rb] = pa ^ pb ^ want;
        if (rankv[ra] == rankv[rb]) ++rankv[ra];
        return true;
    }
};

bool curveIsClosedLoop(const CutCurve& c, const CellularSurface& S) {
    std::map<std::int32_t, int> deg;
    for (auto e : c.edge_sequence) {
        if (e < 0 || e >= (std::int32_t)S.edges.size()) return false;
        ++deg[S.edges[e].first];
        ++deg[S.edges[e].second];
    }
    for (auto& kv : deg) if (kv.second % 2 != 0) return false;
    return !c.edge_sequence.empty();
}

}  // namespace

DiagramValidity validateDiagram(const TrisectionDiagram& diagram) {
    DiagramValidity v;
    const auto& S = diagram.surface;
    v.genus = diagram.genus();
    const std::int32_t nv = S.num_vertices;

    // (1) Connectivity of the 1-skeleton.
    if (nv <= 0) {
        v.messages.push_back("surface has no vertices");
        return v;
    }
    std::vector<std::vector<std::int32_t>> adj(nv);
    for (auto& e : S.edges) {
        if (e.first < nv && e.second < nv) {
            adj[e.first].push_back(e.second);
            adj[e.second].push_back(e.first);
        }
    }
    std::vector<char> seen(nv, 0);
    std::deque<std::int32_t> q{0};
    seen[0] = 1;
    std::int32_t nseen = 1;
    while (!q.empty()) {
        auto x = q.front(); q.pop_front();
        for (auto y : adj[x]) if (!seen[y]) { seen[y] = 1; ++nseen; q.push_back(y); }
    }
    v.surface_connected = (nseen == nv);
    if (!v.surface_connected) v.messages.push_back("surface 1-skeleton disconnected");

    // (2) Closedness + orientability from face boundary traversals.
    // For each undirected edge collect (face, dir) where dir = traversed
    // low->high (true) or high->low (false).
    std::map<std::pair<std::int32_t,std::int32_t>,
             std::vector<std::pair<int,bool>>> inc;
    for (std::size_t fi = 0; fi < S.faces.size(); ++fi) {
        const auto& f = S.faces[fi];
        for (std::size_t i = 0; i < f.size(); ++i) {
            std::int32_t a = f[i], b = f[(i + 1) % f.size()];
            bool low_high = (a < b);
            std::pair<std::int32_t,std::int32_t> key{std::min(a,b), std::max(a,b)};
            inc[key].push_back({(int)fi, low_high});
        }
    }
    bool closed = true;
    bool orientable = true;
    ParityUF uf;
    uf.init((int)S.faces.size());
    for (auto& kv : inc) {
        if (kv.second.size() != 2) { closed = false; continue; }
        auto [f1, d1] = kv.second[0];
        auto [f2, d2] = kv.second[1];
        // consistent orientation requires opposite traversal directions:
        // if d1==d2 the faces must have opposite orientation (want=1).
        int want = (d1 == d2) ? 1 : 0;
        if (!uf.unite(f1, f2, want)) orientable = false;
    }
    v.surface_closed = closed && !S.faces.empty();
    v.surface_orientable = orientable && !S.faces.empty();
    if (!closed) v.messages.push_back("surface has boundary or non-manifold edge");
    if (!v.surface_orientable) v.messages.push_back("surface is non-orientable");

    // (3) Cut systems: each family must have exactly g closed, pairwise
    // edge-disjoint curves.
    auto checkFamily = [&](const std::vector<CutCurve>& fam, const char* name) -> bool {
        if ((std::int32_t)fam.size() != v.genus) {
            v.messages.push_back(std::string(name) + ": expected g curves, got "
                                 + std::to_string(fam.size()));
            return false;
        }
        std::set<std::int32_t> used;
        for (const auto& c : fam) {
            if (!curveIsClosedLoop(c, S)) {
                v.messages.push_back(std::string(name) + ": a curve is not a closed loop");
                return false;
            }
            for (auto e : c.edge_sequence) {
                if (!used.insert(e).second) {
                    v.messages.push_back(std::string(name) + ": curves share an edge");
                    return false;
                }
            }
        }
        return true;
    };
    v.alpha_is_cut_system = checkFamily(diagram.alpha, "alpha");
    v.beta_is_cut_system  = checkFamily(diagram.beta,  "beta");
    v.gamma_is_cut_system = checkFamily(diagram.gamma, "gamma");

    v.ok = v.surface_connected && v.surface_closed && v.surface_orientable
         && v.alpha_is_cut_system && v.beta_is_cut_system && v.gamma_is_cut_system;
    return v;
}

}  // namespace bhrt
