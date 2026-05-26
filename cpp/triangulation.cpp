// triangulation.cpp -- core 4D simplicial-complex implementation.
//
// Provides Triangulation, Pentachoron, skeleton computation (vertices,
// edges, triangles, tetrahedra), and the canonical isomorphism
// signature. No external topology library is required.

#include "bhrt_trisect.hpp"

#include <algorithm>
#include <cstdio>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>

namespace bhrt {

// ---------------------------------------------------------------------- //
// Permutations of {0..4}
// ---------------------------------------------------------------------- //

Perm5 perm_compose(const Perm5& p, const Perm5& q) noexcept {
    return Perm5{p[q[0]], p[q[1]], p[q[2]], p[q[3]], p[q[4]]};
}

Perm5 perm_inverse(const Perm5& p) noexcept {
    Perm5 r{};
    for (std::uint8_t i = 0; i < 5; ++i) r[p[i]] = i;
    return r;
}

int perm_sign(const Perm5& p) noexcept {
    int s = 0;
    for (int i = 0; i < 5; ++i)
        for (int j = i + 1; j < 5; ++j)
            if (p[i] > p[j]) s ^= 1;
    return s ? -1 : 1;
}

// ---------------------------------------------------------------------- //
// Disjoint-set union, parameterised by a hashable key type
// ---------------------------------------------------------------------- //

template <class K>
class UnionFind {
public:
    void touch(const K& x) {
        if (parent_.find(x) == parent_.end()) {
            parent_[x] = x;
            rank_[x] = 0;
        }
    }
    K find(K x) {
        touch(x);
        K root = x;
        while (parent_[root] != root) root = parent_[root];
        while (parent_[x] != root) {
            K next = parent_[x];
            parent_[x] = root;
            x = next;
        }
        return root;
    }
    void unite(K a, K b) {
        K ra = find(a), rb = find(b);
        if (ra == rb) return;
        if (rank_[ra] < rank_[rb]) std::swap(ra, rb);
        parent_[rb] = ra;
        if (rank_[ra] == rank_[rb]) rank_[ra]++;
    }
    std::map<K, std::int32_t> dense_index() {
        std::map<K, K> root_of;
        std::set<K> roots;
        for (auto& kv : parent_) {
            K r = find(kv.first);
            root_of[kv.first] = r;
            roots.insert(r);
        }
        std::map<K, std::int32_t> id_of_root;
        std::int32_t i = 0;
        for (auto& r : roots) id_of_root[r] = i++;
        std::map<K, std::int32_t> out;
        for (auto& kv : root_of) out[kv.first] = id_of_root[kv.second];
        return out;
    }

private:
    std::map<K, K> parent_;
    std::map<K, std::int32_t> rank_;
};

// ---------------------------------------------------------------------- //
// Skeleton cache
// ---------------------------------------------------------------------- //

struct Triangulation::Skeleton {
    // Keys: (pent, local_vertex) -> vertex_id
    std::map<std::pair<std::int32_t, std::uint8_t>, std::int32_t> v_id;
    // Keys: (pent, packed pair of locals) -> edge_id
    std::map<std::pair<std::int32_t, std::uint16_t>, std::int32_t> e_id;
    // Keys: (pent, packed triple of locals) -> tri_id
    std::map<std::pair<std::int32_t, std::uint16_t>, std::int32_t> t_id;
    // Keys: (pent, facet) -> tetra_id
    std::map<std::pair<std::int32_t, std::uint8_t>, std::int32_t> tet_id;

    std::int32_t n_v{0}, n_e{0}, n_t{0}, n_tet{0};

    static std::uint16_t packPair(std::uint8_t a, std::uint8_t b) noexcept {
        if (a > b) std::swap(a, b);
        return static_cast<std::uint16_t>((a << 4) | b);
    }
    static std::uint16_t packTriple(std::uint8_t a, std::uint8_t b,
                                     std::uint8_t c) noexcept {
        std::uint8_t v[3]{a, b, c};
        std::sort(v, v + 3);
        return static_cast<std::uint16_t>((v[0] << 8) | (v[1] << 4) | v[2]);
    }
};

void Triangulation::ensureSkeleton() const {
    if (!dirty_ && skel_) return;
    auto s = std::make_shared<Skeleton>();
    UnionFind<std::pair<std::int32_t, std::uint8_t>>   uf_v;
    UnionFind<std::pair<std::int32_t, std::uint16_t>>  uf_e;
    UnionFind<std::pair<std::int32_t, std::uint16_t>>  uf_t;
    UnionFind<std::pair<std::int32_t, std::uint8_t>>   uf_T;

    const std::int32_t n = static_cast<std::int32_t>(pentachora_.size());
    for (std::int32_t p = 0; p < n; ++p) {
        for (std::uint8_t v = 0; v < 5; ++v) uf_v.touch({p, v});
        for (std::uint8_t a = 0; a < 5; ++a)
            for (std::uint8_t b = a + 1; b < 5; ++b)
                uf_e.touch({p, Skeleton::packPair(a, b)});
        for (std::uint8_t a = 0; a < 5; ++a)
            for (std::uint8_t b = a + 1; b < 5; ++b)
                for (std::uint8_t c = b + 1; c < 5; ++c)
                    uf_t.touch({p, Skeleton::packTriple(a, b, c)});
        for (std::uint8_t f = 0; f < 5; ++f) uf_T.touch({p, f});
    }

    for (std::int32_t p = 0; p < n; ++p) {
        for (std::uint8_t f = 0; f < 5; ++f) {
            const auto& g = pentachora_[p].gluings[f];
            if (!g.has_value()) continue;
            if (g->dst_pent < p) continue;
            // Shared face vertices are {0..4} \ {f}.
            std::uint8_t shared[4];
            std::uint8_t k = 0;
            for (std::uint8_t v = 0; v < 5; ++v) if (v != f) shared[k++] = v;
            for (std::uint8_t v : shared)
                uf_v.unite({p, v}, {g->dst_pent, g->perm[v]});
            for (std::uint8_t i = 0; i < 4; ++i)
                for (std::uint8_t j = i + 1; j < 4; ++j) {
                    auto local_pair = Skeleton::packPair(shared[i], shared[j]);
                    auto dst_pair = Skeleton::packPair(g->perm[shared[i]],
                                                        g->perm[shared[j]]);
                    uf_e.unite({p, local_pair}, {g->dst_pent, dst_pair});
                }
            for (std::uint8_t i = 0; i < 4; ++i)
                for (std::uint8_t j = i + 1; j < 4; ++j)
                    for (std::uint8_t k2 = j + 1; k2 < 4; ++k2) {
                        auto local_tri = Skeleton::packTriple(
                            shared[i], shared[j], shared[k2]);
                        auto dst_tri = Skeleton::packTriple(
                            g->perm[shared[i]], g->perm[shared[j]],
                            g->perm[shared[k2]]);
                        uf_t.unite({p, local_tri}, {g->dst_pent, dst_tri});
                    }
            uf_T.unite({p, f}, {g->dst_pent, g->dst_facet});
        }
    }

    s->v_id = uf_v.dense_index();
    s->e_id = uf_e.dense_index();
    s->t_id = uf_t.dense_index();
    s->tet_id = uf_T.dense_index();
    auto distinct = [](const auto& m) {
        std::set<std::int32_t> u;
        for (auto& kv : m) u.insert(kv.second);
        return static_cast<std::int32_t>(u.size());
    };
    s->n_v = distinct(s->v_id);
    s->n_e = distinct(s->e_id);
    s->n_t = distinct(s->t_id);
    s->n_tet = distinct(s->tet_id);
    skel_ = s;
    dirty_ = false;
}

std::int32_t Triangulation::vertexOf(std::int32_t pent, std::uint8_t local) const {
    ensureSkeleton();
    return skel_->v_id.at({pent, local});
}
std::int32_t Triangulation::edgeOf(std::int32_t pent,
                                    std::uint8_t a, std::uint8_t b) const {
    ensureSkeleton();
    return skel_->e_id.at({pent, Skeleton::packPair(a, b)});
}
std::int32_t Triangulation::triangleOf(std::int32_t pent,
                                        std::uint8_t a,
                                        std::uint8_t b,
                                        std::uint8_t c) const {
    ensureSkeleton();
    return skel_->t_id.at({pent, Skeleton::packTriple(a, b, c)});
}
std::int32_t Triangulation::tetraOf(std::int32_t pent, std::uint8_t facet) const {
    ensureSkeleton();
    return skel_->tet_id.at({pent, facet});
}
std::int32_t Triangulation::nVertices() const { ensureSkeleton(); return skel_->n_v; }
std::int32_t Triangulation::nEdges()    const { ensureSkeleton(); return skel_->n_e; }
std::int32_t Triangulation::nTriangles() const { ensureSkeleton(); return skel_->n_t; }
std::int32_t Triangulation::nTetrahedra() const { ensureSkeleton(); return skel_->n_tet; }

std::int32_t Triangulation::eulerCharacteristic() const {
    return nVertices() - nEdges() + nTriangles() - nTetrahedra()
         + static_cast<std::int32_t>(size());
}

std::vector<std::int32_t> Triangulation::edgeDegrees() const {
    ensureSkeleton();
    std::map<std::int32_t, std::int32_t> counts;
    for (auto& kv : skel_->e_id) ++counts[kv.second];
    std::vector<std::int32_t> out;
    out.reserve(counts.size());
    for (auto& kv : counts) out.push_back(kv.second);
    std::sort(out.begin(), out.end());
    return out;
}

// ---------------------------------------------------------------------- //
// Mutation
// ---------------------------------------------------------------------- //

Pentachoron& Triangulation::newPentachoron() {
    Pentachoron p;
    p.index = static_cast<std::int32_t>(pentachora_.size());
    pentachora_.push_back(p);
    dirty_ = true;
    return pentachora_.back();
}

void Triangulation::addPentachora(std::size_t count) {
    for (std::size_t i = 0; i < count; ++i) (void)newPentachoron();
}

void Triangulation::glue(std::int32_t src, std::uint8_t src_facet,
                          std::int32_t dst, const Perm5& perm) {
    if (perm[src_facet] >= 5)
        throw std::runtime_error("invalid gluing permutation");
    std::uint8_t dst_facet = perm[src_facet];
    pentachora_[src].gluings[src_facet] = Gluing{
        src_facet, dst, dst_facet, perm};
    Perm5 inv = perm_inverse(perm);
    pentachora_[dst].gluings[dst_facet] = Gluing{
        dst_facet, src, src_facet, inv};
    dirty_ = true;
}

void Triangulation::unglue(std::int32_t pent, std::uint8_t facet) {
    auto& g = pentachora_[pent].gluings[facet];
    if (!g.has_value()) return;
    pentachora_[g->dst_pent].gluings[g->dst_facet].reset();
    g.reset();
    dirty_ = true;
}

bool Triangulation::isClosed() const noexcept {
    for (const auto& p : pentachora_)
        for (const auto& g : p.gluings)
            if (!g.has_value()) return false;
    return true;
}

bool Triangulation::isValid() const noexcept {
    for (std::size_t p = 0; p < pentachora_.size(); ++p) {
        for (std::uint8_t f = 0; f < 5; ++f) {
            const auto& g = pentachora_[p].gluings[f];
            if (!g.has_value()) continue;
            const auto& mate = pentachora_[g->dst_pent].gluings[g->dst_facet];
            if (!mate.has_value()
                || mate->dst_pent != static_cast<std::int32_t>(p)
                || mate->dst_facet != f) return false;
            if (perm_compose(g->perm, mate->perm) != perm_identity) return false;
        }
    }
    return true;
}

Triangulation Triangulation::clone() const {
    Triangulation copy(label_);
    copy.pentachora_ = pentachora_;
    copy.dirty_ = true;
    return copy;
}

// ---------------------------------------------------------------------- //
// Quick hash
// ---------------------------------------------------------------------- //

std::string Triangulation::quickHash() const {
    // FNV-1a 64-bit
    std::uint64_t h = 0xcbf29ce484222325ULL;
    auto mix = [&](std::uint8_t b) {
        h ^= b;
        h *= 0x100000001b3ULL;
    };
    for (const auto& p : pentachora_) {
        for (const auto& g : p.gluings) {
            if (!g.has_value()) {
                mix(0xFF);
            } else {
                mix(g->src_facet);
                mix(static_cast<std::uint8_t>(g->dst_pent & 0xFF));
                mix(static_cast<std::uint8_t>((g->dst_pent >> 8) & 0xFF));
                mix(g->dst_facet);
                for (auto v : g->perm) mix(v);
            }
        }
    }
    char buf[17];
    std::snprintf(buf, sizeof(buf), "%016lx",
                   static_cast<unsigned long>(h));
    return std::string(buf);
}

// ---------------------------------------------------------------------- //
// Canonical isomorphism signature
// ---------------------------------------------------------------------- //

namespace {

std::string canonicalRelabel(const Triangulation& T,
                              std::int32_t start_pent,
                              const Perm5& start_perm) {
    const std::int32_t n = static_cast<std::int32_t>(T.size());
    std::map<std::int32_t, std::int32_t> old_to_new;
    std::map<std::int32_t, Perm5>        relabel;
    old_to_new[start_pent] = 0;
    relabel[start_pent] = start_perm;

    std::vector<std::int32_t> queue{start_pent};
    std::vector<std::int32_t> visited{start_pent};
    std::size_t head = 0;
    std::vector<std::uint8_t> code;

    while (head < queue.size()) {
        std::int32_t cur = queue[head++];
        Perm5 cur_perm = relabel[cur];
        for (std::uint8_t new_f = 0; new_f < 5; ++new_f) {
            std::uint8_t old_f = cur_perm[new_f];
            const auto& g = T.pentachora()[cur].gluings[old_f];
            if (!g.has_value()) {
                code.push_back(255);
                continue;
            }
            auto it = old_to_new.find(g->dst_pent);
            if (it == old_to_new.end()) {
                std::int32_t new_idx = static_cast<std::int32_t>(old_to_new.size());
                old_to_new[g->dst_pent] = new_idx;
                queue.push_back(g->dst_pent);
                visited.push_back(g->dst_pent);
                Perm5 neighbour_perm{};
                for (std::uint8_t v = 0; v < 5; ++v) {
                    if (v == new_f) neighbour_perm[v] = g->dst_facet;
                    else            neighbour_perm[v] = g->perm[cur_perm[v]];
                }
                relabel[g->dst_pent] = neighbour_perm;
                code.push_back(0);
                code.push_back(static_cast<std::uint8_t>(new_idx & 0xFF));
                code.push_back(static_cast<std::uint8_t>(new_f));
            } else {
                code.push_back(1);
                code.push_back(static_cast<std::uint8_t>(it->second & 0xFF));
                Perm5 neighbour_perm = relabel[g->dst_pent];
                Perm5 inv = perm_inverse(neighbour_perm);
                code.push_back(inv[g->dst_facet]);
            }
        }
    }
    code.push_back(0xFE);
    code.push_back(static_cast<std::uint8_t>((n - static_cast<std::int32_t>(visited.size())) & 0xFF));

    // Base64 encode
    static const char* alpha =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string out;
    int val = 0, bits = 0;
    for (auto b : code) {
        val = (val << 8) | b;
        bits += 8;
        while (bits >= 6) {
            bits -= 6;
            out.push_back(alpha[(val >> bits) & 0x3F]);
        }
    }
    if (bits > 0) out.push_back(alpha[(val << (6 - bits)) & 0x3F]);
    return out;
}

std::vector<Perm5> all_perms_of_5() {
    std::vector<Perm5> out;
    std::array<std::uint8_t, 5> p{0, 1, 2, 3, 4};
    do { out.push_back(Perm5{p[0], p[1], p[2], p[3], p[4]}); }
    while (std::next_permutation(p.begin(), p.end()));
    return out;
}

}  // namespace

std::string Triangulation::isoSig() const {
    if (pentachora_.empty()) return "EMPTY";
    static const auto perms = all_perms_of_5();
    std::string best;
    for (std::int32_t start = 0; start < static_cast<std::int32_t>(pentachora_.size()); ++start) {
        for (const auto& perm : perms) {
            std::string code = canonicalRelabel(*this, start, perm);
            if (best.empty() || code < best) best = std::move(code);
        }
    }
    return best;
}

std::string Triangulation::edgeDegreeIsoSig() const {
    auto degs = edgeDegrees();
    std::ostringstream oss;
    for (std::size_t i = 0; i < degs.size(); ++i) {
        if (i) oss << ',';
        oss << degs[i];
    }
    oss << '|' << isoSig();
    return oss.str();
}

// ---------------------------------------------------------------------- //
// Canonical examples
// ---------------------------------------------------------------------- //

Triangulation s4_minimal() {
    Triangulation T("S4_2p");
    T.addPentachora(2);
    for (std::uint8_t f = 0; f < 5; ++f)
        T.glue(0, f, 1, perm_identity);
    return T;
}

Triangulation open_pentachoron() {
    Triangulation T("open_pent");
    T.addPentachora(1);
    return T;
}

Triangulation cp2_minimal() {
    // Placeholder: returns an open 4-pentachoron triangulation; the
    // closed Kuhnel CP^2 needs an explicit gluing table that is
    // populated from the Dim4Census .esig record when Regina is
    // available.
    Triangulation T("CP2_4p_placeholder");
    T.addPentachora(4);
    return T;
}

}  // namespace bhrt
