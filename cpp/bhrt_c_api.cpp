// bhrt_c_api.cpp -- C API wrapper implementation.

#include "bhrt_c_api.h"
#include "bhrt_trisect.hpp"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <exception>
#include <new>
#include <vector>

struct bhrt_triangulation {
    bhrt::Triangulation t;
};

extern "C" {

bhrt_triangulation* bhrt_new_triangulation(void) {
    return new (std::nothrow) bhrt_triangulation();
}

void bhrt_free_triangulation(bhrt_triangulation* t) {
    delete t;
}

bhrt_triangulation* bhrt_s4_minimal(void) {
    auto* w = new (std::nothrow) bhrt_triangulation();
    if (!w) return nullptr;
    w->t = bhrt::s4_minimal();
    return w;
}

void bhrt_add_pentachora(bhrt_triangulation* t, size_t count) {
    if (!t) return;
    t->t.addPentachora(count);
}

void bhrt_glue(bhrt_triangulation* t,
               int32_t src, uint8_t src_facet, int32_t dst,
               const uint8_t perm5[5]) {
    if (!t) return;
    bhrt::Perm5 p{perm5[0], perm5[1], perm5[2], perm5[3], perm5[4]};
    t->t.glue(src, src_facet, dst, p);
}

size_t  bhrt_size(const bhrt_triangulation* t)        { return t ? t->t.size() : 0; }
int32_t bhrt_n_vertices(const bhrt_triangulation* t)  { return t ? t->t.nVertices() : 0; }
int32_t bhrt_n_edges(const bhrt_triangulation* t)     { return t ? t->t.nEdges() : 0; }
int32_t bhrt_n_triangles(const bhrt_triangulation* t) { return t ? t->t.nTriangles() : 0; }
int32_t bhrt_n_tetrahedra(const bhrt_triangulation* t){ return t ? t->t.nTetrahedra() : 0; }
int32_t bhrt_euler(const bhrt_triangulation* t)       { return t ? t->t.eulerCharacteristic() : 0; }
int     bhrt_is_closed(const bhrt_triangulation* t)   { return t && t->t.isClosed() ? 1 : 0; }
int     bhrt_is_valid(const bhrt_triangulation* t)    { return t && t->t.isValid() ? 1 : 0; }

static size_t copy_to(const std::string& s, char* out, size_t cap) {
    if (out && cap) {
        size_t n = s.size() < cap - 1 ? s.size() : cap - 1;
        std::memcpy(out, s.data(), n);
        out[n] = '\0';
    }
    return s.size();
}

size_t bhrt_isosig(const bhrt_triangulation* t, char* out, size_t cap) {
    if (!t) { if (out && cap) out[0] = '\0'; return 0; }
    return copy_to(t->t.isoSig(), out, cap);
}

size_t bhrt_edge_degree_isosig(const bhrt_triangulation* t, char* out, size_t cap) {
    if (!t) { if (out && cap) out[0] = '\0'; return 0; }
    return copy_to(t->t.edgeDegreeIsoSig(), out, cap);
}

void bhrt_scan_single(const bhrt_triangulation* t,
                       int collapse_budget,
                       int64_t counters[6]) {
    for (int i = 0; i < 6; ++i) counters[i] = 0;
    if (!t) return;
    std::vector<bhrt::Triangulation> v{t->t};
    auto c = bhrt::countAdmitTrisection(v, collapse_budget);
    counters[0] = c.n_total;
    counters[1] = c.n_pass_precheck;
    counters[2] = c.n_with_tricolouring;
    counters[3] = c.n_with_c_tricolouring;
    counters[4] = c.n_with_ts_tricolouring;
    counters[5] = c.n_ts_tricolourings_total;
}

int bhrt_beam_search_best_pentachora(const bhrt_triangulation* t,
                                      int beam_width,
                                      int excess_height,
                                      int max_iterations,
                                      uint64_t seed) {
    if (!t) return -1;
    bhrt::SearchConfig cfg;
    cfg.beam_width = beam_width;
    cfg.excess_height = excess_height;
    cfg.max_iterations = max_iterations;
    cfg.seed = seed;
    auto pareto = bhrt::beamSearch(t->t, cfg);
    int best = static_cast<int>(t->t.size());
    for (auto& s : pareto.items) best = std::min(best, s.pentachora);
    return best;
}

int bhrt_scan_file(const char* path, int collapse_budget,
                    int64_t counters[6]) {
    for (int i = 0; i < 6; ++i) counters[i] = 0;
    try {
        auto tris = bhrt::loadBhrtFile(path);
        auto c = bhrt::countAdmitTrisection(tris, collapse_budget);
        counters[0] = c.n_total;
        counters[1] = c.n_pass_precheck;
        counters[2] = c.n_with_tricolouring;
        counters[3] = c.n_with_c_tricolouring;
        counters[4] = c.n_with_ts_tricolouring;
        counters[5] = c.n_ts_tricolourings_total;
        return static_cast<int>(tris.size());
    } catch (const std::exception& e) {
        std::fprintf(stderr, "bhrt_scan_file: %s\n", e.what());
        return -1;
    }
}

int bhrt_dump_isosigs(const char* path) {
    try {
        auto tris = bhrt::loadBhrtFile(path);
        for (const auto& T : tris) {
            auto sig = T.edgeDegreeIsoSig();
            std::printf("%s\n", sig.c_str());
        }
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "bhrt_dump_isosigs: %s\n", e.what());
        return -1;
    }
}

int bhrt_write_demo(const char* path) {
    try {
        std::vector<bhrt::Triangulation> tris{bhrt::s4_minimal()};
        bhrt::writeBhrtFile(path, tris);
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "bhrt_write_demo: %s\n", e.what());
        return -1;
    }
}

}  // extern "C"
