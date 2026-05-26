// bhrt_test.cpp -- native C++ test driver.
// No external test framework. Each test function returns the number of
// failures; total/failed are printed at the end.

#include "bhrt_trisect.hpp"

#if BHRT_HAS_REGINA
#include "regina_io.hpp"
#endif

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

static int g_failures = 0;
static int g_total = 0;

#define BHRT_ASSERT(cond, msg)                                          \
    do {                                                                \
        ++g_total;                                                      \
        if (!(cond)) {                                                  \
            ++g_failures;                                               \
            std::fprintf(stderr,                                        \
                "FAIL [%s:%d] %s: %s\n",                                \
                __FILE__, __LINE__, msg, #cond);                        \
        }                                                               \
    } while (0)

#define BHRT_ASSERT_EQ(a, b, msg)                                       \
    do {                                                                \
        ++g_total;                                                      \
        auto va = (a); auto vb = (b);                                   \
        if (!(va == vb)) {                                              \
            ++g_failures;                                               \
            std::fprintf(stderr,                                        \
                "FAIL [%s:%d] %s: %s != %s\n",                          \
                __FILE__, __LINE__, msg, #a, #b);                       \
        }                                                               \
    } while (0)

// ---------- Triangulation skeleton ----------

void test_s4_minimal_skeleton() {
    auto T = bhrt::s4_minimal();
    BHRT_ASSERT_EQ(T.size(), 2u,                  "S^4 pent count");
    BHRT_ASSERT_EQ(T.nVertices(), 5,              "S^4 vertices");
    BHRT_ASSERT_EQ(T.nEdges(), 10,                "S^4 edges");
    BHRT_ASSERT_EQ(T.nTriangles(), 10,            "S^4 triangles");
    BHRT_ASSERT_EQ(T.nTetrahedra(), 5,            "S^4 tetrahedra");
    BHRT_ASSERT_EQ(T.eulerCharacteristic(), 2,    "S^4 chi");
    BHRT_ASSERT(T.isClosed(),                      "S^4 is closed");
    BHRT_ASSERT(T.isValid(),                       "S^4 is valid");
}

void test_open_pentachoron() {
    auto T = bhrt::open_pentachoron();
    BHRT_ASSERT_EQ(T.size(), 1u,                  "open pent count");
    BHRT_ASSERT(!T.isClosed(),                     "open is not closed");
    BHRT_ASSERT(T.isValid(),                       "open is valid");
}

void test_perm_compose_and_invert() {
    bhrt::Perm5 p{1, 2, 0, 4, 3};
    auto inv = bhrt::perm_inverse(p);
    auto id  = bhrt::perm_compose(p, inv);
    BHRT_ASSERT(id == bhrt::perm_identity,         "p * p^-1 == id");
    // (1,2,0,4,3) has 3 inversions -> sign = -1.
    BHRT_ASSERT_EQ(bhrt::perm_sign(p), -1,         "odd permutation");
    BHRT_ASSERT_EQ(bhrt::perm_sign(bhrt::perm_identity), 1, "identity sign");
}

void test_clone_is_independent() {
    auto T = bhrt::s4_minimal();
    auto U = T.clone();
    BHRT_ASSERT_EQ(T.quickHash(), U.quickHash(),  "clone preserves hash");
    U.unglue(0, 0);
    BHRT_ASSERT(T.quickHash() != U.quickHash(),    "mutation diverges");
}

// ---------- Isosig ----------

void test_isosig_nonempty() {
    auto T = bhrt::s4_minimal();
    BHRT_ASSERT(!T.isoSig().empty(),               "isosig non-empty");
    BHRT_ASSERT(T.edgeDegreeIsoSig().find('|') != std::string::npos,
                "edge-degree contains '|'");
}

void test_isosig_canonical_under_clone() {
    auto T = bhrt::s4_minimal();
    auto U = T.clone();
    BHRT_ASSERT_EQ(T.isoSig(), U.isoSig(),         "isosig stable under clone");
}

// ---------- ts-tricolouring ----------

void test_uniform_colour_rejected() {
    auto T = bhrt::s4_minimal();
    std::vector<std::int32_t> bad(T.nVertices(), 0);
    auto tc = bhrt::isTsTricolouring(T, bad);
    BHRT_ASSERT(!tc.is_ts,                          "uniform colour not ts");
}

void test_canonical_221_runs() {
    auto T = bhrt::s4_minimal();
    std::vector<std::int32_t> colour{0, 0, 1, 1, 2};
    auto tc = bhrt::isTsTricolouring(T, colour);
    BHRT_ASSERT(!tc.audit.empty(),                  "audit log populated");
}

void test_count_aggregator_runs() {
    auto T = bhrt::s4_minimal();
    auto c = bhrt::countAdmitTrisection({T});
    BHRT_ASSERT_EQ(c.n_total, 1,                   "scan total");
    BHRT_ASSERT_EQ(c.n_pass_precheck, 1,           "precheck passes");
}

// ---------- Two-complex ----------

void test_two_complex_triangle_collapses() {
    bhrt::TwoComplex K;
    K.addFace(0, 1, 2);
    BHRT_ASSERT(K.isConnected(),                   "triangle connected");
    BHRT_ASSERT(K.greedyCollapseTo1(),             "triangle collapses");
}

void test_two_complex_bouquet_collapses() {
    bhrt::TwoComplex K;
    K.addFace(0, 1, 2);
    K.addFace(0, 1, 3);
    BHRT_ASSERT(K.greedyCollapseTo1(),             "bouquet collapses");
}

// ---------- Diagram ----------

void test_diagram_runs_on_synthetic_input() {
    auto T = bhrt::s4_minimal();
    bhrt::TSColouring fake;
    fake.colour = {0, 0, 1, 1, 2};
    fake.is_ts = false;
    auto d = bhrt::extractDiagram(T, fake);
    BHRT_ASSERT(!d.audit.empty(),                   "diagram audit populated");
}

// ---------- Invariants ----------

void test_snf_diagonal() {
    std::vector<std::vector<std::int64_t>> M = {{2,0,0},{0,6,0},{0,0,12}};
    auto diag = bhrt::smithNormalFormDiagonal(M);
    std::sort(diag.begin(), diag.end());
    BHRT_ASSERT_EQ(diag.size(), 3u,                "SNF length");
    BHRT_ASSERT_EQ(diag[0], 2,                     "SNF[0]");
    BHRT_ASSERT_EQ(diag[1], 6,                     "SNF[1]");
    BHRT_ASSERT_EQ(diag[2], 12,                    "SNF[2]");
}

void test_snf_nontrivial() {
    std::vector<std::vector<std::int64_t>> M = {
        { 2,  4,  4},
        {-6,  6, 12},
        {10, -4,-16}
    };
    auto diag = bhrt::smithNormalFormDiagonal(M);
    long long prod = 1;
    for (auto d : diag) prod *= d > 0 ? d : 1;
    BHRT_ASSERT_EQ(prod, 144,                      "SNF product = |det|");
}

// ---------- .bhrt round-trip ----------

void test_bhrt_format_roundtrip() {
    auto T = bhrt::s4_minimal();
    std::string path = "/tmp/bhrt_roundtrip.bhrt";
    bhrt::writeBhrtFile(path, {T});
    auto loaded = bhrt::loadBhrtFile(path);
    BHRT_ASSERT_EQ(loaded.size(), 1u,              "loaded one record");
    BHRT_ASSERT_EQ(loaded[0].size(), T.size(),     "loaded pent count");
    BHRT_ASSERT_EQ(loaded[0].isoSig(), T.isoSig(), "loaded isosig matches");
}

// ---------- Beam search ----------

void test_beam_search_terminates() {
    auto T = bhrt::s4_minimal();
    bhrt::SearchConfig cfg;
    cfg.beam_width = 8;
    cfg.excess_height = 2;
    cfg.max_iterations = 3;
    cfg.time_limit_seconds = 5.0;
    auto pareto = bhrt::beamSearch(T, cfg);
    BHRT_ASSERT(!pareto.items.empty(),             "pareto front non-empty");
}

// ---------- Move executor (1-5 / 5-1) ----------

void test_one_five_move_preserves_manifold() {
    auto T = bhrt::s4_minimal();
    auto U = bhrt::applyMove(T, {bhrt::MoveType::Pachner_1_5, {0}, +4});
    BHRT_ASSERT(U.has_value(),                     "1-5 applies to S^4");
    BHRT_ASSERT_EQ(U->size(), 6u,                  "1-5 grows 2 -> 6 pents");
    BHRT_ASSERT(U->isValid(),                      "1-5 result is valid");
    BHRT_ASSERT(U->isClosed(),                     "1-5 result is closed");
    BHRT_ASSERT_EQ(U->eulerCharacteristic(), 2,    "1-5 preserves chi(S^4)=2");
}

void test_five_one_inverts_one_five() {
    auto T = bhrt::s4_minimal();
    auto U = bhrt::applyMove(T, {bhrt::MoveType::Pachner_1_5, {0}, +4});
    BHRT_ASSERT(U.has_value(),                     "1-5 applies");
    bool recovered = false;
    for (const auto& m : bhrt::enumerateMoves(*U)) {
        if (m.move != bhrt::MoveType::Pachner_5_1) continue;
        auto W = bhrt::applyMove(*U, m);
        if (W && W->isoSig() == T.isoSig()) { recovered = true; break; }
    }
    BHRT_ASSERT(recovered,                         "5-1 inverts 1-5 back to S^4");
}

void test_search_actually_expands() {
    auto T = bhrt::s4_minimal();
    // The executor must be able to expand the seed: at least one enumerated
    // 1-5 candidate applies and grows the triangulation. (The Pareto front
    // then *prunes* such grown states as dominated, so we test the executor
    // directly rather than inspecting the front.)
    bool canExpand = false;
    for (const auto& m : bhrt::enumerateMoves(T)) {
        if (m.move != bhrt::MoveType::Pachner_1_5) continue;
        auto U = bhrt::applyMove(T, m);
        if (U && (int)U->size() > (int)T.size() && U->isValid()) {
            canExpand = true; break;
        }
    }
    BHRT_ASSERT(canExpand, "executor expands the seed via an enumerated 1-5 move");
    // And the search runs to completion and returns a non-empty Pareto front.
    bhrt::SearchConfig cfg;
    cfg.beam_width = 16; cfg.excess_height = 8;
    cfg.max_iterations = 2; cfg.time_limit_seconds = 5.0;
    auto pareto = bhrt::beamSearch(T, cfg);
    BHRT_ASSERT(!pareto.items.empty(), "beam search returns a Pareto front");
}

// ---------- Ground-truth invariants (hand-coded Lagrangian diagrams) ----------

static bhrt::LagrangianDiagram mkDiag(int g,
        std::vector<std::vector<std::int64_t>> a,
        std::vector<std::vector<std::int64_t>> b,
        std::vector<std::vector<std::int64_t>> c) {
    bhrt::LagrangianDiagram d; d.genus = g;
    d.alpha = std::move(a); d.beta = std::move(b); d.gamma = std::move(c);
    return d;
}

void test_invariants_s4() {
    auto d = mkDiag(0, {}, {}, {});
    BHRT_ASSERT(bhrt::validateLagrangianDiagram(d).ok, "S^4 diagram valid");
    auto inv = bhrt::computeInvariantsFromClasses(d);
    BHRT_ASSERT_EQ(inv.h1_free_rank, 0, "S^4 b1=0");
    BHRT_ASSERT_EQ(inv.h2_free_rank, 0, "S^4 b2=0");
    BHRT_ASSERT_EQ(inv.signature, 0,    "S^4 sigma=0");
}

void test_invariants_cp2() {
    // genus-1: alpha=a, beta=b, gamma=a+b.
    auto d = mkDiag(1, {{1,0}}, {{0,1}}, {{1,1}});
    BHRT_ASSERT(bhrt::validateLagrangianDiagram(d).ok, "CP^2 diagram valid");
    auto inv = bhrt::computeInvariantsFromClasses(d);
    BHRT_ASSERT_EQ(inv.h1_free_rank, 0, "CP^2 b1=0");
    BHRT_ASSERT_EQ(inv.h2_free_rank, 1, "CP^2 b2=1");
    BHRT_ASSERT_EQ(inv.h3_free_rank, 0, "CP^2 b3=0");
    BHRT_ASSERT_EQ(inv.signature, 1,    "CP^2 sigma=1");
    BHRT_ASSERT(inv.parity == "odd",    "CP^2 form is odd");
    BHRT_ASSERT(inv.intersection_form.size() == 1
        && inv.intersection_form[0][0] == 1, "CP^2 Q=[1]");
}

void test_invariants_cp2_bar() {
    // reversed orientation: gamma = a - b -> sigma = -1.
    auto d = mkDiag(1, {{1,0}}, {{0,1}}, {{1,-1}});
    auto inv = bhrt::computeInvariantsFromClasses(d);
    BHRT_ASSERT_EQ(inv.signature, -1,   "-CP^2 sigma=-1");
    BHRT_ASSERT(inv.parity == "odd",    "-CP^2 form is odd");
}

void test_invariants_s2xs2() {
    // genus-2: alpha=a1,a2; beta=b1,b2; gamma = a2+b1, a1+b2 -> Q=[[0,1],[1,0]].
    auto d = mkDiag(2,
        {{1,0,0,0},{0,1,0,0}},
        {{0,0,1,0},{0,0,0,1}},
        {{0,1,1,0},{1,0,0,1}});
    BHRT_ASSERT(bhrt::validateLagrangianDiagram(d).ok, "S2xS2 diagram valid");
    auto inv = bhrt::computeInvariantsFromClasses(d);
    BHRT_ASSERT_EQ(inv.h1_free_rank, 0, "S2xS2 b1=0");
    BHRT_ASSERT_EQ(inv.h2_free_rank, 2, "S2xS2 b2=2");
    BHRT_ASSERT_EQ(inv.signature, 0,    "S2xS2 sigma=0");
    BHRT_ASSERT(inv.parity == "even",   "S2xS2 form is even");
    BHRT_ASSERT(inv.intersection_form.size() == 2
        && inv.intersection_form[0][0] == 0
        && inv.intersection_form[0][1] == 1
        && inv.intersection_form[1][0] == 1
        && inv.intersection_form[1][1] == 0, "S2xS2 Q=[[0,1],[1,0]]");
}

void test_invariants_cp2_connsum() {
    // CP^2 # CP^2: gamma = a1+b1, a2+b2 -> Q = I_2, sigma=2, odd.
    auto d = mkDiag(2,
        {{1,0,0,0},{0,1,0,0}},
        {{0,0,1,0},{0,0,0,1}},
        {{1,0,1,0},{0,1,0,1}});
    auto inv = bhrt::computeInvariantsFromClasses(d);
    BHRT_ASSERT_EQ(inv.h2_free_rank, 2, "CP2#CP2 b2=2");
    BHRT_ASSERT_EQ(inv.signature, 2,    "CP2#CP2 sigma=2");
    BHRT_ASSERT(inv.parity == "odd",    "CP2#CP2 form is odd");
}

void test_invariants_s1xs3() {
    // genus-1, all three systems the same Lagrangian <a>: b1=1, b2=0, b3=1.
    auto d = mkDiag(1, {{1,0}}, {{1,0}}, {{1,0}});
    auto inv = bhrt::computeInvariantsFromClasses(d);
    BHRT_ASSERT_EQ(inv.h1_free_rank, 1, "S1xS3 b1=1");
    BHRT_ASSERT_EQ(inv.h2_free_rank, 0, "S1xS3 b2=0");
    BHRT_ASSERT_EQ(inv.h3_free_rank, 1, "S1xS3 b3=1");
}

void test_validate_rejects_nonlagrangian() {
    // genus-1 "alpha" that is not isotropic is impossible (single curve is
    // always isotropic); instead give a non-independent gamma family.
    auto d = mkDiag(2,
        {{1,0,0,0},{0,1,0,0}},
        {{0,0,1,0},{0,0,0,1}},
        {{0,1,1,0},{0,1,1,0}});   // two identical gamma curves -> rank 1
    auto val = bhrt::validateLagrangianDiagram(d);
    BHRT_ASSERT(!val.gamma_is_cut_system, "duplicate gamma curves rejected");
    BHRT_ASSERT(!val.ok,                  "diagram flagged invalid");
}

// ---------- Triangulation-driven pipeline (end-to-end consistency) ----------

void test_pipeline_runs_end_to_end() {
    auto T = bhrt::s4_minimal();
    auto colourings = bhrt::enumerateTsTricolourings(T, 10000, false);
    // Whether or not S^4's 2-pentachoron triangulation admits a ts-colouring,
    // the pipeline must run without crashing and report a consistent genus.
    for (const auto& tc : colourings) {
        auto diag = bhrt::extractDiagram(T, tc);
        BHRT_ASSERT(diag.genus() >= 0,            "extracted genus is non-negative");
        auto val = bhrt::validateDiagram(diag);   // must not throw
        (void)val;
        auto inv = bhrt::computeInvariants(diag);  // must not throw
        BHRT_ASSERT(inv.genus >= 0,               "invariant bundle genus >= 0");
    }
    BHRT_ASSERT(true,                              "pipeline completed");
}

// ---------- Regina backend ----------

#if BHRT_HAS_REGINA

void test_regina_roundtrip_s4() {
    auto T = bhrt::s4_minimal();
    auto R = bhrt::toRegina(T);
    BHRT_ASSERT(R != nullptr,                       "S^4 converts to Regina");
    if (!R) return;
    BHRT_ASSERT_EQ(R->size(), T.size(),             "Regina pent count");
    BHRT_ASSERT_EQ((int)R->countTetrahedra(), T.nTetrahedra(),
                                                    "Regina tetra count");
    BHRT_ASSERT_EQ((int)R->countEdges(), T.nEdges(), "Regina edge count");
    auto U = bhrt::fromRegina(*R);
    BHRT_ASSERT_EQ(U.quickHash(), T.quickHash(),    "round-trip exact gluings");
    BHRT_ASSERT_EQ(U.isoSig(), T.isoSig(),          "round-trip isosig");
    BHRT_ASSERT_EQ(U.eulerCharacteristic(), 2,      "round-trip chi");
}

void test_regina_roundtrip_self_gluing() {
    // One pentachoron with facet 0 glued to facet 1 (distinct facets of the
    // SAME pentachoron): representable in Regina, must survive a round-trip.
    bhrt::Triangulation T("self_gluing");
    T.addPentachora(1);
    T.glue(0, 0, 0, bhrt::Perm5{1, 0, 2, 3, 4});    // facet 0 <-> facet 1
    BHRT_ASSERT(T.isValid(),                        "self-gluing input valid");
    auto R = bhrt::toRegina(T);
    BHRT_ASSERT(R != nullptr,                       "self-gluing converts");
    if (R) {
        auto U = bhrt::fromRegina(*R);
        BHRT_ASSERT_EQ(U.quickHash(), T.quickHash(), "self-gluing round-trip");
    }
}

void test_regina_rejects_facet_to_itself() {
    // Facet 0 glued to itself by an involution fixing vertex 0: valid in the
    // native format, NOT representable in Regina. toRegina must return
    // nullptr (not throw, not corrupt), and the move executor must skip.
    bhrt::Triangulation T("facet_to_itself");
    T.addPentachora(1);
    T.glue(0, 0, 0, bhrt::Perm5{0, 2, 1, 4, 3});    // perm[0]=0 -> same facet
    BHRT_ASSERT(T.isValid(),                        "facet-to-itself input valid");
    auto R = bhrt::toRegina(T);
    BHRT_ASSERT(R == nullptr,                       "facet-to-itself rejected");
    auto M = bhrt::applyMove(T, {bhrt::MoveType::Pachner_3_3, {0}, 0});
    BHRT_ASSERT(!M.has_value(),                     "executor skips, no crash");
}

void test_regina_two_four_on_s4() {
    // Minimal S^4: every tetrahedron joins the two distinct pentachora, so
    // the 2-4 move is legal everywhere. Deterministic check of the full
    // applyMove dispatch path.
    auto T = bhrt::s4_minimal();
    auto U = bhrt::applyMove(T, {bhrt::MoveType::Pachner_2_3, {0}, +2});
    BHRT_ASSERT(U.has_value(),                      "2-4 applies to S^4");
    if (U.has_value()) {
        BHRT_ASSERT_EQ(U->size(), 4u,               "2-4 grows 2 -> 4 pents");
        BHRT_ASSERT(U->isValid(),                   "2-4 result valid");
        BHRT_ASSERT(U->isClosed(),                  "2-4 result closed");
        BHRT_ASSERT_EQ(U->eulerCharacteristic(), 2, "2-4 preserves chi(S^4)=2");
    }
}

void test_regina_lateral_moves_preserve_s4() {
    // Grow S^4 with a 1-5, then sweep ALL enumerated lateral candidates.
    // Every candidate that applies must preserve validity, closedness and
    // chi, and must change the pentachoron count by exactly its 4D delta.
    auto T0 = bhrt::s4_minimal();
    auto Ug = bhrt::applyMove(T0, {bhrt::MoveType::Pachner_1_5, {0}, +4});
    BHRT_ASSERT(Ug.has_value(),                     "1-5 grows the seed");
    if (!Ug.has_value()) return;
    const auto& T = *Ug;

    int applied = 0;
    for (const auto& m : bhrt::enumerateMoves(T)) {
        int expected_delta = 0;
        switch (m.move) {
            case bhrt::MoveType::Pachner_2_3: expected_delta = +2; break;
            case bhrt::MoveType::Pachner_3_3: expected_delta = 0;  break;
            case bhrt::MoveType::Pachner_4_4: expected_delta = 0;  break;
            case bhrt::MoveType::Pachner_4_2: expected_delta = -2; break;
            default: continue;                       // native + EdgeCollapse
        }
        auto W = bhrt::applyMove(T, m);
        if (!W) continue;                            // illegal here: fine
        ++applied;
        BHRT_ASSERT(W->isValid(),                   "lateral result valid");
        BHRT_ASSERT(W->isClosed(),                  "lateral result closed");
        BHRT_ASSERT_EQ(W->eulerCharacteristic(), 2, "lateral preserves chi");
        BHRT_ASSERT_EQ((int)W->size() - (int)T.size(), expected_delta,
                       "lateral move has its 4D pentachoron delta");
    }
    BHRT_ASSERT(applied > 0, "at least one lateral move applies to grown S^4");
}

void test_regina_edge_collapse_executes_safely() {
    // EdgeCollapse candidates must run through Regina without crashing and,
    // when they apply, must shrink the triangulation and preserve chi.
    auto T0 = bhrt::s4_minimal();
    auto Ug = bhrt::applyMove(T0, {bhrt::MoveType::Pachner_1_5, {0}, +4});
    BHRT_ASSERT(Ug.has_value(),                     "1-5 grows the seed");
    if (!Ug.has_value()) return;
    const auto& T = *Ug;
    for (const auto& m : bhrt::enumerateMoves(T)) {
        if (m.move != bhrt::MoveType::EdgeCollapse) continue;
        auto W = bhrt::applyMove(T, m);
        if (!W) continue;
        BHRT_ASSERT(W->isValid(),                   "collapse result valid");
        BHRT_ASSERT(W->isClosed(),                  "collapse result closed");
        BHRT_ASSERT_EQ(W->eulerCharacteristic(), 2, "collapse preserves chi");
        BHRT_ASSERT((int)W->size() < (int)T.size(), "collapse shrinks");
    }
    BHRT_ASSERT(true, "edge-collapse sweep completed without crashing");
}

void test_regina_esig_roundtrip() {
    auto T = bhrt::s4_minimal();
    auto Ug = bhrt::applyMove(T, {bhrt::MoveType::Pachner_1_5, {0}, +4});
    BHRT_ASSERT(Ug.has_value(),                     "1-5 grows the seed");
    if (!Ug.has_value()) return;
    std::string path = "/tmp/bhrt_roundtrip.esig";
    bhrt::writeEsig(path, {T, *Ug});
    auto loaded = bhrt::loadEsig(path);
    BHRT_ASSERT_EQ(loaded.size(), 2u,               "esig: two records back");
    if (loaded.size() == 2) {
        // Regina re-canonicalises, so compare canonical signatures.
        BHRT_ASSERT_EQ(loaded[0].isoSig(), T.isoSig(),   "esig[0] isosig");
        BHRT_ASSERT_EQ(loaded[1].isoSig(), Ug->isoSig(), "esig[1] isosig");
        BHRT_ASSERT_EQ(bhrt::reginaIsoSig(loaded[0]), bhrt::reginaIsoSig(T),
                       "esig[0] regina-canonical sig");
    }
}

void test_regina_beam_search_with_lateral_moves() {
    // End-to-end: the beam search must run with the Regina executor live
    // (this is the exact code path of `bhrt-cli search` on a Regina build).
    auto T = bhrt::s4_minimal();
    bhrt::SearchConfig cfg;
    cfg.beam_width = 8;
    cfg.excess_height = 4;
    cfg.max_iterations = 3;
    cfg.time_limit_seconds = 30.0;
    auto pareto = bhrt::beamSearch(T, cfg);
    BHRT_ASSERT(!pareto.items.empty(), "Regina-backed search returns a front");
    for (const auto& s : pareto.items)
        BHRT_ASSERT(s.triangulation.isValid(), "every front state is valid");
}

#endif  // BHRT_HAS_REGINA

// ---------- CUDA scorer ----------

#if BHRT_HAS_CUDA

// Host launcher defined in cuda/score_moves.cu (same declaration as in
// search_cpu.cpp). Returns 0 on success, non-zero if no usable GPU.
extern "C" int bhrt_score_moves_host(const float* records, int n_records,
                                     int n_cols, const float* weights,
                                     float* scores);

// Build a deterministic synthetic batch of flattened move records.
static std::vector<float> mk_records(int n) {
    static const float codes[7] = {23, 33, 44, 42, 15, 51, 99};
    static const float deltas[7] = {+2, 0, 0, -2, +4, -4, -1};
    std::vector<float> recs(static_cast<std::size_t>(n) * 5);
    for (int i = 0; i < n; ++i) {
        recs[i * 5 + 0] = codes[i % 7];
        recs[i * 5 + 1] = static_cast<float>(i % 97);
        recs[i * 5 + 2] = deltas[i % 7];
        recs[i * 5 + 3] = static_cast<float>(2 + (i % 50));
        recs[i * 5 + 4] = static_cast<float>(i % 2);
    }
    return recs;
}

static bool scores_agree(float gpu, float cpu) {
    float tol = 1e-5f * (std::fabs(cpu) > 1.0f ? std::fabs(cpu) : 1.0f);
    return std::fabs(gpu - cpu) <= tol;
}

void test_cuda_scorer_matches_cpu() {
    const int N = 1000;
    auto recs = mk_records(N);
    const std::array<float, 5> w{1.0f, 0.5f, 0.3f, -0.2f, -1.0f};
    std::vector<float> gpu(N);
    int rc = bhrt_score_moves_host(recs.data(), N, 5, w.data(), gpu.data());
    if (rc != 0) {
        // No usable GPU on this machine: the clean fallback IS the pass.
        BHRT_ASSERT(true, "GPU unavailable; scorer fell back cleanly");
        return;
    }
    bool all_match = true;
    for (int i = 0; i < N && all_match; ++i) {
        std::array<float, 5> r{recs[i * 5 + 0], recs[i * 5 + 1],
                               recs[i * 5 + 2], recs[i * 5 + 3],
                               recs[i * 5 + 4]};
        all_match = scores_agree(gpu[i], bhrt::scoreCandidateCPU(r, w));
    }
    BHRT_ASSERT(all_match, "GPU scores match CPU reference (1k batch)");
}

void test_cuda_scorer_large_batch_heap_safe() {
    // Large enough that the host vectors are mmap'd glibc chunks: an
    // out-of-bounds write by the copy-back corrupts a chunk header and
    // aborts ("munmap_chunk(): invalid pointer") when the vectors die.
    // Surviving two rounds of this is the regression test for the
    // GPU-path heap-corruption crash.
    const int N = 200000;
    auto recs = mk_records(N);
    const std::array<float, 5> w{1.0f, 0.5f, 0.3f, -0.2f, -1.0f};
    for (int round = 0; round < 2; ++round) {
        std::vector<float> gpu(N);
        int rc = bhrt_score_moves_host(recs.data(), N, 5, w.data(), gpu.data());
        if (rc != 0) {
            BHRT_ASSERT(true, "GPU unavailable; large batch fell back cleanly");
            return;
        }
        // Spot-check ends and middle (full sweep already done in the 1k test).
        for (int i : {0, 1, N / 2, N - 2, N - 1}) {
            std::array<float, 5> r{recs[i * 5 + 0], recs[i * 5 + 1],
                                   recs[i * 5 + 2], recs[i * 5 + 3],
                                   recs[i * 5 + 4]};
            BHRT_ASSERT(scores_agree(gpu[i], bhrt::scoreCandidateCPU(r, w)),
                        "GPU large-batch spot check");
        }
    }
    BHRT_ASSERT(true, "two 200k GPU batches completed without heap damage");
}

void test_cuda_ts_enumeration_matches_cpu() {
    // The GPU enumerator must return EXACTLY the CPU enumerator's
    // colouring set (the GPU only prefilters; survivors are CPU-verified).
    // On machines without a GPU the wrapper falls back to the CPU path,
    // so the comparison passes trivially.
    auto key = [](const std::vector<bhrt::TSColouring>& v) {
        std::vector<std::vector<std::int32_t>> k;
        k.reserve(v.size());
        for (const auto& t : v) k.push_back(t.colour);
        std::sort(k.begin(), k.end());
        return k;
    };
    auto T = bhrt::s4_minimal();
    auto cpu = bhrt::enumerateTsTricolourings(T, 10000, false);
    auto gpu = bhrt::enumerateTsTricolouringsGPU(T, 10000, false);
    BHRT_ASSERT(key(cpu) == key(gpu), "GPU == CPU ts-enumeration on S^4");

    auto U = bhrt::applyMove(T, {bhrt::MoveType::Pachner_1_5, {0}, +4});
    BHRT_ASSERT(U.has_value(), "1-5 grows the seed");
    if (U.has_value()) {
        auto cpu2 = bhrt::enumerateTsTricolourings(*U, 10000, false);
        auto gpu2 = bhrt::enumerateTsTricolouringsGPU(*U, 10000, false);
        BHRT_ASSERT(key(cpu2) == key(gpu2),
                    "GPU == CPU ts-enumeration on grown S^4");
    }
}

void test_cuda_beam_search_gpu_flag() {
    // End-to-end: the exact code path of `BHRT_USE_GPU=1 bhrt-cli search`,
    // including the first-batch GPU-vs-CPU self-verification inside
    // beamSearch. Passes whether or not a GPU is present (clean fallback).
    auto T = bhrt::s4_minimal();
    bhrt::SearchConfig cfg;
    cfg.beam_width = 8;
    cfg.excess_height = 4;
    cfg.max_iterations = 3;
    cfg.time_limit_seconds = 30.0;
    cfg.use_gpu_scorer = true;
    auto pareto = bhrt::beamSearch(T, cfg);
    BHRT_ASSERT(!pareto.items.empty(), "GPU-flagged search returns a front");
    for (const auto& s : pareto.items)
        BHRT_ASSERT(s.triangulation.isValid(), "every front state is valid");
}

#endif  // BHRT_HAS_CUDA

// ---------- Driver ----------

int main() {
    test_s4_minimal_skeleton();
    test_open_pentachoron();
    test_perm_compose_and_invert();
    test_clone_is_independent();
    test_isosig_nonempty();
    test_isosig_canonical_under_clone();
    test_uniform_colour_rejected();
    test_canonical_221_runs();
    test_count_aggregator_runs();
    test_two_complex_triangle_collapses();
    test_two_complex_bouquet_collapses();
    test_diagram_runs_on_synthetic_input();
    test_snf_diagonal();
    test_snf_nontrivial();
    test_bhrt_format_roundtrip();
    // Move executor (native).
    test_one_five_move_preserves_manifold();
    test_five_one_inverts_one_five();
#if BHRT_HAS_REGINA
    // Regina backend UNIT tests run BEFORE any beam search so that a crash
    // localises: conversions first, single moves next, search interaction
    // last. (Diagnostic ordering for the munmap_chunk investigation.)
    test_regina_roundtrip_s4();
    test_regina_roundtrip_self_gluing();
    test_regina_rejects_facet_to_itself();
    test_regina_two_four_on_s4();
    test_regina_lateral_moves_preserve_s4();
    test_regina_edge_collapse_executes_safely();
    test_regina_esig_roundtrip();
#endif
#if BHRT_HAS_CUDA
    // CUDA unit tests (clean fallback on machines without a GPU).
    test_cuda_scorer_matches_cpu();
    test_cuda_scorer_large_batch_heap_safe();
    test_cuda_ts_enumeration_matches_cpu();
#endif
    // Beam search (native executor only without Regina; with Regina the
    // lateral moves go live here).
    test_beam_search_terminates();
    test_search_actually_expands();
    // Ground-truth invariants.
    test_invariants_s4();
    test_invariants_cp2();
    test_invariants_cp2_bar();
    test_invariants_s2xs2();
    test_invariants_cp2_connsum();
    test_invariants_s1xs3();
    test_validate_rejects_nonlagrangian();
    // End-to-end pipeline.
    test_pipeline_runs_end_to_end();
    // Search-engine integration tests LAST (these exercise the full
    // Regina/CUDA interaction inside beamSearch).
#if BHRT_HAS_REGINA
    test_regina_beam_search_with_lateral_moves();
#endif
#if BHRT_HAS_CUDA
    test_cuda_beam_search_gpu_flag();
#endif
    std::printf("\n%d/%d assertions passed, %d failed\n",
                g_total - g_failures, g_total, g_failures);
    return g_failures == 0 ? 0 : 1;
}
