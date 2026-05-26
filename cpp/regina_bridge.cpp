// regina_bridge.cpp -- Regina-backed lateral moves + shared conversions.
//
// Compiled only when BHRT_HAS_REGINA is defined. This file owns the single
// copy of the bhrt <-> Regina conversion code (declared in regina_io.hpp,
// also used by tri_io.cpp and the test driver) and provides
// applyMoveRegina(), the executor for every move the built-in engine does
// not implement natively:
//
//   MoveType::Pachner_2_3  ->  Regina pachner on a tetrahedron  (2-4,  +2 pents)
//   MoveType::Pachner_3_3  ->  Regina pachner on a triangle     (3-3,   0 pents)
//   MoveType::Pachner_4_4  ->  Regina move44 on an edge         (4-4,   0 pents)
//   MoveType::Pachner_4_2  ->  Regina pachner on an edge        (4-2,  -2 pents)
//   MoveType::EdgeCollapse ->  Regina collapseEdge              (-deg(e) pents)
//
// The native 1-5 / 5-1 moves are handled by search_cpu.cpp and never reach
// here. All Regina move functions check legality and return true iff the
// move was performed; an illegal locator yields nullopt. We operate on a
// fresh copy and self-verify the result (validity, closedness, Euler
// characteristic) before returning, so a surprise in either direction of
// the conversion degrades to "move skipped", never to a corrupted state.
//
// Locator semantics: enumerateMoves() emits NATIVE skeleton face ids
// (tetrahedron / triangle / edge equivalence classes). Native and Regina
// face *counts* agree for a valid triangulation but the *orderings* need
// not, so the locator is reduced mod the Regina face count and used as a
// deterministic selector. Any legal Pachner / collapse move preserves the
// PL type, so this affects only which move is tried, never correctness.

#include "bhrt_trisect.hpp"

#if BHRT_HAS_REGINA

#include "regina_io.hpp"

#include <memory>
#include <stdexcept>
#include <vector>

namespace bhrt {

// ---------------------------------------------------------------------- //
// Shared conversions (declared in regina_io.hpp)
// ---------------------------------------------------------------------- //

std::unique_ptr<regina::Triangulation<4>> toRegina(const Triangulation& T) {
    auto R = std::make_unique<regina::Triangulation<4>>();
    std::vector<regina::Pentachoron<4>*> pents;
    pents.reserve(T.size());
    for (std::size_t i = 0; i < T.size(); ++i)
        pents.push_back(R->newPentachoron());
    for (std::size_t p = 0; p < T.size(); ++p) {
        const auto& src = T.pentachora()[p];
        for (std::uint8_t f = 0; f < 5; ++f) {
            const auto& g = src.gluings[f];
            if (!g.has_value()) continue;
            // Realise each gluing exactly once (join() glues both sides and
            // throws if a facet is already joined):
            //   * cross-pentachoron pairs from the smaller pentachoron;
            //   * self-gluings (dst_pent == p, distinct facets) from the
            //     smaller facet;
            //   * a facet glued to ITSELF has no Regina representation.
            if (g->dst_pent < static_cast<std::int32_t>(p)) continue;
            if (g->dst_pent == static_cast<std::int32_t>(p)) {
                if (g->dst_facet == f) return nullptr;   // unrepresentable
                if (g->dst_facet < f) continue;          // already realised
            }
            regina::Perm<5> perm(
                g->perm[0], g->perm[1], g->perm[2], g->perm[3], g->perm[4]);
            pents[p]->join(f, pents[g->dst_pent], perm);
        }
    }
    return R;
}

Triangulation fromRegina(const regina::Triangulation<4>& R) {
    Triangulation T;
    for (std::size_t i = 0; i < R.size(); ++i) T.newPentachoron();
    for (std::size_t p = 0; p < R.size(); ++p) {
        const auto* P = R.pentachoron(p);
        for (std::uint8_t f = 0; f < 5; ++f) {
            const auto* adj = P->adjacentPentachoron(f);
            if (!adj) continue;                          // boundary facet
            auto perm = P->adjacentGluing(f);
            // Mirror toRegina: realise each gluing once. (Triangulation::glue
            // writes both sides; the second visit would be a harmless
            // identical overwrite, but skip it for symmetry.)
            if (adj->index() < p) continue;
            if (adj->index() == p && perm[f] < static_cast<int>(f)) continue;
            Perm5 p5{
                static_cast<std::uint8_t>(perm[0]),
                static_cast<std::uint8_t>(perm[1]),
                static_cast<std::uint8_t>(perm[2]),
                static_cast<std::uint8_t>(perm[3]),
                static_cast<std::uint8_t>(perm[4]),
            };
            T.glue(p, f, adj->index(), p5);
        }
    }
    return T;
}

std::string reginaIsoSig(const Triangulation& T) {
    auto R = toRegina(T);
    if (!R)
        throw std::runtime_error(
            "reginaIsoSig: facet glued to itself is not representable in Regina");
    return R->isoSig();
}

// ---------------------------------------------------------------------- //
// Move executor
// ---------------------------------------------------------------------- //

std::optional<Triangulation> applyMoveRegina(const Triangulation& T,
                                             const MoveCandidate& c) {
  try {
    auto R = toRegina(T);                       // fresh copy; safe on failure
    if (!R) return std::nullopt;                // facet-to-itself gluing
    const long idx = c.locator.empty() ? 0 : static_cast<long>(c.locator[0]);
    bool ok = false;
    switch (c.move) {
        case MoveType::Pachner_2_3: {            // 2-4 move, about a tetrahedron
            const long n = static_cast<long>(R->countTetrahedra());
            if (n == 0) return std::nullopt;
            ok = R->pachner(R->tetrahedron(idx % n));
            break;
        }
        case MoveType::Pachner_3_3: {            // 3-3 move, about a triangle
            const long n = static_cast<long>(R->countTriangles());
            if (n == 0) return std::nullopt;
            ok = R->pachner(R->triangle(idx % n));
            break;
        }
        case MoveType::Pachner_4_4: {            // 4-4 move, about an edge
            const long n = static_cast<long>(R->countEdges());
            if (n == 0) return std::nullopt;
            ok = R->move44(R->edge(idx % n));
            break;
        }
        case MoveType::Pachner_4_2: {            // 4-2 move, about an edge
            const long n = static_cast<long>(R->countEdges());
            if (n == 0) return std::nullopt;
            ok = R->pachner(R->edge(idx % n));
            break;
        }
        case MoveType::EdgeCollapse: {           // collapse an edge entirely
            const long n = static_cast<long>(R->countEdges());
            if (n == 0) return std::nullopt;
            ok = R->collapseEdge(R->edge(idx % n));
            break;
        }
        default:
            return std::nullopt;                 // 1-5 / 5-1 handled natively
    }
    if (!ok) return std::nullopt;

    Triangulation U = fromRegina(*R);
    U.setLabel(T.label());

    // Self-verify (mirrors the 5-1 executor's philosophy): every legal move
    // above is a PL homeomorphism, so the result must be a valid
    // triangulation with the same closedness and Euler characteristic.
    // Anything else means a bridge bug -> skip the move, never commit it.
    if (!U.isValid()) return std::nullopt;
    if (U.isClosed() != T.isClosed()) return std::nullopt;
    if (U.eulerCharacteristic() != T.eulerCharacteristic()) return std::nullopt;
    return U;
  } catch (...) {
    // Regina signals precondition violations (e.g. join on an already-glued
    // facet, bad isosig) with exceptions; treat all of them as "move not
    // available" rather than letting them unwind through the beam search.
    return std::nullopt;
  }
}

}  // namespace bhrt

#endif  // BHRT_HAS_REGINA
