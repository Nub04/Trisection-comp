// regina_io.hpp -- internal shared declarations for the Regina backend.
//
// Include ONLY from translation units compiled with BHRT_HAS_REGINA
// (regina_bridge.cpp, tri_io.cpp, and the Regina-gated tests). The
// definitions live in regina_bridge.cpp so there is exactly one copy of
// the bhrt::Triangulation <-> regina::Triangulation<4> conversion code.

#pragma once

#include "bhrt_trisect.hpp"

#if BHRT_HAS_REGINA

#include <regina/triangulation/dim4.h>

#include <memory>
#include <string>

namespace bhrt {

// Convert a native triangulation to Regina. Pentachoron indices and facet
// permutations are preserved verbatim. Self-gluings (two DISTINCT facets
// of the same pentachoron) are supported; a facet glued to ITSELF cannot
// be represented in Regina, in which case nullptr is returned.
// Precondition: T.isValid().
std::unique_ptr<regina::Triangulation<4>> toRegina(const Triangulation& T);

// Convert back. Index- and permutation-preserving inverse of toRegina:
// fromRegina(*toRegina(T)) has a gluing table identical to T's.
Triangulation fromRegina(const regina::Triangulation<4>& R);

// Regina's canonical isomorphism signature of T (independent cross-check
// of the native isoSig). Throws std::runtime_error if T has a facet glued
// to itself.
std::string reginaIsoSig(const Triangulation& T);

}  // namespace bhrt

#endif  // BHRT_HAS_REGINA
