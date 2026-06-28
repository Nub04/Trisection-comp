// tri_io.cpp — Regina-backed .esig loading and writing.
//
// This compilation unit links against libregina and provides the file-level
// I/O on top of the shared conversions in regina_io.hpp (defined once in
// regina_bridge.cpp). It is only compiled when BHRT_HAS_REGINA is set, so
// the #else branches below exist purely for documentation symmetry.
//
// .esig format: one Regina isomorphism signature per line, optionally
// prefixed with "<edge-degrees>|" (the Dim4Census record style); blank
// lines and '#' comments are skipped.

#include "bhrt_trisect.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

#if BHRT_HAS_REGINA
#include "regina_io.hpp"
#endif

namespace bhrt {

// Public API: load .esig file ------------------------------------------

std::vector<Triangulation> loadEsig(const std::string& path) {
    std::vector<Triangulation> out;
#if BHRT_HAS_REGINA
    std::ifstream fh(path);
    if (!fh) throw std::runtime_error("cannot open " + path);
    std::string line;
    while (std::getline(fh, line)) {
        // Tolerate CRLF line endings (.esig files authored on Windows).
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        // Optional "<prefix>|" before the signature. The prefix is taken as a
        // human-readable label (e.g. "CP2|iLvL..."), which lets downstream
        // tools (e.g. bhrt-bench --known) key on a manifold name; legacy files
        // that put an edge-degree string there simply get it as the label.
        auto bar = line.find('|');
        std::string label = (bar == std::string::npos) ? std::string()
                                                        : line.substr(0, bar);
        std::string sig   = (bar == std::string::npos) ? line
                                                        : line.substr(bar + 1);
        // Regina 7.x: fromIsoSig returns a Triangulation<4> by value and
        // throws on an invalid signature.
        try {
            auto R = regina::Triangulation<4>::fromIsoSig(sig);
            Triangulation T = fromRegina(R);
            if (!label.empty()) T.setLabel(label);
            out.push_back(std::move(T));
        } catch (const std::exception&) {
            throw std::runtime_error("bad isosig: " + sig);
        }
    }
#else
    (void)path;
    throw std::runtime_error(
        "Regina was not compiled in; rebuild with BHRT_HAS_REGINA=1");
#endif
    return out;
}

// Public API: load one triangulation from a single isosig ---------------

Triangulation fromReginaIsoSig(const std::string& sig) {
#if BHRT_HAS_REGINA
    try {
        auto R = regina::Triangulation<4>::fromIsoSig(sig);
        return fromRegina(R);
    } catch (const std::exception&) {
        throw std::runtime_error("bad isosig: " + sig);
    }
#else
    (void)sig;
    throw std::runtime_error(
        "Regina was not compiled in; rebuild with BHRT_HAS_REGINA=1");
#endif
}

// Public API: write .esig file ------------------------------------------

void writeEsig(const std::string& path,
               const std::vector<Triangulation>& tris) {
#if BHRT_HAS_REGINA
    std::ofstream fh(path);
    if (!fh) throw std::runtime_error("cannot open " + path + " for writing");
    fh << "# .esig written by bhrt_trisect (one Regina isosig per line)\n";
    for (const auto& T : tris) {
        auto R = toRegina(T);
        if (!R)
            throw std::runtime_error(
                "writeEsig: triangulation '" + T.label() +
                "' has a facet glued to itself; not representable in Regina");
        fh << R->isoSig() << '\n';
    }
    if (!fh) throw std::runtime_error("write failed: " + path);
#else
    (void)path; (void)tris;
    throw std::runtime_error(
        "Regina was not compiled in; rebuild with BHRT_HAS_REGINA=1");
#endif
}

}  // namespace bhrt
