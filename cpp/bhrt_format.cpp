// bhrt_format.cpp -- .bhrt text-format triangulation reader/writer.
//
// Format (line-based, ASCII):
//
//   # bhrt-format v1
//   # label: <optional label>
//   N <pentachora_count>
//   G <src_pent> <src_facet> <dst_pent> <p0> <p1> <p2> <p3> <p4>
//   G ...
//   END
//
// Where each G line defines a face gluing: facet src_facet of pent
// src_pent is glued to pent dst_pent under the 5-permutation (p0..p4).
// Only one direction of each gluing needs to appear; the loader writes
// the mate automatically. Lines beginning with `#` are comments.
//
// A whole .bhrt file may contain multiple consecutive blocks (each
// ending with END), which makes it the C++ equivalent of an .esig list.
//
// This file format is bridged from Regina's .esig format by the small
// Python helper at python/bhrt_trisect/regina_bridge.py.

#include "bhrt_trisect.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace bhrt {

namespace {

bool startsWith(const std::string& s, const std::string& pfx) {
    return s.size() >= pfx.size()
        && s.compare(0, pfx.size(), pfx) == 0;
}

}  // namespace

// ---------------------------------------------------------------------- //
// Reader
// ---------------------------------------------------------------------- //

std::vector<Triangulation> loadBhrtFile(const std::string& path) {
    std::ifstream fh(path);
    if (!fh) throw std::runtime_error("cannot open " + path);
    std::vector<Triangulation> out;

    Triangulation current;
    bool in_block = false;
    int pending_pents = -1;
    std::string line;
    int line_no = 0;

    auto flush = [&] {
        if (in_block) {
            if (!current.isValid())
                throw std::runtime_error(
                    "block at line " + std::to_string(line_no)
                    + " is not a valid triangulation");
            out.push_back(std::move(current));
            current = Triangulation();
            in_block = false;
            pending_pents = -1;
        }
    };

    while (std::getline(fh, line)) {
        ++line_no;
        // Strip trailing whitespace and CR (for Windows line endings).
        while (!line.empty() && (line.back() == '\r' || line.back() == ' '
                                  || line.back() == '\t')) line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '#') {
            if (startsWith(line, "# label:") && !in_block) {
                current.setLabel(line.substr(8));
                while (!current.label().empty()
                       && current.label().front() == ' ')
                    current.setLabel(current.label().substr(1));
            }
            continue;
        }

        std::istringstream iss(line);
        std::string tag;
        iss >> tag;
        if (tag == "N") {
            if (in_block)
                throw std::runtime_error("nested N at line "
                                          + std::to_string(line_no));
            int n;
            if (!(iss >> n) || n < 0)
                throw std::runtime_error("bad N at line "
                                          + std::to_string(line_no));
            current.addPentachora(n);
            pending_pents = n;
            in_block = true;
        } else if (tag == "G") {
            if (!in_block)
                throw std::runtime_error("G before N at line "
                                          + std::to_string(line_no));
            int src, src_facet, dst;
            int p0, p1, p2, p3, p4;
            if (!(iss >> src >> src_facet >> dst >> p0 >> p1 >> p2 >> p3 >> p4))
                throw std::runtime_error("malformed G at line "
                                          + std::to_string(line_no));
            if (src < 0 || dst < 0 || src_facet < 0 || src_facet > 4
                || src >= pending_pents || dst >= pending_pents)
                throw std::runtime_error("bad pent or facet at line "
                                          + std::to_string(line_no));
            Perm5 perm{
                static_cast<std::uint8_t>(p0),
                static_cast<std::uint8_t>(p1),
                static_cast<std::uint8_t>(p2),
                static_cast<std::uint8_t>(p3),
                static_cast<std::uint8_t>(p4),
            };
            if (perm[src_facet] != static_cast<std::uint8_t>(
                    (int)((uint8_t)perm[src_facet])))
                ; // suppress warning
            // Only apply the gluing if not already set (handles
            // double-listed mate lines gracefully).
            const auto& g = current.pentachoron(src).gluings[src_facet];
            if (!g.has_value())
                current.glue(src, static_cast<std::uint8_t>(src_facet),
                              dst, perm);
        } else if (tag == "END") {
            flush();
        } else {
            throw std::runtime_error("unknown tag '" + tag + "' at line "
                                      + std::to_string(line_no));
        }
    }
    flush();
    return out;
}

// ---------------------------------------------------------------------- //
// Writer
// ---------------------------------------------------------------------- //

void writeBhrtFile(const std::string& path,
                   const std::vector<Triangulation>& tris) {
    std::ofstream fh(path);
    if (!fh) throw std::runtime_error("cannot write " + path);
    fh << "# bhrt-format v1\n";
    for (const auto& T : tris) {
        if (!T.label().empty()) fh << "# label: " << T.label() << "\n";
        fh << "N " << T.size() << "\n";
        for (std::size_t p = 0; p < T.size(); ++p) {
            const auto& pent = T.pentachoron(p);
            for (std::uint8_t f = 0; f < 5; ++f) {
                const auto& g = pent.gluings[f];
                if (!g.has_value()) continue;
                if (g->dst_pent < static_cast<std::int32_t>(p)) continue;
                fh << "G " << p << " " << (int)f << " " << g->dst_pent
                   << " " << (int)g->perm[0] << " " << (int)g->perm[1]
                   << " " << (int)g->perm[2] << " " << (int)g->perm[3]
                   << " " << (int)g->perm[4] << "\n";
            }
        }
        fh << "END\n";
    }
}

}  // namespace bhrt
