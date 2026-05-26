// bhrt_bench.cpp -- standalone C++ benchmark driver.
//
// Replaces the Python scripts/reproduce_paper.py for the C++-only build
// path. Iterates over a .bhrt file, runs the ts-tricolouring counter,
// (optionally) extracts diagrams and computes invariants for each
// supported triangulation, and runs the beam-search per record.
//
// Output is one JSONL row per triangulation, written to stdout.
//
// Usage:
//   bhrt-bench <input.bhrt> [--invariants] [--search]
//              [--beam B] [--height H] [--iters I] [--seed S]

#include "bhrt_trisect.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {

struct Args {
    std::string path;
    bool invariants = false;
    bool search = false;
    int  beam   = 64;
    int  height = 4;
    int  iters  = 50;
    std::uint64_t seed = 0;
};

Args parse_args(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr,
            "usage: %s <input.bhrt> [--invariants] [--search]\n"
            "       [--beam B] [--height H] [--iters I] [--seed S]\n",
            argv[0]);
        std::exit(2);
    }
    Args a;
    a.path = argv[1];
    for (int i = 2; i < argc; ++i) {
        std::string opt = argv[i];
        if (opt == "--invariants") a.invariants = true;
        else if (opt == "--search") a.search = true;
        else if (opt == "--beam"   && i + 1 < argc) a.beam   = std::atoi(argv[++i]);
        else if (opt == "--height" && i + 1 < argc) a.height = std::atoi(argv[++i]);
        else if (opt == "--iters"  && i + 1 < argc) a.iters  = std::atoi(argv[++i]);
        else if (opt == "--seed"   && i + 1 < argc) a.seed   = std::strtoull(argv[++i], nullptr, 10);
        else {
            std::fprintf(stderr, "unknown option '%s'\n", opt.c_str());
            std::exit(2);
        }
    }
    return a;
}

std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 2);
    out.push_back('"');
    for (char c : s) {
        if (c == '"' || c == '\\') { out.push_back('\\'); out.push_back(c); }
        else if (c == '\n') out += "\\n";
        else if (c == '\r') out += "\\r";
        else if (c == '\t') out += "\\t";
        else out.push_back(c);
    }
    out.push_back('"');
    return out;
}

void emit_record(int idx, const bhrt::Triangulation& T,
                  const Args& args) {
    auto t0 = std::chrono::steady_clock::now();
    auto ts_list = bhrt::enumerateTsTricolourings(T);
    auto t1 = std::chrono::steady_clock::now();

    std::ostringstream js;
    js << "{";
    js << "\"idx\":" << idx;
    js << ",\"label\":" << json_escape(T.label());
    js << ",\"pent\":" << T.size();
    js << ",\"vertices\":" << T.nVertices();
    js << ",\"edges\":" << T.nEdges();
    js << ",\"euler\":" << T.eulerCharacteristic();
    js << ",\"isosig\":" << json_escape(T.edgeDegreeIsoSig());
    js << ",\"ts_count\":" << ts_list.size();
    js << ",\"ts_seconds\":"
       << std::chrono::duration<double>(t1 - t0).count();

    if (args.invariants && !ts_list.empty()) {
        auto diag = bhrt::extractDiagram(T, ts_list.front());
        auto inv = bhrt::computeInvariants(diag);
        js << ",\"genus\":" << diag.genus();
        js << ",\"h1_free\":" << inv.h1_free_rank;
        js << ",\"h2_free\":" << inv.h2_free_rank;
        js << ",\"h3_free\":" << inv.h3_free_rank;
        js << ",\"signature\":" << inv.signature;
        js << ",\"parity\":" << json_escape(inv.parity);
    }
    if (args.search) {
        bhrt::SearchConfig cfg;
        cfg.beam_width = args.beam;
        cfg.excess_height = args.height;
        cfg.max_iterations = args.iters;
        cfg.seed = args.seed;
        auto pareto = bhrt::beamSearch(T, cfg);
        int best_pent = static_cast<int>(T.size());
        for (auto& s : pareto.items) best_pent = std::min(best_pent, s.pentachora);
        js << ",\"search_best_pent\":" << best_pent
           << ",\"pareto_size\":" << pareto.items.size();
    }
    js << "}";
    std::cout << js.str() << "\n";
}

}  // namespace

int main(int argc, char** argv) {
    Args args = parse_args(argc, argv);
    std::vector<bhrt::Triangulation> tris;
    try {
        tris = bhrt::loadBhrtFile(args.path);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "load error: %s\n", e.what());
        return 1;
    }
    for (std::size_t i = 0; i < tris.size(); ++i)
        emit_record(static_cast<int>(i), tris[i], args);
    return 0;
}
