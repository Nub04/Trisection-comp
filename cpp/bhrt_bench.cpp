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

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

struct Args {
    std::string path;
    bool invariants = false;
    bool search = false;
    bool gpu    = false;   // route the beam-search scorer through CUDA
    int  beam   = 64;
    int  height = 4;
    int  iters  = 50;
    std::uint64_t seed = 0;
    bool scan_bench = false;   // GPU-vs-CPU exhaustive ts-colouring scan
    int  scan_vertices = 14;   // grow S^4 to this many vertices (3^(n-1) space)

    // Census mode (Pachner-graph trisection-genus census, Burton 1110.6080).
    bool        census          = false; // run a census instead of the per-record loop
    std::string census_path;             // .bhrt (native) or .esig (Regina) node list
    bool        census_generate = false; // self-generate a Pachner level set from S^4
    int         census_max_pent = 6;     // pentachora cap for --census-generate
    int         census_max_nodes = 5000; // hard node cap (dedup'd) for generation safety
    std::string known_path;              // optional "<label-or-isosig> <expected_genus>" file
};

Args parse_args(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr,
            "usage: %s <input.bhrt> [--invariants] [--search] [--gpu]\n"
            "       [--beam B] [--height H] [--iters I] [--seed S]\n"
            "  --gpu             route the beam-search scorer through CUDA\n"
            "  --scan-bench [N]  GPU-vs-CPU exhaustive ts-colouring scan; grows the\n"
            "                    S^4 seed to N vertices (3^(N-1) space), default 14\n"
            "  --census <file>   trisection-genus census over a node list (.bhrt or .esig):\n"
            "                    dedup by isoSig, GPU/CPU ts-scan each, report min genus\n"
            "  --census-generate [N]  self-generate a Pachner level set from S^4 up to N\n"
            "                    pentachora (default 6) and census it (machinery check)\n"
            "  --max-nodes M     hard cap on generated nodes (default 5000)\n"
            "  --known <file>    validate min genus against \"<label-or-isosig> <g>\" lines\n",
            argv[0]);
        std::exit(2);
    }
    Args a;
    int i = 1;
    if (argv[1][0] != '-') { a.path = argv[1]; i = 2; }  // first arg = file unless an option
    for (; i < argc; ++i) {
        std::string opt = argv[i];
        if (opt == "--invariants") a.invariants = true;
        else if (opt == "--search") a.search = true;
        else if (opt == "--gpu")    a.gpu = true;
        else if (opt == "--scan-bench") {
            a.scan_bench = true;
            if (i + 1 < argc && std::isdigit((unsigned char)argv[i + 1][0]))
                a.scan_vertices = std::atoi(argv[++i]);
        }
        else if (opt == "--census" && i + 1 < argc) {
            a.census = true;
            a.census_path = argv[++i];
        }
        else if (opt == "--census-generate") {
            a.census = true;
            a.census_generate = true;
            if (i + 1 < argc && std::isdigit((unsigned char)argv[i + 1][0]))
                a.census_max_pent = std::atoi(argv[++i]);
        }
        else if (opt == "--max-nodes" && i + 1 < argc) a.census_max_nodes = std::atoi(argv[++i]);
        else if (opt == "--known"     && i + 1 < argc) a.known_path = argv[++i];
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
        js << ",\"form_computed\":"
           << (inv.intersection_form.empty() ? "false" : "true");
        if (!inv.intersection_form.empty()) {
            js << ",\"form\":[";
            for (std::size_t r = 0; r < inv.intersection_form.size(); ++r) {
                if (r) js << ",";
                js << "[";
                for (std::size_t c = 0; c < inv.intersection_form[r].size(); ++c) {
                    if (c) js << ",";
                    js << inv.intersection_form[r][c];
                }
                js << "]";
            }
            js << "]";
        }
        // The audit carries the dev signal for cut-curve extraction: whether the
        // curves are null-homologous (degenerate extractDiagram) and how/why the
        // intersection form was computed or deferred.
        js << ",\"audit\":[";
        for (std::size_t k = 0; k < inv.audit.size(); ++k) {
            if (k) js << ",";
            js << json_escape(inv.audit[k]);
        }
        js << "]";
    }
    if (args.search) {
        bhrt::SearchConfig cfg;
        cfg.beam_width = args.beam;
        cfg.excess_height = args.height;
        cfg.max_iterations = args.iters;
        cfg.seed = args.seed;
        cfg.use_gpu_scorer = args.gpu;
        auto s0 = std::chrono::steady_clock::now();
        auto pareto = bhrt::beamSearch(T, cfg);
        auto s1 = std::chrono::steady_clock::now();
        int best_pent = static_cast<int>(T.size());
        for (auto& s : pareto.items) best_pent = std::min(best_pent, s.pentachora);
        js << ",\"search_best_pent\":" << best_pent
           << ",\"pareto_size\":" << pareto.items.size()
           << ",\"search_gpu\":" << (args.gpu ? "true" : "false")
           << ",\"search_seconds\":"
           << std::chrono::duration<double>(s1 - s0).count();
    }
    js << "}";
    // Flush each record immediately: on CUDA builds the process may be torn
    // down abnormally at exit (see main()), and we must not lose output that
    // is already computed.
    std::cout << js.str() << std::endl;
}

// GPU-vs-CPU exhaustive ts-tricolouring scan. Grows the minimal S^4 to N
// vertices via 1-5 Pachner moves, then races the GPU on-device precheck over
// the whole 3^(N-1) colouring space against the CPU enumerator. This measures
// the parallelism that actually matters (the colouring scan), not the tiny
// beam-search scorer.
void run_scan_bench(const Args& args) {
    int target = args.scan_vertices;
    if (target < 5) target = 5;
    if (target > 31) {
        std::fprintf(stderr,
            "[scan-bench] capping at 31 vertices (GPU scan supports <= 32)\n");
        target = 31;
    }
    bhrt::Triangulation T = bhrt::s4_minimal();
    while (T.nVertices() < target) {
        bool grew = false;
        for (auto& m : bhrt::enumerateMoves(T)) {
            if (m.move == bhrt::MoveType::Pachner_1_5) {
                auto next = bhrt::applyMove(T, m);
                if (next) { T = std::move(*next); grew = true; break; }
            }
        }
        if (!grew) break;
    }
    const int n = T.nVertices();
    const double total = std::pow(3.0, n - 1);   // colouring space (vertex 0 pinned)
    const bool run_cpu = total <= 5.0e7;         // CPU enumerate is ~O(3^(n-1))
    std::fprintf(stderr,
        "[scan-bench] %zu pentachora, %d vertices -> 3^%d = %.4g colourings%s\n",
        T.size(), n, n - 1, total,
        run_cpu ? "" : "  (CPU side skipped: too large)");

    double cpu_s = -1.0;
    std::size_t cpu_found = 0;
    if (run_cpu) {
        auto c0 = std::chrono::steady_clock::now();
        auto cpu = bhrt::enumerateTsTricolourings(T);
        auto c1 = std::chrono::steady_clock::now();
        cpu_s = std::chrono::duration<double>(c1 - c0).count();
        cpu_found = cpu.size();
    }

    std::ostringstream js;
    js << "{\"pentachora\":" << T.size()
       << ",\"vertices\":" << n
       << ",\"colourings\":" << static_cast<unsigned long long>(total);
    if (run_cpu)
        js << ",\"cpu_seconds\":" << cpu_s
           << ",\"cpu_colourings_per_sec\":" << (cpu_s > 0 ? total / cpu_s : 0.0)
           << ",\"cpu_ts_found\":" << cpu_found;
    else
        js << ",\"cpu_seconds\":null";
#if BHRT_HAS_CUDA
    auto g0 = std::chrono::steady_clock::now();
    auto gpu = bhrt::enumerateTsTricolouringsGPU(T);
    auto g1 = std::chrono::steady_clock::now();
    const double gpu_s = std::chrono::duration<double>(g1 - g0).count();
    js << ",\"gpu_seconds\":" << gpu_s
       << ",\"gpu_colourings_per_sec\":" << (gpu_s > 0 ? total / gpu_s : 0.0)
       << ",\"gpu_ts_found\":" << gpu.size();
    if (run_cpu)
        js << ",\"speedup\":" << (gpu_s > 0 ? cpu_s / gpu_s : 0.0)
           << ",\"results_match\":" << (gpu.size() == cpu_found ? "true" : "false");
#else
    js << ",\"gpu\":\"not built (needs a CUDA build)\"";
#endif
    js << "}";
    std::cout << js.str() << std::endl;
}

// ---------------------------------------------------------------------- //
// Census mode: Pachner-graph census of trisection genus.
//
// Methodology follows Burton, "Simplification paths in the Pachner graphs of
// closed orientable 3-manifold triangulations" (arXiv:1110.6080), lifted to
// 4-manifold trisections:
//   * nodes are 4-manifold triangulations, deduplicated by the COMPLETE
//     isomorphism signature isoSig() (Burton's sigma, a complete invariant;
//     edgeDegreeIsoSig is only an edge-degree invariant and is NOT used here);
//   * --census-generate walks the Pachner graph outward from S^4 via the
//     native 1-5/5-1 moves (the full 2-4/3-3/4-4/4-2 catalogue is available in
//     Regina builds), forming a finite level set up to a pentachora cap;
//   * each node is GPU- (or CPU-) scanned for every ts-tricolouring, and the
//     minimum central-surface genus g(Sigma) is the per-triangulation
//     trisection genus.  g(Sigma) uses only the colouring and the central
//     surface -- both of which extractDiagram builds correctly -- and does NOT
//     touch the (still-degenerate) cut curves.
//
// Self-generated S^4 nodes are all S^4, so generation validates the machinery
// (dedup, GPU==CPU agreement, genus tabulation), not interesting genera.  Real
// genera come from a census of distinct manifolds fed via --census <file>.
// ---------------------------------------------------------------------- //

std::vector<bhrt::TSColouring> enumerate_ts(const bhrt::Triangulation& T, bool gpu) {
#if BHRT_HAS_CUDA
    if (gpu) return bhrt::enumerateTsTricolouringsGPU(T);
#else
    (void)gpu;
#endif
    return bhrt::enumerateTsTricolourings(T);
}

bool ends_with(const std::string& s, const char* suffix) {
    const std::size_t n = std::strlen(suffix);
    return s.size() >= n && s.compare(s.size() - n, n, suffix) == 0;
}

std::vector<bhrt::Triangulation> load_census_file(const std::string& path) {
    if (ends_with(path, ".esig")) {
#if BHRT_HAS_REGINA
        return bhrt::loadEsig(path);
#else
        throw std::runtime_error(
            ".esig input requires a Regina build (BHRT_HAS_REGINA); "
            "convert it to .bhrt or rebuild with Regina");
#endif
    }
    return bhrt::loadBhrtFile(path);
}

// Breadth-first walk of the Pachner graph from the minimal S^4, deduplicated
// by the complete isoSig, bounded by a pentachora cap and a hard node cap.
std::vector<bhrt::Triangulation> generate_level_set(int max_pent, int max_nodes) {
    std::vector<bhrt::Triangulation> out;
    std::unordered_set<std::string>  seen;
    std::deque<bhrt::Triangulation>  queue;

    bhrt::Triangulation seed = bhrt::s4_minimal();
    seen.insert(seed.isoSig());
    queue.push_back(std::move(seed));

    while (!queue.empty() && static_cast<int>(out.size()) < max_nodes) {
        bhrt::Triangulation T = std::move(queue.front());
        queue.pop_front();
        out.push_back(T);  // keep a copy; T stays valid for move enumeration

        for (const auto& m : bhrt::enumerateMoves(T)) {
            auto next = bhrt::applyMove(T, m);
            if (!next) continue;                                  // illegal / no native executor
            if (static_cast<int>(next->size()) > max_pent) continue;
            std::string sig = next->isoSig();
            if (!seen.insert(sig).second) continue;               // already seen
            queue.push_back(std::move(*next));
        }
    }
    return out;
}

struct NodeResult {
    std::string   label;
    std::string   isosig;
    int           pent      = 0;
    int           vertices  = 0;
    std::size_t   ts_count  = 0;
    bool          admits    = false;
    int           min_genus = -1;     // valid iff admits
    std::set<int> genera;
    double        seconds   = 0.0;
};

NodeResult census_node(const bhrt::Triangulation& T,
                       const std::string& sig, bool gpu) {
    NodeResult r;
    r.label    = T.label();
    r.isosig   = sig;
    r.pent     = static_cast<int>(T.size());
    r.vertices = T.nVertices();

    auto t0 = std::chrono::steady_clock::now();
    std::vector<bhrt::TSColouring> ts = enumerate_ts(T, gpu);
    r.ts_count = ts.size();
    r.admits   = !ts.empty();
    for (const auto& c : ts) {
        bhrt::TrisectionDiagram diag = bhrt::extractDiagram(T, c);
        int g = diag.genus();
        r.genera.insert(g);
        if (r.min_genus < 0 || g < r.min_genus) r.min_genus = g;
    }
    auto t1 = std::chrono::steady_clock::now();
    r.seconds = std::chrono::duration<double>(t1 - t0).count();
    return r;
}

std::map<std::string, int> load_known(const std::string& path) {
    std::map<std::string, int> known;
    if (path.empty()) return known;
    std::ifstream in(path);
    if (!in) {
        std::fprintf(stderr, "[census] warning: cannot open --known file '%s'\n",
                     path.c_str());
        return known;
    }
    std::string line;
    while (std::getline(in, line)) {
        std::size_t h = line.find('#');
        if (h != std::string::npos) line.erase(h);
        std::istringstream ls(line);
        std::string key;
        int g;
        if (ls >> key >> g) known[key] = g;
    }
    return known;
}

void run_census(const Args& args) {
    std::vector<bhrt::Triangulation> nodes;
    try {
        if (args.census_generate)
            nodes = generate_level_set(args.census_max_pent, args.census_max_nodes);
        else
            nodes = load_census_file(args.census_path);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "[census] load error: %s\n", e.what());
        std::exit(1);
    }

    std::map<std::string, int> known = load_known(args.known_path);

    std::fprintf(stderr,
        "[census] %zu input node(s); mode=%s, scan=%s%s\n",
        nodes.size(),
        args.census_generate ? "generate" : "file",
        args.gpu ? "gpu" : "cpu",
        known.empty() ? "" : ", validating against --known");

    std::unordered_set<std::string> seen;
    std::map<int, int> genus_hist;      // min_genus -> #nodes
    int n_unique = 0, n_admit = 0, n_dup = 0;
    int n_checked = 0, n_pass = 0, n_fail = 0;
    int overall_min = -1;

    for (const bhrt::Triangulation& T : nodes) {
        std::string sig = T.isoSig();
        if (!seen.insert(sig).second) { ++n_dup; continue; }
        ++n_unique;

        NodeResult r = census_node(T, sig, args.gpu);
        if (r.admits) {
            ++n_admit;
            ++genus_hist[r.min_genus];
            if (overall_min < 0 || r.min_genus < overall_min) overall_min = r.min_genus;
        }

        bool have_expect = false;
        int  expect = -1;
        bool match = false;
        if (!known.empty()) {
            auto it = known.find(r.label);
            if (it == known.end()) it = known.find(r.isosig);
            if (it != known.end()) {
                have_expect = true;
                expect = it->second;
                match = r.admits && (r.min_genus == expect);
                ++n_checked;
                if (match) ++n_pass; else ++n_fail;
            }
        }

        std::ostringstream js;
        js << "{\"label\":" << json_escape(r.label)
           << ",\"isosig\":" << json_escape(r.isosig)
           << ",\"pent\":" << r.pent
           << ",\"vertices\":" << r.vertices
           << ",\"ts_count\":" << r.ts_count
           << ",\"admits_trisection\":" << (r.admits ? "true" : "false");
        if (r.admits) {
            js << ",\"min_genus\":" << r.min_genus << ",\"genera\":[";
            bool first = true;
            for (int g : r.genera) { if (!first) js << ","; js << g; first = false; }
            js << "]";
        } else {
            js << ",\"min_genus\":null";
        }
        js << ",\"scan\":" << (args.gpu ? "\"gpu\"" : "\"cpu\"")
           << ",\"seconds\":" << r.seconds;
        if (have_expect)
            js << ",\"expected_genus\":" << expect
               << ",\"match\":" << (match ? "true" : "false");
        js << "}";
        std::cout << js.str() << std::endl;
    }

    std::ostringstream sm;
    sm << "{\"summary\":true"
       << ",\"nodes_input\":" << nodes.size()
       << ",\"nodes_unique\":" << n_unique
       << ",\"nodes_duplicate\":" << n_dup
       << ",\"nodes_admitting_trisection\":" << n_admit
       << ",\"scan\":" << (args.gpu ? "\"gpu\"" : "\"cpu\"");
    if (overall_min >= 0) sm << ",\"min_genus_overall\":" << overall_min;
    else                  sm << ",\"min_genus_overall\":null";
    sm << ",\"genus_histogram\":{";
    bool first = true;
    for (const auto& kv : genus_hist) {
        if (!first) sm << ",";
        sm << "\"" << kv.first << "\":" << kv.second;
        first = false;
    }
    sm << "}";
    if (!known.empty())
        sm << ",\"validated\":" << n_checked
           << ",\"validation_pass\":" << n_pass
           << ",\"validation_fail\":" << n_fail;
    sm << "}";
    std::cout << sm.str() << std::endl;
}

}  // namespace

int main(int argc, char** argv) {
    Args args = parse_args(argc, argv);
    if (args.census) {
        run_census(args);
    } else if (args.scan_bench) {
        run_scan_bench(args);
    } else {
        std::vector<bhrt::Triangulation> tris;
        try {
            tris = bhrt::loadBhrtFile(args.path);
        } catch (const std::exception& e) {
            std::fprintf(stderr, "load error: %s\n", e.what());
            return 1;
        }
        for (std::size_t i = 0; i < tris.size(); ++i)
            emit_record(static_cast<int>(i), tris[i], args);
    }

    std::cout.flush();
#if BHRT_HAS_CUDA
    // On some WSL2 + CUDA driver stacks the CUDA runtime's own process-exit
    // teardown corrupts the glibc heap ("munmap_chunk(): invalid pointer"),
    // which aborts the process AFTER all work is done and would otherwise
    // swallow buffered output. The computation is complete at this point, so
    // terminate immediately and let the OS reclaim resources, bypassing the
    // fragile teardown. (Plain non-CUDA builds run normal destructors.)
    std::fflush(nullptr);
    std::_Exit(0);
#endif
    return 0;
}
