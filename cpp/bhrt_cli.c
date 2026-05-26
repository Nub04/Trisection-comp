/* bhrt_cli.c -- C-language command-line driver for the BHRT pipeline.
 *
 * Subcommands:
 *   info                 print skeleton statistics of a built-in example
 *   isosig               print canonical isomorphism signature
 *   scan                 ts-tricolouring counters on the built-in example
 *   search [B H I S]     beam search (beam, height, iters, seed)
 *   scan-file <path>     ts-tricolouring counters on every triangulation
 *                        in a .bhrt file
 *   dump-isosigs <path>  print isosig for every triangulation in a .bhrt file
 *   write-demo <path>    write a small demo .bhrt file
 *   reduce-file <path> [B H I S]  beam-search every triangulation in a
 *                        .bhrt file, print best-pent per record
 *   help                 print this message
 */

#include "bhrt_c_api.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void usage(const char* prog) {
    fprintf(stderr,
        "usage: %s <command> [options]\n\n"
        "commands:\n"
        "  info                       built-in S^4 skeleton statistics\n"
        "  isosig                     canonical isomorphism signature of S^4\n"
        "  scan                       ts-tricolouring counters on S^4\n"
        "  search [B H I S]           beam search on S^4\n"
        "  scan-file <path>           scan every triangulation in a .bhrt file\n"
        "  dump-isosigs <path>        emit isosig for every record\n"
        "  write-demo <path>          write a small demo .bhrt file\n"
        "  reduce-file <p> [B H I S]  beam-search every record in a .bhrt file\n"
        "  help                       this message\n",
        prog);
}

static int cmd_info(void) {
    bhrt_triangulation* t = bhrt_s4_minimal();
    if (!t) return 1;
    printf("label:        S4_2p\n");
    printf("pentachora:   %zu\n", bhrt_size(t));
    printf("vertices:     %d\n",  bhrt_n_vertices(t));
    printf("edges:        %d\n",  bhrt_n_edges(t));
    printf("triangles:    %d\n",  bhrt_n_triangles(t));
    printf("tetrahedra:   %d\n",  bhrt_n_tetrahedra(t));
    printf("euler chi:    %d\n",  bhrt_euler(t));
    printf("is_closed:    %d\n",  bhrt_is_closed(t));
    printf("is_valid:     %d\n",  bhrt_is_valid(t));
    bhrt_free_triangulation(t);
    return 0;
}

static int cmd_isosig(void) {
    bhrt_triangulation* t = bhrt_s4_minimal();
    if (!t) return 1;
    char buf[4096];
    bhrt_edge_degree_isosig(t, buf, sizeof(buf));
    printf("%s\n", buf);
    bhrt_free_triangulation(t);
    return 0;
}

static const char* counter_names[6] = {
    "n_total",
    "n_pass_precheck",
    "n_with_tricolouring",
    "n_with_c_tricolouring",
    "n_with_ts_tricolouring",
    "n_ts_tricolourings_total"
};

static void print_counters(const int64_t counters[6]) {
    for (int i = 0; i < 6; ++i)
        printf("%-26s = %lld\n", counter_names[i], (long long)counters[i]);
}

static int cmd_scan(void) {
    bhrt_triangulation* t = bhrt_s4_minimal();
    if (!t) return 1;
    int64_t c[6];
    bhrt_scan_single(t, 10000, c);
    print_counters(c);
    bhrt_free_triangulation(t);
    return 0;
}

static int cmd_search(int argc, char** argv) {
    int beam = argc > 0 ? atoi(argv[0]) : 64;
    int height = argc > 1 ? atoi(argv[1]) : 4;
    int iters = argc > 2 ? atoi(argv[2]) : 50;
    unsigned long seed = argc > 3 ? strtoul(argv[3], NULL, 10) : 0;
    bhrt_triangulation* t = bhrt_s4_minimal();
    if (!t) return 1;
    int best = bhrt_beam_search_best_pentachora(
        t, beam, height, iters, (uint64_t)seed);
    printf("best pentachora found: %d\n", best);
    bhrt_free_triangulation(t);
    return 0;
}

static int cmd_scan_file(int argc, char** argv) {
    if (argc < 1) { fprintf(stderr, "scan-file requires a path\n"); return 2; }
    int64_t c[6];
    int n = bhrt_scan_file(argv[0], 10000, c);
    if (n < 0) return 1;
    printf("# %d triangulations read from %s\n", n, argv[0]);
    print_counters(c);
    return 0;
}

static int cmd_dump_isosigs(int argc, char** argv) {
    if (argc < 1) { fprintf(stderr, "dump-isosigs requires a path\n"); return 2; }
    return bhrt_dump_isosigs(argv[0]) == 0 ? 0 : 1;
}

static int cmd_write_demo(int argc, char** argv) {
    if (argc < 1) { fprintf(stderr, "write-demo requires a path\n"); return 2; }
    if (bhrt_write_demo(argv[0]) != 0) return 1;
    printf("wrote demo to %s\n", argv[0]);
    return 0;
}

/* reduce-file walks the file via dump-isosigs as a quick round-trip
 * proxy; the full per-record beam search uses the C++ core directly
 * (see cpp/bhrt_bench.cpp). */
static int cmd_reduce_file(int argc, char** argv) {
    if (argc < 1) { fprintf(stderr, "reduce-file requires a path\n"); return 2; }
    /* For now the binary returns the per-file scan counters plus the
     * isosigs; the full beam-driven reduction over a corpus is invoked
     * via the bhrt-bench binary (built by the same Makefile rule). */
    int64_t c[6];
    int n = bhrt_scan_file(argv[0], 10000, c);
    if (n < 0) return 1;
    printf("# scanned %d records; beam search per record disabled in the\n"
           "# default build; use bhrt-bench for the full sweep.\n", n);
    print_counters(c);
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(argv[0]); return 2; }
    const char* cmd = argv[1];
    if (strcmp(cmd, "info") == 0)         return cmd_info();
    if (strcmp(cmd, "isosig") == 0)       return cmd_isosig();
    if (strcmp(cmd, "scan") == 0)         return cmd_scan();
    if (strcmp(cmd, "search") == 0)       return cmd_search(argc - 2, argv + 2);
    if (strcmp(cmd, "scan-file") == 0)    return cmd_scan_file(argc - 2, argv + 2);
    if (strcmp(cmd, "dump-isosigs") == 0) return cmd_dump_isosigs(argc - 2, argv + 2);
    if (strcmp(cmd, "write-demo") == 0)   return cmd_write_demo(argc - 2, argv + 2);
    if (strcmp(cmd, "reduce-file") == 0)  return cmd_reduce_file(argc - 2, argv + 2);
    if (strcmp(cmd, "help") == 0)         { usage(argv[0]); return 0; }
    usage(argv[0]);
    return 2;
}
