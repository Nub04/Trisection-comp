/* bhrt_c_api.h -- C-language wrapper around the C++ core.
 *
 * Exposes a stable C ABI for the BHRT pipeline so that it can be linked
 * from a C-language CLI, from Python via ctypes/cffi, or from any other
 * language with a C FFI.
 */

#ifndef BHRT_C_API_H
#define BHRT_C_API_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdint.h>

typedef struct bhrt_triangulation bhrt_triangulation;

/* Lifecycle */
bhrt_triangulation* bhrt_new_triangulation(void);
void                bhrt_free_triangulation(bhrt_triangulation* t);
bhrt_triangulation* bhrt_s4_minimal(void);

/* Construction */
void bhrt_add_pentachora(bhrt_triangulation* t, size_t count);
void bhrt_glue(bhrt_triangulation* t,
               int32_t src, uint8_t src_facet, int32_t dst,
               const uint8_t perm5[5]);

/* Skeleton */
size_t  bhrt_size(const bhrt_triangulation* t);
int32_t bhrt_n_vertices(const bhrt_triangulation* t);
int32_t bhrt_n_edges(const bhrt_triangulation* t);
int32_t bhrt_n_triangles(const bhrt_triangulation* t);
int32_t bhrt_n_tetrahedra(const bhrt_triangulation* t);
int32_t bhrt_euler(const bhrt_triangulation* t);
int     bhrt_is_closed(const bhrt_triangulation* t);
int     bhrt_is_valid(const bhrt_triangulation* t);

/* Canonical signature. Writes up to `cap` bytes to `out` and returns the
 * total length needed (excluding the trailing NUL). */
size_t bhrt_isosig(const bhrt_triangulation* t, char* out, size_t cap);
size_t bhrt_edge_degree_isosig(const bhrt_triangulation* t, char* out, size_t cap);

/* ts-tricolouring scan: writes the six counters into `counters[6]`:
 *   counters[0] = n_total
 *   counters[1] = n_pass_precheck
 *   counters[2] = n_with_tricolouring
 *   counters[3] = n_with_c_tricolouring
 *   counters[4] = n_with_ts_tricolouring
 *   counters[5] = n_ts_tricolourings_total
 */
void bhrt_scan_single(const bhrt_triangulation* t,
                       int collapse_budget,
                       int64_t counters[6]);

/* Beam search. Returns the best (lowest pentachora, lowest genus) state's
 * pentachora count; the full pareto front is exposed via the C++ layer. */
int bhrt_beam_search_best_pentachora(const bhrt_triangulation* t,
                                      int beam_width,
                                      int excess_height,
                                      int max_iterations,
                                      uint64_t seed);

/* .bhrt file I/O. The scan-file entry point reads every triangulation
 * from a .bhrt file, runs the ts-tricolouring counter, and writes the
 * aggregate result to `counters[6]` (same layout as bhrt_scan_single).
 * Returns the number of triangulations read, or -1 on error. */
int bhrt_scan_file(const char* path, int collapse_budget,
                   int64_t counters[6]);

/* Dump the canonical edge-degree isosig of every triangulation in a
 * .bhrt file to stdout, one per line. Returns 0 on success, -1 on error. */
int bhrt_dump_isosigs(const char* path);

/* Write a .bhrt file containing the built-in S^4 example (for testing
 * round-trip). Returns 0 on success, -1 on error. */
int bhrt_write_demo(const char* path);

#ifdef __cplusplus
}
#endif

#endif /* BHRT_C_API_H */
