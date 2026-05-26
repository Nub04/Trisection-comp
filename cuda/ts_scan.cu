// ts_scan.cu -- GPU-exhaustive ts-tricolouring precheck filter.
//
// The exponential bottleneck of the BHRT pipeline is enumerating the
// 3^(n_vertices - 1) candidate vertex tricolourings of a triangulation
// (vertex 0 pinned to colour 0). Each candidate is an INDEPENDENT,
// fixed-size combinatorial check, so the enumeration parallelises
// perfectly: one GPU thread takes one colouring and applies, in order,
//
//   (1) all three colour classes non-empty,
//   (2) every pentachoron is type (2,2,1)            [Spreer-Tillmann],
//   (3) no monochromatic triangle,
//   (4) each monochromatic graph Gamma_k is connected,
//
// which is exactly the cheap prefix of bhrt::isTsTricolouring. Survivor
// indices are appended to a device list with an atomic counter; the host
// then re-runs the FULL certified check (including the gamma 2-complex
// collapsibility step, which needs heap data structures) on the survivors
// only. The GPU is therefore a *filter*, never an oracle: every colouring
// the pipeline reports has been verified by the CPU reference code, and a
// GPU failure of any kind degrades to the CPU path.
//
// Input packing (built by bhrt::enumerateTsTricolouringsGPU):
//   pent_verts : n_pents * 5  global vertex ids (one per pent corner)
//   tri_verts  : n_tris  * 3  global vertex ids (one UNIQUE triangle class
//                             per row, from its first-seen representative)
//   edge_verts : n_edges * 2  global vertex ids (one unique edge class per
//                             row; loop edges a==b are permitted)
// Limits: 3 <= n_vertices <= 32 (colour classes fit a 32-bit mask).

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>

// Weak LeakSanitizer hooks: resolve to null in non-sanitized builds.
extern "C" {
    void __lsan_disable(void) __attribute__((weak));
    void __lsan_enable(void)  __attribute__((weak));
}
namespace {
struct LsanGuard {
    LsanGuard()  { if (__lsan_disable) __lsan_disable(); }
    ~LsanGuard() { if (__lsan_enable)  __lsan_enable();  }
};
}  // namespace

extern "C" {

__global__ void bhrt_ts_precheck(
    const unsigned char* __restrict__ pent_verts, int n_pents,
    const unsigned char* __restrict__ tri_verts,  int n_tris,
    const unsigned char* __restrict__ edge_verts, int n_edges,
    int n_vertices,
    unsigned long long total,
    unsigned long long* __restrict__ out,
    unsigned long long cap,
    unsigned long long* __restrict__ n_out)
{
    const unsigned long long stride =
        static_cast<unsigned long long>(blockDim.x) * gridDim.x;
    for (unsigned long long idx =
             static_cast<unsigned long long>(blockIdx.x) * blockDim.x
             + threadIdx.x;
         idx < total; idx += stride) {

        // ---- decode base-3 colouring (vertex 0 pinned to colour 0) ----
        unsigned char colour[32];
        colour[0] = 0;
        int cnt[3] = {1, 0, 0};
        unsigned long long x = idx;
        for (int v = 1; v < n_vertices; ++v) {
            const unsigned char c = static_cast<unsigned char>(x % 3ULL);
            x /= 3ULL;
            colour[v] = c;
            ++cnt[c];
        }
        if (cnt[1] == 0 || cnt[2] == 0) continue;        // (1)

        // ---- (2) every pentachoron type (2,2,1) ----
        // Three counts summing to 5 with each <= 2 forces {2,2,1}.
        bool ok = true;
        for (int p = 0; p < n_pents && ok; ++p) {
            int pc[3] = {0, 0, 0};
            for (int k = 0; k < 5; ++k)
                ++pc[colour[pent_verts[p * 5 + k]]];
            if (pc[0] > 2 || pc[1] > 2 || pc[2] > 2) ok = false;
        }
        if (!ok) continue;

        // ---- (3) no monochromatic triangle ----
        for (int t = 0; t < n_tris && ok; ++t) {
            const unsigned char a = colour[tri_verts[t * 3 + 0]];
            const unsigned char b = colour[tri_verts[t * 3 + 1]];
            const unsigned char c = colour[tri_verts[t * 3 + 2]];
            if (a == b && b == c) ok = false;
        }
        if (!ok) continue;

        // ---- (4) Gamma_k connected for k = 0,1,2 (bitmask flood fill) ----
        for (int k = 0; k < 3 && ok; ++k) {
            unsigned int verts_k = 0u;
            for (int v = 0; v < n_vertices; ++v)
                if (colour[v] == k) verts_k |= (1u << v);
            unsigned int seen = verts_k & (~verts_k + 1u);   // lowest bit
            for (int round = 0; round < n_vertices; ++round) {
                const unsigned int before = seen;
                for (int e = 0; e < n_edges; ++e) {
                    const unsigned int ea = 1u << edge_verts[e * 2 + 0];
                    const unsigned int eb = 1u << edge_verts[e * 2 + 1];
                    if ((verts_k & ea) && (verts_k & eb)) {
                        if (seen & ea) seen |= eb;
                        if (seen & eb) seen |= ea;
                    }
                }
                if (seen == before) break;
            }
            if (seen != verts_k) ok = false;
        }
        if (!ok) continue;

        // ---- survivor ----
        const unsigned long long slot = atomicAdd(n_out, 1ULL);
        if (slot < cap) out[slot] = idx;
    }
}

// Host launcher. Returns:
//   0  success (survivors + count written)
//   1  no usable GPU / CUDA error          -> caller uses the CPU path
//   2  survivor buffer overflow            -> caller uses the CPU path
//   3  input outside supported limits      -> caller uses the CPU path
int bhrt_ts_scan_host(const unsigned char* pent_verts, int n_pents,
                      const unsigned char* tri_verts,  int n_tris,
                      const unsigned char* edge_verts, int n_edges,
                      int n_vertices,
                      unsigned long long total,
                      unsigned long long* survivors,
                      unsigned long long cap,
                      unsigned long long* n_survivors)
{
    *n_survivors = 0;
    if (n_vertices < 3 || n_vertices > 32) return 3;
    if (n_pents <= 0 || total == 0 || cap == 0) return 3;

    LsanGuard lsan_guard;

    int dev = 0;
    cudaError_t e = cudaGetDeviceCount(&dev);
    if (e != cudaSuccess || dev == 0) return 1;

    const size_t pent_bytes = static_cast<size_t>(n_pents) * 5;
    const size_t tri_bytes  = static_cast<size_t>(n_tris) * 3;
    const size_t edge_bytes = static_cast<size_t>(n_edges) * 2;
    const size_t out_bytes  = static_cast<size_t>(cap) * sizeof(unsigned long long);

    unsigned char* d_pent = nullptr;
    unsigned char* d_tri  = nullptr;
    unsigned char* d_edge = nullptr;
    unsigned long long* d_out = nullptr;
    unsigned long long* d_cnt = nullptr;
    int rc = 1;
    cudaError_t err = cudaSuccess;

    // cudaMalloc(0) is implementation-defined; always allocate >= 1 byte.
    if ((err = cudaMalloc(&d_pent, pent_bytes ? pent_bytes : 1)) == cudaSuccess &&
        (err = cudaMalloc(&d_tri,  tri_bytes  ? tri_bytes  : 1)) == cudaSuccess &&
        (err = cudaMalloc(&d_edge, edge_bytes ? edge_bytes : 1)) == cudaSuccess &&
        (err = cudaMalloc(&d_out,  out_bytes)) == cudaSuccess &&
        (err = cudaMalloc(&d_cnt,  sizeof(unsigned long long))) == cudaSuccess &&
        (err = cudaMemset(d_cnt, 0, sizeof(unsigned long long))) == cudaSuccess &&
        (pent_bytes == 0 || (err = cudaMemcpy(d_pent, pent_verts, pent_bytes,
                                cudaMemcpyHostToDevice)) == cudaSuccess) &&
        (tri_bytes  == 0 || (err = cudaMemcpy(d_tri, tri_verts, tri_bytes,
                                cudaMemcpyHostToDevice)) == cudaSuccess) &&
        (edge_bytes == 0 || (err = cudaMemcpy(d_edge, edge_verts, edge_bytes,
                                cudaMemcpyHostToDevice)) == cudaSuccess)) {
        const int threads = 256;
        unsigned long long want =
            (total + threads - 1) / static_cast<unsigned long long>(threads);
        const int blocks = static_cast<int>(want < 4096ULL ? want : 4096ULL);
        bhrt_ts_precheck<<<blocks, threads>>>(
            d_pent, n_pents, d_tri, n_tris, d_edge, n_edges,
            n_vertices, total, d_out, cap, d_cnt);
        err = cudaGetLastError();
        if (err == cudaSuccess) err = cudaDeviceSynchronize();
        unsigned long long count = 0;
        if (err == cudaSuccess)
            err = cudaMemcpy(&count, d_cnt, sizeof(count),
                             cudaMemcpyDeviceToHost);
        if (err == cudaSuccess) {
            const unsigned long long n_copy = count < cap ? count : cap;
            if (n_copy > 0)
                err = cudaMemcpy(survivors, d_out,
                                 n_copy * sizeof(unsigned long long),
                                 cudaMemcpyDeviceToHost);
            if (err == cudaSuccess) {
                *n_survivors = n_copy;
                rc = (count > cap) ? 2 : 0;   // overflow -> incomplete list
            }
        }
    }
    if (d_pent) cudaFree(d_pent);
    if (d_tri)  cudaFree(d_tri);
    if (d_edge) cudaFree(d_edge);
    if (d_out)  cudaFree(d_out);
    if (d_cnt)  cudaFree(d_cnt);
    return rc;
}

}  // extern "C"
