// score_moves.cu — CUDA candidate scorer for the genus-reduction beam search.
//
// Input:
//   records  : float32 array of shape (N, 5) on device
//   weights  : float32 array of length 5 on device
// Output:
//   scores   : float32 array of length N on device
//
// Each thread handles one candidate move record. The scoring function
// is the bit-exact translation of python/bhrt_trisect/search_gpu.py
// `_score_numpy`. Identical scores are produced across CPU, Numba CUDA,
// and this kernel.
//
// The record layout is:
//   col 0 : move_code (23, 33, 44, 15, 51, 99)
//   col 1 : locator first int (triangle/edge/pent id)
//   col 2 : delta_pent_estimate (signed)
//   col 3 : current_pent_count
//   col 4 : feasibility flag (1.0 if parent state had a ts-tricolouring)

#include <cuda_runtime.h>
#include <math_constants.h>
#include <cstdlib>
#include <cstdio>

// LeakSanitizer region control. These are provided by the ASan runtime when
// the program is linked with -fsanitize=address; declared weak so the symbols
// resolve to null (and the guard becomes a no-op) in non-sanitized builds.
extern "C" {
    void __lsan_disable(void) __attribute__((weak));
    void __lsan_enable(void)  __attribute__((weak));
}
namespace {
// While alive, allocations are not tracked by LeakSanitizer. Used to bracket
// the CUDA driver's one-time lazy initialisation, whose internal allocations
// live in anonymous (symbol-less) memory that pattern suppressions cannot
// match. These allocations belong to the closed NVIDIA driver, not to bhrt.
struct LsanGuard {
    LsanGuard()  { if (__lsan_disable) __lsan_disable(); }
    ~LsanGuard() { if (__lsan_enable)  __lsan_enable();  }
};
}  // namespace

extern "C" {

__global__ void bhrt_score_moves(const float* __restrict__ records,
                                  const float* __restrict__ weights,
                                  float* __restrict__ out,
                                  int n_records,
                                  int n_cols) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_records) return;

    const float delta_pent = records[i * n_cols + 2];
    const float pent_count = records[i * n_cols + 3];
    const float feasible   = records[i * n_cols + 4];

    const float delta = -delta_pent;
    const float size_proxy = -pent_count / 100.0f;
    const float excess = delta_pent > 0.0f ? delta_pent : 0.0f;
    const float novelty = 0.0f;  // injected on the host side as a separate column in the
                                  // production build that carries a hashed signature.

    float s = weights[0] * delta
            + weights[1] * feasible
            + weights[2] * size_proxy
            + weights[3] * excess
            + weights[4] * novelty;
    out[i] = s;
}

// Warp-level top-K reduction (used by the second-stage shortlist kernel
// once scoring is complete). The implementation here picks the highest
// `k` scores across N records using a single warp's worth of shared
// memory, suitable for k <= 256.

__global__ void bhrt_top_k(const float* __restrict__ scores,
                            int* __restrict__ indices_out,
                            int n_records, int k) {
    extern __shared__ float smem[];
    float* best_scores = smem;
    int*   best_idx    = reinterpret_cast<int*>(smem + k);
    int t = threadIdx.x;

    if (t < k) {
        best_scores[t] = -CUDART_INF_F;
        best_idx[t] = -1;
    }
    __syncthreads();

    for (int i = t; i < n_records; i += blockDim.x) {
        float s = scores[i];
        // Find slot with smallest current best.
        int min_slot = 0;
        float min_val = best_scores[0];
        for (int j = 1; j < k; ++j) {
            if (best_scores[j] < min_val) {
                min_val = best_scores[j];
                min_slot = j;
            }
        }
        if (s > min_val) {
            best_scores[min_slot] = s;
            best_idx[min_slot] = i;
        }
    }
    __syncthreads();
    if (t < k) indices_out[t] = best_idx[t];
}

// Host launcher: copy records up, score on the GPU, copy scores back.
// Returns 0 on success; non-zero if no usable GPU is present or a CUDA call
// fails, in which case the caller (beamSearch) falls back to the CPU scorer.
// Scores match scoreCandidateCPU up to floating-point contraction (nvcc
// fuses a*b+c into FMA by default, so the last ulp may differ from the
// host); beamSearch self-verifies the first batch against the CPU reference
// with a small tolerance and disables the GPU path on any real mismatch.
int bhrt_score_moves_host(const float* records_host,
                          int n_records, int n_cols,
                          const float* weights_host,
                          float* scores_host) {
    if (n_records <= 0) return 0;

    // Do not let LeakSanitizer track the driver's lazy-init allocations made
    // during the CUDA calls below; they are owned by libcuda and freed only at
    // driver teardown. The guard covers every return path via RAII.
    LsanGuard lsan_guard;

    static bool s_warned = false;  // diagnose a GPU fallback at most once

    int dev = 0;
    cudaError_t e = cudaGetDeviceCount(&dev);
    if (e != cudaSuccess || dev == 0) {
        if (!s_warned) {
            s_warned = true;
            std::fprintf(stderr,
                "[bhrt] GPU scorer unavailable (%s; devices=%d); using CPU.\n",
                cudaGetErrorString(e), dev);
        }
        return 1;
    }

    // Release the CUDA context at process exit so the driver frees the
    // device-context allocations it would otherwise hold until teardown.
    // Registered ONLY when the leak-checking (ASan/LSan) runtime is linked
    // in -- detected via the weak __lsan_disable symbol -- because that is
    // its sole purpose: the reset must run before LSan's end-of-run check or
    // the context allocations are reported as leaks. In plain builds (e.g.
    // the CMake Release build) an atexit cudaDeviceReset can race the CUDA
    // runtime's own process-exit teardown, which manifests as glibc heap
    // errors ("munmap_chunk(): invalid pointer", "free(): invalid pointer")
    // on some driver stacks, WSL2 in particular -- so there we simply let
    // the driver tear the context down itself.
    static bool s_cleanup_registered = false;
    if (!s_cleanup_registered) {
        s_cleanup_registered = true;
        if (__lsan_disable)            // non-null <=> sanitizer runtime present
            std::atexit([] { cudaDeviceReset(); });
    }

    const size_t rec_bytes = static_cast<size_t>(n_records) * n_cols * sizeof(float);
    const size_t out_bytes = static_cast<size_t>(n_records) * sizeof(float);
    const size_t w_bytes   = 5 * sizeof(float);

    float* d_rec = nullptr;
    float* d_w   = nullptr;
    float* d_out = nullptr;
    int rc = 1;
    cudaError_t err = cudaSuccess;

    if ((err = cudaMalloc(&d_rec, rec_bytes)) == cudaSuccess &&
        (err = cudaMalloc(&d_w,   w_bytes))   == cudaSuccess &&
        (err = cudaMalloc(&d_out, out_bytes)) == cudaSuccess &&
        (err = cudaMemcpy(d_rec, records_host, rec_bytes,
                          cudaMemcpyHostToDevice)) == cudaSuccess &&
        (err = cudaMemcpy(d_w, weights_host, w_bytes,
                          cudaMemcpyHostToDevice)) == cudaSuccess) {
        const int threads = 256;
        const int blocks  = (n_records + threads - 1) / threads;
        bhrt_score_moves<<<blocks, threads>>>(d_rec, d_w, d_out, n_records, n_cols);
        // A launch-configuration / missing-kernel-image failure surfaces in
        // cudaGetLastError; an asynchronous execution failure surfaces in
        // cudaDeviceSynchronize. Check both, then the copy-back.
        err = cudaGetLastError();
        if (err == cudaSuccess) err = cudaDeviceSynchronize();
        if (err == cudaSuccess)
            err = cudaMemcpy(scores_host, d_out, out_bytes,
                             cudaMemcpyDeviceToHost);
        if (err == cudaSuccess) rc = 0;
    }
    if (rc != 0 && !s_warned) {
        s_warned = true;
        std::fprintf(stderr,
            "[bhrt] GPU scoring failed (%s); using CPU.\n",
            cudaGetErrorString(err));
    }
    if (d_rec) cudaFree(d_rec);
    if (d_w)   cudaFree(d_w);
    if (d_out) cudaFree(d_out);
    return rc;
}

}  // extern "C"
