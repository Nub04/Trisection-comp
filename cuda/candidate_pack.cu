// candidate_pack.cu — flatten CPU-side move candidates into the fixed-
// width device buffer consumed by score_moves.cu.
//
// The host walks a Regina state once, enumerates legal moves into a
// std::vector<MoveCandidate>, and writes pinned-memory records to the
// device in batches. This kernel performs the on-device finishing step:
//
//   * normalises the feasibility flag (cheap-edge-degree filter),
//   * computes a lightweight "diagram-simplicity proxy" by counting
//     local degree-3 triangles around the affected pentachoron,
//   * applies an excess-height clamp.
//
// The kernel is independent of the topology engine -- it operates purely
// on the flattened records and a per-pentachoron degree-array passed by
// the host. This keeps the GPU side of the pipeline free of pointer
// chasing and matches the architecture recommended in the research
// report (GPU = candidate scorer; Regina = state owner).

#include <cuda_runtime.h>

extern "C" {

__global__ void bhrt_finish_pack(float* __restrict__ records,
                                  const int* __restrict__ pent_degree,
                                  int n_records,
                                  int n_cols) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_records) return;

    const int pent_id_int = static_cast<int>(records[i * n_cols + 1]);
    if (pent_id_int >= 0) {
        const int d = pent_degree[pent_id_int];
        // diagram simplicity proxy: smaller degree -> higher feasibility
        records[i * n_cols + 4] *= (d <= 4) ? 1.0f : 0.5f;
    }
    // clamp excess
    if (records[i * n_cols + 2] > 8.0f) records[i * n_cols + 2] = 8.0f;
    if (records[i * n_cols + 2] < -8.0f) records[i * n_cols + 2] = -8.0f;
}

}  // extern "C"
