// smoke.cu -- minimal CUDA sanity check, independent of the bhrt code.
// Build:  nvcc -arch=native cuda/smoke.cu -o build/smoke
// Run:    ./build/smoke
// If this crashes (munmap/abort), the problem is the CUDA runtime/driver on
// this machine, not the bhrt scorer. If it prints "OK", the issue is in our
// launcher and I'll dig in.

#include <cstdio>
#include <cuda_runtime.h>

__global__ void add_one(float* a, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) a[i] += 1.0f;
}

int main() {
    int n = 0;
    cudaError_t e = cudaGetDeviceCount(&n);
    std::printf("cudaGetDeviceCount: %s (devices=%d)\n", cudaGetErrorString(e), n);
    if (e != cudaSuccess || n < 1) return 1;

    cudaDeviceProp p{};
    cudaGetDeviceProperties(&p, 0);
    std::printf("device 0: %s, sm_%d%d\n", p.name, p.major, p.minor);

    const int N = 256;
    float h[N];
    for (int i = 0; i < N; ++i) h[i] = (float)i;

    float* d = nullptr;
    e = cudaMalloc(&d, N * sizeof(float));
    std::printf("cudaMalloc:        %s\n", cudaGetErrorString(e));
    if (e != cudaSuccess) return 2;

    e = cudaMemcpy(d, h, N * sizeof(float), cudaMemcpyHostToDevice);
    std::printf("cudaMemcpy H2D:    %s\n", cudaGetErrorString(e));

    add_one<<<(N + 255) / 256, 256>>>(d, N);
    e = cudaGetLastError();
    std::printf("kernel launch:     %s\n", cudaGetErrorString(e));
    e = cudaDeviceSynchronize();
    std::printf("synchronize:       %s\n", cudaGetErrorString(e));

    e = cudaMemcpy(h, d, N * sizeof(float), cudaMemcpyDeviceToHost);
    std::printf("cudaMemcpy D2H:    %s\n", cudaGetErrorString(e));

    e = cudaFree(d);
    std::printf("cudaFree:          %s\n", cudaGetErrorString(e));

    std::printf("result h[10] = %.1f (expected 11.0)\nOK\n", h[10]);
    return 0;
}
