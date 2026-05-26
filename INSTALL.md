# Installation and Usage Guide

This is the full, end-to-end guide for building and running the
**bhrt_trisect** BHRT trisection pipeline plus GPU genus-reduction engine.

The pipeline is designed so the **C++ core works on its own** with
nothing but a C++20 compiler -- you can scan triangulations, certify
ts-tricolourings, extract trisection diagrams, compute untwisted
invariants, and run the beam search without installing anything else.

Three optional backends unlock additional capability:

| Backend  | Purpose                                                                |
|----------|------------------------------------------------------------------------|
| Regina   | Load `.esig` / `.rga` census files, apply topology-aware Pachner moves |
| SageMath | Exact integer linear algebra for the invariants engine (Z, Z[t,t^-1]) |
| CUDA     | GPU-batched move-candidate scoring for the beam search                 |

You can install them in any order, or skip them entirely.

---

## Table of contents

1. [System requirements](#1-system-requirements)
2. [Quick install (no dependencies)](#2-quick-install-no-dependencies)
3. [Install Regina](#3-install-regina)
4. [Install SageMath](#4-install-sagemath)
5. [Install CUDA](#5-install-cuda)
6. [Build with all backends](#6-build-with-all-backends)
7. [Downloading the Dim4Census data](#7-downloading-the-dim4census-data)
8. [Running the C++ CLI](#8-running-the-c-cli)
9. [Running the Python orchestration layer](#9-running-the-python-orchestration-layer)
10. [End-to-end recipe: census -> diagram -> invariants](#10-end-to-end-recipe-census---diagram---invariants)
11. [Reproducing the paper](#11-reproducing-the-paper)
12. [Troubleshooting](#12-troubleshooting)

---

## 1. System requirements

| Platform                 | Status                                                       |
|--------------------------|--------------------------------------------------------------|
| Linux (Ubuntu 22.04 LTS) | Recommended; all backends tested                              |
| Linux (Debian/Fedora)    | Supported; Regina/Sage have official repos                    |
| macOS (12+)              | Supported; Regina has signed installer, Sage via brew or DMG  |
| Windows 10/11            | Supported via WSL2 (strongly recommended); native build works for the C++ core but Regina is fiddly |

Minimum hardware:

* CPU with x86-64 or arm64
* 4 GB RAM (for the 2p/4p census)
* 16 GB RAM (for the 6p census)
* 64-128 GB RAM (for paper-grade Pachner-graph search)
* (optional) NVIDIA GPU with 8 GB+ VRAM for the CUDA scorer

Required tools:

| Tool      | Version | Notes                              |
|-----------|---------|------------------------------------|
| C++ compiler | g++ 11+, clang 14+, or MSVC 19.30+ | C++20 features used |
| C compiler   | gcc 11+ or clang 14+               | for the CLI driver  |
| GNU make     | 4.0+                               | or CMake 3.20+      |
| Python       | 3.10+                              | orchestration only  |

Optional tools:

| Tool      | Version | Purpose                          |
|-----------|---------|----------------------------------|
| CMake     | 3.20+   | alternative to plain Makefile    |
| CUDA      | 12.0+   | GPU candidate scoring            |
| pybind11  | 2.10+   | Python bindings to the C++ core  |

---

## 2. Quick install (no dependencies)

This is the **fastest path** to a working binary. No Regina, no Sage,
no CUDA needed.

### 2.1 Get the code

```bash
git clone https://github.com/<your-fork>/bhrt_trisect.git
cd bhrt_trisect
```

(If you have it as a tarball, untar it and `cd` into the directory.)

### 2.2 Build the C++ CLI

```bash
make                # ~5 seconds on a modern laptop
```

This produces `build/bhrt-cli`, a stand-alone binary that contains the
full mathematical pipeline. No installation step.

If `make` is not available, the manual build is two commands:

```bash
g++ -O2 -std=c++20 -c cpp/triangulation.cpp cpp/ts_color.cpp \
                      cpp/diagram.cpp cpp/invariants.cpp \
                      cpp/search_cpu.cpp cpp/bhrt_c_api.cpp
gcc -O2 -c cpp/bhrt_cli.c -Icpp
g++ -O2 *.o -o bhrt-cli
```

### 2.3 Smoke test

```bash
./build/bhrt-cli info       # prints the 2-pentachoron S^4 skeleton
./build/bhrt-cli scan       # runs the ts-tricolouring counters
./build/bhrt-cli search 32 3 10 0   # beam search demo
```

Expected output of `info`:

```
label:        S4_2p
pentachora:   2
vertices:     5
edges:        10
triangles:    10
tetrahedra:   5
euler chi:    2
is_closed:    1
is_valid:     1
```

If you see those numbers, the C++ core is working. Skip to
[Section 8](#8-running-the-c-cli) for usage.

### 2.4 Optional: run the Python test suite

```bash
PYTHONPATH=python python3 -m pytest tests/ -q
```

You should see `34 passed, 2 skipped` (the two skipped tests need
Regina; see the next section).

---

## 3. Install Regina

Regina is the production-grade 4-manifold triangulation library
(<https://regina-normal.github.io>). bhrt_trisect uses it for two
things:

1. Loading `.esig` (isomorphism-signature) and `.rga` (Regina XML)
   census files.
2. Applying topology-aware Pachner moves and edge collapses in the
   beam search.

You need Regina **7.4 or later**.

### 3.1 Linux (Ubuntu / Debian)

Add Regina's official APT repository:

```bash
# Ubuntu 22.04 / 24.04
sudo bash -c 'echo "deb [signed-by=/etc/apt/keyrings/regina.gpg] \
    https://people.debian.org/~bcasselm/regina/ ./" \
    > /etc/apt/sources.list.d/regina.list'

sudo mkdir -p /etc/apt/keyrings
sudo wget -O /etc/apt/keyrings/regina.gpg \
    https://people.debian.org/~bcasselm/regina/regina.gpg

sudo apt update
sudo apt install -y regina-normal regina-normal-doc \
                     python3-regina
```

Verify:

```bash
regina-python -c "import regina; print(regina.versionString())"
```

You should see `7.4` or higher.

### 3.2 Linux (Fedora / RHEL / Arch)

Fedora packages Regina natively from version 38 onward:

```bash
sudo dnf install regina-normal python3-regina
```

Arch users: it is in the AUR as `regina-normal`:

```bash
yay -S regina-normal
```

### 3.3 macOS

Regina ships a signed `.dmg` for macOS at:

> <https://regina-normal.github.io/install.html#macos>

Download the 7.4+ DMG, drag `Regina.app` to `/Applications`, then add
the embedded Python bindings to your shell:

```bash
echo 'export PYTHONPATH="/Applications/Regina.app/Contents/MacOS/python:$PYTHONPATH"' \
    >> ~/.zshrc
source ~/.zshrc
```

Verify:

```bash
python3 -c "import regina; print(regina.versionString())"
```

Alternative: install via Homebrew:

```bash
brew install --cask regina-normal
```

### 3.4 Windows (via WSL2)

Strongly recommended over native Windows. Install WSL2 with Ubuntu:

```powershell
wsl --install -d Ubuntu-22.04
```

Reboot, log into Ubuntu, then follow [Section 3.1](#31-linux-ubuntu--debian).

### 3.5 Building Regina from source

If your distro has Regina < 7.4, build from source:

```bash
git clone https://github.com/regina-normal/regina.git
cd regina
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
sudo make install
```

Dependencies that Regina needs first (Ubuntu):

```bash
sudo apt install -y build-essential cmake libgmp-dev libgmpxx4ldbl \
                    libpopt-dev libxml2-dev libsource-highlight-dev \
                    libtokyocabinet-dev libjansson-dev libtool \
                    libqt5svg5-dev qtbase5-dev qttools5-dev-tools \
                    python3-dev pybind11-dev
```

### 3.6 Linking bhrt_trisect against Regina

After Regina is installed, rebuild bhrt_trisect with Regina enabled:

```bash
cd bhrt_trisect
cmake -S cpp -B build -DBHRT_HAS_REGINA=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

This produces `build/bhrt-cli` *and* `build/libbhrt_core.a`, both with
the `BHRT_HAS_REGINA` compile-time flag set. The Python loader
auto-detects Regina; no further configuration needed.

Verify:

```bash
PYTHONPATH=python python3 -c \
    "import regina; from bhrt_trisect.tri_io import load_esig_to_triangulations; \
     print('Regina', regina.versionString(), 'loaded')"
```

---

## 4. Install SageMath

SageMath (<https://www.sagemath.org>) is used for **exact integer linear
algebra** in the invariants engine: Smith normal form over Z, symplectic
basis changes, twisted homology over Z[t, t^-1]. The Python fallback
(SymPy) is correct but ~10x slower on larger matrices.

You need **SageMath 10.0 or later**.

### 4.1 Linux (Ubuntu / Debian)

```bash
sudo apt install -y sagemath sagemath-doc
```

(Note: Ubuntu's package is typically 1-2 minor versions behind. For
the latest, use Sage's prebuilt tarball -- see 4.4.)

### 4.2 Linux (Fedora / Arch)

```bash
sudo dnf install sagemath              # Fedora
yay -S sagemath                         # Arch
```

### 4.3 macOS

Download the macOS DMG from <https://www.sagemath.org/download-mac.html>
and drag SageMath.app to /Applications.

To use it from the terminal:

```bash
echo 'alias sage="/Applications/SageMath.app/Contents/Resources/sage/sage"' \
    >> ~/.zshrc
source ~/.zshrc
sage --version
```

Or install via Homebrew (slower but cleaner):

```bash
brew install --cask sage
```

### 4.4 Prebuilt binary tarball (any Linux)

For the latest Sage on any distro:

```bash
cd /opt
sudo wget https://mirrors.mit.edu/sage/linux/64bit/sage-10.4-Linux.tar.bz2
sudo tar xjf sage-10.4-Linux.tar.bz2
sudo ln -s /opt/SageMath/sage /usr/local/bin/sage
sage --version
```

### 4.5 Building Sage from source

This is the canonical way to get an optimised, integrated Sage build,
but it takes 2-4 hours and 10 GB of disk:

```bash
git clone https://github.com/sagemath/sage.git
cd sage
make configure
./configure
MAKE='make -j$(nproc)' make
sudo ln -s $PWD/sage /usr/local/bin/sage
```

### 4.6 Wiring Sage into bhrt_trisect

The Sage backend is opt-in. Tell the invariants engine to use it:

```python
from bhrt_trisect.invariants import set_algebra_backend
set_algebra_backend("sage")
```

Or, from the CLI, run the invariants script under `sage -python`:

```bash
sage -python -m bhrt_trisect.cli invariants \
    --input results/best_diagrams.jsonl \
    --backend sage \
    --out results/invariants.jsonl
```

The Sage helpers live in `sage/invariants.py` and
`sage/twisted_invariants.py`; both are imported automatically when
SageMath is available on `PYTHONPATH`.

---

## 5. Install CUDA

The CUDA candidate scorer is **optional**. The Python NumPy backend
and the C++ CPU scorer both produce identical scores; CUDA just runs
them ~10-50x faster on batches of >10000 candidates.

You need **CUDA 12.0 or later** and an NVIDIA GPU with compute
capability >= 7.0 (V100, T4, A100, RTX 20-series, RTX 30-series, RTX
40-series, H100).

### 5.1 Linux (Ubuntu / Debian)

NVIDIA's official CUDA installer:

```bash
# Step 1: install the NVIDIA driver
sudo apt install -y nvidia-driver-550

# Step 2: install the CUDA toolkit
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt update
sudo apt install -y cuda-toolkit-12-4

# Step 3: add CUDA to PATH
echo 'export PATH=/usr/local/cuda-12.4/bin:$PATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=/usr/local/cuda-12.4/lib64:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc

# Step 4: verify
nvidia-smi      # should list your GPU
nvcc --version  # should print CUDA 12.4
```

### 5.2 macOS

CUDA is not supported on Apple Silicon Macs. If you have an older
Intel Mac with an NVIDIA GPU, follow NVIDIA's legacy installer; we
recommend using a Linux box or cloud GPU instead.

### 5.3 Windows

Download the installer from <https://developer.nvidia.com/cuda-downloads>
and follow the prompts. The bhrt_trisect build supports MSVC 19.30+
and the CUDA-MSVC toolchain.

### 5.4 Python GPU backends

Two Python paths are supported:

```bash
# Option A: PyCUDA (preferred; lower overhead)
pip install pycuda

# Option B: Numba CUDA (easier to install)
pip install numba
```

The bhrt_trisect package picks PyCUDA if available, then Numba, then
NumPy. Force a specific backend with:

```bash
export BHRT_GPU_BACKEND=cuda    # PyCUDA
export BHRT_GPU_BACKEND=numba   # Numba
export BHRT_GPU_BACKEND=cpu     # NumPy reference
```

### 5.5 Building the C++ binary with CUDA

```bash
cmake -S cpp -B build -DBHRT_HAS_CUDA=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

Now `bhrt-cli search --device 0` dispatches scoring to the GPU.

---

## 6. Build with all backends

The recommended production build:

```bash
cd bhrt_trisect

cmake -S cpp -B build \
      -DBHRT_HAS_REGINA=ON \
      -DBHRT_HAS_CUDA=ON \
      -DBHRT_BUILD_PYBIND=ON \
      -DCMAKE_BUILD_TYPE=Release

cmake --build build -j$(nproc)

# Install the Python package in editable mode
cd python && pip install -e . && cd ..
```

You now have:

* `build/bhrt-cli` -- standalone C-language CLI
* `build/libbhrt_core.a` -- the C++ core library
* `build/bhrt_trisect_core*.so` -- pybind11 Python extension
* `python/bhrt_trisect/` -- Python package (orchestration + CLI)

Verify everything is wired:

```bash
./build/bhrt-cli info
python3 -c "import bhrt_trisect, regina; print('all good')"
sage -python -c "from bhrt_trisect.invariants import set_algebra_backend; set_algebra_backend('sage')"
```

---

## 7. Downloading the Dim4Census data

The Dim4Census repository provides the public closed orientable 4D
triangulation censuses used by every paper in this area. We benchmark
against the documented counts:

* 2 pentachora orientable census: 8 triangulations, 2 admit a ts-tricolouring
* 4 pentachora orientable census: 784 triangulations, 15 admit
* 6 pentachora orientable census: 440,495 triangulations, 445 admit

### 7.1 Clone the data

```bash
cd bhrt_trisect/data
git clone https://github.com/regina-normal/Dim4Census.git Census
ls Census/
```

You should see `2p-closedOrientable.esig`, `4p-closedOrientable.esig`,
`6p-closedOrientable.esig`, and several `.rga` files.

### 7.2 Validate against the manifest

```bash
cat data/manifests/census.json
```

The manifest pins exact expected counts for each file. Run the
regression test:

```bash
PYTHONPATH=python python3 -m pytest tests/test_ts_color.py -k spreer -q
```

(This test is skipped automatically if `data/Census/` is not present.)

---

## 8. Running the C++ CLI

The C-language CLI driver `build/bhrt-cli` is the fastest way to
exercise the C++ core.

### 8.1 Subcommands

| Subcommand                  | Purpose                                           |
|-----------------------------|---------------------------------------------------|
| `bhrt-cli info`             | print 2p S^4 skeleton statistics                  |
| `bhrt-cli isosig`           | print canonical edge-degree isomorphism signature |
| `bhrt-cli scan`             | run the six ts-tricolouring counters              |
| `bhrt-cli search B H I S`   | beam search (beam=B, height=H, iter=I, seed=S)    |

### 8.2 Examples

```bash
$ ./build/bhrt-cli info
label:        S4_2p
pentachora:   2
vertices:     5
edges:        10
triangles:    10
tetrahedra:   5
euler chi:    2
is_closed:    1
is_valid:     1

$ ./build/bhrt-cli isosig
2,2,2,2,2,2,2,2,2,2|AAEAAQEBAQECAQEDAQEEAQAAAQABAQACAQADAQAE/gA

$ ./build/bhrt-cli scan
n_total                    = 1
n_pass_precheck            = 1
n_with_tricolouring        = 1
n_with_c_tricolouring      = 1
n_with_ts_tricolouring     = 0
n_ts_tricolourings_total   = 0

$ ./build/bhrt-cli search 64 4 20 42
best pentachora found: 2
```

### 8.3 Embedding the C ABI in your own program

Link against the C++ core and include `bhrt_c_api.h`:

```c
#include "bhrt_c_api.h"
#include <stdio.h>

int main(void) {
    bhrt_triangulation* T = bhrt_s4_minimal();
    char sig[2048];
    bhrt_edge_degree_isosig(T, sig, sizeof(sig));
    printf("isosig = %s\n", sig);
    bhrt_free_triangulation(T);
    return 0;
}
```

Compile:

```bash
gcc -O2 my_program.c -Icpp -Lbuild -lbhrt_core -lstdc++ -o my_program
```

---

## 9. Running the Python orchestration layer

The Python layer drives census scans, Sage-backed invariants, and the
search engine. Install it editable:

```bash
cd python && pip install -e . && cd ..
```

### 9.1 Subcommands

| Subcommand                          | Purpose                                |
|--------------------------------------|----------------------------------------|
| `python -m bhrt_trisect.cli scan-census <input>` | ts-tricolouring scan of an .esig file |
| `python -m bhrt_trisect.cli reduce-genus`         | beam search for lower trisection genus |
| `python -m bhrt_trisect.cli invariants`           | compute FKSZ untwisted invariants      |

Use `--input demo` to run any subcommand on the built-in 2-pentachoron
S^4 triangulation, requiring no external data.

### 9.2 Examples

```bash
# scan the 6p closed orientable census (requires Regina)
python3 -m bhrt_trisect.cli scan-census \
    data/Census/6p-closedOrientable.esig \
    --out results/scan_6p.json

# beam search with GPU scorer
python3 -m bhrt_trisect.cli reduce-genus \
    --input data/Census/6p-closedOrientable.esig \
    --beam 256 --height 6 --device 0 --seed 42 \
    --checkpoint-dir results/search_6p/

# invariants via Sage (exact algebra)
sage -python -m bhrt_trisect.cli invariants \
    --input results/best_diagrams.jsonl \
    --backend sage \
    --out results/invariants.jsonl
```

### 9.3 Programmatic API

```python
from bhrt_trisect.triangulation import s4_minimal
from bhrt_trisect.ts_color import enumerate_ts_tricolourings, count_admit_trisection
from bhrt_trisect.diagram import extract_diagram
from bhrt_trisect.invariants import compute_invariants
from bhrt_trisect.search_cpu import SearchConfig, beam_search

T = s4_minimal()

# 1. ts-tricolouring scan
counters = count_admit_trisection([T])
print(counters)

# 2. enumerate ts-tricolourings
ts_list = enumerate_ts_tricolourings(T)
print(f"{len(ts_list)} ts-tricolourings found")

# 3. extract a trisection diagram
if ts_list:
    diag = extract_diagram(T, ts_list[0])
    print(f"genus = {diag.genus}")

# 4. compute untwisted invariants
    inv = compute_invariants(diag)
    print(inv.as_dict())

# 5. beam search
config = SearchConfig(beam_width=64, excess_height=4,
                     time_limit_seconds=60.0, max_iterations=100, seed=0)
pareto, beam = beam_search(T, config)
for s in pareto.items:
    print(s.pentachora, s.best_genus)
```

---

## 10. End-to-end recipe: census -> diagram -> invariants

This recipe reproduces the four-stage pipeline on the 4-pentachoron
closed orientable census. Requires Regina to be installed.

### 10.1 Stage 1 -- scan the census

```bash
python3 -m bhrt_trisect.cli scan-census \
    data/Census/4p-closedOrientable.esig \
    --out results/scan_4p.json
cat results/scan_4p.json
```

Expected output (matches the Spreer-Tillmann published count):

```json
{
  "corpus_path": "data/Census/4p-closedOrientable.esig",
  "n_total": 784,
  "n_pass_precheck": 784,
  "n_with_tricolouring": 30,
  "n_with_c_tricolouring": 22,
  "n_with_ts_tricolouring": 15,
  "n_ts_tricolourings_total": 28,
  "wall_seconds": 1.4
}
```

### 10.2 Stage 2 -- extract diagrams for the 15 positive cases

```bash
python3 -m bhrt_trisect.cli invariants \
    --input data/Census/4p-closedOrientable.esig \
    --out results/invariants_4p.jsonl
wc -l results/invariants_4p.jsonl
```

Each line is a JSON record:

```json
{
  "triangulation_label": "...",
  "pentachora": 4,
  "genus": 3,
  "h1": {"free": 0, "torsion": []},
  "h2": {"free": 1, "torsion": []},
  "h3": {"free": 0, "torsion": []},
  "signature": 1,
  "parity": "odd"
}
```

### 10.3 Stage 3 -- beam-search any high-genus case

```bash
python3 -m bhrt_trisect.cli reduce-genus \
    --input data/Census/4p-closedOrientable.esig \
    --beam 128 --height 4 --device -1 \
    --checkpoint-dir results/search_4p/
```

The search writes `results/search_4p/beam.jsonl` and
`results/search_4p/pareto.jsonl` as it progresses. You can interrupt
with Ctrl-C; restart resumes from the last checkpoint.

### 10.4 Stage 4 -- twisted invariants (Sage)

```bash
sage -python -m bhrt_trisect.cli invariants \
    --input results/best_diagrams.jsonl \
    --backend sage \
    --out results/twisted_invariants.jsonl
```

---

## 11. Reproducing the paper

The single command that rebuilds every table and figure:

```bash
python3 scripts/reproduce_paper.py \
    --out-dir results/paper \
    --manifest data/manifests/census.json \
    --seeds 8
```

Outputs:

* `results/paper/census_table.json` -- Spreer-Tillmann counter table
* `results/paper/genus_histogram.png` -- distribution of trisection genus
* `results/paper/search_pareto.png` -- (pentachora, genus) Pareto front
* `results/paper/timing_table.csv` -- per-corpus wall time

See `docs/reproducibility.md` for the full negative-result protocol and
deterministic-seed setup.

---

## 12. Troubleshooting

### 12.1 `make` succeeds but `./build/bhrt-cli info` segfaults

You probably linked against a different C++ standard library than was
used to compile the objects. Clean rebuild:

```bash
make clean && make
```

### 12.2 `import regina` works but bhrt_trisect can't load .esig

Check the Regina version:

```bash
python3 -c "import regina; print(regina.versionString())"
```

If it prints `< 7.4`, upgrade Regina; the .esig parser changed in 7.4.

### 12.3 Sage backend rejects matrices over Z

Sage requires `Matrix(ZZ, ...)` to construct an integer matrix. The
shim `sage/invariants.py` does this for you; if you bypass it and pass
floats, Sage will fall back to the rational backend. Always pass
`numpy.int64` arrays into bhrt_trisect.

### 12.4 CUDA kernel fails with "no kernel image is available"

The kernels are compiled for compute capability 7.0, 8.0, and 8.9.
For older GPUs (Maxwell/Pascal), rebuild with:

```bash
cmake -S cpp -B build -DBHRT_HAS_CUDA=ON \
      -DCMAKE_CUDA_ARCHITECTURES="60;61;70;80;89"
```

### 12.5 Beam search "runs forever"

Set tighter budgets:

```bash
./build/bhrt-cli search 32 2 50 0    # smaller beam, lower height
```

or use the Python CLI with `--time` cap:

```bash
python3 -m bhrt_trisect.cli reduce-genus --input demo \
    --beam 64 --height 3 --time 60
```

Regina's documentation also warns explicitly that retriangulation
searches can be superexponential; cap aggressively for large inputs.

### 12.6 Pytest reports `source code string cannot contain null bytes`

Bytecode cache corruption from a parent process. Force pytest to use a
fresh cache directory:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/tmp/bhrt_pyc \
    python3 -B -m pytest tests/ -p no:cacheprovider \
    --import-mode=importlib
```

### 12.7 macOS: "library not loaded: libRegina.dylib"

Add the Regina library path to your environment:

```bash
echo 'export DYLD_LIBRARY_PATH="/Applications/Regina.app/Contents/Frameworks:$DYLD_LIBRARY_PATH"' \
    >> ~/.zshrc
source ~/.zshrc
```

### 12.8 Windows native build: "C++20 not supported"

Use MSVC 19.30 (Visual Studio 2022) or later, or switch to WSL2.

### 12.9 Where do I file bugs?

* C++ core bugs: open an issue in this repo with the output of
  `./build/bhrt-cli info` and the relevant input.
* Regina bugs: <https://github.com/regina-normal/regina/issues>.
* Sage bugs: <https://github.com/sagemath/sage/issues>.

---

## Appendix A: One-shot install scripts

### A.1 Ubuntu 22.04, no GPU, no Sage

```bash
sudo apt update
sudo apt install -y build-essential gcc g++ make python3 python3-pip \
                    python3-numpy python3-sympy python3-pytest

# Regina
sudo apt install -y regina-normal python3-regina

# bhrt_trisect
git clone https://github.com/<your-fork>/bhrt_trisect.git
cd bhrt_trisect
make
./build/bhrt-cli info
```

### A.2 Ubuntu 22.04, NVIDIA GPU, all backends

```bash
# C++, Regina, Sage, CUDA
sudo apt update
sudo apt install -y build-essential gcc g++ make cmake python3 python3-pip \
                    regina-normal python3-regina sagemath
sudo apt install -y nvidia-driver-550

wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt update
sudo apt install -y cuda-toolkit-12-4
pip install numba pycuda numpy sympy pytest

# bhrt_trisect
git clone https://github.com/<your-fork>/bhrt_trisect.git
cd bhrt_trisect
cmake -S cpp -B build \
      -DBHRT_HAS_REGINA=ON -DBHRT_HAS_CUDA=ON \
      -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
cd python && pip install -e . && cd ..

# verify
./build/bhrt-cli scan
python3 -c "import bhrt_trisect, regina; print('OK')"
sage -python -c "from bhrt_trisect.invariants import set_algebra_backend; set_algebra_backend('sage'); print('OK')"
```

### A.3 macOS Apple Silicon, no GPU, all CPU backends

```bash
brew install gcc make cmake python@3.11
brew install --cask regina-normal sage

git clone https://github.com/<your-fork>/bhrt_trisect.git
cd bhrt_trisect
make
./build/bhrt-cli info

cd python && pip3 install -e . && cd ..
PYTHONPATH=python python3 -m pytest tests/ -q
```

---

## Appendix B: Directory map

```
bhrt_trisect/
  README.md
  INSTALL.md          (this file)
  Makefile
  build/              (compiled artefacts; created by `make`)
  cpp/                (C++ core + C CLI driver)
    bhrt_trisect.hpp
    bhrt_c_api.h
    triangulation.cpp
    ts_color.cpp
    diagram.cpp
    invariants.cpp
    search_cpu.cpp
    bhrt_c_api.cpp
    bhrt_cli.c
    CMakeLists.txt
  cuda/               (CUDA scoring kernels)
    score_moves.cu
    candidate_pack.cu
  python/             (orchestration + executable spec)
    bhrt_trisect/
      __init__.py
      triangulation.py
      isosig.py
      tri_io.py
      ts_color.py
      bhrt.py
      diagram.py
      invariants.py
      moves.py
      search_cpu.py
      search_gpu.py
      bench.py
      cli.py
    pyproject.toml
  sage/               (SageMath exact-algebra helpers)
    invariants.py
    twisted_invariants.py
  data/
    manifests/census.json
    Census/           (Dim4Census files go here)
  scripts/            (user-facing entry points)
    scan_census.py
    reduce_genus.py
    compute_invariants.py
    reproduce_paper.py
  tests/              (pytest regression suite)
  docs/
    architecture.md
    reproducibility.md
```

---

Last updated: 2026-05.
