# Top-level Makefile for the BHRT trisection pipeline.
#
# Produces three binaries under build/:
#
#   bhrt-cli    -- C-language CLI driver (info, scan, search, file I/O)
#   bhrt-bench  -- C++ benchmark driver (corpus scan + invariants + search)
#   bhrt-test   -- C++ native test suite (no Python required)
#
# Optional flags:
#   BHRT_HAS_REGINA=1   link against libregina (requires Regina 7.4+)
#   BHRT_HAS_CUDA=1     build the CUDA scoring kernels into bhrt-bench

# Use whatever g++/gcc is available (respects the environment and make's
# defaults; not pinned to any version).
# NOTE: if you ever see "g++-11: No such file", that name is coming from your
# shell environment (often conda exports CXX=g++-11), not from this file.
# Clear it with:  unset CXX CC
CXX      ?= g++
CC       ?= gcc
CXXFLAGS ?= -g -fsanitize=address -std=c++20 -Wall -Wno-unused-variable
CFLAGS   ?= -g -fsanitize=address -Wall
LDFLAGS  ?=

CORE_CPP := \
    cpp/triangulation.cpp \
    cpp/ts_color.cpp      \
    cpp/diagram.cpp       \
    cpp/invariants.cpp    \
    cpp/search_cpu.cpp    \
    cpp/bhrt_format.cpp   \
    cpp/bhrt_c_api.cpp

# --- optional Regina backend: lateral 2-4/3-3/4-2 moves + .esig loading ---
# Enable with:  make BHRT_HAS_REGINA=1 ...
# If Regina is INSTALLED (e.g. under /usr/local), the default header search
# already finds it; usually you only need to name the library:
#   make BHRT_HAS_REGINA=1 REGINA_CFLAGS="" REGINA_LIBS="-lregina-engine"
# (add -L/usr/local/lib to REGINA_LIBS if the linker can't find it, and set
#  LD_LIBRARY_PATH=/usr/local/lib at run time if needed).
# If a pkg-config file is present the defaults below pick it up automatically:
#   make BHRT_HAS_REGINA=1
# Do NOT add -I pointing at the Regina *source* tree (regina-7.4.1/engine):
# those headers need a generated regina-config.h from the build dir and will
# fail with "regina-config.h: No such file or directory".
ifeq ($(BHRT_HAS_REGINA),1)
  CORE_CPP      += cpp/tri_io.cpp cpp/regina_bridge.cpp
  REGINA_CFLAGS ?= $(shell pkg-config --cflags regina-engine 2>/dev/null)
  REGINA_LIBS   ?= $(shell pkg-config --libs regina-engine 2>/dev/null)
  # pkg-config is often absent for a /usr/local install; fall back to the
  # standard locations. Regina's headers include relative to the regina/
  # dir ("triangulation/dim4/triangulation4.h"), so BOTH include dirs are
  # required, and Regina's inline GMP usage requires linking gmp/gmpxx
  # directly (ld will not resolve them from an indirect DSO).
  ifeq ($(strip $(REGINA_CFLAGS)),)
    REGINA_CFLAGS := -I/usr/local/include -I/usr/local/include/regina
  endif
  ifeq ($(strip $(REGINA_LIBS)),)
    REGINA_LIBS := -L/usr/local/lib -lregina-engine -lgmp -lgmpxx
  endif
  CXXFLAGS      += -DBHRT_HAS_REGINA=1 $(REGINA_CFLAGS)
  LDFLAGS       += $(REGINA_LIBS)
endif

CLI_C    := cpp/bhrt_cli.c
BENCH_CPP := cpp/bhrt_bench.cpp
TEST_CPP  := cpp/bhrt_test.cpp
DIAG_CPP  := cpp/bhrt_diagram.cpp

BUILD := build
CORE_OBJ  := $(patsubst cpp/%.cpp,$(BUILD)/%.o,$(CORE_CPP))
CLI_OBJ   := $(patsubst cpp/%.c,$(BUILD)/%.o,$(CLI_C))
BENCH_OBJ := $(patsubst cpp/%.cpp,$(BUILD)/%.o,$(BENCH_CPP))
TEST_OBJ  := $(patsubst cpp/%.cpp,$(BUILD)/%.o,$(TEST_CPP))
DIAG_OBJ  := $(patsubst cpp/%.cpp,$(BUILD)/%.o,$(DIAG_CPP))

# --- optional CUDA backend: GPU batch scorer for the beam search ---
# Enable with:  make BHRT_HAS_CUDA=1 ...    (requires the CUDA toolkit + a GPU)
# The GPU path runs only when SearchConfig::use_gpu_scorer is set, and falls
# back to the CPU scorer if no GPU is present. Scores are bit-identical, so
# results never depend on whether CUDA is used. Override the arch if -arch=native
# is unsupported, e.g. CUDA_ARCH="-arch=sm_86".
#
# If your GPU is NEWER than your CUDA toolkit (e.g. a compute-capability 12.0
# Blackwell card with CUDA 12.0, which errors "Unsupported gpu architecture
# 'compute_120'"), compile to forward-compatible PTX for an arch the toolkit
# supports and let the driver JIT it at run time:
#   make BHRT_HAS_CUDA=1 CUDA_ARCH="-arch=compute_90"
ifeq ($(BHRT_HAS_CUDA),1)
  NVCC        ?= nvcc
  CUDA_ARCH   ?= -arch=native
  CUDA_LIBDIR ?= /usr/local/cuda/lib64
  CUDA_SRC    := cuda/score_moves.cu cuda/candidate_pack.cu cuda/ts_scan.cu
  CUDA_OBJ    := $(patsubst cuda/%.cu,$(BUILD)/%.o,$(CUDA_SRC))
  CXXFLAGS    += -DBHRT_HAS_CUDA=1
  LDFLAGS     += -L$(CUDA_LIBDIR) -lcudart
  CORE_OBJ    += $(CUDA_OBJ)
endif

CLI_BIN   := $(BUILD)/bhrt-cli
BENCH_BIN := $(BUILD)/bhrt-bench
TEST_BIN  := $(BUILD)/bhrt-test
DIAG_BIN  := $(BUILD)/bhrt-diagram

.PHONY: all cli bench test diag test-run smoke clean sage-check sage-demo

all: cli bench test diag

cli:   $(CLI_BIN)
bench: $(BENCH_BIN)
test:  $(TEST_BIN)
diag:  $(DIAG_BIN)

$(CLI_BIN):   $(CORE_OBJ) $(CLI_OBJ)   | $(BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(BENCH_BIN): $(CORE_OBJ) $(BENCH_OBJ) | $(BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(TEST_BIN):  $(CORE_OBJ) $(TEST_OBJ)  | $(BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(DIAG_BIN):  $(CORE_OBJ) $(DIAG_OBJ)  | $(BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(BUILD)/%.o: cpp/%.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -Icpp -c $< -o $@

$(BUILD)/%.o: cpp/%.c | $(BUILD)
	$(CC) $(CFLAGS) -Icpp -c $< -o $@

$(BUILD)/%.o: cuda/%.cu | $(BUILD)
	$(NVCC) $(CUDA_ARCH) -std=c++17 -DBHRT_HAS_CUDA=1 -Icpp -c $< -o $@

$(BUILD):
	mkdir -p $@

test-run: $(TEST_BIN)
	$(TEST_BIN)

smoke: $(CLI_BIN)
	@echo "--- info ---";   $(CLI_BIN) info
	@echo "--- isosig ---"; $(CLI_BIN) isosig
	@echo "--- scan ---";   $(CLI_BIN) scan
	@echo "--- write-demo ---"; $(CLI_BIN) write-demo /tmp/bhrt_smoke.bhrt
	@echo "--- scan-file ---";  $(CLI_BIN) scan-file /tmp/bhrt_smoke.bhrt

# --- SageMath (exact-arithmetic) cross-checks ---
sage-check:
	sage -python -c "import sys; sys.path.insert(0,'sage'); import invariants as I; I.smoke()"

sage-demo:
	sage -python sage/trisection_invariants.py --demo cp2

clean:
	rm -rf $(BUILD)
