# CMakeLists.txt -- build the full BHRT trisection pipeline.
#
# Builds four executables (bhrt-cli, bhrt-bench, bhrt-test, bhrt-diagram) from a
# shared static core library, with optional Regina, CUDA and AddressSanitizer.
# Unlike the hand-written Makefile, CMake tracks option changes and header
# dependencies, so toggling BHRT_HAS_CUDA / BHRT_HAS_REGINA rebuilds correctly
# without a manual clean.
#
# Quick start (CPU only):
#   cmake -S cpp -B build-cmake -DCMAKE_BUILD_TYPE=Release
#   cmake --build build-cmake -j
#   ctest --test-dir build-cmake            # runs the native test suite
#   ./build-cmake/bhrt-cli info
#
# With CUDA and Regina (Regina auto-detected from /usr/local etc.):
#   cmake -S cpp -B build-cmake -DBHRT_HAS_CUDA=ON -DBHRT_HAS_REGINA=ON
#   cmake --build build-cmake -j
#
# AddressSanitizer (use with CUDA OFF for clean leak reports):
#   cmake -S cpp -B build-asan -DBHRT_ASAN=ON
#
# Use a build directory SEPARATE from the Makefile's build/ (e.g. build-cmake).

cmake_minimum_required(VERSION 3.24)   # 3.24+ for CUDA_ARCHITECTURES "native"
project(bhrt_trisect LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(BHRT_HAS_REGINA   "Build with the Regina backend (lateral moves, .esig)" OFF)
option(BHRT_HAS_CUDA     "Build with the CUDA GPU move scorer"                   OFF)
option(BHRT_ASAN         "Build with AddressSanitizer"                           OFF)
option(BHRT_BUILD_PYBIND "Build the pybind11 module"                            OFF)

# ---- AddressSanitizer (C/C++ only; never passed to nvcc) -------------------
if(BHRT_ASAN)
    add_compile_options(
        $<$<COMPILE_LANGUAGE:C>:-g>   $<$<COMPILE_LANGUAGE:C>:-fsanitize=address>
        $<$<COMPILE_LANGUAGE:CXX>:-g> $<$<COMPILE_LANGUAGE:CXX>:-fsanitize=address>)
    add_link_options(-fsanitize=address)
endif()

# ---- Core library ----------------------------------------------------------
set(BHRT_CORE_SRC
    triangulation.cpp
    ts_color.cpp
    diagram.cpp
    invariants.cpp
    search_cpu.cpp
    bhrt_format.cpp
    bhrt_c_api.cpp)

add_library(bhrt_core STATIC ${BHRT_CORE_SRC})
target_include_directories(bhrt_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# ---- Optional Regina backend ----------------------------------------------
if(BHRT_HAS_REGINA)
    # Prefer a Regina CMake package if present; otherwise locate the installed
    # headers/library directly (covers a /usr/local install with
    # libregina-engine). Do NOT use the Regina *source* tree -- its headers need
    # a generated regina-config.h that lives only in Regina's own build dir.
    find_package(Regina QUIET)
    if(Regina_FOUND)
        target_link_libraries(bhrt_core PUBLIC Regina::engine)
        target_sources(bhrt_core PRIVATE tri_io.cpp regina_bridge.cpp)
        target_compile_definitions(bhrt_core PUBLIC BHRT_HAS_REGINA=1)
    else()
        find_path(REGINA_INCLUDE_DIR NAMES regina/triangulation/dim4.h)
        find_library(REGINA_LIBRARY  NAMES regina-engine regina)
        if(REGINA_INCLUDE_DIR AND REGINA_LIBRARY)
            message(STATUS "Regina: ${REGINA_LIBRARY} (headers ${REGINA_INCLUDE_DIR})")
            target_sources(bhrt_core PRIVATE tri_io.cpp regina_bridge.cpp)
            # Regina's own headers include relative to the regina/ dir (e.g.
            # "triangulation/dim4/triangulation4.h"), so put BOTH the prefix
            # (for <regina/...>) and the regina/ subdir on the include path.
            target_include_directories(bhrt_core PUBLIC
                ${REGINA_INCLUDE_DIR} ${REGINA_INCLUDE_DIR}/regina)
            # Regina's Integer uses GMP inline, so callers must link gmp/gmpxx
            # explicitly (ld won't use symbols from an indirect DSO).
            target_link_libraries(bhrt_core PUBLIC ${REGINA_LIBRARY} gmp gmpxx)
            target_compile_definitions(bhrt_core PUBLIC BHRT_HAS_REGINA=1)
        else()
            message(WARNING "Regina requested but not found; building without it")
        endif()
    endif()
endif()

target_link_libraries(bhrt_core
    PUBLIC
    gmp
    gmpxx
)
# ---- Optional CUDA backend -------------------------------------------------
if(BHRT_HAS_CUDA)
    # Prefer a /usr/local/cuda toolkit over an older apt nvcc (the
    # nvidia-cuda-toolkit package), whose ptxas may not support your GPU and
    # causes "Value 'sm_XX' is not defined" mismatches. Override explicitly
    # with -DCMAKE_CUDA_COMPILER=$(which nvcc) if your nvcc is elsewhere.
    if(NOT CMAKE_CUDA_COMPILER AND EXISTS /usr/local/cuda/bin/nvcc)
        set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
    endif()
    # MUST be set before enable_language(CUDA). Otherwise CMake tries to
    # auto-detect a default architecture and fails on some WSL/driver setups
    # ("Failed to detect a default CUDA architecture"). 'native' targets the
    # installed GPU; if it's rejected, override with a portable PTX value:
    #   -DCMAKE_CUDA_ARCHITECTURES=90-virtual   (compiles to PTX, JITs at runtime)
    # or a specific number (e.g. 120 for a Blackwell card on CUDA 13).
    if(NOT CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES native)
    endif()
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    set(CMAKE_CUDA_STANDARD 17)
    target_sources(bhrt_core PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../cuda/score_moves.cu
        ${CMAKE_CURRENT_SOURCE_DIR}/../cuda/candidate_pack.cu
        ${CMAKE_CURRENT_SOURCE_DIR}/../cuda/ts_scan.cu)
    set_target_properties(bhrt_core PROPERTIES
        CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
    target_compile_definitions(bhrt_core PUBLIC BHRT_HAS_CUDA=1)
    target_link_libraries(bhrt_core PUBLIC CUDA::cudart)
endif()

# ---- Executables -----------------------------------------------------------
add_executable(bhrt-cli     bhrt_cli.c)
add_executable(bhrt-bench   bhrt_bench.cpp)
add_executable(bhrt-test    bhrt_test.cpp)
add_executable(bhrt-diagram bhrt_diagram.cpp)
foreach(exe bhrt-cli bhrt-bench bhrt-test bhrt-diagram)
    target_link_libraries(${exe} PRIVATE bhrt_core)
endforeach()
# bhrt-cli is a C source but links the C++/CUDA core: link it as C++ so the
# C++ runtime (and CUDA runtime) are pulled in correctly.
set_target_properties(bhrt-cli PROPERTIES LINKER_LANGUAGE CXX)

# ---- Tests -----------------------------------------------------------------
enable_testing()
add_test(NAME unit COMMAND bhrt-test)

# ---- Optional pybind11 module ----------------------------------------------
if(BHRT_BUILD_PYBIND)
    find_package(pybind11 QUIET)
    if(pybind11_FOUND)
        pybind11_add_module(bhrt_trisect_core bindings.cpp)
        target_link_libraries(bhrt_trisect_core PRIVATE bhrt_core)
    else()
        message(WARNING "pybind11 not found; skipping python bindings")
    endif()
endif()
