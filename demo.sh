#!/usr/bin/env bash
#
# demo.sh -- guided demonstration of the bhrt_trisect pipeline.
#
# Usage:
#   bash demo.sh              # interactive: pauses between sections
#   NOPAUSE=1 bash demo.sh    # run straight through, no pauses
#
# Uses only the CPU build (no Regina / CUDA needed) -- the validated core.

set -u
cd "$(dirname "$0")"

# The conda 'base' env often exports a non-existent CXX=g++-11; clear it so the
# build uses the system g++/gcc.
unset CXX CC

bold()  { printf '\n\033[1;36m========== %s ==========\033[0m\n' "$1"; }
note()  { printf '\033[0;33m%s\033[0m\n' "$1"; }
pause() { [ "${NOPAUSE:-0}" = "1" ] && return 0
          printf '\n\033[2m(press Enter to continue)\033[0m '; read -r _ || true; }

# ---------------------------------------------------------------------------
bold "0. Build (CPU, no external dependencies)"
note "Compiling the C++ core and four binaries with the available g++..."
make -s clean
if ! make -s; then
    echo "Build failed. Run 'unset CXX CC' and check your compiler, then retry."
    exit 1
fi
note "Build complete: build/bhrt-test, bhrt-cli, bhrt-bench, bhrt-diagram"
pause

# ---------------------------------------------------------------------------
bold "1. Correctness -- native test suite"
note "73 assertions, including GROUND-TRUTH invariants of S^4, CP^2, S^2xS^2,"
note "CP^2#CP^2 and S^1xS^3: signatures, intersection forms, homology."
./build/bhrt-test
pause

# ---------------------------------------------------------------------------
bold "2. Invariants of standard trisection diagrams"
note "The engine takes a trisection diagram and returns the 4-manifold's"
note "homology, intersection form, signature and parity."
for m in s4 cp2 s2xs2 s1xs3; do
    note "----- demo: $m -----"
    ./build/bhrt-diagram --demo "$m"
    echo
done
pause

# ---------------------------------------------------------------------------
bold "3. Custom trisection input"
note "Any diagram can be supplied as a small text file (genus + alpha/beta/gamma"
note "curve classes). Here is examples/cp2.tri:"
echo
cat examples/cp2.tri
echo
note "Computing its invariants:"
./build/bhrt-diagram examples/cp2.tri
pause

# ---------------------------------------------------------------------------
bold "4. Triangulation side -- the built-in S^4 example"
note "Skeleton statistics:"
./build/bhrt-cli info
echo
note "ts-tricolouring scan counters:"
./build/bhrt-cli scan
pause

# ---------------------------------------------------------------------------
bold "Summary"
note "Pipeline:  triangulation -> ts-tricolouring -> trisection diagram -> invariants"
note ""
note "Rigorous + validated today: the Spreer-Tillmann ts-tricolouring detector"
note "and the homology / intersection-form engine (checked against the manifolds"
note "above)."
note "In progress: turning a raw triangulation into genuine cut systems, and the"
note "general FKSZ intersection form for non-simply-connected diagrams."
note "Optional backends: Regina (lateral Pachner moves), CUDA (GPU move scorer)."
note ""
note "See docs/PRESENTATION.md for the full mathematical walkthrough."
