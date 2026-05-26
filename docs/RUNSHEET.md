# Presentation runsheet — bhrt_trisect

## Before they arrive (2 min)
- Open a WSL terminal in the project: `cd ~/"Math rs"/bhrt_trisect`
- Clear the broken env var (this has broken every build): `unset CXX CC`
- Dry run so nothing surprises you live: `NOPAUSE=1 bash demo.sh`
  (should end with everything printed and no errors)
- Have open on the side: `docs/HANDOUT.md`, `docs/figures/cp2_diagram.svg`

## One-sentence opener
"It's a pipeline that takes a triangulated 4-manifold, finds a trisection of
it combinatorially, and computes the manifold's homology and intersection
form — implementing the Bell–Hass–Rubinstein–Tillmann / Spreer–Tillmann
approach, with a search that tries to lower the trisection genus."

## Live demo — just run `bash demo.sh` and narrate each section
| Section | Command (auto in demo.sh) | Say |
|---|---|---|
| Build | `make` | "Builds from clean with no external dependencies." |
| **Correctness** | `./build/bhrt-test` | "73 assertions, including ground-truth invariants of S⁴, ℂP², S²×S², ℂP²#ℂP², S¹×S³ — right signatures and intersection forms. This is the validation that the math is correct." |
| **Invariants** | `./build/bhrt-diagram --demo cp2` (and s2xs2, s1xs3) | "Give it any trisection diagram, it returns homology, intersection form, signature, parity. ℂP² → (1), σ=+1, odd; S²×S² → hyperbolic, σ=0, even." |
| **Custom input** | `./build/bhrt-diagram examples/cp2.tri` | "I can hand-write any diagram in a small text file and compute its invariants." |
| Triangulation side | `./build/bhrt-cli info` / `scan` | "The other end: skeleton stats and the ts-tricolouring scan on a triangulation." |

## Show the figures
Open `docs/figures/cp2_diagram.svg` and `s1xs3_diagram.svg`:
"These are the standard genus-1 diagrams the engine consumes — a torus with the
three cut systems. ℂP²'s curves are in general position (b₂=1, σ=+1); S¹×S³'s
are all parallel (b₁=1, b₂=0)."

## The pipeline (one slide of intuition)
triangulation → ts-tricolouring → trisection diagram → invariants, plus a
genus-reduction search. Point at `docs/PRESENTATION.md` for the full math.

## Be honest about scope (this earns trust)
- **Rigorous + validated today:** the Spreer–Tillmann ts-tricolouring detector,
  and the homology / intersection-form engine (checked against known manifolds).
- **In progress:** turning a raw triangulation into *genuine* cut systems, and
  the general FKSZ intersection form for non-simply-connected diagrams.
- **Optional backends:** Regina (lateral Pachner moves) and a CUDA GPU move
  scorer — both wired in behind build flags; CPU path does the actual math.

## If they ask about Regina / a specific manifold
Load an isosig live in the Regina GUI (or `regina-python`):
```python
t = Triangulation4.fromIsoSig("<isosig, no spaces>")
print(t.size(), t.isValid(), t.isClosed()); print(t.homology())
```
Then it can be pushed through the pipeline via `regina_bridge.py` → `.bhrt`.

## Likely questions + honest answers
- *"How do you get the intersection form?"* → "Exactly, via the linking matrix
  for simply-connected (g;0,0,0) diagrams; the general FKSZ case needs the
  geometric crossing data, which is the next step."
- *"Does the search reduce genus?"* → "The 1↔5 moves are in and the search
  explores; the lateral 2-4/3-3/4-2 moves come from the Regina backend."
- *"Is it GPU-accelerated?"* → "There's a CUDA scorer wired in behind a flag,
  bit-identical to the CPU path; it's a throughput option for large searches,
  not a result. CPU does the mathematics."
- *"Have you reproduced known trisection genera?"* → point to the validated
  invariants table; full census reproduction is future work.

## Golden rules
- Start the session with `unset CXX CC`.
- Demo from the **CPU build** (`demo.sh`) — it cannot crash on Regina/CUDA.
- Don't run `BHRT_USE_GPU=1` live (GPU path has an open driver/integration issue).
- Claim only what's validated; frame the rest as "the active frontier."
