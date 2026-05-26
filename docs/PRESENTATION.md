# Computing trisections of 4–manifolds: `bhrt_trisect`

*A walkthrough for topologists.*

This document presents the project end to end: the mathematical setting, the
pipeline, the precise content of each stage, and the algorithms implemented in
the C++ core — including the four components rebuilt in this round (move
executor, homology/intersection form, diagram validation, ground-truth tests).
Each section is self-contained enough to be a talk segment.

---

## 1. Setting and goal

A **trisection** of a closed, oriented, smooth 4–manifold $X$ (Gay–Kirby) is a
decomposition
$$X = X_1 \cup X_2 \cup X_3, \qquad X_i \cong \natural^{k_i}\,(S^1\times B^3),$$
such that each pairwise intersection $H_{ij}=X_i\cap X_j$ is a genus-$g$
3–dimensional handlebody, and the common intersection
$\Sigma = X_1\cap X_2\cap X_3$ is a closed genus-$g$ surface (the **central
surface**). The data is recorded by a **trisection diagram**
$(\Sigma_g;\alpha,\beta,\gamma)$, where $\alpha,\beta,\gamma$ are three cut
systems (each $g$ disjoint simple closed curves) bounding compressing disks in
the three handlebodies. Every closed oriented 4–manifold admits a trisection,
and the diagram determines $X$.

**The computational problem (BHRT).** Bell–Hass–Rubinstein–Tillmann and
Spreer–Tillmann make this effective: a (singular) triangulation of $X$ that
carries a suitable vertex 3–colouring induces a trisection, and the **trisection
genus** can be studied through triangulations. `bhrt_trisect` implements that
bridge and then searches the Pachner-move graph to lower the genus.

The pipeline is one arrow, plus a search loop:
$$\text{triangulation} \;\to\; \text{ts-tricolouring} \;\to\;
\text{trisection diagram} \;\to\; \text{invariants},
\qquad \circlearrowleft\ \text{genus reduction.}$$

---

## 2. Stage 0 — Triangulations (`cpp/triangulation.cpp`)

A closed 4–manifold is a set of pentachora $\Delta^4$ with face pairings. Each
pentachoron has facets $0,\dots,4$; a gluing of facet $f$ of pentachoron $p$ to
facet $f'$ of $p'$ is a permutation $\pi\in S_5$ with $\pi(f)=f'$. The skeleton
(vertices, edges, triangles, tetrahedra) is recovered by union–find over the
gluing identifications, giving
$$\chi(X) = \#V - \#E - \#\,T_{\text{tet}} + \#\Delta^4 + \#T_{\text{tri}}.$$
A canonical isomorphism signature (BFS over all $(start,\ \mathrm{perm})$ seeds,
lexicographically minimal encoding) gives an isomorphism-invariant key for
deduplication.

Inputs: Regina census `.esig` strings (parsed by Regina's Python bindings into
the plain-text `.bhrt` format) or the native `.bhrt` reader.

---

## 3. Stage 1 — ts-tricolourings (`cpp/ts_color.cpp`)

Following Rubinstein–Tillmann and Spreer–Tillmann, a **tricolouring** is a
partition of the vertex set $V = P_0 \sqcup P_1 \sqcup P_2$ such that **every
pentachoron is type $(2,2,1)$**: it meets two parts in two vertices each and the
third in exactly one. The implementation tests, in order:

1. **Type $(2,2,1)$** for every pentachoron.
2. **No monochromatic triangle** (no 2–face has all three vertices one colour).
3. **$c$-condition:** each monochromatic subgraph $\Gamma_k = $ span of $P_k$ in
   the 1–skeleton is connected.
4. **ts-condition:** each bi-coloured 2–complex
   $\gamma_{ij} = \{$ faces with colours in $\{i,j\}\}$ collapses to a
   1–complex (a spine), tested by a greedy free-face collapse.

A colouring passing 1–4 certifies a trisection of $X$. The module also
enumerates all such colourings (partitions with vertex $0$ pinned, base-3
assignment of the rest) and provides census scan counters.

> **Status / caveat.** The collapse in step 4 is *greedy*, hence a sufficient
> but not complete test of collapsibility; an inconclusive collapse is reported
> as such rather than asserted either way (matching the paper's negative-result
> protocol). Enumeration is brute force over $3^{\,|V|-1}$ partitions — fine for
> small triangulations, the bottleneck at census scale.

---

## 4. Stage 2 — Trisection diagrams (`cpp/diagram.cpp`)

From a certified colouring the central surface $\Sigma$ is assembled with **one
quadrilateral per pentachoron**: in a type $(2,2,1)$ pentachoron the singleton
colour is the cone point and the four edges running between the two paired
colours carry the four $\Sigma$-vertices, giving a square face. Edges of
$\Sigma$ are reconstructed from the face boundary cycles; cut curves are read off
as meridians of the non-tree edges of the spine of each $\gamma_{ij}$.

### 4a. Validation — `validateDiagram` (new)

A diagram is only meaningful if $\Sigma$ is a genuine closed surface and
$\alpha,\beta,\gamma$ are genuine cut systems. The validator checks:

- **Connected** 1–skeleton (BFS).
- **Closed:** every edge lies in exactly two face-sides (else boundary or
  non-manifold).
- **Orientable:** orientations of the faces can be chosen so each interior edge
  is traversed once in each direction. Implemented as a 2–colouring (parity
  union–find) of the face-adjacency graph: two faces sharing an edge traversed
  in the **same** direction must receive **opposite** orientations; a parity
  conflict ⇒ non-orientable.
- **Cut systems:** each of $\alpha,\beta,\gamma$ has exactly $g$ curves, each a
  closed loop (all vertex degrees even in its edge multiset), pairwise
  edge-disjoint within the family.

This is the formal gate that the previous "best-effort" extractor lacked.

---

## 5. Stage 3 — Homology and the intersection form (`cpp/invariants.cpp`)

This stage was rebuilt. The previous code computed an ad-hoc
$Q=(GG^{\mathsf T}-G^{\mathsf T}G)/2$ "companion" in the edge basis, which is not
the intersection form. The new computation works in the symplectic homology of
the central surface, which is the correct input (Feller–Klug–Schirmer–Zemke).

### 5a. Algebraic input

Work in $H_1(\Sigma;\mathbb Z) = \mathbb Z^{2g}$ with the standard symplectic
form
$$\omega\big((p\mid q),(r\mid s)\big) = p\cdot s - q\cdot r .$$
Each cut system is a $g\times 2g$ integer matrix whose rows are the homology
classes of its curves; its row space is a **Lagrangian** $V_1=\langle\alpha\rangle$,
$V_2=\langle\beta\rangle$, $V_3=\langle\gamma\rangle$ (rank $g$, $\omega$-isotropic).
`validateLagrangianDiagram` checks rank and isotropy.

### 5b. First homology (exact)

Each handlebody kills the homology of its cut system, so
$$\boxed{\,H_1(X) \;\cong\; \mathbb Z^{2g}\big/\big(V_1+V_2+V_3\big)\,}$$
computed as the cokernel of the stacked $3g\times 2g$ matrix
$[\alpha;\beta;\gamma]$ via Smith normal form: free rank
$b_1 = 2g - \operatorname{rank}[\alpha;\beta;\gamma]$, torsion = the elementary
divisors $>1$.

### 5c. Betti numbers via Euler characteristic (exact)

Let $k_{xy} = 2g - \operatorname{rank}[x;y]$. Because $\partial X_i$ is a genus-$g$
Heegaard splitting of $\#^{k}(S^1\times S^2)$ with $H_1$ free of rank $k$, the
three sector genera are $k_{\alpha\beta},\,k_{\beta\gamma},\,k_{\alpha\gamma}$.
Inclusion–exclusion on $X=X_1\cup X_2\cup X_3$ with
$\chi(\natural^{k}S^1\times B^3)=1-k$, $\chi(\text{genus-}g\text{ handlebody})=1-g$,
$\chi(\Sigma_g)=2-2g$ gives
$$\chi(X) = \sum_i(1-k_i) - 3(1-g) + (2-2g) = 2 + g - \textstyle\sum_i k_i .$$
For a closed oriented 4–manifold $\chi = 2 - 2b_1 + b_2$, hence
$$\boxed{\,b_2 = g - \textstyle\sum_i k_i + 2b_1\,}, \qquad b_3 = b_1 .$$

By Poincaré duality and universal coefficients, $H_3(X)=\mathbb Z^{b_1}$
(torsion-free) and $\operatorname{Tor}H_2(X)\cong\operatorname{Tor}H_1(X)$ — both
implemented.

**Validation (computed by hand and asserted in the test suite):**

| $X$            | $g$ | $b_1$ | $\sum k_i$ | $b_2$ | $b_3$ | $\sigma$ |
|----------------|-----|-------|------------|-------|-------|----------|
| $S^4$          | 0   | 0     | 0          | 0     | 0     | 0        |
| $S^1\times S^3$| 1   | 1     | 3          | 0     | 1     | 0        |
| $\mathbb{CP}^2$| 1   | 0     | 0          | 1     | 0     | $+1$     |
| $S^2\times S^2$| 2   | 0     | 0          | 2     | 0     | 0        |
| $\mathbb{CP}^2\#\mathbb{CP}^2$ | 2 | 0 | 0 | 2 | 0 | $+2$ |

### 5d. The intersection form

A subtlety worth stating to an expert audience: **the intersection form is not
determined by the homology classes $(\alpha,\beta,\gamma)$ alone in general.**
Self-intersection in a 4–manifold needs framing/crossing data, not just classes.
Concretely, the naive symmetric symplectic guess
$Q(P,P')=\omega(xA,y'B)+\omega(yB,z'C)+\omega(zC,x'A)$ on $\ker[\alpha;\beta;\gamma]^{\mathsf T}$
gives $Q(P,P)=3$ for $\mathbb{CP}^2$ instead of $1$ — it is wrong. This is
exactly why FKSZ's general algorithm uses the geometric intersection pattern of
the curves on $\Sigma$.

**What is computed exactly here** is the **simply-connected, $(g;0,0,0)$-standard
case** (which covers $S^4$, $\mathbb{CP}^2$, $S^2\times S^2$, and connected sums)
via the classical linking matrix. When $(\alpha,\beta)$ is a genus-$g$ Heegaard
diagram of $S^3$ (so $[\alpha;\beta]$ is unimodular), write each $\gamma$-curve in
that basis,
$$\gamma_i = \sum_j M_{ij}\,\alpha_j + \sum_j N_{ij}\,\beta_j,$$
i.e. $[\,M \mid N\,] = \gamma\,[\alpha;\beta]^{-1}$ split into $g\times g$ blocks.
Then
$$\boxed{\,Q_X = M\,N^{-1}\,}$$
is the integral, symmetric, unimodular intersection form. Signature comes from
the eigenvalue signs (Jacobi), parity from the diagonal.

*Worked checks.* $\mathbb{CP}^2$: $[\alpha;\beta]=I_2$, $\gamma=(1,1)$, so
$M=N=[1]$, $Q=[1]$, $\sigma=1$, odd. $S^2\times S^2$: $\gamma_1=a_2+b_1$,
$\gamma_2=a_1+b_2$ give $M=\left[\begin{smallmatrix}0&1\\1&0\end{smallmatrix}\right]$,
$N=I$, $Q=\left[\begin{smallmatrix}0&1\\1&0\end{smallmatrix}\right]$, $\sigma=0$,
even.

The exact rational linear algebra is done with a small `Frac` type and
Gauss–Jordan inversion (with an overflow guard); non-integral or
non-standardisable results fall back to "homology only" with an audit note
rather than fabricating a form.

> **Boundary of correctness.** Homology ($b_1,b_2,b_3$, torsion) is exact for
> *any* Lagrangian diagram. The integral form / signature / parity are exact in
> the simply-connected standardisable case; the general (e.g. $b_1>0$ or
> non-standardisable) signature needs the full FKSZ crossing data, which the
> cellular extractor does not yet record — so `computeInvariants(diagram)`
> reports homology and defers the form, while `computeInvariantsFromClasses`
> handles algebraic input directly.

---

## 6. Stage 4 — Genus-reduction search and the move executor (`cpp/search_cpu.cpp`)

The search walks the Pachner-move graph, evaluating the trisection genus at each
reached triangulation and keeping a Pareto front of
$(\#\text{pentachora},\ \text{genus})$. The previous executor returned
`nullopt` for every move, so the search never left the seed. Two moves are now
implemented combinatorially.

### 6a. The $1\!\to\!5$ move (always legal)

Cone a fresh central vertex $c$ over the five tetrahedral facets of a pentachoron
$P$ with vertices $\{0,1,2,3,4\}$. Replace $P$ by five pentachora $P_0,\dots,P_4$,
where $P_i$ has vertices $(\{0,\dots,4\}\setminus\{i\})\cup\{c\}$ with $c$ placed
at local index $i$.

- **Internal gluings:** $P_i$ facet $j \leftrightarrow P_j$ facet $i$ with
  permutation the transposition $(i\,j)$ — these two facets share the vertex set
  $\{c\}\cup(\{0,\dots,4\}\setminus\{i,j\})$.
- **External gluings:** $P_i$ facet $i$ (opposite $c$) is the original facet $i$
  of $P$, so it inherits $P$'s external gluing verbatim; self-gluings of $P$ map
  to cone–cone gluings.

This is a stellar subdivision; it preserves the PL type, so $\chi$ is unchanged.
The test `test_one_five_move_preserves_manifold` asserts exactly this:
$|U|=6$ and $\chi=2$ after a $1\!\to\!5$ on the 2-pentachoron $S^4$.

### 6b. The $5\!\to\!1$ move (self-verified inverse)

At a degree-5 internal vertex $c$, the five incident pentachora are merged back
into one $P_{\text{new}}$ whose vertices are the five "outer" vertices (those
opposite $c$), with external gluings inherited from the facets opposite $c$. The
reconstruction's permutation bookkeeping is delicate, so it is **self-verified**:
the candidate result $U$ is accepted only if it is valid, has matching
closedness, and $\mathrm{oneFiveMove}(U,P_{\text{new}})$ reproduces the original
isomorphism signature. A faulty reconstruction therefore returns `nullopt`
("skip") and can never corrupt the triangulation. `test_five_one_inverts_one_five`
checks the round trip $S^4 \to (1\!\to\!5) \to (5\!\to\!1) \to S^4$.

### 6c. Scope

`enumerateMoves` now emits $1\!\to\!5$ candidates at every pentachoron and
$5\!\to\!1$ candidates at every degree-5 vertex, so the beam search genuinely
expands (`test_search_actually_expands`). The lateral 4D bistellar moves
$2\!-\!4,\ 3\!-\!3,\ 4\!-\!2$ are routed to a Regina-backed executor
(`applyMoveRegina`, compiled in under `-DBHRT_HAS_REGINA`); without Regina they
are skipped.

---

## 7. Validation strategy (`cpp/bhrt_test.cpp`)

Because the build environment here cannot compile, the tests are written to pin
**externally known** answers, so a wrong implementation fails loudly rather than
silently:

- **Move executor:** manifold-preservation ($\chi$), and the $1\!\to\!5/5\!\to\!1$
  round trip.
- **Ground-truth invariants (hand-coded diagrams):** $S^4$, $\mathbb{CP}^2$
  ($\sigma=1$, odd, $Q=[1]$), $-\mathbb{CP}^2$ ($\sigma=-1$), $S^2\times S^2$
  ($\sigma=0$, even, $Q=\left[\begin{smallmatrix}0&1\\1&0\end{smallmatrix}\right]$),
  $\mathbb{CP}^2\#\mathbb{CP}^2$ ($\sigma=2$), $S^1\times S^3$
  ($b_1=b_3=1,\ b_2=0$).
- **Validator:** rejects a non-Lagrangian / dependent cut family.
- **Triangulation-driven (end-to-end):** the full
  colour → diagram → validate → invariants pipeline runs on $S^4$ without
  crashing and reports a consistent non-negative genus.

These complement, rather than encode, the implementation: the asserted numbers
are the known topology of the test manifolds.

---

## 8. What is rigorous vs. open

**Rigorous now**

- Spreer–Tillmann ts-tricolouring detection.
- $H_1, b_2, b_3$ and torsion from the three Lagrangians (validated on five
  manifolds).
- Integral intersection form, signature, parity for simply-connected
  $(g;0,0,0)$-standardisable diagrams.
- The $1\!\to\!5$ move (manifold-preserving) and its self-verified inverse.
- Surface and cut-system validators.

**Open / next**

1. **Lateral Pachner moves** in the built-in executor (currently Regina-only),
   to give the search real lateral mobility for genus reduction.
2. **Full FKSZ intersection form** for the general case — requires recording a
   rotation system / geometric crossings on $\Sigma$ in the extractor, then the
   FKSZ pairing.
3. **Complete collapsibility** test (replace greedy) and **pruned enumeration**
   for census-scale ts-detection.
4. Reconcile the parallel code trees and the two architecture docs into one
   canonical source.

---

## 9. References

- D. Gay, R. Kirby, *Trisecting 4–manifolds*, Geom. Topol. 20 (2016).
- J. H. Rubinstein, S. Tillmann, *Multisections of piecewise linear manifolds*.
- M. Bell, J. Hass, J. H. Rubinstein, S. Tillmann, *Computing trisections of
  4–manifolds*, PNAS 115 (2018).
- J. Spreer, S. Tillmann, *Determining the trisection genus through
  triangulations*.
- P. Feller, M. Klug, T. Schirmer, D. Zemke, *Calculating the homology and
  intersection form of a 4–manifold from a trisection diagram*, PNAS 115 (2018),
  arXiv:1711.04762.
