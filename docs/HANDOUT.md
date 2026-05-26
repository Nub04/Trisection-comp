# `bhrt_trisect` — one-page overview

**What it is.** An end-to-end computational pipeline that takes a triangulated
closed oriented 4–manifold, produces a **trisection** of it combinatorially,
draws the **trisection diagram**, and computes the manifold's **homology and
intersection form** — with a Pachner-move search to lower the trisection genus.
It realizes the Bell–Hass–Rubinstein–Tillmann program via Spreer–Tillmann
colourings.

$$\text{triangulation} \to \text{ts-tricolouring} \to
(\Sigma_g;\alpha,\beta,\gamma) \to H_*(X),\,Q_X \qquad
\circlearrowleft\ \text{genus reduction}$$

---

### Rigorous and validated

- **ts-tricolouring detector** — exact Spreer–Tillmann test: every pentachoron
  type $(2,2,1)$, no monochromatic triangle, each $\Gamma_k$ connected, each
  $\gamma_{ij}$ collapses to a spine. A pass certifies a trisection.
- **Homology, exactly, from the three Lagrangians** $V_i\subset H_1(\Sigma)=\mathbb Z^{2g}$:
  $$H_1(X)=\mathbb Z^{2g}/(V_1{+}V_2{+}V_3),\quad
  \chi(X)=2+g-\textstyle\sum_i k_i,\ \ k_{xy}=2g-\operatorname{rank}[x;y],$$
  giving $b_2=g-\sum_i k_i+2b_1$, $\,b_3=b_1$, with correct torsion
  ($\operatorname{Tor}H_2\cong\operatorname{Tor}H_1$).
- **Intersection form** via the classical linking matrix $Q_X=MN^{-1}$ for
  simply-connected $(g;0,0,0)$-standard diagrams ($\gamma$ written in the
  $[\alpha;\beta]$ basis). Signature and parity follow.
- **Move executor**: $1\!\to\!5$ Pachner move (proved manifold-preserving:
  $\chi$ invariant) and a **self-verified** $5\!\to\!1$ inverse.

### Validation (asserted against known topology in the test suite)

| $X$ | $g$ | $b_1$ | $b_2$ | $b_3$ | $\sigma$ | $Q_X$ |
|---|---|---|---|---|---|---|
| $S^4$ | 0 | 0 | 0 | 0 | 0 | — |
| $S^1\times S^3$ | 1 | 1 | 0 | 1 | 0 | — |
| $\mathbb{CP}^2$ | 1 | 0 | 1 | 0 | $+1$ | $(1)$, odd |
| $S^2\times S^2$ | 2 | 0 | 2 | 0 | 0 | $\left[\begin{smallmatrix}0&1\\1&0\end{smallmatrix}\right]$, even |
| $\mathbb{CP}^2\#\mathbb{CP}^2$ | 2 | 0 | 2 | 0 | $+2$ | $I_2$, odd |

---

### Open frontier (state honestly)

1. **Triangulation $\to$ genuine cut systems.** Diagram extraction is currently
   heuristic; a validator now gates whether $\Sigma$ is a closed orientable
   surface and $\alpha,\beta,\gamma$ are honest cut systems.
2. **General intersection form.** The full Feller–Klug–Schirmer–Zemke
   computation needs the *geometric crossing data* of the curves on $\Sigma$,
   not just their homology classes — not yet recorded. We do the simply-connected
   case exactly and defer the rest.
3. **Lateral moves.** $2\text{-}4,\,3\text{-}3,\,4\text{-}2$ run only through the
   Regina backend; adding them to the built-in executor is what gives the search
   real genus-reduction power.
4. **Scalability.** Greedy collapsibility (sufficient, not complete) and
   brute-force colour enumeration ($3^{|V|-1}$) are the census-scale bottlenecks.

---

**Framing:** *"A working pipeline from triangulation to certified trisection to
homology and intersection form. The ts-detection and homology are rigorous and
validated against standard 4–manifolds; getting from a triangulation to genuine
cut systems, and the general FKSZ intersection form, are the active frontier."*

**Refs:** Gay–Kirby; Rubinstein–Tillmann; Bell–Hass–Rubinstein–Tillmann;
Spreer–Tillmann; Feller–Klug–Schirmer–Zemke (arXiv:1711.04762).
