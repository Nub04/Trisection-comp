// bhrt_trisect.hpp -- public C++ API for the BHRT trisection pipeline.
//
// This header declares a fully self-contained C++ implementation of
//
//     triangulation -> ts-tricolouring -> trisection diagram -> invariants
//
// plus a CPU beam-search engine that drives the genus-reduction loop.
// The implementation does NOT require Regina: every datum lives inside
// this translation unit's own simplicial-complex types, and the standard
// pipeline runs end-to-end through pure C++ (with optional CUDA scoring).
//
// Regina/Sage are still useful as optional adapters: the bindings layer
// (cpp/bindings.cpp) can import/export Regina triangulations when
// BHRT_HAS_REGINA is defined, and the Sage exact-algebra layer is wired
// into the invariants engine via Python orchestration.

#pragma once

#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace bhrt {

// ---------------------------------------------------------------------- //
// Permutations of {0,1,2,3,4}
// ---------------------------------------------------------------------- //

using Perm5 = std::array<std::uint8_t, 5>;

constexpr Perm5 perm_identity{0, 1, 2, 3, 4};

Perm5 perm_compose(const Perm5& p, const Perm5& q) noexcept;
Perm5 perm_inverse(const Perm5& p) noexcept;
int   perm_sign(const Perm5& p) noexcept;

// ---------------------------------------------------------------------- //
// Pentachoron + Triangulation
// ---------------------------------------------------------------------- //

struct Gluing {
    std::uint8_t  src_facet;
    std::int32_t  dst_pent;
    std::uint8_t  dst_facet;
    Perm5         perm;
};

struct Pentachoron {
    std::int32_t                                  index{0};
    std::array<std::optional<Gluing>, 5>          gluings{};
};

class Triangulation {
public:
    Triangulation() = default;
    explicit Triangulation(std::string label) : label_(std::move(label)) {}

    Pentachoron&       newPentachoron();
    void               addPentachora(std::size_t count);

    void               glue(std::int32_t src, std::uint8_t src_facet,
                            std::int32_t dst, const Perm5& perm);
    void               unglue(std::int32_t pent, std::uint8_t facet);

    // Sizes ------------------------------------------------------------
    std::size_t        size() const noexcept { return pentachora_.size(); }
    bool               isClosed() const noexcept;
    bool               isValid() const noexcept;

    // Skeleton (computed lazily) --------------------------------------
    std::int32_t       vertexOf(std::int32_t pent, std::uint8_t local) const;
    std::int32_t       edgeOf(std::int32_t pent,
                              std::uint8_t a, std::uint8_t b) const;
    std::int32_t       triangleOf(std::int32_t pent,
                                  std::uint8_t a,
                                  std::uint8_t b,
                                  std::uint8_t c) const;
    std::int32_t       tetraOf(std::int32_t pent, std::uint8_t facet) const;

    std::int32_t       nVertices() const;
    std::int32_t       nEdges() const;
    std::int32_t       nTriangles() const;
    std::int32_t       nTetrahedra() const;
    std::int32_t       eulerCharacteristic() const;
    std::vector<std::int32_t> edgeDegrees() const;

    // Canonical signature (BFS over (start_pent, start_perm) seeds, take
    // lexicographically smallest encoding). Matches the python reference.
    std::string        isoSig() const;
    std::string        edgeDegreeIsoSig() const;

    // Cheap content-addressed hash (NOT canonical).
    std::string        quickHash() const;

    // Access ----------------------------------------------------------
    const std::vector<Pentachoron>& pentachora() const noexcept { return pentachora_; }
    Pentachoron&                    pentachoron(std::size_t i) { return pentachora_[i]; }
    const Pentachoron&              pentachoron(std::size_t i) const { return pentachora_[i]; }
    const std::string&              label() const noexcept { return label_; }
    void                            setLabel(std::string l) { label_ = std::move(l); }

    Triangulation                   clone() const;

private:
    struct Skeleton;
    void  ensureSkeleton() const;

    std::vector<Pentachoron>          pentachora_;
    std::string                       label_;
    mutable bool                      dirty_ = true;
    mutable std::shared_ptr<Skeleton> skel_;
};

// Convenience builders --------------------------------------------------

Triangulation s4_minimal();           // 2-pentachoron S^4
Triangulation open_pentachoron();
Triangulation cp2_minimal();          // 4-pentachoron CP^2 (Kuhnel; placeholder)

// ---------------------------------------------------------------------- //
// ts-tricolouring (Spreer-Tillmann)
// ---------------------------------------------------------------------- //

struct TwoComplex {
    std::unordered_set<std::int64_t>                                          faces;
    std::unordered_set<std::int64_t>                                          edges;
    std::unordered_set<std::int32_t>                                          vertices;
    std::unordered_map<std::int64_t, std::unordered_set<std::int64_t>>        edge_to_faces;

    static std::int64_t packEdge(std::int32_t a, std::int32_t b) noexcept;
    static std::int64_t packTri(std::int32_t a, std::int32_t b,
                                std::int32_t c) noexcept;

    void addEdge(std::int32_t a, std::int32_t b);
    void addFace(std::int32_t a, std::int32_t b, std::int32_t c);
    bool isConnected() const;
    bool greedyCollapseTo1(int budget = 10000);
    std::size_t nFaces() const noexcept { return faces.size(); }
};

struct TSColouring {
    std::vector<std::int32_t>                                                 colour;
    bool                                                                      is_c{false};
    bool                                                                      is_ts{false};
    std::map<int, std::vector<std::set<std::int32_t>>>                        monochromatic_components;
    std::map<std::pair<int, int>, TwoComplex>                                 gamma_2complexes;
    std::vector<std::string>                                                  audit;
};

TSColouring isTsTricolouring(const Triangulation& T,
                              const std::vector<std::int32_t>& colour,
                              int collapse_budget = 10000);

std::vector<TSColouring> enumerateTsTricolourings(
    const Triangulation& T,
    int  collapse_budget = 10000,
    bool stop_at_first   = false);

#if BHRT_HAS_CUDA
// GPU-exhaustive variant (CUDA builds only; cpp/ts_color.cpp +
// cuda/ts_scan.cu): the 3^(n-1) candidate colourings are filtered on the
// GPU by the cheap prechecks ((2,2,1) type, monochromatic triangles,
// Gamma_k connectivity); only the survivors run the full certified CPU
// check, so the result set is IDENTICAL to enumerateTsTricolourings.
// Falls back to the CPU enumerator when no usable GPU is present or the
// input is outside the supported limits (3 <= n_vertices <= 32).
std::vector<TSColouring> enumerateTsTricolouringsGPU(
    const Triangulation& T,
    int  collapse_budget = 10000,
    bool stop_at_first   = false);
#endif  // BHRT_HAS_CUDA

struct ScanCounters {
    std::int64_t n_total                  = 0;
    std::int64_t n_pass_precheck          = 0;
    std::int64_t n_with_tricolouring      = 0;
    std::int64_t n_with_c_tricolouring    = 0;
    std::int64_t n_with_ts_tricolouring   = 0;
    std::int64_t n_ts_tricolourings_total = 0;
};

ScanCounters countAdmitTrisection(const std::vector<Triangulation>& tris,
                                   int collapse_budget = 10000);

// ---------------------------------------------------------------------- //
// Trisection diagram
// ---------------------------------------------------------------------- //

struct CellularSurface {
    std::int32_t                                       num_vertices{0};
    std::vector<std::pair<std::int32_t, std::int32_t>> edges;
    std::vector<std::vector<std::int32_t>>             faces;
    std::map<std::pair<std::int32_t, std::int32_t>, std::int32_t> edge_index;

    std::int32_t euler() const noexcept;
    std::int32_t genus() const noexcept;
};

struct CutCurve {
    std::vector<std::int32_t> edge_sequence;
    std::string               role;
};

struct TrisectionDiagram {
    CellularSurface          surface;
    std::vector<CutCurve>    alpha;
    std::vector<CutCurve>    beta;
    std::vector<CutCurve>    gamma;
    std::vector<std::string> audit;
    std::int32_t             genus() const noexcept { return surface.genus(); }
};

TrisectionDiagram extractDiagram(const Triangulation& T,
                                  const TSColouring&  colouring);

// ---------------------------------------------------------------------- //
// Invariants
// ---------------------------------------------------------------------- //

struct InvariantBundle {
    std::int32_t              h1_free_rank{0};
    std::vector<std::int64_t> h1_torsion;
    std::int32_t              h2_free_rank{0};
    std::vector<std::int64_t> h2_torsion;
    std::int32_t              h3_free_rank{0};
    std::vector<std::int64_t> h3_torsion;
    std::vector<std::vector<std::int64_t>> intersection_form;
    std::int32_t              signature{0};
    std::string               parity = "even";
    std::int32_t              genus{0};
    std::vector<std::string>  audit;
};

InvariantBundle computeInvariants(const TrisectionDiagram& diagram);

// ---------------------------------------------------------------------- //
// Homology / intersection form from the algebraic (Lagrangian) data of a
// trisection diagram: the homology classes of the alpha/beta/gamma cut
// curves in H_1(Sigma; Z) = Z^{2g} with the standard symplectic form
// omega((p|q),(r|s)) = p.s - q.r.  This is the mathematically precise
// input for the Feller-Klug-Schirmer-Zemke homology calculation and is
// what the ground-truth test suite exercises.
//
// Each of alpha/beta/gamma is a list of g rows, every row a length-2g
// integer vector; the three row-spaces must be Lagrangian (rank g,
// isotropic).  Use validateLagrangianDiagram() to check.
struct LagrangianDiagram {
    std::int32_t                           genus{0};
    std::vector<std::vector<std::int64_t>> alpha;
    std::vector<std::vector<std::int64_t>> beta;
    std::vector<std::vector<std::int64_t>> gamma;
};

struct DiagramValidity {
    bool                     ok{false};
    bool                     surface_closed{false};
    bool                     surface_orientable{false};
    bool                     surface_connected{false};
    bool                     alpha_is_cut_system{false};
    bool                     beta_is_cut_system{false};
    bool                     gamma_is_cut_system{false};
    std::int32_t             genus{0};
    std::vector<std::string> messages;
};

// Validate a Lagrangian diagram (rank g, isotropic, each pair a Heegaard
// diagram of #S^1xS^2). Pure linear algebra over Z.
DiagramValidity validateLagrangianDiagram(const LagrangianDiagram& d);

// H_*, intersection form, signature and parity of the closed oriented X.
// Homology is exact for any input; the integral intersection form (hence
// signature/parity) is computed exactly when X is simply connected and
// the diagram is (g;0,0,0)-standardisable (the classical linking-matrix
// case, which covers S^4, CP^2, S^2xS^2, connected sums, ...).  For other
// inputs the form is left empty with an audit note, since the full FKSZ
// signature needs the geometric crossing data of the curves on Sigma.
InvariantBundle computeInvariantsFromClasses(const LagrangianDiagram& d);

// Standard symplectic intersection pairing on H_1(Sigma;Z)=Z^{2g}.
std::int64_t symplecticOmega(const std::vector<std::int64_t>& u,
                             const std::vector<std::int64_t>& v);

// Validate a cellular trisection diagram: surface closed/orientable/
// connected and each curve family a genuine cut system of g curves.
DiagramValidity validateDiagram(const TrisectionDiagram& diagram);

// Smith normal form for an integer matrix (returns elementary divisors,
// padded with zeros to min(rows, cols)).
std::vector<std::int64_t> smithNormalFormDiagonal(
    const std::vector<std::vector<std::int64_t>>& M);

// Integer rank (number of non-zero elementary divisors).
std::int32_t integerRank(const std::vector<std::vector<std::int64_t>>& M);

// ---------------------------------------------------------------------- //
// CPU beam search (genus reduction)
// ---------------------------------------------------------------------- //

// 4D move catalogue. The enum values are historical (3D-flavoured names);
// the 4-dimensional semantics and pentachoron deltas are:
//
//   Pachner_1_5   1-5 about a pentachoron   (+4)  native executor
//   Pachner_5_1   5-1 about a vertex        (-4)  native executor
//   Pachner_2_3   2-4 about a tetrahedron   (+2)  Regina executor
//   Pachner_3_3   3-3 about a triangle      ( 0)  Regina executor
//   Pachner_4_4   4-4 about an edge         ( 0)  Regina executor
//   Pachner_4_2   4-2 about an edge         (-2)  Regina executor
//   EdgeCollapse  collapse an edge          (-deg(e)) Regina executor
enum class MoveType : std::uint8_t {
    Pachner_2_3 = 23,
    Pachner_3_3 = 33,
    Pachner_4_4 = 44,
    Pachner_4_2 = 42,
    Pachner_1_5 = 15,
    Pachner_5_1 = 51,
    EdgeCollapse = 99,
};

struct MoveCandidate {
    MoveType                  move;
    std::vector<std::int32_t> locator;
    int                       delta_pent_estimate;
};

std::vector<MoveCandidate> enumerateMoves(const Triangulation& T);

// Apply a move using the conservative built-in executor. Returns nullopt
// when the move's preconditions are not satisfied; the safe Regina path
// is in cpp/regina_bridge.cpp when BHRT_HAS_REGINA is defined.
std::optional<Triangulation> applyMove(const Triangulation& T,
                                        const MoveCandidate& c);

struct SearchConfig {
    int          beam_width         = 64;
    int          excess_height      = 4;
    double       time_limit_seconds = 600.0;
    int          max_iterations     = 1000;
    std::uint64_t seed              = 0;
    bool         use_gpu_scorer     = false;
    std::string  checkpoint_dir;
};

struct SearchState {
    Triangulation             triangulation;
    int                       pentachora{0};
    std::string               isosig;
    bool                      ts_feasible{false};
    std::optional<int>        best_genus;
    int                       excess_height{0};
    double                    score{0.0};
    std::vector<std::string>  move_history;
};

struct ParetoFront {
    std::vector<SearchState> items;
    void consider(const SearchState& s);
};

ParetoFront beamSearch(const Triangulation& start, const SearchConfig& cfg);

// ---------------------------------------------------------------------- //
// Flatten a candidate to a fixed-width record for the CUDA scorer.
// ---------------------------------------------------------------------- //

std::array<float, 5> flattenCandidate(const SearchState& parent,
                                       const MoveCandidate& c);

float scoreCandidateCPU(const std::array<float, 5>& record,
                         const std::array<float, 5>& weights);

// ---------------------------------------------------------------------- //
// .bhrt text format I/O
// ---------------------------------------------------------------------- //

std::vector<Triangulation> loadBhrtFile(const std::string& path);
void writeBhrtFile(const std::string& path,
                   const std::vector<Triangulation>& tris);

#if BHRT_HAS_REGINA
// ---------------------------------------------------------------------- //
// .esig I/O (Regina builds only; defined in cpp/tri_io.cpp).
// One Regina isomorphism signature per line, optional "<degrees>|" prefix
// on input, '#' comments skipped.
// ---------------------------------------------------------------------- //

std::vector<Triangulation> loadEsig(const std::string& path);
void writeEsig(const std::string& path,
               const std::vector<Triangulation>& tris);
#endif  // BHRT_HAS_REGINA

}  // namespace bhrt
