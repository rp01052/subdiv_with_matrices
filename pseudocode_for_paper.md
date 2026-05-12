# Pseudo-code for paper — extended context

This document accompanies `pseudocode_for_paper_trimmed.md`. The trimmed file
contains paper-ready pseudo-code; this file adds the surrounding context,
motivation, and implementation notes that would otherwise need to be
re-derived. It is intended as background material for the LaTeX translation
of the trimmed file and as a reading reference, not for inclusion in the
paper body.

Sources: `src/doerfler.{h,cpp}`, `src/activeDoFs.{h,cpp}`,
`src/residual_estimator.h`, `adaptive_vector_laplacian/{main,solvers,
residual_estimator}.h`, weak forms in `weak_forms/laplace/`.

---

## 0. What makes the adaptive scheme non-standard

A few features distinguish the adaptive procedure here from textbook AFEM
on standard finite elements, and motivate the algorithms that follow.

1. **Subdivision basis functions are not element-local.** A subdivision
   basis function on level `l` has its support spread over more than one
   incident face: under repeated subdivision its limit support extends one
   or more rings of fine-mesh faces beyond the simplex it nominally lives
   on. Consequently, the set of coarse-level basis functions whose support
   intersects a finer region — the **interface DoFs** at level boundaries
   — must be enumerated explicitly. In standard FE this set is empty
   because basis functions are element-local; here it is large enough to
   need its own construction (Algorithm 5).

2. **The active subdivision space is a subspace of a single
   embedding-level FE space.** Bilinear forms are assembled once on the
   finest mesh `T_L` against standard FE spaces (CG1 for 0-forms, N1curl
   for 1-forms). The active subdivision space at any adaptive step is
   then expressed by a sparse matrix `A^k` whose columns are level-`L`
   coefficient expansions of the active basis functions, obtained by
   composing per-level subdivision operators (Algorithm 6).

3. **Refinement is recorded as a partition, not a tree of cells.** The
   *subdivision partition* `P` records the level at which each region of
   the domain is "resolved": `f ∈ P[l]` means face `f` on level `l` is
   active. A parent and its children are never both active. This is
   different from many `hp`-AFEM implementations where the active leaves
   of a refinement tree may coexist with their parents in the basis.

4. **Local enrichment can break exactness of the discrete de Rham
   complex.** Standard FEEC adaptivity preserves exactness automatically
   because element-by-element refinement keeps the differential operators
   well-defined. With subdivision elements, the discrete gradient and
   curl–curl couplings `B`, `C` are obtained by Galerkin restriction of
   the level-`L` operators to the active space via `A^k`. Certain local
   refinement patterns produce active 1-form spaces that are *no longer
   the image of any active 0-form space*, opening a hole in the exact
   sequence. Algorithm 4 patches the known geometric origin of this
   defect at domain-boundary vertices.

5. **The marked-vs-refined relation is not one-to-one.** A marked face is
   replaced by its four children only after a padding step that enforces
   an admissibility condition on the resulting patches (Algorithm 3).
   Padding can promote additional faces to the marked set on coarser
   levels.

## Notation

- `P` — *subdivision partition*: a per-level set of active faces, stored as
  `P[l] ⊂ F_l` for `l = 0..L`. Faces in `P[l]` are not present in `P[l+1]` (no
  parent/child duplication). Implemented as `DoFs` (`vector<unordered_set<int>>`).
- `D[k][l]` — active `k`-form DoFs on level `l` (`k = 0,1,2`), with
  `D[2] = P`. Implemented as `KFormDoFs`.
- `T_L`, `F_L`, `E_L`, `V_L` — the embedding mesh at the finest level and
  its faces, edges, vertices.
- For a face `f` on level `l`, its four children on level `l+1` are
  `{4f, 4f+1, 4f+2, 4f+3}`; for an edge `e`, the two children are
  `{2e, 2e+1}`. The integer-division ancestor walk `f → f/4`, `e → e/2`
  is therefore the parent map.
- `N_e(f)` — set of edge-neighbours of face `f` on the same level.
- `vert_faces(v), vert_edges(v)` — fan of faces / edges incident to a
  vertex `v` on its own level.
- `S^k_l` — per-level subdivision matrix for `k`-forms (level `l` → level
  `l + 1`). `S^k_{l → L}` is the accumulated composition
  `S^k_{L−1} · S^k_{L−2} · … · S^k_l` mapping level-`l` coefficients to
  level-`L` coefficients.

---

## 1. Outer adaptive loop

The driver follows the standard `SOLVE → ESTIMATE → MARK → REFINE` shape of
AFEM. What differs from a textbook scheme:

- The `SOLVE` step takes place in the active subdivision space, which is a
  *non-nested* subspace of the level-`L` embedding FE space — non-nested in
  the sense that two adaptive steps `n` and `n+1` produce active spaces
  that are subspaces of the same level-`L` FE space but generally do not
  embed into one another. The previous step's solution can still be tested
  against the new active space via the level-`L` representation (see
  `checkNestedness` calls in `main.cpp:626-628`).
- The `REFINE` step is decomposed into three sub-steps: Dörfler marking
  (§2), admissibility padding (§3), and topological repair at the boundary
  (§4). All three operate on the partition `P`; only at the end of the
  step are the active DoFs reconstructed.

```text
Input: initial partition P, max steps S
for step = 0, 1, ..., S-1:
    D                 ← BuildActiveDoFs(P)              # §5
    (A0, A1)          ← BuildKFormSubdivMatrices(D)     # §6
    (M0, C0, B0)      ← (A0ᵀ M A0, A1ᵀ C A1, A0ᵀ B A1)
    (σ, u)            ← MINRES on K(σ, u) = (0, A1ᵀ rhs) using §7
    η_f               ← residual estimators per face f ∈ P    # §8
    P'                ← DoerflerMarking(P, η, θ)        # §2
    P'                ← FixIsolatedBoundaryVertices(P') # §4
    P                 ← P'
```

A few implementation comments:

- The right-hand side is assembled against the level-`L` N1curl space and
  must be projected to the active space by `A1ᵀ`. A swap matrix
  `P_swap`, recording the DoF re-ordering between SEC's enumeration and
  DOLFIN's, is applied between assembly and projection (see CLAUDE.md).
- `M, C, B` are the level-`L` mass, curl-curl, and gradient-coupling
  matrices respectively. The Galerkin restriction
  `(A0ᵀ M A0, A1ᵀ C A1, A0ᵀ B A1)` does not require re-assembly — the
  embedding-level operators are assembled once before the loop.
- `EnforceGrading` is defined in `src/doerfler.cpp` but is **not called**
  by the adaptive driver. The one-level admissibility behaviour is
  delivered instead by `PadMarkedFaces` (Algorithm 3). The function is
  flagged at the bottom as currently dead code.

---

## 2. Dörfler marking

The marking step is a standard bulk-chasing greedy selection: order the
indicators in decreasing magnitude and accumulate until the squared sum
reaches `θ` times the total.

One detail is forced by the subdivision setting: the partition `P` admits a
maximum level `L`, and faces already on `P[L]` cannot be refined further
(refining them would require a level `L+1` mesh that the hierarchy does not
provide). They are therefore **excluded from the threshold total**
(`doerfler.cpp:396-398`). The level-`L` η² is still reported and is still
part of the global error estimate; it simply cannot be reduced by the
present round of marking. This means that as the partition matures, the
effective Dörfler parameter (fraction of *refineable* `η²` selected) is
slightly higher than the user-set `θ` would suggest if one accumulated
against the full total — the difference is small whenever the level-`L`
share of the total error is small.

```text
Algorithm: DoerflerMarking(P, η, θ)

Input:  partition P, estimators η[l][f], parameter θ ∈ (0, 1]
Output: refined partition P_new

1.  total ← Σ { η[l][f]²  :  f ∈ P[l],  l < L }            # exclude l = L
    if total = 0: return P

2.  sort faces { (l, f) : f ∈ P[l], l < L } by descending η[l][f]²

3.  accumulated ← 0;  marked[l] ← ∅ for all l
    for each (l, f) in sorted order while accumulated < θ · total:
        marked[l] ← marked[l] ∪ {f}
        accumulated ← accumulated + η[l][f]²

4.  if padForAdmissibility:
        PadMarkedFaces(marked, P)                         # §3

5.  P_new[l] ← ∅
    for each f ∈ P[l]:
        if f ∈ marked[l]: insert {4f, 4f+1, 4f+2, 4f+3} into P_new[l+1]
        else:             insert f into P_new[l]

6.  if padForAdmissibility and not ValidateSubdivisionPartition(P_new):
        # fallback: every patch failing Def. 38 has its faces' coarse-level
        # parents identified, and the edge-neighbours of those parents are
        # added to `marked` on the coarser level, then step 5 is redone.
        ... (see doerfler.cpp:504-583)
```

Note: The code uses a priority-queue formulation with an α-boost on edge
neighbours of selected faces (a "cluster-aware" bias). Experiments showed it
does not change the marked-vs-padded ratio meaningfully, so it is omitted
above. With `α = 1` (default off) the priority-queue version reduces to a
plain greedy sort.

Flag: The fallback in step 6 (`doerfler.cpp:504-583`) reports the partition as
invalid via `ValidateSubdivisionPartition`, which checks **both** Def. 38
conditions (I) and (II). However, `PadMarkedFaces` only enforces (I) — the
condition-(II) check is commented out (lines 265-284). It is therefore
possible for the fallback to fire because the patch satisfies (I) but not (II);
I have not verified whether this actually occurs in practice. Worth confirming.

---

## 3. Padding for subdivision-partition admissibility

The subdivision-partition condition (Def. 38, stated in the paper text)
requires every face in a level-`l` patch to have at least one *interior
vertex* (a vertex all of whose adjacent faces lie in the patch) or to be
edge-adjacent to an interior face (a face all of whose edges are shared
with other faces in the patch). The condition is what guarantees that the
local active basis on each patch is the full subdivision space on that
patch — without it, the patch is too thin (e.g. a single strip of faces)
to support all of the level-`l` basis functions that would nominally live
on it, and the active 2-form DoFs become linearly dependent on finer
DoFs after one or two further subdivisions.

In practice the partition produced by raw Dörfler marking can violate the
condition in two ways: a marked face's children form a single isolated
patch on `P[l+1]` whose neighbours on level `l` are all still coarse, or
patches are too elongated to satisfy the interior-vertex requirement.
`PadMarkedFaces` repairs this by iteratively marking additional ancestors
until every face of every patch in the resulting partition has an
interior vertex.

The local repair rule is: for an offending face `f` on level `l`, choose
the vertex `v` of `f` whose fan is *cheapest* to complete (fewest missing
faces) and mark the missing faces' currently-active ancestors. Walking up
to the active ancestor is necessary because a missing face may already be
deeper in the hierarchy — marking it directly would have no effect
because the marked set only records faces that are currently in `P`. The
walk uses the SEC parent map `f → f/4`.

```text
Algorithm: PadMarkedFaces(marked, P)

Input/Output: marked[l] (in/out), partition P
Constant:     MAX_ITER = 50

repeat up to MAX_ITER times:
    # Build the tentative partition that would result from `marked`
    tent[l] ← ∅ for all l
    for each f ∈ P[l]:
        if f ∈ marked[l] and l < L: insert children into tent[l+1]
        else:                       insert f into tent[l]

    new_marks[l] ← ∅ for all l
    for each level l with tent[l] ≠ ∅:
        # Interior vertices: every adjacent face is in tent[l]
        I_v ← { v : ∀ g ∈ vert_faces(v), g ∈ tent[l] }

        for each f ∈ tent[l]:
            if any vertex of f is in I_v: continue          # cond. (I) holds

            # cond. (I) fails — find the cheapest vertex v of f and try to
            # complete the fan around it by marking missing faces' ancestors
            v ← arg min_{v ∈ verts(f)} #{ g ∈ vert_faces(v) : g ∉ tent[l] }

            fixable ← false
            for each g ∈ vert_faces(v) with g ∉ tent[l]:
                # walk up: find ancestor of g that is currently in P
                a ← g;  for a_l = l down to 0:
                    if a ∈ P[a_l]:
                        if a_l < L and a ∉ marked[a_l]:
                            new_marks[a_l] ← new_marks[a_l] ∪ {a}
                            fixable ← true
                        break
                    a ← a / 4
            if not fixable:
                # fall back: try to mark the ancestor of f itself
                (same walk-up applied to f)

    if new_marks is empty: break
    marked ← marked ∪ new_marks
```

Flag: The condition-(II) check exists in the source as a commented-out block
(lines 265-284 of `doerfler.cpp`); only condition (I) is currently enforced
here. Worth flagging in the paper if Def. 38 is stated with both conditions.

Flag: The inner loop walks ancestors up to level 0 by integer division `a/4`
(faces) — this exploits the SEC convention `children(f) = {4f, ..., 4f+3}`. If
the paper uses different parent/child notation this should be matched up.

---

## 4. Restoring exactness at isolated boundary vertices

This step has no analogue in standard AFEM and is specific to FEEC over
subdivision spaces.

After marking and padding the partition is admissible in the
patch-interior sense of §3, but a different failure mode can appear at
the *domain boundary*. Consider a boundary vertex `v` of the level-`l`
mesh with two boundary edges `e1, e2` and a fan of interior faces between
them. If both faces adjacent to `e1` and `e2` along the boundary remain
coarse while some face in the interior of the fan has been refined, then
the local 1-form basis around `v` admits a curl-free element supported on
the boundary-adjacent edges that is **not** the gradient of any active
0-form — geometrically, the refined patch separates the two boundary
edges and the discrete `∇` cannot connect them. The discrete de Rham
sequence

```math
0 → V^0 ─ ∇ → V^1 ─ curl → V^2 → 0
```

then ceases to be exact at `V^1`, and the Hodge-Laplacian block system
`K` becomes rank-deficient on its (2,2) block.

The fix is geometric: walk the fan around `v` from one boundary face
`bF1` to the other `bF2`, find the side along which the refined region is
*closer*, and refine all coarse faces on that side. This guarantees that
the refined patch reaches the boundary on at least one side and the
topological loop closes. The "shorter side" rule keeps the additional
refinement minimal.

The boundary fan around a vertex `v` with exactly two boundary edges is
well-defined as an ordered list `bF1 → … → bF2` because the interior
edges incident to `v` give a unique walk in the dual graph from one
boundary face to the other (`traverseFan`, lines 19-39 of
`doerfler.cpp`).

```text
Algorithm: FixIsolatedBoundaryVertices(P)

Output: number of new child faces inserted

count ← 0
for l = 0 .. L:
    for each f ∈ P[l]:
        for each boundary vertex v of f:
            (e1, e2)        ← the two boundary edges incident to v
            (bF1, bF2)      ← unique adjacent faces of e1, e2
            if bF1 ∉ P[l] or bF2 ∉ P[l]: continue          # already fixed

            if no face in vert_faces(v) is missing from P[l]: continue
                                                            # not isolated

            fan ← traverseFan(v, bF1, bF2)                  # ordered list
            if fan invalid or does not terminate at bF2: continue

            leftDist  ← index of first refined face from the bF1 end
            rightDist ← index of first refined face from the bF2 end

            if leftDist ≤ rightDist:
                toMark ← fan[0 .. leftDist - 1]
            else:
                toMark ← fan[|fan| - rightDist .. |fan| - 1]

            if l + 1 > L: continue                          # cannot refine

            for each f' ∈ toMark:
                P[l]   ← P[l]   \ {f'}
                P[l+1] ← P[l+1] ∪ {4f', 4f'+1, 4f'+2, 4f'+3}
                count ← count + insertions
return count
```

`traverseFan(v, start, end)` (lines 19-39 of `doerfler.cpp`) walks the fan of
faces around `v` from `start` to `end` by hopping across interior edges
incident to `v`; on a boundary vertex with exactly two boundary edges this is
well-defined.

Flag: I describe the *motivation* (de Rham exactness) from the comment block
in `doerfler.h:78-86`. The code itself only implements the geometric fix; the
analytical claim is not separately verified here.

---

## 5. Building the active-DoF set from a subdivision partition

The active DoFs are needed to assemble the Galerkin projection of the
level-`L` system to the active space (§6). Each active DoF corresponds to
one subdivision basis function whose level-`L` coefficient expansion is a
column of `A^k`.

There are two kinds of contributions per level:

- **Native DoFs**: vertices and edges of every active face on that level.
  These are the basis functions whose nominal home simplex lies in `P[l]`.
- **Interface DoFs**: extra basis functions at the boundary between
  level-`l` active faces and a level-`(l+1)`-or-deeper refined region.
  These are subdivision basis functions whose support is *split* by the
  level boundary — on the coarse side they are level-`l` basis functions
  carried into the refined region by repeated subdivision, and on the
  fine side they are level-`(l+1)` basis functions whose support reaches
  back into the coarse region.

If the interface DoFs were not added, the active basis would fail to
represent functions that are non-zero across the level boundary, and
projection of the level-`L` operators to the active space would have a
non-trivial kernel.

The interface DoFs are enumerated by a single sweep that maintains a
**frontier** — the set of edges that separate a resolved (coarse-side)
region from an unresolved (fine-side) region. The frontier is split each
level, and the new interface vertices it produces (one *edge-vertex* at
each split midpoint, plus *vertex-vertices* at the endpoints) emit the
DoFs introduced by their level. A frontier edge is *consumed* when the
fine side it borders has been resolved by the partition; this terminates
propagation downward.

The five sub-steps in Algorithm 5 (consume → spawn → emit → split →
collect endpoints) execute in a fixed order. Two ordering invariants
matter:
- Spawning new frontier edges on level `l` must happen **before** vertex
  emission for level `l`, because emission writes to `D[2][l]` and would
  otherwise turn the level-`l` view used by `FindInterfaceEdgesOnLevel`
  into something inconsistent.
- An edge can only enter the frontier as either a *newly spawned* edge or
  a *split child* of a level-`(l−1)` frontier edge. Edges that were just
  consumed on level `l` must be excluded from spawning, because they
  border a resolved coarser region (not a finer one) and so are not
  interface edges in the sense the sweep cares about.

### 5.1 Native DoFs

```text
Algorithm: BuildNativeDoFs(P)

D[0][l] ← ∅;  D[1][l] ← ∅;  D[2][l] ← P[l]  for all l
for l = 0 .. L:
    for each f ∈ P[l]:
        D[0][l] ← D[0][l] ∪ verts(f)
        D[1][l] ← D[1][l] ∪ edges(f)
return D
```

### 5.2 Interface DoF frontier sweep

```text
Algorithm: GetInterfaceDoFs(D)

frontier             ← ∅   # set of InterfaceEdge
verticesToProcess    ← ∅   # spawned by edges split on level l-1

for l = 0 .. L:
    # (1) consume frontier edges whose fine (inactive) side is now resolved
    consumed ← ∅
    for each e ∈ frontier:
        if all faces on e's inactive side are absent from P[l]:
            # this edge no longer borders an unresolved finer region
            (frontier retains e)
        else if all faces on e's inactive side are now present in P[l]:
            consumed ← consumed ∪ {index(e)};  remove e from frontier

    # (2) spawn new frontier edges from the active faces on level l
    #     (must run BEFORE step 3, otherwise process() pollutes P[l] view)
    new ← FindInterfaceEdgesOnLevel(l, D)        # see helper below
    new ← new \ { e : index(e) ∈ consumed }      # those border coarser, not finer
    frontier ← frontier ∪ new

    # (3) emit interface DoFs from vertices spawned on level l - 1
    for each v ∈ verticesToProcess:
        D[0][l] ← D[0][l] ∪ {index(v)}
        D[1][l] ← D[1][l] ∪ activeEdges(v)
        D[2][l] ← D[2][l] ∪ activeFaces(v)
    verticesToProcess ← ∅

    # (4) split frontier edges into 2 child edges + 1 midpoint vertex
    if l < L:
        refined ← ∅;  edgeVerts ← ∅
        for each e ∈ frontier:
            (e1, e2, midVertex) ← SplitInterfaceEdge(e)   # §5.3
            refined ← refined ∪ {e1, e2}
            edgeVerts ← edgeVerts ∪ {midVertex}

        # (5) endpoint (vertex-vertex) collection: each endpoint may be
        #     shared by several frontier edges; collect it once.
        endpointsSeen ← ∅
        for each e ∈ frontier:
            for each v ∈ verts(e):
                if v ∈ endpointsSeen: continue
                endpointsSeen ← endpointsSeen ∪ {v}
                (activeF, inactiveF) ← split vert_faces(v) by D[2][l]
                (activeE, inactiveE) ← split vert_edges(v) by D[1][l]
                if inactiveF = ∅: continue       # not an interface vertex
                vv ← BuildVertexVertex(v, l, activeF, activeE)   # §5.3
                verticesToProcess ← verticesToProcess ∪ {vv}

        verticesToProcess ← verticesToProcess ∪ edgeVerts
        frontier ← refined

# final emission of leftover verticesToProcess
```

**Helper: FindInterfaceEdgesOnLevel(l, D).** An edge `e ∈ E_l` is an
*interface edge* iff exactly one of its two adjacent faces lies in `D[2][l]`
and `e` is not on the domain boundary. (Domain-boundary edges have only one
adjacent face and are skipped.) Lines `346-385` of `activeDoFs.cpp`.

### 5.3 Splitting a frontier edge / building interface vertices

The split rule encodes how a parent interface basis function's support is
inherited by its children. Each child vertex's active faces are the
children of the parent's active fan that are adjacent to that vertex; the
child's active edges are the edges of those active faces incident to the
vertex. Concretely:

When an interface edge `e` on level `l` is split:
- Its active side (the level-`l` active face on one side of `e`) is refined,
  giving four child faces on level `l+1`.
- The two child edges `2·index(e)` and `2·index(e)+1` of `e` each inherit the
  active/inactive face classification of `e`.
- The midpoint vertex of `e` becomes an *edge-vertex* on level `l+1`; its
  active faces are the children of the parent's active side that are incident
  to the midpoint, and its active edges are those edges around the midpoint
  that belong to at least one of these active faces.
- The two endpoint vertices of `e` are *vertex-vertices*: they exist at the
  same index on level `l+1` (subdivision keeps coarse vertices) and their
  active faces are computed by refining the coarse active fan and intersecting
  with the fine vertex's face fan.

```text
Algorithm: SplitInterfaceEdge(e)            # e on level l, returns level l+1

e1 ← 2 · index(e);  e2 ← 2 · index(e) + 1
refinedActive ← refineFaces(activeFaces(e))

for ec in (e1, e2):
    for each face f adjacent to ec on level l+1:
        if f ∈ refinedActive: ec.activeFaces.insert(f)
        else:                 ec.inactiveFaces.insert(f)

midpoint ← the unique vertex shared by e1 and e2
return (e1, e2, BuildEdgeVertex(midpoint, e))


Algorithm: BuildEdgeVertex(v, parent e)     # parent on level l, result on l+1

refinedActive ← refineFaces(activeFaces(e))
for each f ∈ vert_faces(v) on level l+1:
    if f ∈ refinedActive: activeFaces.insert(f)
    else:                 inactiveFaces.insert(f)
edgesOfActive ← ⋃_{f ∈ activeFaces} edges(f)
for each ε ∈ vert_edges(v):
    if ε ∈ edgesOfActive: activeEdges.insert(ε)
    else:                 inactiveEdges.insert(ε)


Algorithm: BuildVertexVertex(v, l, activeF_coarse, activeE_coarse)
            # v is a coarse vertex on level l, kept at same index on level l+1

refinedActive ← refineFaces(activeF_coarse)
for each f ∈ vert_faces(v) on level l+1:
    if f ∈ refinedActive: activeFaces.insert(f)
    else:                 inactiveFaces.insert(f)
edgesOfActive ← ⋃_{f ∈ activeFaces} edges(f)
for each ε ∈ vert_edges(v) on level l+1:
    if ε ∈ edgesOfActive: activeEdges.insert(ε)
    else:                 inactiveEdges.insert(ε)
```

Flag: The base-class comment in `activeDoFs.cpp:98-114` says the vertex
itself contributes a 0-form DoF, *its active edges* contribute 1-form DoFs,
and *its active faces* contribute 2-form DoFs — i.e. step (3) of §5.2 above.
The 2-form contribution is interesting: native 2-form DoFs are already added
on the parent's active side via `BuildNativeDoFs`, so the vertex-emitted
2-form DoFs are precisely the child faces of the *coarser-side* active fan
that meet the interface. Worth double-checking the framing here against the
intent of the paper — the code does what is written, but the interpretation
should be confirmed.

---

## 6. Assembly of the active subdivision matrices

The Galerkin projection of the level-`L` operators to the active space is
realised by sparse matrices `A^0, A^1, A^2` (one per form degree), each of
shape `K_L × K_D` with `K_L = |V_L|, |E_L|, |F_L|` respectively and `K_D`
the total number of active DoFs for that form degree. Each column of `A^k`
is the level-`L` coefficient expansion of one active basis function.

For an active basis function indexed by simplex `s ∈ D[k][l]`, this
expansion is the `s`-th column of `S^k_{l → L}`, the accumulated
subdivision matrix from level `l` to level `L`. `S^k_{l → L}` is the
product `S^k_{L-1} … S^k_{l+1} S^k_l` of the per-level subdivision matrices
(Loop for 0-forms, Wang's scheme for 1-forms, bilinear for 2-forms in the
present implementation).

Two structural conventions matter downstream:

- **Columns are ordered level-by-level.** Within `A^k`, the first
  `|D[k][0]|` columns correspond to level-0 active DoFs, the next
  `|D[k][1]|` to level-1, and so on. This ordering is what the
  preconditioner (§7) exploits when it identifies the top-left
  `n0_coarse × n0_coarse` block of `M0` with the level-0 sub-problem;
  the two-level multigrid variant uses the same convention to construct
  its coarse space.
- **Rows are in DOLFIN's level-`L` ordering.** When the right-hand side
  is assembled by DOLFIN on the embedding mesh, the swap matrix
  `P_swap` (level-`L`) converts between SEC's enumeration of edges /
  vertices and DOLFIN's; `A^kᵀ · P_swapᵀ · rhs_dolfin` gives the active
  RHS.

```text
Algorithm: BuildKFormSubdivMatrix(D, k)

A ← (K_L × K_D) sparse matrix
col ← 0
for l = 0, …, L:
    S ← S^k_{l → L}              # accumulated subdivision matrix
    for each simplex s ∈ D[k][l]:
        A[ : , col ] ← S[ : , s ]
        col ← col + 1
return A
```

Once `A^0, A^1` are built, the active system reads
```math
M0 = (A^0)ᵀ M A^0,   C0 = (A^1)ᵀ C A^1,   B0 = (A^0)ᵀ B A^1,
```
with `M, C, B` the level-`L` mass, curl-curl, and gradient-coupling
matrices respectively. The block system is

```math
K = [[ −M0   B0 ],
     [  B0ᵀ  C0 ]],
```

symmetric indefinite with zero block in the (1,1) entry of the dual
formulation. (`-M0` rather than `+M0` is the convention here so that `K`
has the form of a saddle-point system.)

The driver builds only `A^0` and `A^1` (`main.cpp:622-623`); a 2-form
active matrix `A^2` is not needed because the residual estimator (§8) is
aggregated from per-element DG0 contributions on `T_L` directly, without
projecting through any active 2-form space.

---

## 7. Preconditioner

`K` is symmetric indefinite and is solved with MINRES. Two preconditioner
variants are implemented in `solvers.h`; the one used in all reported
experiments is **block-diagonal direct**:

```math
P = [[ M0   0 ],
     [  0   S ]],   S = C0 + B0ᵀ · diag(M0 · 1)⁻¹ · B0,
```

where `1` is the all-ones vector of length `dim D[0]`. Both blocks of `P`
are SPD and are factorised once per adaptive step by sparse Cholesky
(LDLT); each MINRES iteration costs one back-solve per block.

A few design choices worth recording:

- **`S` is a row-sum-lumped Schur approximation, not bare diagonal
  lumping.** For piecewise-linear CG1 elements on triangles the row sums
  of the mass matrix are uniformly equivalent to the local element
  volumes, while bare diagonal entries are not. On the irregular meshes
  produced by adaptive refinement the row-sum approximation has
  noticeably better mesh-independence as a Schur approximation; this is
  the empirical justification for the choice, mentioned in the
  `main.cpp:651-654` comment block.
- **Direct factorisation is appropriate for the problem sizes
  considered**: the bottleneck on the configurations reported is the
  adaptive cycle itself (refinement + estimator + DoF construction),
  not the solve. The factorisation is dropped at each adaptive step
  because the active spaces are non-nested.
- A **two-level V-cycle variant** (`VCycleBlockDiagPreconditioner`)
  exists in the source but is not used. It would replace each direct
  block solve by a damped-Jacobi smoothed V-cycle with the level-0
  block as the coarse space. It currently does not converge: the
  zero-padding prolongation used by the V-cycle is not compatible with
  the globally-supported subdivision basis (level-0 corrections leak
  into all finer DoFs through the basis functions' support, not only
  the level-0 portion of the active vector). Fixing this would require
  a prolongation that respects `A^k`.

---

## 8. Residual estimator

Two estimator variants are implemented; the one used in all reported
results is the standard residual element estimator for the 1-form
Hodge Laplacian.

For each triangle `T` of the embedding mesh `T_L`, the per-element
indicator is

```math
η_T² = h_T² · ‖f − ∇σ_h‖²_{L²(T)}
       +  Σ_{e ∈ ∂T \ ∂Ω}  h_e · ‖[[ curl u_h ]]_e‖²_{L²(e)},
```

where `[[·]]_e` is the jump across the interior facet `e` and the sum
runs over edges of `T` not on the domain boundary. The two contributions
correspond to the two parts of the strong-form residual of the mixed
Hodge-Laplacian:

- **Interior element residual.** For degree-1 N1curl on triangles,
  `curl(u_h)` is piecewise constant per element, so `curl curl u_h ≡ 0`
  inside each triangle. The strong-form residual `f − curl curl u_h −
  ∇σ_h` therefore reduces to `f − ∇σ_h` in the element interior.
- **Edge jump residual.** All of the `curl curl` information is carried
  by the discontinuity of `curl u_h` across interior facets, expressed
  as the jump term.

Both contributions are scaled with the usual mesh-size factors `h_T²`,
`h_e` so that `η_T` is dimensionally consistent and comparable across
levels — i.e. the per-level scaling does **not** drift as the partition
becomes deeper, which matters because Dörfler accumulates contributions
from many levels at once. (`weak_forms/laplace/eta_form_Hcurl.ufl` line
22-23.)

The UFL implementation has two minor wrinkles relative to the formula
above:

- The edge-size `h_e` is approximated by `avg(h)` — the mean of the cell
  diameters of the two adjacent triangles. For shape-regular meshes this
  is comparable to the true edge length.
- The interior-facet contribution `avg(h) · ‖[[curl u_h]]‖²_{L²(e)}` is
  distributed to both adjacent elements via `avg(q)` (a DG0 test
  function). Each adjacent element therefore receives **half** of the
  edge contribution. The global sum `Σ_T η_T²` then contains each
  interior-edge contribution exactly once, which matches the standard
  AFEM upper-bound statement.

A second variant is implemented as `EstimatorType::DiagM1`: solve `M1
r_h = r1` for the dual residual representation `r_h` (with `r1` the
N1curl-tested residual functional), then evaluate
`η_T² = ∫_T |r_h|²` on each element. The diagonal version replaces `M1`
by its diagonal. This variant is included for comparison but not used
in the reported experiments because the `M1⁻¹` step is sensitive to
high-frequency modes that the subdivision space cannot represent.

The element indicators `{η_T : T ∈ T_L}` are computed as a `DG0`
assembly on the embedding mesh (one quadrature contribution per face,
one jump contribution per interior facet, distributed to the two
adjacent faces via `avg(q)`). They are then aggregated onto the
**active** faces of the partition by walking each level-`L` face up its
quad-tree until an ancestor in `P` is found:

```text
Algorithm: ProjectToActiveDoFs(η²_f over F_L, P)
            # η²_f is the level-L assembled value at face f ∈ F_L

ζ[l][f] ← 0  for all l, f ∈ P[l]
for each fine face f ∈ F_L:
    g ← f
    for l = L, L−1, …, 0:
        if g ∈ P[l]:
            ζ[l][g] ← ζ[l][g] + η²_f       # accumulate the original η²_f
            break
        g ← g / 4                          # parent face on level l − 1

η[l][f] ← √ ζ[l][f]                        # consumed by Alg. 2
```

The resulting `η[l][f]` is the square-root of the summed `η²_T` over all
descendants of `f` that live on `T_L`. The square-root step exists
because `DoerflerMarking` squares its inputs internally — keeping `η`
in linear (not squared) units lets the same indicator be reported in
solution log files in the natural error-norm scale.

Faces on level `L` receive their own `η_T²` (no descendants). Recall
that they are excluded from the threshold total in Dörfler (§2) because
they cannot be refined, even though they remain part of the reported
global error.

---

---

## Things flagged for re-checking

These are items I noticed while extracting the algorithms that may or may
not be deliberate; please confirm before publication.

1. **`EnforceGrading` is dead code in the adaptive driver.** It is defined in
   `doerfler.cpp:593-684` but only `DoerflerMarking` + `FixIsolatedBoundary
   Vertices` are called by `adaptive_vector_laplacian/main.cpp`. The
   one-level admissibility constraint that the function would enforce is
   delivered instead by the interior-vertex completion in `PadMarkedFaces`.
2. **`UpdateActiveDoFs` is declared but undefined.** `activeDoFs.h:204`
   declares `void UpdateActiveDoFs(DoFs&)`; no definition exists in
   `activeDoFs.cpp`. `ActiveDoFs` is currently always reconstructed from
   scratch each adaptive step. Incremental update would be a sensible
   optimisation but is not implemented.
3. **Condition (II) of Def. 38 is enforced by validation but not by padding.**
   `PadMarkedFaces` fixes only condition (I); the condition-(II) branch is
   commented out. `ValidateSubdivisionPartition` and `ValidatePatch` still
   test both. The fallback patch-expansion step in `DoerflerMarking` is
   triggered by `ValidateSubdivisionPartition` failures, so it can in
   principle fire on a (I)-satisfying but (II)-violating patch — I did not
   run a check on whether this happens in practice.
4. **Finest-level faces dropped from the Dörfler total.** The total against
   which `θ · total` is compared excludes faces on level `L`. This changes
   the meaning of `θ` slightly relative to the textbook Dörfler statement;
   already discussed in §2 of this document.
5. **The cluster-aware `α`-boost.** The priority-queue form of Dörfler in
   `doerfler.cpp` boosts edge-neighbours of selected faces by `α`. Experiments
   showed it does not meaningfully change the marked-vs-padded ratio, so the
   trimmed pseudo-code drops it. The function-signature default is `α = 2`,
   but the adaptive driver explicitly passes `alpha = 1.0`
   (`main.cpp:451`), so cluster-awareness is **off by default** in the
   reported runs. With `α = 1` the priority-queue formulation is equivalent
   to a plain descending sort by raw `η²`.
6. **Two-level multigrid preconditioner.** `VCycleBlockDiagPreconditioner`
   exists but does not converge for the subdivision basis (the prolongation
   does not respect the globally-supported basis functions). Reported as
   future work in §7 above.

These items are easy to fix in the prose; flagging them here so we can decide
explicitly whether they are described in the paper as-is or whether the code
should change first.
