#include "calculate.hpp"
#include "variant.hpp"
#include "operators.h"
#include "tree.hpp"
#include "utils.hpp"
#include <iostream>
#include <omp.h>

// A cell's z centre for a real 3D tree, or 0.0 for a planar (D=2) one, where
// Cell<D>::centre only has 2 components. Every cell-cell displacement in this
// file goes through this rather than storing a z=0 shell per cell.
template <int D>
static inline double zc(const std::array<double, D> &centre) {
  if constexpr (D == 3) return centre[2];
  else return 0.0;
}

// Dispatches to the genuinely 2-argument-position P2P_batchxy for a planar
// (D=2) tree, or the ordinary 3-argument P2P_batch otherwise. Guarded by
// FMMGEN_PLANAR so example driver builds that never set planar=True still
// compile for D=3 -- the fallback (#else branch) pads with the tz/bz the
// caller already has, exactly as every other operator call in this file
// does via zc<D>()/Tree::body_z_zero, and is what every call site used
// unconditionally before this existed.
#ifdef FMMGEN_PLANAR
template <int D>
static inline void p2p_batch_dispatch(double tx, double ty, double tz,
                                      const double *bx, const double *by, const double *bz,
                                      const double *bS, size_t begin, size_t end, double *F) {
  if constexpr (D == 2) {
    P2P_batchxy(tx, ty, bx, by, bS, begin, end, F);
  } else {
    P2P_batch(tx, ty, tz, bx, by, bz, bS, begin, end, F);
  }
}
#else
template <int D>
static inline void p2p_batch_dispatch(double tx, double ty, double tz,
                                      const double *bx, const double *by, const double *bz,
                                      const double *bS, size_t begin, size_t end, double *F) {
  P2P_batch(tx, ty, tz, bx, by, bz, bS, begin, end, F);
}
#endif

template <int D>
void M_sanity_check(const std::vector<Cell<D>> &cells) {
	double M0 = 0;
	for(size_t c = 1; c < cells.size(); c++) {
      if (cells[c].nchild == 0) {
		    M0 += cells[c].M[0];
	    }
  }
	std::cout << "Cell 0 has M0 = " << cells[0].M[0] << std::endl;
	std::cout << "Other cells   = " << M0 << std::endl;
  if (std::abs((cells[0].M[0] - M0)/M0) > 10e-10) {
    throw std::runtime_error("M0 sanity check failed");
  }
}


// interact_dehnen (the non-lazy, immediate-evaluation traversal) and P2P_Cells
// were removed here. interact_dehnen was only ever called from itself, so it
// was unreachable; P2P_Cells lost its last caller when evaluate_P2P_lazy was
// restructured to group by target below.

// One step of the dual-tree traversal: either classify (A,B) as a terminal
// interaction, or emit the child pairs it decomposes into.
//
// Factored out so the traversal can be driven two ways: breadth-first to carve
// out independent work, then depth-first inside each piece.
template <int D>
static inline void dehnen_step(const size_t A, const size_t B,
                               const std::vector<Cell<D>> &cells,
                               const double theta, const size_t ncrit,
                               std::vector<Interaction> &M2L_list,
                               std::vector<Interaction> &P2P_list,
                               std::vector<Interaction> *children) {
  const double dx = cells[A].centre[0] - cells[B].centre[0];
  const double dy = cells[A].centre[1] - cells[B].centre[1];
  const double dz = zc<D>(cells[A].centre) - zc<D>(cells[B].centre);

  // Squared multipole acceptance criterion; see the note in the header.
  const double R2 = dx*dx + dy*dy + dz*dz;
  const double rsum = cells[A].rmax + cells[B].rmax;

  if (R2 * theta * theta > rsum * rsum) {
    M2L_list.push_back(Interaction{(uint32_t)A, (uint32_t)B});
  }
  else if (cells[A].nchild == 0 && cells[B].nchild == 0) {
    // The previous version had an `if (cells[B].nleaf >= ncrit)` branch here
    // that pushed an M2L pair and also evaluated M2L immediately. It was
    // unreachable: nchild == 0 implies nleaf < ncrit by construction, so the
    // condition can never hold. Worse, it ran during build_tree, before
    // cells[].M and cells[].L had been pointed at their storage. Dropping it
    // leaves the interaction lists byte-identical.
    P2P_list.push_back(Interaction{(uint32_t)A, (uint32_t)B});
  }
  else if (cells[B].nchild == 0 || (cells[A].rmax >= cells[B].rmax && cells[A].nchild != 0)) {
    for (size_t oa = 0; oa < Cell<D>::NCHILD; oa++) {
      if (cells[A].nchild & (1 << oa))
        children->push_back(Interaction{(uint32_t)cells[A].child[oa], (uint32_t)B});
    }
  }
  else {
    for (size_t ob = 0; ob < Cell<D>::NCHILD; ob++) {
      if (cells[B].nchild & (1 << ob))
        children->push_back(Interaction{(uint32_t)A, (uint32_t)cells[B].child[ob]});
    }
  }
}

// Depth-first traversal of one subtree pair, into private lists.
template <int D>
static void dehnen_dfs(const size_t A, const size_t B,
                       const std::vector<Cell<D>> &cells,
                       const double theta, const size_t ncrit,
                       std::vector<Interaction> &M2L_list,
                       std::vector<Interaction> &P2P_list) {
  std::vector<Interaction> kids;
  dehnen_step<D>(A, B, cells, theta, ncrit, M2L_list, P2P_list, &kids);
  for (size_t k = 0; k < kids.size(); k++) {
    dehnen_dfs<D>(kids[k].target, kids[k].source, cells, theta, ncrit, M2L_list, P2P_list);
  }
}

template <int D>
void interact_dehnen_lazy(const size_t A, const size_t B,
                          const std::vector<Cell<D>> &cells,
                          const std::vector<Particle<D>> &particles,
                          const double theta, const size_t order,
                          const size_t ncrit,
                          std::vector<Interaction> &M2L_list,
                          std::vector<Interaction> &P2P_list) {
  (void)particles; (void)order;

  // The traversal used to be a serial recursion, and was one of the larger
  // remaining costs as well as a hard Amdahl ceiling.
  //
  // exafmm parallelises it with `#pragma omp task untied`, but that makes the
  // order in which interactions are appended depend on the scheduler, which
  // would cost us the bit-reproducibility that removing the atomics bought.
  //
  // Instead: expand breadth-first until there are enough independent subtree
  // pairs to keep the threads busy, then run each depth-first into its OWN
  // buffer and concatenate the buffers in index order. The work split and the
  // concatenation order are both fixed, so the result is identical to the
  // serial traversal regardless of thread count or scheduling.
  const size_t target_tasks = 8 * (size_t)omp_get_max_threads();

  std::vector<Interaction> frontier{Interaction{(uint32_t)A, (uint32_t)B}};
  std::vector<Interaction> next;
  while (frontier.size() < target_tasks) {
    next.clear();
    for (size_t f = 0; f < frontier.size(); f++) {
      dehnen_step<D>(frontier[f].target, frontier[f].source, cells, theta, ncrit,
                     M2L_list, P2P_list, &next);
    }
    // Swap before checking emptiness: on break, frontier must become the
    // (now-empty) set of pairs still needing depth-first work, not the set
    // just fully classified by dehnen_step above -- otherwise the depth-first
    // phase below re-classifies and double-counts every pair from the last
    // breadth-first round. Only reachable for trees small enough to fully
    // classify before target_tasks pairs accumulate.
    frontier.swap(next);
    if (frontier.empty()) break;
  }

  std::vector<std::vector<Interaction>> m2l_buf(frontier.size());
  std::vector<std::vector<Interaction>> p2p_buf(frontier.size());

  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t f = 0; f < frontier.size(); f++) {
    dehnen_dfs<D>(frontier[f].target, frontier[f].source, cells, theta, ncrit,
                 m2l_buf[f], p2p_buf[f]);
  }

  for (size_t f = 0; f < frontier.size(); f++) {
    M2L_list.insert(M2L_list.end(), m2l_buf[f].begin(), m2l_buf[f].end());
    P2P_list.insert(P2P_list.end(), p2p_buf[f].begin(), p2p_buf[f].end());
  }
}

template <int D>
void evaluate_P2M(std::vector<Cell<D>> &cells,
                  const double *bx, const double *by, const double *bz,
                  const double *bS, size_t ncrit, size_t exporder) {
  // Uses the generated S2M rather than seeding a scratch multipole with the
  // source's moments and calling the general M2M. The old form encoded an
  // undocumented basis assumption in the driver -- that a point source's
  // multipole coefficients ARE its moments, true for the Cartesian basis but
  // not in general -- and paid for a full M2M when all but FMMGEN_SOURCESIZE
  // input coefficients are known to be zero. S2M is that specialisation:
  // 16.9x fewer multiplies at p=11, 4.1x at p=5.
  #pragma omp for
  for(size_t c = 0; c < cells.size(); c++) {
    if (cells[c].nleaf < ncrit) {
      const size_t o = cells[c].body_offset;
      for(size_t m = o; m < o + cells[c].nleaf; m++) {
        fmm->s2m(cells[c].centre[0] - bx[m], cells[c].centre[1] - by[m],
            zc<D>(cells[c].centre) - bz[m],
            const_cast<double*>(&bS[FMMGEN_SOURCESIZE*m]), cells[c].M, exporder);
      }
   }
  }
}

template <int D>
void evaluate_M2M(std::vector<Cell<D>> &cells,
                  const std::vector<std::vector<size_t>> &levels, size_t exporder) {
  /*
  evaluate_M2M with level-synchronous parallel traversal.

  Process tree bottom-up, parallelizing within each level.
  Cells at the same level can be processed in parallel since
  they write to different parent cells.
  */
  // Bottom-up: start from deepest level and work up to level 1 (skip root at level 0).
  //
  // Parallelise over the PARENT cells at level l-1, not over the children at
  // level l. Siblings share a parent, so distributing children across threads
  // means up to 2^D threads accumulate into the same cells[p].M concurrently.
  // Iterating parents gives each thread exclusive ownership of the cell it
  // writes, which is the same pattern evaluate_L2L already uses.
  for (int l = levels.size() - 1; l > 0; l--) {
    #pragma omp for schedule(static)
    for (size_t i = 0; i < levels[l-1].size(); i++) {
      size_t p = levels[l-1][i];
      for (size_t octant = 0; octant < Cell<D>::NCHILD; octant++) {
        if (cells[p].nchild & (1 << octant)) {
          size_t c = cells[p].child[octant];
          double dx = (cells[p].centre[0] - cells[c].centre[0]);
          double dy = (cells[p].centre[1] - cells[c].centre[1]);
          double dz = zc<D>(cells[p].centre) - zc<D>(cells[c].centre);
          fmm->m2m(dx, dy, dz, cells[c].M, cells[p].M, exporder);
        }
      }
    }
  }
}


template <int D>
void evaluate_M2L_lazy(std::vector<Cell<D>> &cells,
                       std::vector<Interaction> &M2L_list,
                       std::vector<size_t> &M2L_group, size_t order) {
    // M2L_list is grouped by target, so one thread can own an entire target.
    // That thread is the only writer of cells[A].L, which lets M2L accumulate
    // straight into it. Compared with the previous version this removes, per
    // interaction: a heap allocation, the zeroing of the temporary, and an
    // lsize-long loop of atomic adds (364 atomics per interaction at p=11,
    // roughly 660M over a typical run).
    //
    // Entries for target A live in [M2L_group[A], M2L_group[A+1]); empty
    // targets are skipped. Groups vary widely in size, hence dynamic
    // scheduling.
    const size_t ngroups = M2L_group.empty() ? 0 : M2L_group.size() - 1;
    #pragma omp for schedule(dynamic, 16)
    for (size_t A = 0; A < ngroups; A++) {
        const size_t begin = M2L_group[A], end = M2L_group[A+1];
        if (begin == end) continue;
        double *const L = cells[A].L;
        const double ax = cells[A].centre[0], ay = cells[A].centre[1], az = zc<D>(cells[A].centre);
        for (size_t i = begin; i < end; i++) {
            const size_t B = M2L_list[i].source;
            fmm->m2l(ax - cells[B].centre[0], ay - cells[B].centre[1], az - zc<D>(cells[B].centre),
                cells[B].M, L, order);
        }
    }

}

template <int D>
void evaluate_P2P_lazy(std::vector<Cell<D>> &cells,
                       const double *bx, const double *by, const double *bz,
                       const double *body_S,
                       std::vector<Interaction> &P2P_list,
                       std::vector<size_t> &P2P_group, double *F) {
   // Grouped by target, so the owning thread is the sole writer of F for that
   // target's particles: no atomics, and the result is bit-reproducible.
   //
   // F is in Morton order here too, so the target write is contiguous and the
   // body_perm indirection is gone entirely; the caller scatters once at the end.
   //
   // Sources are handed to the GENERATED P2P_batch kernel a contiguous run at
   // a time. The per-pair P2P() lives in another translation unit, so calling
   // it per interaction cost a real function call and made vectorisation
   // impossible -- the object code had zero vector instructions here.
   const size_t ngroups = P2P_group.empty() ? 0 : P2P_group.size() - 1;
   #pragma omp for schedule(dynamic, 16)
   for (size_t A = 0; A < ngroups; A++) {
       const size_t begin = P2P_group[A], end = P2P_group[A+1];
       if (begin == end) continue;
       const size_t ao = cells[A].body_offset, an = cells[A].nleaf;
       for (size_t i = begin; i < end; i++) {
           const size_t B = P2P_list[i].source;
           const size_t bo = cells[B].body_offset, bn = cells[B].nleaf;
           for (size_t t = ao; t < ao + an; t++) {
               double *const Fl = &F[FMMGEN_OUTPUTSIZE * t];
               // Self-interaction is excluded by SPLITTING the source range
               // around t rather than testing inside the loop: the kernel stays
               // branch-free, and index-based exclusion preserves the original
               // semantics exactly. Masking on r2 > 0 instead would also drop
               // coincident-but-distinct particles, a silent behaviour change.
               if (t >= bo && t < bo + bn) {
                   p2p_batch_dispatch<D>(bx[t], by[t], bz[t], bx, by, bz, body_S, bo, t, Fl);
                   p2p_batch_dispatch<D>(bx[t], by[t], bz[t], bx, by, bz, body_S, t+1, bo+bn, Fl);
               } else {
                   p2p_batch_dispatch<D>(bx[t], by[t], bz[t], bx, by, bz, body_S, bo, bo+bn, Fl);
               }
           }
       }
   }
}




template <int D>
void evaluate_L2L(std::vector<Cell<D>> &cells, const std::vector<std::vector<size_t>> &levels,
                  size_t exporder) {
  /*
  evaluate_L2L with level-synchronous parallel traversal.

  Process tree top-down, parallelizing within each level.
  Cells at the same level can be processed in parallel.
  */
  // Top-down: start from root level and work down
  for (size_t l = 0; l < levels.size() - 1; l++) {
    #pragma omp for schedule(static)
    for (size_t i = 0; i < levels[l].size(); i++) {
      size_t p = levels[l][i];
      for (size_t octant = 0; octant < Cell<D>::NCHILD; octant++) {
        if (cells[p].nchild & (1 << octant)) {
          size_t c = cells[p].child[octant];
          double dx = cells[c].centre[0] - cells[p].centre[0];
          double dy = cells[c].centre[1] - cells[p].centre[1];
          double dz = zc<D>(cells[c].centre) - zc<D>(cells[p].centre);
          // Each thread processes different parent cells, so no race condition
          fmm->l2l(dx, dy, dz, cells[p].L, cells[c].L, exporder);
        }
      }
    }
  }
}

template <int D>
void evaluate_L2P(std::vector<Cell<D>> &cells,
                  const double *bx, const double *by, const double *bz,
                  double *F, size_t ncrit, size_t exporder) {
  #pragma omp for schedule(runtime)
  for (size_t i = 0; i < cells.size(); i++) {
    if (cells[i].nleaf < ncrit) {
      const size_t o = cells[i].body_offset;
      for (size_t m = o; m < o + cells[i].nleaf; m++) {
        fmm->l2p(bx[m] - cells[i].centre[0], by[m] - cells[i].centre[1], bz[m] - zc<D>(cells[i].centre),
            cells[i].L, &F[FMMGEN_OUTPUTSIZE*m], exporder);
      }
    }
  }
}

template <int D>
void evaluate_direct(const double *bx, const double *by, const double *bz,
                     const double *bS, double *F, size_t n) {
  #pragma omp parallel for schedule(runtime)
  for (size_t i = 0; i < n; i++) {
    double *const Fl = &F[FMMGEN_OUTPUTSIZE*i];
    // Split around i so the batched kernel stays branch-free.
    if (i > 0)     p2p_batch_dispatch<D>(bx[i], by[i], bz[i], bx, by, bz, bS, 0, i, Fl);
    if (i + 1 < n) p2p_batch_dispatch<D>(bx[i], by[i], bz[i], bx, by, bz, bS, i+1, n, Fl);
  }
}

template <int D>
void evaluate_M2P_and_P2P(const double *bx, const double *by, const double *bz,
  const double *bS, unsigned int p, size_t m,
  std::vector<Cell<D>> &cells, double *F, unsigned int n_crit, double theta,
  unsigned int exporder) {
  // For Morton-ordered particle m, start at cell p.
  double dx, dy, dz;
  if (cells[p].nleaf >= n_crit) {
    for (size_t octant = 0; octant < Cell<D>::NCHILD; octant++) {
      if (cells[p].nchild & (1 << octant)) {
        const size_t c = cells[p].child[octant];
        dx = bx[m] - cells[c].centre[0];
        dy = by[m] - cells[c].centre[1];
        dz = bz[m] - zc<D>(cells[c].centre);
        // Squared form of  cells[c].r > theta * r,  as in the FMM criterion.
        const double r2 = dx*dx + dy*dy + dz*dz;
        if (cells[c].r * cells[c].r > theta * theta * r2) {
            evaluate_M2P_and_P2P<D>(bx, by, bz, bS, c, m, cells, F, n_crit, theta, exporder);
        }
        else {
            fmm->m2p(dx, dy, dz, cells[c].M, &F[FMMGEN_OUTPUTSIZE*m], exporder);
        }
      }
    }
  }
  else {
    // Leaf: its particles are the contiguous run [o, o + nleaf). Split around
    // m so the batched kernel stays branch-free, as elsewhere.
    const size_t o = cells[p].body_offset, n = cells[p].nleaf;
    double *const Fl = &F[FMMGEN_OUTPUTSIZE*m];
    if (m >= o && m < o + n) {
      if (m > o)         p2p_batch_dispatch<D>(bx[m], by[m], bz[m], bx, by, bz, bS, o, m, Fl);
      if (m + 1 < o + n) p2p_batch_dispatch<D>(bx[m], by[m], bz[m], bx, by, bz, bS, m+1, o+n, Fl);
    } else {
      p2p_batch_dispatch<D>(bx[m], by[m], bz[m], bx, by, bz, bS, o, o+n, Fl);
    }
  }
}

// D = 2 (quadtree, planar) and D = 3 (octree) are the only supported
// instantiations -- see the static_assert in Particle<D>.
template void M_sanity_check<2>(const std::vector<Cell<2>> &);
template void M_sanity_check<3>(const std::vector<Cell<3>> &);
template void evaluate_P2M<2>(std::vector<Cell<2>> &, const double *, const double *, const double *, const double *, size_t, size_t);
template void evaluate_P2M<3>(std::vector<Cell<3>> &, const double *, const double *, const double *, const double *, size_t, size_t);
template void evaluate_M2M<2>(std::vector<Cell<2>> &, const std::vector<std::vector<size_t>> &, size_t);
template void evaluate_M2M<3>(std::vector<Cell<3>> &, const std::vector<std::vector<size_t>> &, size_t);
template void evaluate_L2L<2>(std::vector<Cell<2>> &, const std::vector<std::vector<size_t>> &, size_t);
template void evaluate_L2L<3>(std::vector<Cell<3>> &, const std::vector<std::vector<size_t>> &, size_t);
template void evaluate_L2P<2>(std::vector<Cell<2>> &, const double *, const double *, const double *, double *, size_t, size_t);
template void evaluate_L2P<3>(std::vector<Cell<3>> &, const double *, const double *, const double *, double *, size_t, size_t);
template void evaluate_direct<2>(const double *, const double *, const double *, const double *, double *, size_t);
template void evaluate_direct<3>(const double *, const double *, const double *, const double *, double *, size_t);
template void interact_dehnen_lazy<2>(const size_t, const size_t, const std::vector<Cell<2>> &, const std::vector<Particle<2>> &, const double, const size_t, const size_t, std::vector<Interaction> &, std::vector<Interaction> &);
template void interact_dehnen_lazy<3>(const size_t, const size_t, const std::vector<Cell<3>> &, const std::vector<Particle<3>> &, const double, const size_t, const size_t, std::vector<Interaction> &, std::vector<Interaction> &);
template void evaluate_M2L_lazy<2>(std::vector<Cell<2>> &, std::vector<Interaction> &, std::vector<size_t> &, size_t);
template void evaluate_M2L_lazy<3>(std::vector<Cell<3>> &, std::vector<Interaction> &, std::vector<size_t> &, size_t);
template void evaluate_P2P_lazy<2>(std::vector<Cell<2>> &, const double *, const double *, const double *, const double *, std::vector<Interaction> &, std::vector<size_t> &, double *);
template void evaluate_P2P_lazy<3>(std::vector<Cell<3>> &, const double *, const double *, const double *, const double *, std::vector<Interaction> &, std::vector<size_t> &, double *);
template void evaluate_M2P_and_P2P<2>(const double *, const double *, const double *, const double *, unsigned int, size_t, std::vector<Cell<2>> &, double *, unsigned int, double, unsigned int);
template void evaluate_M2P_and_P2P<3>(const double *, const double *, const double *, const double *, unsigned int, size_t, std::vector<Cell<3>> &, double *, unsigned int, double, unsigned int);
