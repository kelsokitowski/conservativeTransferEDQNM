# Cyclic Symmetry Implementation for EDQNM

## Overview

This implementation achieves energy conservation through **geometric cyclic symmetry** of the triad domain integration. Instead of relying solely on algebraic delta corrections, we sum contributions from all three cyclic permutations (k,p,q) → (p,q,k) → (q,k,p).

## Mathematical Foundation

The triad domain |p-q| < k < p+q is **invariant under cyclic permutations**:
- If (p,q,k) is a valid triad, then so are (q,k,p) and (k,p,q)
- The volume element is also invariant: dv(i,j,k) = dv(j,k,i) = dv(k,i,j)
- And symmetric under reflection: dv(i,j,k) = dv(j,i,k)

This means each triad can be integrated in three equivalent ways, and summing all three ensures symmetric treatment of all modes.

## Key Changes

### 1. buildTriadWeightsCentroidsExact.m

**Core Insight**: Fill `dv` array with all 6 permutations (3 cyclic + 3 p↔q mirrors), then compute weights.

The `dv` array satisfies:
- **Cyclic symmetry**: `dv(i,j,k) = dv(j,k,i) = dv(k,i,j)`
- **Reflection symmetry**: `dv(i,j,k) = dv(j,i,k)`

Then compute three weight arrays:
- **k-slices**: `weight_k(pj,qj,kj) = dv(pj,qj,kj) / dk(kj)`
  - Integrate over (p,q) for fixed k
  - Centroids: **Geometric centroids (p*, q*)** from polygon moments
  - These are exact centroids of the triad domain slice

- **p-slices**: `weight_p(qj,kj,pj) = dv(qj,kj,pj) / dk(pj)`
  - Integrate over (q,k) for fixed p
  - Centroids: **Bin centers** for both q and k coordinates
  - Storage uses cyclically permuted indices

- **q-slices**: `weight_q(kj,pj,qj) = dv(kj,pj,qj) / dk(qj)`
  - Integrate over (k,p) for fixed q
  - Centroids: **Bin centers** for both k and p coordinates
  - Storage uses cyclically permuted indices

**Important**:
1. Due to cyclic symmetry of `dv` and the fact that all three divide by the 3rd index, the three weight arrays are **identical**
2. Only k-slice uses geometric centroids; p-slice and q-slice use bin centers (standard FV approximation)
3. Computing exact centroids for all three cyclic permutations would require triple the computational cost

### 2. transfer_scatter_add_kernelE0.m

Now sums contributions from **three separate integration loops**:

1. **k-slice loop**: Integrates over (p,q) slices for each k
2. **p-slice loop**: Integrates over (q,k) slices for each p
3. **q-slice loop**: Integrates over (k,p) slices for each q

Each loop:
- Retrieves appropriate weights and centroids
- Evaluates kernel at (kstar, pstar, qstar)
- Applies delta correction for local energy conservation
- Scatter-adds to bins (kj, pj, qj)
- Handles p↔q mirroring

## Energy Conservation

Energy conservation is achieved through **two mechanisms**:

1. **Local (triad-level)**: Delta correction ensures each triad conserves energy
   ```matlab
   delta = (Jk*Sk_raw + Jp*Sp_raw + Jq*Sq_raw) / (Jk + Jp + Jq)
   ```

2. **Global (geometric)**: Cyclic symmetry ensures balanced treatment of all modes
   - Each mode k receives contributions where it plays k-leg, p-leg, and q-leg roles
   - Symmetric integration measure prevents systematic bias

## Notes on Centroids

**K-slice (exact geometry):**
- Full geometric centroids (p*, q*) computed via polygon moments (Sutherland-Hodgman clipping + Green's theorem)
- These are the exact centroids of the triad domain |p-q| < k < p+q intersected with the (p,q) bin
- Guaranteed to lie inside the triad domain

**P-slice and Q-slice (bin center approximation):**
- Use bin centers for both coordinates: `kVals(i)` and `kVals(j)`
- Standard FV approximation - avoids computing full 3D geometric centroids for all permutations
- Reduces computational cost by 2/3 while maintaining conservation properties
- The cyclic symmetry of the integration domain ensures balanced treatment regardless

## Expected Behavior

When running main.m, expect:
- `Cyclic symmetry errors in dv`: Should be O(1e-15) or exactly 0
- `p<->q symmetry error in dv`: Should be O(1e-15) or exactly 0
- `p<->q symmetry errors in weights`: Should be O(1e-15) or exactly 0
- `Weight consistency`: All three weight arrays should be identical (error ~ 0)
- `FV total energy transfer`: Should be ≈0 (limited by kernel symmetry and numerical precision)
- `max triad energy residual`: Should be O(1e-15) (machine precision)

## Future Work

When replacing toy kernel with actual EDQNM kernel:
- Ensure kernel has form S_k(k,p,q) = f(E(k), E(p), E(q))
- Verify kernel satisfies S_k + S_p + S_q = 0 algebraically
- Check that cyclic symmetry improves energy conservation vs. single-slice integration
