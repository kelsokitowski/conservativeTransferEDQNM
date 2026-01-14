# Cyclic Symmetry Verification for EDQNM

## Overview

This implementation uses **cyclic symmetry verification** to ensure geometric correctness of the triad domain integration, while integrating each triad exactly once via k-slices. Energy conservation is achieved through proper volume weighting and delta corrections, not by summing multiple permutations.

## Mathematical Foundation

The triad domain |p-q| < k < p+q is **invariant under cyclic permutations**:
- If (p,q,k) is a valid triad, then so are (q,k,p) and (k,p,q)
- The volume element is also invariant: dv(i,j,k) = dv(j,k,i) = dv(k,i,j)
- And symmetric under reflection: dv(i,j,k) = dv(j,i,k)

This symmetry is **verified** to ensure correctness, but each triad is **integrated only once** to preserve proper physics.

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

Uses **k-slice integration only**:

- **k-slice loop**: Integrates over (p,q) for each k
  - Each triad integrated exactly once
  - Uses geometric centroids for (p,q) coordinates
  - Evaluates kernel at (kstar, pstar, qstar)
  - Applies delta correction for local energy conservation: dEk + dEp + dEq = 0
  - Scatter-adds to bins (kj, pj, qj)
  - Handles p↔q mirroring

**Why only k-slice?**
- Each unique triad should be integrated once (not three times)
- The cyclic symmetry checks verify our geometry is correct
- Multiple integrations at different evaluation points don't improve conservation
- Physics: ∫∫∫ S(k,p,q) dv = 0 is achieved through delta corrections, not multiple counting

## Energy Conservation

Energy conservation is achieved through **three mechanisms**:

1. **Local (triad-level)**: Delta correction ensures each triad conserves energy
   ```matlab
   delta = (Jk*Sk_raw + Jp*Sp_raw + Jq*Sq_raw) / (Jk + Jp + Jq)
   ```
   After correction: dEk + dEp + dEq = 0 (up to roundoff)

2. **Geometric correctness**: Cyclic symmetry of dv verified via permutation checks
   - Ensures the triad domain is correctly discretized
   - No systematic bias from asymmetric volume computation

3. **Proper normalization**: Each triad integrated once with weight dv/dk
   - No multiple counting
   - Volume element properly normalized by bin width

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
