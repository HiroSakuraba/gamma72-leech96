# Leech96: globally valid exact shortest-vector enumeration

**Research date:** October 13, 2026  
**Status:** Exact global search framework implemented and tested; **global minimum remains unproved**.

## Why the algorithm changes the investigation

Previous searches enumerated all vectors with specified *small support in the first 48 integer coordinates*. That was valuable but left high-support vectors untouched. This implementation searches the **entire integer coefficient lattice in an exactly verified, unimodularly reduced 96-dimensional basis**. It does not restrict support or make a random-search assumption.

A completed exhaustive run up to squared norm 8 with no nonzero witness would, together with evenness, determinant one, and the rank-96 extremal bound, prove that this explicit candidate has minimum 10. That run **has not completed**; the present reports make no claim to the contrary.

## Data and exact basis validation

Source: original `candidate_3cycle.json` from the published Leech96 work; corresponding 96×96 LLL-reduced Gram matrix and integer basis transformation from the previous research package. The matrix identities were **rechecked afresh** using exact arithmetic.

- Candidate `A` SHA-256 (signed 64-bit entries, little endian): `929e71c3a3f54f3361d9a8c318786b1da726e141250a93969e9b2cec8cd3c656`.
- Reduced Gram JSON SHA-256 (canonical compact JSON): `05aa5b6137cec45bf7fb1524d9fed0847878639df00ebfd63dc261df15452526`.
- The original Gram is `G0=[[ (I+AA^T)/16, A ], [A^T, 16I]]`; congruences guarantee its integer even entries.
- Its determinant is exactly 1 by Schur complement, and the exact transform `U` satisfies `det(U)=+1` and `G=U G0 U^T` entry for entry.
- Fraction-free `LDL^T` decomposition of `G` also gives determinant one and was independently reconstructed in rational arithmetic.
- The last Gram–Schmidt squared length of the published reduced basis is `1/50`, so the initial norm-8 search branches heavily. This is a major computational difficulty, not a defect in the proof method.

## Exact enumeration lemma

For the integer Gram matrix `G`, set `p[-1]=1` and let `p[i]` be the successive leading principal minors. Fraction-free symmetric elimination produces integers `h[j,i]` for `j>i` with

\[
xGx^{T}=\sum_{i=0}^{95}
 \frac{\left(p_i x_i+\sum_{j=i+1}^{95}h_{j,i}x_j\right)^2}
      {p_{i-1}p_i}.
\]

All pivots are strictly positive, all divisions in elimination are exact, and the identity is independently checked against every Gram entry using rational arithmetic.

For a partially assigned suffix `x[95],...,x[i+1]` and remaining squared-norm budget `R=a/b`, put `s=sum(h[j,i]*x[j])`. The bound `xGx^T <= 8` implies

\[
 -H\le p_i x_i+s\le H,
 \qquad H=\left\lfloor\sqrt{\left\lfloor R p_{i-1}p_i\right\rfloor}\right\rfloor.
\]

The implementation computes `H` using `math.isqrt`, both coefficient endpoints using exact integer floor/ceiling division, and the remaining norm using Python's exact `Fraction`. **No floating-point comparison can discard a short vector.**

It enumerates modulo the symmetry `x~-x` by requiring the highest nonzero coefficient in the reduced basis to be positive. The search is equivalent to unrestricted enumeration for the question “does any nonzero vector have squared norm <=8?”

## Findings of the current full-dimensional prototype

- At radius 8, the sign-reduced top-coordinate split `x[95]` has precisely **21 admissible choices**: `0,...,20`. They partition every nonzero vector up to sign.
- **The entire `x[95]=20` partition is excluded.** The `i=95` squared-norm term is `20²/50=8`, leaving zero for all previous terms. At `i=94`, the factor gives `(50*x[94]-460)^2/(2471*50)`; no integer `x[94]` makes this term zero. Thus this disjoint family contains no vector of squared norm <=8, irrespective of the original coordinate support.
- A 110,000-node full search used an exact checkpoint and reported **INCOMPLETE**, one zero leaf reached, **no nonzero witness**. It does **not** certify any of the other 20 partitions.
- All 96-dimensional basis, `E8` roots, complete suffix coverage on an 8-dimensional test, checkpoint/resumption, and intentional-corruption tests passed.

This constitutes a new complete exclusion for one *global reduced-coordinate* partition, not a full lattice minimum proof.

## Files and commands

Python 3.9+ (no optional libraries for the enumerator itself); SymPy is needed only for the independent 96×96 unimodularity test.

```bash
python exact_svp.py selftest
python check_basis_96.py z16_route/candidate_3cycle.json \
  z16_route/candidate_3cycle_LLL_transform.json \
  z16_route/candidate_3cycle_LLL_gram.json basis_check.json
python test_global.py
python exact_svp.py inspect z16_route/candidate_3cycle_LLL_gram.json
python exact_svp.py prefixes z16_route/candidate_3cycle_LLL_gram.json \
  --depth 1 --bound 8 --output global_prefixes_depth1.json
python exact_svp.py search z16_route/candidate_3cycle_LLL_gram.json \
  --bound 8 --suffix 20 --max-nodes 1000 --state partition_20.json
# Genuine 96D but INCOMPLETE until every branch has been exhausted:
python exact_svp.py search z16_route/candidate_3cycle_LLL_gram.json \
  --bound 8 --max-nodes 10000 --state partial.json
python exact_svp.py search z16_route/candidate_3cycle_LLL_gram.json \
  --resume partial.json --max-nodes 10000 --state partial.json
```

Use `global_aggregate.py` **only after** each of the 21 partition reports has status `COMPLETE_NO_SHORT_VECTOR`. It independently regenerates all mathematical suffix partitions from the supplied Gram matrix, then checks SHA, radius, exact partition coverage and completeness, but is only a *coverage audit*: it cannot establish that external JSON log claims are truthful without rerunning their enumerations. Partial states and timeouts are *not* certificates.

`global_witness.py` independently maps any found reduced-lattice integer coefficient vector back into the original 96 integer coordinates `(u,v)` and checks the exact modular code constraint and squared norm.

## Practical implications

The initial enumeration is expensive even after the previous LLL basis transformation. Full proof probably requires a substantially better reduced basis (for example a mathematically auditable Block Korkine–Zolotarev (BKZ) basis transformation), high-performance integer-enumeration code, and additional mathematically **necessary** pruning bounds. Merely increasing a randomized BKZ search effort cannot certify a lower bound. The present implementation gives a rigorous ground truth for checking more aggressive proof candidates.

We deliberately do **not** create an unbounded GitHub Actions job: the full search's resource requirements are unknown and could be enormous. Automated checks should run only the finite small-rank regressions, the exact basis validation, the `x[95]=20` exclusion, and a bounded resume smoke test.
