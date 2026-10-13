# Leech96: exact support-five exclusion by modulo-8 lower bounds

**13 October 2026 · Status: exact finite computation for a restricted family, not a minimum-norm proof**

This update is for the explicit three-cycle Construction-A lattice recorded in
[`z16_candidate_checks/candidate_3cycle.json`](z16_candidate_checks/candidate_3cycle.json).
The matrix `A` is unchanged; SHA-256 (48×48 signed 64-bit little-endian integer entries):
`929e71c3a3f54f3361d9a8c318786b1da726e141250a93969e9b2cec8cd3c656`.

## Main new result

No vector of squared norm at most 8 can have **exactly five nonzero
first-half integer coordinates, after centering modulo 16**. This includes all
odd/even coefficient mixtures.

An additional finite check excludes vectors of squared norm at most 8 with
**exactly six nonzero first-half coordinates, all six odd**.

The problem of excluding short vectors with six mixed-parity first-half
coordinates, or with at least seven first-half coordinates and some odd
coordinates, remains completely open. The earlier even-coordinate certificate
already covers *all* first-half-even vectors. Neither new finite search
proves the 96-dimensional lattice is extremal.

## Mathematical statement and proof of the search reduction

The lattice consists of

```text
(u,v)/4, where u,v in Z^48, v == u A (mod 16),
squared norm = (sum(u_i^2) + sum(v_i^2))/16.
```

We search for vectors whose numerator is at most `16*8=128`.
The independent structural checker verifies `A A^T == -I (mod 16)`,
the Type II row divisibility by 32, the precise matrix hash, and
`A mod 2 = J - I`.

**Centering.** Adding any multiple of 16 to an individual coordinate of
`u` or `v` preserves lattice membership. For a norm-at-most-8 vector,
replacing coordinates by their least-absolute-value representatives modulo 16
never increases the norm. Consequently every nonzero centered `u_i`
is in `{-8,...,-1,1,...,8}`; each `v_i` can be independently centered
modulo 16. Negating the vector permits the *first* nonzero `u` coefficient
to be positive.

**Mod-8 lower bound.** For a residue `r=0,...,7` put

```text
c8(r) = min{|x|^2 : x in Z, x == r (mod 8)}
      = [0,1,4,9,16,9,4,1][r].
```

For any fixed support `S` and nonzero mod-8 coefficients `r`
on the support, the actual norm numerator satisfies

```text
||u||² + ||v||² >= sum_i m_i(r_i)
                    + sum_{j=0}^{47} c8((r A)_j mod 8),
```

where the first supported coordinate has `m_0(r)=r²` if `r=1..7`,
or 64 for residue zero (`u=+8`), while every other supported coordinate
has `m_i(r)=c8(r)` for `r!=0`, or 64 when `r=0`.
The exceptional residue-zero bound is essential: the support requires
`u_i != 0`, so `u_i=+8` or `-8`, not zero.

Any residue class with this exact integer lower bound greater than 128 is
impossible. This rejection cannot eliminate a genuine short vector.

**Exhaustive lifts.** For each surviving mod-8 pattern, the first
coefficient has its unique positive representative `1..8`.
For each remaining coefficient with residue `r=1..7`,
the two centered mod-16 representatives are `r` and `r-8`;
for residue `r=0`, they are `+8` and `-8`.
All `2^(support_size-1)` lifts are considered.

The code uses a second safe lower bound before the full calculation.
Let `k` be the number of odd `u_i`. Since `A mod 2=J-I`,
if `k` is even then the number of odd entries of `v` is `k`;
if `k` is odd that number is `48-k`. Hence

```text
||v||² >= (k if k is even else 48-k).
```

Every lift surviving this test is evaluated with the **exact** 48-coordinate
modulo-16 residue vector `v=uA (mod 16)` and the squared-norm table
`[0,1,4,9,16,25,36,49,64,49,36,25,16,9,4,1]`.
The search stops only when the running nonnegative sum exceeds 128.

The code's AVX2 byte-vector operations only compute the low three bits
of sums and lookup eight-entry squared-distance tables. They do not
introduce floating-point arithmetic. A scalar Python reference separately
tests three support patterns, including one with two mod-8 survivors.

## Exact computed results

| Scope | Supports | Residue classes | Passed mod-8 bound | Checked lifts | Norm ≤ 8 |
|---|---:|---:|---:|---:|---:|
| All 5-coordinate patterns | 1,712,304 | 56,108,777,472 | 15,278 | 127,019 | **0** |
| All 6-coordinate patterns with all 6 odd | 12,271,512 | 50,264,113,152 | 14,898 | 252,059 | **0** |

The five-coordinate survivor distribution by number of odd supported
coordinates is `[0, 5714, 635, 3895, 1498, 3536]`.
The sum is 15,278.

Both computations have been run to completion on the published candidate.
The reported class counts count mod-8 *coefficient patterns after the
global-sign convention*, not complete integer lattice vectors.

## How to reproduce

On x86-64 with AVX2, C++17, OpenMP, and Python 3.9+:

```bash
python3 z16_candidate_checks/support_mod8_reference.py \
  z16_candidate_checks/candidate_3cycle.json

OMP_NUM_THREADS=4 python3 z16_candidate_checks/support_mod8_verify.py \
  z16_candidate_checks/candidate_3cycle.json five /tmp/support5.json

OMP_NUM_THREADS=4 python3 z16_candidate_checks/support_mod8_verify.py \
  z16_candidate_checks/candidate_3cycle.json six /tmp/support6.json
```

The wrapper generates a temporary exact candidate matrix header, compiles
the scanner, asserts all documented full-search counters, rejects any
norm-8 witness, and outputs JSON. GitHub Actions reruns the checks in CI.
The scalar reference is a spot-check and must not be represented as
an independent 56-billion-case implementation.

## Next feasible mathematical move

A complete six-coordinate all-parities run would require substantially more
work (naively `choose(48,6)*8^6` mod-8 classes). Prefer a cheaper
support-level necessary-condition filter, then shard the remaining
six-coordinate mixed-parity families. An exact short vector would immediately
kill this candidate; an exhaustive global exclusion is still the objective.

The theta-series and automorphism results from the earlier work are
unaffected. No Hermitian obstruction or global rank-96 existence claim
is implied by the new sparse-support certificates.
