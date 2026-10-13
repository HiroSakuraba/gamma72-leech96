# Leech96: support-seven odd-coordinate modular-syndrome search

**October 13, 2026 · Exact computation with explicit scope, not a proof of extremality**

The explicit candidate is the rank-96 even unimodular lattice defined by
\[
L_A=\{(u,v)/4:u,v\in\mathbb Z^{48},\quad v\equiv uA\pmod{16}\},\qquad
\|(u,v)/4\|^2=(\|u\|^2+\|v\|^2)/16.
\]
The original candidate matrix in
`z16_candidate_checks/candidate_3cycle.json` is unchanged. The new
verifiers check its SHA-256 fingerprint and integral structural congruences.
There is a verified norm-10 lattice vector, but the global absence of
norm-8-or-shorter vectors is not yet proved.

## First new exact result: seven ±1 coordinates

A complete *Gray-code* sign enumeration over all
\(\binom{48}{7}=73,629,072\) seven-element supports and all
\(2^6=64\) sign vectors, with the first sign fixed by global negation,
covered **4,712,260,608** assignments.

- **95,377** assignments passed the *necessary* modulo-8 norm bound.
- **Zero** had squared norm at most 8 under the full exact modulo-16 test.
- This is a complete exclusion of support-seven vectors with all
  nonzero first-half integer coordinates equal to ±1.
- It does **not** exclude support-seven vectors with other magnitudes.

The independent C++ checker and script are
`support7_pm1.cpp` and `support7_pm1_verify.py`.

## Broader exact family: all seven odd

An all-odd support-seven vector can have any of the centered odd
coefficients \(\pm1,\pm3,\pm5,\pm7\). Fix a seven-element support
\(S=\{i_0<\dots<i_6\}\). Center the first-half integer coordinates
modulo 16, which cannot increase their norm, and negate the whole
vector so \(u_{i_0}>0\). Then every vector is represented uniquely by

- \(r_0,\ldots,r_6\in\{1,3,5,7\}\), the residues modulo 8;
- the first centered coefficient \(u_{i_0}=r_0\);
- for \(j=1,\ldots,6\), exactly one of the centered lifts
  \(u_{i_j}=r_j\) or \(r_j-8\).

Hence for every support the algorithm must examine
\(4^7=16,384\) modulo-8 patterns and, for each one that passes
the necessary lower bound, its at most \(2^6=64\) centered
modulo-16 lifts. The total *modulo-8 pattern count* across all
supports is **1,206,338,715,648**.

### Exact pruning lemma

Let
\[
c_8(r)=\min\{z^2:z\in\mathbb Z,\ z\equiv r\pmod8\}
=(0,1,4,9,16,9,4,1)[r].
\]
For a fixed support and residues \(r\), define the lower bound
\[
B(r)=r_0^2+\sum_{j=1}^6c_8(r_j)
 +\sum_{\ell=0}^{47}c_8\left(
    \sum_{j=0}^{6}r_jA_{i_j,\ell}\pmod8
 \right).
\]
Because the actual integer first-half coefficients reduce to
\(r_j\) modulo 8 and \(v\equiv uA\pmod{16}\), **every corresponding
lattice vector has squared-norm numerator at least \(B(r)\)**.
Therefore \(B(r)>128\) excludes the entire residue pattern.

The matrix satisfies \(A\bmod2=J-I\). With seven odd entries
in \(u\), precisely 41 entries of \(v\) are odd; consequently
\(\|v\|^2\ge41\), an additional safe lower bound used before exact
lift enumeration. Each surviving lift is independently checked using
48 exact squared norms of its centered residues modulo 16.

The implementation splits the seven row contributions into blocks
of three and four, caching their modulo-8 syndrome vectors and
partial squared-cost bounds. AVX2 only performs **integer** byte
additions and table lookups; no floating-point operations occur.

### Reproducibility

Requires Python 3.9+, a C++17/OpenMP compiler, and x86-64 AVX2.

```bash
python3 z16_candidate_checks/support7_scalar_reference.py \
  z16_candidate_checks/candidate_3cycle.json

python3 z16_candidate_checks/support7_pm1_verify.py \
  z16_candidate_checks/candidate_3cycle.json /tmp/pm1.json

# Scan one disjoint first-index partition; repeat FIRST=0,...,41.
OMP_NUM_THREADS=4 python3 z16_candidate_checks/support7_allodd_verify.py \
  z16_candidate_checks/candidate_3cycle.json 39 /tmp/odd7_first39.json
```

The fully completed first-index 39 partition has 28 supports,
458,752 modulo-8 assignments, one survivor, seven exact lifts,
and zero norm-8 vectors. The nontrivial support
\((39,41,43,44,45,46,47)\) independently has exactly one
modulo-8 survivor and seven admissible lifts in a portable
Python scalar implementation.

The `leech96-support7-odd.yml` GitHub Actions workflow runs all
42 disjoint partitions and writes a separate JSON certificate per
partition, along with the independent scalar and full ±1 checks.
**Do not claim exclusion of every all-odd support-seven vector until
all 42 partitions pass**. Any short-vector witness is to be
preserved and independently checked, not discarded.

## Remaining mathematical work

Even when all 42 partitions pass, this only excludes seven
nonzero **odd** first-half coordinates. Vectors with seven
mixed-parity supported coordinates, and vectors with eight or more
supported coordinates containing at least one odd entry, remain
unresolved. The previously certified all-even first-half sector
is separate.

The natural next improvement is a binary-code/syndrome
support-level prefilter (particularly for masks with even coefficients),
rather than naively multiplying the \(4^7\) residue search by all
128 parity masks.

**No extremal 96-dimensional lattice has been proved.**
