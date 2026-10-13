# Leech96: exhaustive six-coordinate parity-mask test

**October 13, 2026 — exact finite computation, not a proof of extremality.**

This update targets the explicit rank-96 even unimodular lattice supplied as
`z16_candidate_checks/candidate_3cycle.json`. The construction is

\[
L_A = \{(u,v)/4 : u,v\in\mathbb Z^{48},\ v\equiv uA\pmod{16}\},
\qquad \|(u,v)/4\|^2=(\|u\|^2+\|v\|^2)/16.
\]

The remaining extremality question is whether \(\|(u,v)/4\|^2\le8\) has a
nonzero solution. The existing exact certificate rules out **all-even \(u\)**.
Other completed finite searches rule out \(u\) of support at most five,
after centering modulo 16.

## Main result: every support-six parity mask is excluded

**Complete exact computation:** all 64 parity masks passed in
[GitHub Actions run #38424589542](https://github.com/HiroSakuraba/gamma72-leech96/actions/runs/38424589542),
together with a separate scalar reference-check job (**65/65 success**).
Thus **no lattice vector of squared norm at most 8 has exactly six nonzero
centered first-half coordinates**. This is an exact finite
computer-assisted exclusion, not a global shortest-vector proof.

For each of the 64 masks the checker exhausted all
`12,271,512 × 4096 = 50,264,113,152` modulo-8 residue patterns.
Across the full experiment, it covered **3,216,903,241,728** patterns,
tested every lifted candidate that survived the necessary lower bounds,
and found **zero** short vectors. Each mask produces a separate
machine-readable artifact in the Actions run.

Combined with the existing even-coordinate and supports-one-through-five
certificates, any as-yet-undetected vector of squared norm at most 8 has
**at least seven nonzero centered first-half coordinates**.

## New completed subfamilies

The October 13 finite calculations have now completed all configurations
of **exactly six nonzero centered coordinates in \(u\)** with:

- exactly **one odd** and five even coefficients, and
- exactly **five odd** and one even coefficient.

Each family has 6 possibilities for which ordered support position has
the minority parity. Every one of these **12** cases examined
\(\binom{48}{6}4^6=50,264,113,152\) modulo-8 coefficient classes.
The combined 12-case search examined **603,169,357,824** such classes,
retained **60,006** modulo-8 survivors, checked **469,248** centered
modulo-16 lifts exactly, and found **zero** vectors of squared norm at
most 8. See `z16_candidate_checks/support6_mixed_regression.json`
for the 12 separate counters.

The earlier **six-odd** family also passed. The two-, three-, and four-odd families were subsequently **completed**
in the full 64-mask run cited above.

## Generalized exact verifier

The added `support6_mask_exact.cpp` enumerates one **oddmask** in
`0,...,63` per invocation. Bit \(i\) of the mask records whether
the \(i\)-th coordinate of the increasing six-index support is odd.

For each of the \(\binom{48}{6}=12,271,512\) supports, each coordinate
has four possible residues modulo 8 with the required parity.
There are therefore exactly \(4^6=4096\) residue vectors per support.

For a vector of squared norm at most 8, the numerator must be at most
\(128\). We use the following safe rejection tests.

1. **Centering and global sign:** Each coordinate can be centered
   modulo 16 without increasing its squared norm. A nonzero supported
   coordinate congruent to zero modulo 8 must equal \(+8\) or \(-8\).
   Negate the entire vector if necessary so its first supported
   coordinate is **positive**. All remaining coordinates have exactly
   two centered lifts for each fixed residue modulo 8.

2. **Necessary modulo-8 bound:** Put
   \[
   c_8(r)=\min_{x\equiv r\pmod8}x^2
   = (0,1,4,9,16,9,4,1)[r].
   \]
   The first positive supported coefficient \(u_0\) is either
   \(r\) when \(r=1,\dots,7\), or \(8\) when \(r=0\). The other
   supported coefficients have minimum squared cost \(c_8(r)\) if
   \(r\ne0\), and **64** when \(r=0\). Since
   \(v\equiv uA\pmod8\), add the 48 values
   \(c_8((uA)_j)\). Reject only when the resulting bound is \(>128\).

3. **Parity lower bound:** The candidate satisfies
   \(A\bmod2=J-I\). If \(k\) coordinates of \(u\) are odd,
   then \(v\) has \(k\) odd coordinates for even \(k\), and
   \(48-k\) odd coordinates for odd \(k\). Thus
   \(\|v\|^2\ge (k\text{ if }k\text{ even else }48-k)\).
   This rejects further impossible lifted vectors.

4. **Exact lift check:** For every remaining residue assignment,
   exhaust all \(2^5=32\) centered lifts after fixing the first sign.
   Compute \(v=uA\bmod16\) using ordinary integer sums, center
   its 48 residues coordinatewise, and calculate its exact squared norm.
   No floating-point calculation, lattice-reduction heuristic, or
   statistical inference enters the exclusion.

The AVX2 instructions are used **only** to sum byte residues and
lookup exact squared-norm values modulo 8. They cannot round or
approximate anything. The Python verifier checks the exact candidate
matrix fingerprint and Type II congruences before compiling the C++ code,
and checks \(\binom{48}{6}4^6\) classes for the selected mask.

The separate `support6_mask_reference.py` computes five nontrivial
small subproblems **without vector instructions**, including masked
cases with one, two, three, four, and five odd coordinates. It agrees
with the optimized scanner on the number of modulo-8 survivors, number
of admissible lifts, and zero short vectors. This is a spot-check,
**not** a second independent exhaustive scan.

## Reproduction and full 64-mask check

```bash
python3 z16_candidate_checks/support6_mask_reference.py \
  z16_candidate_checks/candidate_3cycle.json

OMP_NUM_THREADS=4 python3 z16_candidate_checks/support6_mask_verify.py \
  z16_candidate_checks/candidate_3cycle.json 7 /tmp/mask7.json
```

Requirements: Python 3.9+, C++17, OpenMP, x86-64 CPU with AVX2.
Set the oddmask from 0 to 63 to cover every support-six parity type.
The GitHub Actions workflow `leech96-support6-all-masks.yml`
runs all 64 masks separately and uploads exact JSON reports.

**All 64 mask jobs succeeded** in the run cited above. A future
revision that reports a norm-8 witness must instead check membership
exactly and withdraw this candidate's extremality prospects.

## What remains even after all 64 masks pass

The all-even \(u\) case and small supports have separate exact
certificates. Excluding all six-support vectors would still leave
first-half support **seven or greater with at least one odd
coordinate**. No extremal 96-dimensional lattice is established
by this work. The general 96-dimensional existence problem stays open.
