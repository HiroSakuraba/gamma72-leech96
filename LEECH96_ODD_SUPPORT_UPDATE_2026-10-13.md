# Leech96 odd-coordinate exact search — October 13, 2026

The published deformed Z/16 candidate has
[
L_A=\{(u,v)/4 : u,v\in\mathbb Z^{48},\ v\equiv uA\pmod{16}\},
\qquad \|(u,v)/4\|^2=(\|u\|^2+\|v\|^2)/16.
]

The earlier [even-u certificate](LEECH96_RESEARCH_UPDATE_2026-10-13.md)
excludes norm 2, 4, 6, or 8 whenever every coordinate of \(u\) is even.
The new test excludes these norms for **any u with exactly four nonzero
coordinates and at least one odd coordinate**. That is an exact additional
subcase, **not a global minimum-norm proof**.

## Why the scan is complete for support four

1. Replace each integer \(u_i\) by the least-absolute-value representative
   modulo 16, hence \(-8\le u_i\le8\). This preserves lattice membership
   after replacing \(v\) by its centered representative, without increasing
   norm. For a nonzero first coordinate, negate the whole vector if
   necessary, so the first nonzero \(u_i\) lies in \(1,\dots,8\).
2. Choose **all** four-element supports from the 48 first-half coordinates.
   Enumerate the first coefficient in \(1,\dots,8\) and each remaining
   coefficient in \(\{-8,\dots,-1,1,\dots,8\}\).
3. For \(k=\#\{i:u_i\text{ odd}\}\), the candidate satisfies
   \(A\bmod2=J-I\). Therefore if \(k\) is even,
   \(v\bmod2=u\bmod2\); if \(k\) is odd,
   \(v\bmod2=\mathbf1+u\bmod2\). Consequently
   \(\|v\|^2\ge k\) for even \(k\), and
   \(\|v\|^2\ge48-k\) for odd \(k\).
   These are safe lower bounds; also discard partial coefficient choices
   whose accumulated \(\|u\|^2\) is already at least 128.
4. For everything remaining, compute all 48 centered residues of
   \(v=uA\pmod {16}\) in exact integer arithmetic, stop only
   if the running numerator exceeds 128.

The first-half u parity-zero sector is covered independently by the older
certificate, and skipped by the new test. Hence all support-four inputs are
excluded when the results are combined.

## Reproduction

```bash
python3 z16_candidate_checks/odd_support4.py \
 z16_candidate_checks/candidate_3cycle.json /tmp/odd_support4.json
```

Requirements: Python 3.9+, C++17 compiler, OpenMP. The script independently
checks the candidate matrix's Type II congruences and SHA-256 fingerprint,
compiles the C++ enumerator, and checks the output against regression
counts. Counts depend on when safe arithmetic pruning is applied;
they are not the number of unpruned formal assignments.

| Metric | Count |
|---|---:|
| Reached signed coefficient tuples after safe prefix pruning | 5,753,341,440 |
| Have at least one odd u coordinate | 5,435,786,880 |
| Pass the parity-based cost lower bound | 3,226,914,720 |
| Short vectors with norm <=8 | **0** |

An additional completed *bounded* exploration searched exactly five
nonzero first-half coordinates where the number of odd coordinates is
odd and every coefficient has absolute value at most three:
6,657,437,952 tuples (3,342,417,408 with odd parity), no norm <=8 witness.
Unlike the four-coordinate result, this five-coordinate exploration is
**not** an exhaustive five-coordinate exclusion. The fully general
odd-u case starts at first-half support five and is still open.

## Next research direction

The raw five-coordinate enumeration becomes expensive very quickly.
Rather than claim infeasibility from further small-support scans, apply
mod-4 and mod-8 necessary-code constraints at the *support level*,
or enumerate a well-reduced 96-dimensional Gram matrix by a certified
short-vector routine (and audit its completeness). Existing BKZ searches
that returned a norm-10 vector are evidence only.

**Status:** The candidate may have minimum 8 or 10. Neither existence
nor nonexistence of an extremal rank-96 even unimodular lattice has been
established here.
