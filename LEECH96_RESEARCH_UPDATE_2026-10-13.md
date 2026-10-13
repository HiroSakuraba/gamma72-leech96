# Leech96 exact-certificate update — October 13, 2026

This update covers all three research priorities, using the published
\`z16_candidate_checks/candidate_3cycle.json\` as input. **It does not prove
an extremal rank-96 lattice exists.** The full shortest-vector question
remains open for first-half vectors with odd integer coordinates.

## 1. Even-coordinate short-vector exclusion

The construction is

\[
L_A=\{(u,v)/4:\ u,v\in\mathbb Z^{48},\
             v\equiv uA\pmod {16}\},\quad
\|(u,v)/4\|^2=(\|u\|^2+\|v\|^2)/16.
\]

The existing \`layer4.py\` independently excludes norm ≤8 for
\`u ≡ 0 (mod 4)\`. The new code considers \`u=2w\`, \`v=2t\`:
\`t ≡ wA (mod 8)\` and norm \`(||w||²+||t||²)/4\`.

### Binary reduction and its exact enumeration

Define the binary matrix \`M\` by \`M_ij=1\` if
\`A_ij ≡ -1 (mod 4)\`, and zero otherwise. It satisfies
\`M+M^T=J-I (mod 2)\`. For \`x=w mod 2\`, let its support be \`S\`
and its Hamming weight \`k\`.

* For odd \`k\`, the parity supports of \`w\` and \`t\` partition
  48 coordinates. Hence \`||w||²+||t||²≥48>32\`: impossible.
* For even \`k\`, the parity supports coincide, implying a
  baseline squared-coordinate cost of at least \`2k\`. Thus
  \`k≤16\` for norm ≤8.
* Set \`y=xM\`, \`z=xM^T\`; for even \`k\`, \`y+z=x\`.
  Outside \`S\`, a direct modulo-4 computation yields a **single**
  \`b∈{0,1}\` such that
  \[
  \|w\|^2+\|t\|^2 \ge 2k+4m,\qquad
  m=\#\{j\notin S:y_j\ne b\}.
  \]
  Therefore \`k+2m≤16\` is necessary.

The program enumerates **all 2²⁴ words** in each of the row-space codes
\`im(M)\`, \`im(M^T)\` and the affine coset \`1+im(M)\` (the last
shares the first code's basis). Their least **positive even**
weights are each 12. The all-ones word belongs to \`im(M^T)\`,
but not to \`im(M)\`.

If \`b=0\`, \`wt(y)+wt(z)=k+2m\`. If both are nonzero,
the sum is at least 24, impossible. Hence \`x\` must lie in
\`im(M)\` or \`im(M^T)\`. If \`b=1\`, apply the same
argument to \`1+y\` and \`1+z\`; then \`x\`
must lie in \`1+im(M)\`. Both alternatives require an even code
weight of **12 or 16** (none of the three has even weight 14).

| Family | Weight 12 | Weight 16 |
|---|---:|---:|
| im(M) | 12,972 | 356,730 |
| im(M^T) | 17,296 | 535,095 |
| 1 + im(M) | 4,324 | 178,365 |

### Finite modulo-8 checks

With 16 odd positions, the baseline cost is 32; \`w\` must consist
entirely of ±1 on those positions and zeros elsewhere. Across the
1,070,190 supports, a necessary modular consistency system leaves
only **26** cases (7, 16, 3 in the three families). Exhaustive
\`2^16\` sign search produces **no norm ≤8 vector**.

With 12 odd positions, the baseline cost is 24, leaving a budget
of 8. Every admissible \`w\` is one of:

1. ±1 on its 12 support coordinates and zero elsewhere;
2. one additional ±2 at a support coordinate (making it ±3);
3. one additional ±2 outside the support.

Two outside ±2 modifications are impossible: on the two modified
coordinates, the modulo-4 equations force an additional total
squared norm cost of at least 8 in \`t\`, atop their cost 8 in
\`w\`. Other modifications already exceed the budget.

All **3,355,424** candidate modifications, including the 34,592
unmodified cases, fail a *necessary* modulo-8 linear system. The
linear system is solved over F₂ for each possible sign sum modulo 4.
It deliberately relaxes the actual sign-sum consistency condition:
this can create false positives but cannot reject a real solution.
Therefore no additional sign testing is needed for weight 12.

**Result:** with the preexisting mod-4 layer certificate, no vector
of norm 2, 4, 6 or 8 with even \`u\` exists. This is an
exact computer-assisted exclusion for one explicit candidate \`A\`,
not a full shortest-vector certificate.

Reproduce in a checkout with Python 3.9+ and a C++17 compiler:

\`\`\`bash
python3 z16_candidate_checks/even_layer_certificate.py \
  z16_candidate_checks/candidate_3cycle.json /tmp/even_layer.json
python3 z16_candidate_checks/layer4.py \
  z16_candidate_checks/candidate_3cycle.json /tmp/layer4.json
\`\`\`

The new verifier checks the Type II matrix congruences before
enumeration, uses Python integer arithmetic to enumerate the binary
codes, generates an integer matrix header, and compiles the standalone
C++ exact scanner. It never modifies the candidate.

## 2. 7-modular J24: a cyclotomic family excluded without GRH

Quebbemann's standard modular-form description implies that any
hypothetical *extremal* even 7-modular rank-24 lattice has a
norm-10 shell of **84,672** vectors. Since

\[
84672\equiv2\pmod5,\qquad84672\equiv3\pmod {13},
\]

such a lattice admits **no fixed-point-free isometry of order 5 or
13**: a prime-order group acting freely would partition the
norm-10 vectors into orbits of that length.

In particular, the rank-one cyclotomic Hermitian ideal-lattice
families over \`Q(zeta_35)\` and \`Q(zeta_70)\` are excluded,
because multiplication by \`zeta_5\` is a fixed-point-free
unitary isometry. Unlike earlier ideal class calculations, this
result does **not** rely on a generalized Riemann hypothesis
assumption. It does not exclude 7-modular lattices with an
order-5 isometry having a nontrivial fixed lattice, nor the
rank-one families over \`Q(zeta_56)\` and \`Q(zeta_84)\`.

The remaining plausible investigations are rank-greater-than-one
Hermitian constructions over smaller cyclotomic fields, and
new searches within the 7-modular genus, with explicit candidate
verification. No extremal \`J24\` has been found.

## 3. Stronger rank-96 automorphism screen

The conditional extremal theta series in dimension 96 gives
**565,866,362,880** norm-10 vectors, and

\[
565866362880\equiv12\pmod {97}.
\]

An order-97 isometry in dimension 96 is necessarily
fixed-point-free: its nontrivial rational irreducible
representation has dimension \`phi(97)=96\`. Its orbits on the
nonzero norm-10 shell would all have 97 elements, a contradiction.

**Therefore no extremal rank-96 even unimodular lattice admits
an automorphism of order 97.** This closes the last 97 case
that the preexisting fixed-lattice Hermite-bound screen
necessarily left open. It does not independently verify all
other exclusion claims in \`aut_screen.py\`; those still rely
on the prime-order discriminant lemma and the Hermite/Blichfeldt
constants used there.

Reproduce *both* new modular-form computations using only exact
Python integers:

\`\`\`bash
python3 leech96_symmetry_screen/theta_orbit_cert.py /tmp/orbit.json
\`\`\`

## Mathematical boundaries

* The theta coefficients are conditional on **existence** of the
  proposed extremal lattice. Their computation does not prove
  existence.
* The computer-aided even-\`u\` exclusion assumes the existing
  separate \`layer4.py\` finite certificate for \`u ≡ 0 mod 4\`.
* The shorter odd-\`u\` cases are **not excluded**.
* These are independently rerunnable exact finite calculations;
  publication should include a second implementation or formal
  checking of the proof of the finite case split.

**Modular-form sources:** H.-G. Quebbemann, *Modular lattices in
Euclidean spaces*, J. Number Theory 54 (1995), 190–202,
[DOI 10.1006/jnth.1995.1111](https://doi.org/10.1006/jnth.1995.1111);
*Atkin-Lehner eigenforms and strongly modular lattices*,
L'Enseignement Mathématique 43 (1997), 55–65.
