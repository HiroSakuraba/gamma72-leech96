# E8 ⊗ J24 over Z[(1+√−7)/2] — first pass (Oct 12, 2026)

Target: an extremal 96-dim lattice as J ⊗ E8, where J is an extremal 7-modular 24-dim lattice
(minimum 10) with Hermitian Z[α]-structure, α = (1+√−7)/2. Class number 1, split norm 5 × 2 = 10.

- `e8_seven.py` → `e8_seven.json`: two Z[α]-structures on E8 (cyclic Hamming code coordinates; √−7 =
  Gauss sum over the order-7 shift + rotation on the 2-dim fixed space). All identities checked exactly.
- `leech_pb_e8.py` → `leech_o7.json`: Barnes P_b ⊗ E8 over Z[α] is the Leech lattice (det 1, even,
  minimum 4, 196,560 minimal vectors), giving Leech with an explicit √−7.
- `j24_search.py`: J candidates M_C = π⁻¹(C) ⊂ Leech for Lagrangian C in Leech/√−7·Leech ≅ F7^12
  (symplectic form b(zS, m) mod 7). Each M_C is 7-modular of det 7^12 (checked). Score = exact count of
  vector pairs of norm ≤ 8 (0 means extremal).
- Results: random Lagrangians have ≈ 1,700 short pairs; hill-climbing with symplectic transvections
  (≈ 4,000 moves over two runs) only reached 1,705 (`run11.log`, `run22.log`). Fully symmetric choices
  split as P_b ⊗ (7-modular sublattice of E8) and have minimum ≤ 8 by Hermite.
- Extremal theta series of J checked to q^199: 1 + 84672 q^5 + 91728 q^6 + …, all positive integers.

Open: whether an extremal 7-modular J exists at all (Nebe's table entry for dimension 24 did not
extract cleanly). Next ideas: Lagrangians invariant under large non-split subgroups of the unitary
group of P_b ⊗ E8; or J outside Leech via Kneser neighbours in the genus; then the rank-2–4 tensor
checks for E8 ⊗ J, which the Coulangeon inequality alone does not settle.
