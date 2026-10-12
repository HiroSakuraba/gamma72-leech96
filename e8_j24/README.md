# E8 ⊗ J24 over Z[(1+√−7)/2]

Target: an extremal 96-dim lattice as E8 ⊗ J, where J is an extremal 7-modular 24-dim lattice
(minimum 10) with Hermitian Z[α]-structure, α = (1+√−7)/2. Class number 1, split norm 5 × 2 = 10.
Status (Oct 12, 2026): **no J found.** Four search families come up empty; the route stays open.

## Building blocks (exact)

- `e8_seven.py` → `e8_seven.json`: two Z[α]-structures on E8 (cyclic Hamming code coordinates; √−7 =
  Gauss sum over the order-7 shift + rotation on the 2-dim fixed space). All identities checked exactly.
- `leech_pb_e8.py` → `leech_o7.json`: Barnes P_b ⊗ E8 over Z[α] is the Leech lattice (det 1, even,
  minimum 4, 196,560 minimal vectors), giving Leech with an explicit √−7.
- Extremal theta series of J checked to q^199: 1 + 84672 q^5 + 91728 q^6 + …, all positive integers
  (`../leech96_symmetry_screen/theta7.py`).

## Searches for J

Score = exact count of vector pairs of norm ≤ 8 (0 would mean extremal).

| Family | Script / output | Result |
|---|---|---|
| M_C = π⁻¹(C) ⊂ Leech, C random Lagrangian in Leech/√−7·Leech ≅ F7^12 | `j24_search.py`, `run11.log`, `run22.log` | ≈ 1,700 short pairs; transvection hill-climb (≈ 4,000 moves) best 1,705 |
| Same, C invariant under an order-3 symmetry (9 usable ones from 48 signed perms of P_b times x ↦ 2x on E8) | `j24_symmetric.py`, `j24_sym_best.json` | best 1,711, minimum 6 |
| Kneser 2-neighbours in the Z[α]-genus, starting from the best M_C | `j24_neighbours.py` | neighbours valid (even, det 7^12, α-stable) but scores 1,711–1,720; no progress |
| Ideal lattices (I, Tr(λ x ȳ)) in Q(ζ_m), φ(m) = 24, 7 \| m, with λ·𝔇·I·Ī = (√−7), λ ≫ 0 — **exhaustive** over ideal classes and totally positive units mod squares | `j24_ideal.py`, `ideal{35,56,84}.log/json` | m = 35: no lattice exists (no totally positive λ; real-unit signs have F2-rank 11 of 12). m = 56: 4 lattices, min 4 (5,656 pairs) or min 8 (2,016 pairs). m = 84: 2 lattices, min 8 (2,016 pairs). None extremal. |

So an extremal J cannot carry a fixed-point-free automorphism of order 35, 56, 70 or 84 compatible with
its √−7 (Q(ζ70) = Q(ζ35)), for this rank-1 ideal-lattice family. Computed with PARI (cypari2);
class groups and units come from `bnfinit(…, 1)`, which assumes GRH (`bnfcertify` was not run).

## Next ideas

- Hermitian lattices over Z[ζ_m] of rank > 1 for smaller m with 7 | m (m = 7, 14, 21, 28, 42): a J with a
  fixed-point-free automorphism of order 7 or 21 is still possible.
- Mass formula / full classification of the Z[α]-genus (rank 12, det 7^12) is likely too big, but a
  neighbour walk that maximises a different score (e.g. kissing number at norm 8 rather than ≤ 8) may
  move further than the current hill-climb.
- If a J appears: build E8 ⊗ J and test rank-2 to rank-4 tensors (Coulangeon's inequality alone does
  not settle the minimum).

Requires `pip install cypari2` for `j24_ideal.py`. Run: `python3 j24_ideal.py 84 out.json`.
