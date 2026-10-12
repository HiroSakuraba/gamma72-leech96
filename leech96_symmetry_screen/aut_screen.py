"""Which prime-order automorphism types are impossible for an extremal (min 10) even unimodular 96-lattice?
sigma of prime order p, fixed lattice F of dim dF = 96 - k(p-1), k = number of Z[zeta_p]-blocks.
Standard lemma (unimodular L): L/(F + Z) = F_p^s, det F = p^s, s <= min(dF, k).  F is even with min >= 10.
Hermite: min(F) <= gamma_dF * det(F)^(1/dF).  gamma exact for dF <= 8 and dF = 24; Blichfeldt bound otherwise:
gamma_n <= (2/pi) * Gamma(2 + n/2)^(2/n).  Type excluded if the bound is < 10 for every admissible s."""
import math
exact = {1: 1, 2: (4/3)**0.5, 3: 2**(1/3), 4: 2**0.5, 5: 8**0.2, 6: (64/3)**(1/6), 7: 64**(1/7), 8: 2.0, 24: 4.0}
def gamma_ub(n):
    if n in exact: return exact[n]
    return (2 / math.pi) * math.exp(math.lgamma(2 + n / 2) * 2 / n)
def primes(n):
    return [p for p in range(2, n + 1) if all(p % d for d in range(2, int(p**0.5) + 1))]
rows = []
for p in primes(97):
    for k in range(1, 96 // (p - 1) + 1):
        dF = 96 - k * (p - 1)
        if dF == 0:
            rows.append((p, k, dF, "fixed-point-free: no constraint from F")); continue
        smax = min(dF, k)
        worst = max(gamma_ub(dF) * p ** (s / dF) for s in range(0, smax + 1)
                    if not (s == 0 and dF % 8 != 0))
        status = "EXCLUDED" if worst < 10 else "open"
        rows.append((p, k, dF, f"max Hermite bound on min(F) = {worst:.2f} -> {status}"))
for r in rows:
    if r[0] >= 11:
        print(f"p={r[0]:2d}  blocks k={r[1]:2d}  dim F={r[2]:2d}  {r[3]}")
