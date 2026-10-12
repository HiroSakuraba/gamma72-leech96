"""Exact certificate: no vector of norm <= 8 in the Z16 candidate with code coordinates u = 0 mod 4.

L = {(u,v)/4 : v = uA mod 16},  norm = (|u|^2+|v|^2)/16.
If u = 4w then v = 4t with t = wA (mod 4), norm = |w|^2 + |t|^2, minimised coordinatewise by centred
representatives: cost(0)=0, cost(+-1)=1, cost(2)=4.  Write x = w mod 2.
 * x = 0  (u = 0 mod 8): w = 2y, t = 2(yA mod 2): cost = 4(|y| + |yA mod 2|)       [checked below]
 * x != 0: odd entries of w cost |x|; odd entries of t are x*(A mod 2).
Everything is enumerated exactly from A itself; the structure of A mod 2 is checked, not assumed."""
import sys, json, itertools
import numpy as np

cand = json.load(open(sys.argv[1]))
A = np.array(cand["A"], dtype=np.int64)
n = 48
A2 = A % 2
A4 = A % 4
J = np.ones((n, n), dtype=np.int64); I = np.eye(n, dtype=np.int64)
print("A mod 2 == J - I :", (A2 == (J - I)).all())
assert (A2 == (J - I)).all()

def cost4(vec):
    r = vec % 4
    return int((r == 1).sum() + (r == 3).sum() + 4 * (r == 2).sum())

# x = 0 case: binary y; cost = 4(|y| + wt(y (J-I) mod 2)); y(J-I) = y if |y| even else complement
best0 = min(4 * (k + (k if k % 2 == 0 else n - k)) for k in range(1, n + 1))
print("layer u=0 mod 8 (and u=0 with v=0 mod 16 gives 16): minimum cost", best0)

# x != 0: parity of t = x (J-I) mod 2 ; |x| odd -> t odd off supp(x): cost >= 48
# |x| even -> t odd exactly on supp(x): cost >= 2|x|; need <= 8 -> |x| in {2,4}
found = []
best = 10 ** 9
# |x| = 4 : w entries on support in {+-1} (mod 4 reps), no entries = 2 allowed (budget exhausted)
for S in itertools.combinations(range(n), 4):
    for signs in itertools.product([1, -1], repeat=3):
        s = (1,) + signs
        t = sum(si * A4[i] for si, i in zip(s, S)) % 4
        w = np.zeros(n, dtype=np.int64); w[list(S)] = s
        c = cost4(w) + cost4(t)
        best = min(best, c)
        if c <= 8:
            found.append((list(map(int, w)), c))
# |x| = 2 : w = signed pair + optional single 2 elsewhere
for S in itertools.combinations(range(n), 2):
    for sgn in (1, -1):
        w0 = np.zeros(n, dtype=np.int64); w0[S[0]] = 1; w0[S[1]] = sgn
        base_t = (A4[S[0]] + sgn * A4[S[1]]) % 4
        for k in [None] + [k for k in range(n) if k not in S]:
            w = w0.copy(); t = base_t.copy()
            if k is not None:
                w[k] = 2; t = (t + 2 * A4[k]) % 4
            c = cost4(w) + cost4(t)
            best = min(best, c)
            if c <= 8:
                found.append((list(map(int, w)), c))
print("layer u=0 mod 4, u!=0 mod 8: minimum cost over the exhaustive (necessary) cases with norm<=8 possible:",
      best, "; vectors of norm <= 8 found:", len(found))
json.dump(dict(A_mod2_is_J_minus_I=True, layer8_min=best0, layer4_min_in_enumeration=best,
               layer4_short_vectors=found[:20], certified_no_norm_le_8_in_layers=len(found) == 0 and best0 > 8),
          open(sys.argv[2], "w"), indent=1)
