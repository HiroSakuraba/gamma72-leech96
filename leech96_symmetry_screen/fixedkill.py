"""Explicit short vectors in the fixed sublattice of the order-47 shift for every symmetric Paley lattice.
Coordinates: index 0 = infinity, 1..47 = F_47 (as in paley_sym.py). Fixed vectors:
  u = (b, a,...,a),  v = (d, e,...,e).  With A = mC + cI:  uA = (-47 m a + c b, (m b + c a) 1_47).
Membership: d = -47 m a + c b (mod N), e = m b + c a (mod N).  Norm = (b^2 + 47 a^2 + d^2 + 47 e^2)/N.
Each witness is re-checked against the full 96-dim congruence v = uA (mod N)."""
import sys, json, itertools
import numpy as np
sys.argv_saved = list(sys.argv)

q = 47
QR = {(x * x) % q for x in range(1, q)}
chi = lambda t: 0 if t % q == 0 else (1 if t % q in QR else -1)
Qm = np.array([[chi(j - i) for j in range(q)] for i in range(q)], dtype=np.int64)
C = np.zeros((48, 48), dtype=np.int64); C[0, 1:] = 1; C[1:, 0] = -1; C[1:, 1:] = Qm
I = np.eye(48, dtype=np.int64)
res = json.load(open(sys.argv[1]))
out = []
for r in res:
    N, m, c = r["N"], r["m"], r["c"]
    A = m * C + c * I
    best = None
    R = 20
    for a in range(-3, 4):
        for b in range(-R, R + 1):
            if a == 0 and b == 0:
                pass
            d0 = (-47 * m * a + c * b) % N
            e0 = (m * b + c * a) % N
            for d in (d0, d0 - N):
                for e in (e0, e0 - N):
                    if a == b == d == e == 0:
                        continue
                    n2 = b * b + 47 * a * a + d * d + 47 * e * e
                    if n2 % N:
                        continue
                    if best is None or n2 < best[0]:
                        best = (n2, a, b, d, e)
    n2, a, b, d, e = best
    u = np.array([b] + [a] * 47); v = np.array([d] + [e] * 47)
    member = ((v - u @ A) % N == 0).all()
    norm = (u @ u + v @ v) / N
    out.append(dict(N=N, m=m, c=c, bkz_shortest=r["shortest"], fixed_witness_u=u.tolist(), fixed_witness_v=v.tolist(),
                    norm=norm, member=bool(member)))
    print(f"N={N:2d} m={m:3d} c={c:3d}: BKZ said {r['shortest']:5}; fixed-sublattice vector norm {norm} "
          f"(a={a}, b={b}, d={d}, e={e}) member={member}")
print("all symmetric lattices killed:", all(o["member"] and o["norm"] <= 8 for o in out))
json.dump(out, open(sys.argv[2], "w"))
