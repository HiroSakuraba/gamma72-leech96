"""Realise target Gram patterns inside Gamma72 and produce checked certificates.
Usage: realize.py WORKDIR gamma72.json targets.json out.json
targets.json: list of [k, K] pairs (tried in order per k; first realised + saturated wins)."""
import sys, json, time
import numpy as np
sys.path.insert(0, sys.argv[1])
from antipode import load, gram_of, is_saturated, best_anticodes, exact_delta

G = load(sys.argv[2])
targets = json.load(open(sys.argv[3]))
out_path = sys.argv[4]
n = 72
E = np.eye(n, dtype=np.int64)

# pool: basis vectors, 2-term and selected 3-term combinations of norm <= 12, both signs
vecs = {}
def add(v):
    nm = int(v @ G @ v)
    if nm <= 12:
        key = v.tobytes()
        if key not in vecs:
            vecs[key] = v
for i in range(n):
    add(E[i])
for i in range(n):
    for j in range(i + 1, n):
        for s in (1, -1):
            add(E[i] + s * E[j])
two = list(vecs.values())
# 3-term: basis vector + norm-8 two-term vector where it lowers the norm
base8 = [v for v in two if int(v @ G @ v) == 8 and np.abs(v).sum() == 2]
for v in base8:
    gv = G @ v
    for i in range(n):
        if v[i] != 0:
            continue
        for s in (1, -1):
            w = v + s * E[i]
            add(w)
pool = np.array(list(vecs.values()))
pool = np.concatenate([pool, -pool])
norms = np.einsum("ij,jk,ik->i", pool, G, pool)
IP = pool @ G @ pool.T
print("pool size", len(pool), {int(m): int((norms == m).sum()) for m in sorted(set(norms.tolist()))}, flush=True)

def realise(T, max_nodes=200000):
    k = len(T)
    T = np.array(T)
    nodes = 0
    def rec(chosen):
        nonlocal nodes
        nodes += 1
        if nodes > max_nodes:
            return None
        i = len(chosen)
        if i == k:
            return chosen
        mask = norms == T[i, i]
        for j, c in enumerate(chosen):
            mask &= IP[c] == T[i, j]
        cands = np.where(mask)[0]
        if i == 0:
            cands = cands[: 50]
        for c in cands[:400]:
            if c in chosen:
                continue
            r = rec(chosen + [int(c)])
            if r is not None:
                return r
        return None
    return rec([])

results = {}
for k, T in targets:
    if str(k) in results:
        continue
    t0 = time.time()
    idx = realise(T)
    if idx is None:
        print(f"k={k} pattern not realised in pool: {T}", flush=True)
        continue
    C = pool[idx]
    Kg = gram_of(G, C)
    assert (Kg == np.array(T)).all()
    sat, ed = is_saturated(C)
    if not sat:
        print(f"k={k} realised but NOT saturated (elementary divisors {ed}); skipping", flush=True)
        continue
    res, Kinv = best_anticodes(Kg, beta_max=0.45, max_pts=800)
    val, s, beta, S = res[0]
    bF, s2, detK = exact_delta(Kg, S)
    results[str(k)] = dict(dim=72 - k, log2delta=val, s=s, beta=str(bF), detK=detK, K=Kg.tolist(),
                           C=C.tolist(), S=[list(map(int, p)) for p in S])
    print(f"k={k} dim={72-k}: realised & saturated. log2delta={val:.5f} s={s} beta={bF} detK={detK} "
          f"({time.time()-t0:.1f}s)", flush=True)
    json.dump(results, open(out_path, "w"))
print("done")
