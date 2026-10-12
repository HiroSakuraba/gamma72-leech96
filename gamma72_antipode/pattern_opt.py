"""Abstract optimisation of the cross-section Gram pattern K (rank k), ignoring Gamma72 for now:
maximise the antipode center density  s * detK^{-1/2} * ((8-beta)/4)^{(72-k)/2}
over even integer Gram matrices with diagonal >= 8 (needed: K sits in a min-8 lattice).
Note: K must also have min >= 8, checked with Fincke-Pohst."""
import sys, json, math, random, time
import numpy as np
sys.path.insert(0, sys.argv[1])
from antipode import best_anticodes, fincke_pohst

k = int(sys.argv[2]); iters = int(sys.argv[3]); seed = int(sys.argv[4]); out = sys.argv[5]
rng = random.Random(seed)

def minnorm_ok(K):
    pts = fincke_pohst(K.astype(float), 7.5)
    return len(pts) == 0

cache = {}
def score(K):
    key = K.tobytes()
    if key in cache:
        return cache[key]
    if np.linalg.eigvalsh(K.astype(float)).min() <= 1e-9 or not minnorm_ok(K):
        cache[key] = (-1e9, None)
        return cache[key]
    res, _ = best_anticodes(K, beta_max=0.45, max_pts=500)
    cache[key] = (res[0][0], res[0][:3]) if res else (-1e9, None)
    return cache[key]

def random_K():
    K = np.zeros((k, k), dtype=np.int64)
    for i in range(k):
        K[i, i] = rng.choice([8, 8, 8, 10])
        for j in range(i + 1, k):
            K[i, j] = K[j, i] = rng.choice([-3, -2, -2, -1, -1, 0])
    return K

best_all = []
t0 = time.time()
for restart in range(iters):
    K = random_K()
    cur, info = score(K)
    tries = 0
    while tries < 60:
        K2 = K.copy()
        i, j = rng.sample(range(k), 2)
        if rng.random() < 0.15:
            K2[i, i] = rng.choice([8, 10, 12])
        else:
            K2[i, j] = K2[j, i] = K[i, j] + rng.choice([-1, 1])
            if abs(K2[i, j]) > 4:
                continue
        s2, info2 = score(K2)
        if s2 > cur:
            K, cur, info = K2, s2, info2
            tries = 0
        else:
            tries += 1
    best_all.append((cur, K.tolist(), [float(x) if not isinstance(x, (int,)) else x for x in info] if info else None))
    best_all.sort(key=lambda t: -t[0])
    best_all = best_all[:15]
    if restart % 5 == 0:
        print(f"k={k} restart={restart} best={best_all[0][0]:.5f} t={time.time()-t0:.0f}s", flush=True)
        json.dump(best_all, open(out, "w"))
json.dump(best_all, open(out, "w"))
print("final", best_all[0], flush=True)
