"""Anatomy of a short vector y in L = P (x) Leech, given in ambient coords (4 blocks of 24).
Prints block norms (Leech norms) and the O-rank of the span of the 4 blocks (= tensor rank)."""
import sys, json
import numpy as np
from sympy import Matrix

D = sys.argv[1]
leech = {s["sign"]: s for s in json.load(open(f"{D}/leech23.json"))}

def analyse(witness, sign, p):
    st = leech[sign]
    H = Matrix(st["H"]); A = np.array(st["A"], dtype=np.int64)
    y = np.array(witness, dtype=np.int64).reshape(4, 24)
    blocks = [Matrix([list(map(int, b))]) * H.inv() for b in y]   # Leech coordinates
    blocks = [np.array(b.tolist()[0], dtype=object) for b in blocks]
    norms = [int(sum(int(v) ** 2 for v in b)) // 8 for b in y]
    # O-rank: Q-rank of {x_i, x_i A} divided by 2
    rows = []
    for b in blocks:
        b = np.array([int(v) for v in b], dtype=np.int64)
        rows.append(b); rows.append(b @ A)
    r = Matrix(rows).rank()
    return norms, r // 2, sum(norms) / p

if __name__ == "__main__":
    files = sys.argv[2:]
    for f in files:
        for rec in json.load(open(f)):
            w = rec.get("witness") or rec.get("witness_ambient")
            norms, rank, nm = analyse(w, rec.get("sign", rec.get("leech_sign")), rec["p"])
            print(f"P#{rec['P_index']} p={rec['p']}: norm {nm}  block Leech-norms {norms}  tensor rank {rank}")
