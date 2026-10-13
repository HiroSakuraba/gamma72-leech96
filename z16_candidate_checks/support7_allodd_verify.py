#!/usr/bin/env python3
"""Complete exact support-seven all-odd scan for one disjoint first-index partition.

Usage:
  python3 support7_allodd_verify.py candidate_3cycle.json FIRST output.json

Run FIRST=0,...,41 to cover every support of size 7 with ALL SEVEN
first-half entries odd. Their magnitudes are unrestricted except by
centering modulo 16. Mixed parity and support >=8 remain unresolved.
Requires Python 3.9+, C++17, OpenMP, and x86-64 AVX2.
"""
import hashlib, json, os, platform, re, subprocess, sys, tempfile
from math import comb
from pathlib import Path

MATRIX_HASH="929e71c3a3f54f3361d9a8c318786b1da726e141250a93969e9b2cec8cd3c656"

def candidate_matrix(doc):
    A=doc["A"]
    assert len(A)==48 and all(len(row)==48 for row in A)
    assert all(type(x) is int and -128<=x<=127 for row in A for x in row)
    blob=b"".join(x.to_bytes(8,"little",signed=True) for row in A for x in row)
    sha=hashlib.sha256(blob).hexdigest()
    assert sha==MATRIX_HASH==doc["A_sha256"]
    for i in range(48):
        assert (1+sum(x*x for x in A[i]))%32==0
        for j in range(48):
            assert A[i][j]%2==int(i!=j)
        for j in range(i+1):
            assert (sum(A[i][t]*A[j][t] for t in range(48))+int(i==j))%16==0
    return A,sha

def main():
    if len(sys.argv)!=4: sys.exit(__doc__)
    candidate=json.loads(Path(sys.argv[1]).read_text())
    A,sha=candidate_matrix(candidate)
    first=int(sys.argv[2])
    if first<0 or first>41:
        raise ValueError("FIRST must be 0..41")
    if platform.machine().lower() not in ("x86_64","amd64"):
        raise RuntimeError("This optimized checker requires x86-64 AVX2")
    source=Path(__file__).with_name("support7_allodd.cpp")
    with tempfile.TemporaryDirectory(prefix="leech96_odd7_") as t:
        tmp=Path(t)
        header="#include <cstdint>\nstatic const int8_t A_pad[48][64]={\n"
        header+=",\n".join("{"+",".join(map(str,row+[0]*16))+"}" for row in A)
        header+="\n};\n"
        (tmp/"matrix48.h").write_text(header)
        exe=tmp/"scan"
        subprocess.run([os.environ.get("CXX","g++"),"-O3","-std=c++17",
                        "-mavx2","-fopenmp","-I",str(tmp),
                        str(source),"-o",str(exe)],check=True)
        p=subprocess.run([str(exe),str(first),str(first+1)],
                         capture_output=True,text=True,
                         env={**os.environ,"OMP_NUM_THREADS":
                              os.environ.get("OMP_NUM_THREADS","4")})
    print(p.stdout,flush=True)
    m=re.search(r"SUMMARY first=\[(\d+),(\d+)\) supports=(\d+) "
                r"mod8_patterns=(\d+) passed_mod8=(\d+) lifts=(\d+) "
                r"short=(\d+) best_numerator=(\d+)",p.stdout)
    assert m,(p.returncode,p.stdout,p.stderr)
    start,end,supports,patterns,passed,lifts,short,best=map(int,m.groups())
    assert (start,end)==(first,first+1)
    assert supports==comb(47-first,6)
    assert patterns==supports*(4**7)
    assert lifts<=64*passed
    assert p.returncode==(3 if short else 0),(p.returncode,p.stderr)
    if first==39:
        assert passed==1 and lifts==7  # independently scalar checked
    report={
        "candidate_A_sha256":sha,
        "first_index":first,
        "scope":"seven nonzero centred first-half integer coordinates, all seven odd",
        "supports":supports,
        "mod8_residue_patterns":patterns,
        "mod8_survivors":passed,
        "centered_mod16_lifts_checked":lifts,
        "short_vectors_norm_sq_le_8":short,
        "complete_for_first_index":True,
        "full_minimum_norm_proved":False,
    }
    Path(sys.argv[3]).write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2),flush=True)
    if short:
        raise SystemExit("Explicit norm-8-or-shorter witness found: inspect stdout")

if __name__=="__main__":main()
