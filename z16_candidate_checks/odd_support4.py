#!/usr/bin/env python3
"""Exact sparse odd-u exclusion for the published rank-96 Z/16 candidate.

Usage:
  python3 odd_support4.py candidate_3cycle.json report.json [maxabs]

The default maxabs=8 checks every possible vector with exactly four nonzero
centred first-half coordinates and at least one odd coordinate, subject to
squared lattice norm <=8. No claim is made about support >=5.
Python 3.9+ and a C++17 compiler with OpenMP are required.
"""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

def verify_matrix(A):
    n = 48
    assert len(A)==n and all(len(row)==n for row in A)
    assert all(type(x)==int and -8<=x<=8 for row in A for x in row)
    assert all(A[i][j]%2==int(i!=j) for i in range(n) for j in range(n))
    for i in range(n):
        assert (1+sum(x*x for x in A[i]))%32==0
        for j in range(i+1):
            assert (sum(A[i][t]*A[j][t] for t in range(n))+int(i==j))%16==0

def main():
    if len(sys.argv) not in (3,4):
        sys.exit(__doc__)
    candidate=json.loads(Path(sys.argv[1]).read_text())
    A=candidate["A"]
    verify_matrix(A)
    limit=int(sys.argv[3]) if len(sys.argv)==4 else 8
    assert 1<=limit<=8
    digest=hashlib.sha256(
        b"".join(int(x).to_bytes(8,"little",signed=True) for row in A for x in row)
    ).hexdigest()
    if "A_sha256" in candidate:
        assert digest==candidate["A_sha256"]
    header=("#include <cstdint>\nstatic const int8_t A[48][48]={\n"
        + ",\n".join("{"+",".join(map(str,row))+"}" for row in A)
        + "\n};\n")
    cpp=Path(__file__).with_name("odd_support4.cpp")
    with tempfile.TemporaryDirectory(prefix="leech96_odd4_") as tmp:
        tmp=Path(tmp)
        (tmp/"matrix48.h").write_text(header)
        exe=tmp/"scan"
        subprocess.run([os.environ.get("CXX","g++"),"-O3","-std=c++17",
            "-fopenmp",str(cpp),"-I",str(tmp),"-o",str(exe)],check=True)
        run=subprocess.run([str(exe),"0","45",str(limit)],text=True,
            capture_output=True,env={**os.environ,"OMP_NUM_THREADS":os.environ.get("ODD_THREADS","4")})
    print(run.stdout,flush=True)
    match=re.search(
        r"SUMMARY i=\[(\d+),(\d+)\) maxabs=(\d+) tested=(\d+) u_odd=(\d+) normbound=(\d+) short=(\d+) min_numerator=(\d+)",
        run.stdout)
    assert match,(run.returncode,run.stderr,run.stdout)
    start,end,used,tested,odd,admitted,short,minimum=map(int,match.groups())
    assert (start,end,used)==(0,45,limit)
    assert run.returncode==(3 if short else 0)
    if limit==8:
        assert (tested,odd,admitted)==(5753341440,5435786880,3226914720)
    report={
        "candidate_matrix_sha256":digest,
        "construction":"L={(u,v)/4 : v == u*A (mod 16)}",
        "complete_scope":"exactly 4 nonzero first-half coordinates, at least one odd",
        "coefficient_bound":limit,
        "exhaustive_for_scope":limit==8,
        "pruned_before_full_candidate_count":"nonzero coefficients with prefix squares >=128",
        "candidate_count_after_safe_prefix_pruning":tested,
        "candidates_with_odd_first_half":odd,
        "passed_parity_cost_lower_bound":admitted,
        "short_vectors_found":short,
        "smallest_reported_numerator":minimum if short else None,
        "full_lattice_minimum_certified":False,
        "next_open_case":"u has >=5 nonzero coordinates and at least one odd"}
    Path(sys.argv[2]).write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
    if short:sys.exit(3)

if __name__=="__main__":
    main()
