#!/usr/bin/env python3
"""Run exact all-support-seven +/-1 sign scan; independent Gray-code regression.
Usage: python3 support7_pm1_verify.py candidate_3cycle.json result.json
"""
import json,os,re,subprocess,sys,tempfile
from math import comb
from pathlib import Path
from support7_allodd_verify import candidate_matrix

def main():
    if len(sys.argv)!=3:sys.exit(__doc__)
    A,sha=candidate_matrix(json.loads(Path(sys.argv[1]).read_text()))
    with tempfile.TemporaryDirectory(prefix="leech96_pm1_") as t:
        tmp=Path(t)
        header="#include <cstdint>\nstatic const int8_t A_pad[48][64]={\n"
        header+=",\n".join("{"+",".join(map(str,row+[0]*16))+"}" for row in A)
        header+="\n};\n"
        (tmp/"matrix48.h").write_text(header)
        exe=tmp/"scan"
        source=Path(__file__).with_name("support7_pm1.cpp")
        subprocess.run([os.environ.get("CXX","g++"),"-O3","-std=c++17","-mavx2",
                        "-fopenmp","-I",str(tmp),str(source),"-o",str(exe)],check=True)
        proc=subprocess.run([str(exe),"0","42"],capture_output=True,text=True,
                            env={**os.environ,"OMP_NUM_THREADS":
                                 os.environ.get("OMP_NUM_THREADS","4")})
    print(proc.stdout,flush=True)
    m=re.search(r"SUMMARY first=\[(\d+),(\d+)\) supports=(\d+) "
                r"signed_patterns=(\d+) passed_mod8=(\d+) short=(\d+) "
                r"best_numerator=(\d+)",proc.stdout)
    assert m,(proc.returncode,proc.stdout,proc.stderr)
    first,end,supports,signs,passed,short,best=map(int,m.groups())
    assert (first,end)==(0,42)
    assert supports==comb(48,7)==73629072
    assert signs==supports*64==4712260608
    assert passed==95377
    assert proc.returncode==(3 if short else 0)
    report={
        "candidate_A_sha256":sha,
        "scope":"exactly seven nonzero centered u coordinates, all +/-1",
        "supports":supports,"sign_patterns":signs,"passed_mod8":passed,
        "norm_sq_le_8_vectors":short,"scope_exhaustive":True,
        "global_minimum_certified":False,
    }
    Path(sys.argv[2]).write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
    if short:raise SystemExit("Explicit <=8 short vector found")

if __name__=="__main__":main()
