#!/usr/bin/env python3
"""Audit the complete collection of 42 exact all-odd support-seven certificates.
Usage: python3 odd7_aggregate.py DIRECTORY combined.json
DIRECTORY may contain nested GitHub Actions artifact folders.
Refuses missing/duplicate partitions, incorrect matrix hash or counts,
or any explicit <=8 short vector. Does not repeat the C++ enumerations.
"""
import json,sys
from pathlib import Path
from math import comb
SHA="929e71c3a3f54f3361d9a8c318786b1da726e141250a93969e9b2cec8cd3c656"
def main():
    if len(sys.argv)!=3:sys.exit(__doc__)
    certs={}
    for path in sorted(Path(sys.argv[1]).rglob("*.json")):
        r=json.loads(path.read_text())
        if "first_index" not in r:continue
        i=r["first_index"]
        assert type(i) is int and 0<=i<=41
        assert i not in certs,("duplicate index",i,path)
        assert r["candidate_A_sha256"]==SHA
        assert r["complete_for_first_index"] is True
        assert r["full_minimum_norm_proved"] is False
        assert r["short_vectors_norm_sq_le_8"]==0
        assert r["supports"]==comb(47-i,6)
        assert r["mod8_residue_patterns"]==r["supports"]*4**7
        assert 0<=r["mod8_survivors"]<=r["mod8_residue_patterns"]
        assert 0<=r["centered_mod16_lifts_checked"]<=r["mod8_survivors"]*64
        certs[i]=r
    assert sorted(certs)==list(range(42)),("missing",sorted(set(range(42))-set(certs)))
    supports=sum(r["supports"] for r in certs.values())
    patterns=sum(r["mod8_residue_patterns"] for r in certs.values())
    assert supports==comb(48,7)==73629072
    assert patterns==supports*4**7==1206338715648
    report={
      "candidate_A_sha256":SHA,
      "scope":"first-half vector has exactly seven nonzero entries, all odd",
      "independent_first_index_partitions":42,
      "supports":supports,
      "mod8_residue_patterns":patterns,
      "mod8_survivors":sum(r["mod8_survivors"] for r in certs.values()),
      "centered_mod16_lifts":sum(r["centered_mod16_lifts_checked"] for r in certs.values()),
      "short_vectors_norm_sq_le_8":0,
      "complete_for_stated_scope":True,
      "global_minimum_certified":False
    }
    Path(sys.argv[2]).write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
if __name__=="__main__":main()
