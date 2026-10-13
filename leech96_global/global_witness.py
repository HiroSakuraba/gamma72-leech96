#!/usr/bin/env python3
"""Independently verify a short vector in the exact 96D Z/16 candidate.

Usage: python global_witness.py candidate.json U.json reduced_gram.json witness.json
Witness file may be an exact_svp checkpoint containing `witness`, or a
standalone JSON dict containing `coefficients` and `norm`.
"""
import json, sys
from pathlib import Path
from check_basis_96 import gram
from exact_svp import digest

def check(candidate_path, transform_path, gram_path, witness_path):
    candidate=json.loads(Path(candidate_path).read_text())
    A=candidate['A']; U=json.loads(Path(transform_path).read_text())
    G=json.loads(Path(gram_path).read_text())
    data=json.loads(Path(witness_path).read_text())
    witness=data.get('witness',data)
    if not isinstance(witness,dict):raise ValueError('no witness in result')
    coeff=witness['coefficients']
    if not len(coeff)==96 or not all(type(x) is int for x in coeff):
        raise ValueError('invalid integer coefficients')
    if data.get('gram_sha256',digest(G))!=digest(G):
        raise ValueError('wrong Gram SHA-256')
    y=[sum(coeff[i]*U[i][j] for i in range(96)) for j in range(96)]
    u=y[:48]; v=[sum(u[i]*A[i][j] for i in range(48))+16*y[48+j] for j in range(48)]
    numerator=sum(x*x for x in u)+sum(x*x for x in v)
    reduced_norm=sum(coeff[i]*sum(G[i][j]*coeff[j] for j in range(96))
                     for i in range(96))
    if numerator!=16*reduced_norm or numerator%16:
        raise AssertionError('basis norm disagreement')
    if not all((v[j]-sum(u[i]*A[i][j] for i in range(48)))%16==0 for j in range(48)):
        raise AssertionError('not a code lattice vector')
    if witness['norm']!=reduced_norm:raise AssertionError('false norm')
    if reduced_norm<=0:raise AssertionError('zero witness')
    bound=data.get('bound',8)
    if reduced_norm>bound:raise AssertionError('witness outside requested radius')
    result={'verified_membership':True,'norm':reduced_norm,
            'numerator':numerator,'reduced_gram_sha256':digest(G),
            'u':u,'v':v,'candidate_killed':reduced_norm<=8}
    print(json.dumps(result,indent=2))
    return result

if __name__=='__main__':
    if len(sys.argv)!=5:raise SystemExit(__doc__)
    check(*sys.argv[1:])
