#!/usr/bin/env python3
"""Exact proof-method regression tests, including 96D checkpoint safety."""
import json,sys
from pathlib import Path
from tempfile import TemporaryDirectory
from exact_svp import (Enumerator,exact_ldlt,verify_factor,read_gram,
                       prefixes,digest,selftest)
from check_basis_96 import verify
from global_aggregate import main as aggregate_main

HERE=Path(__file__).resolve().parent
DATA=HERE/'z16_route'
CAND=DATA/'candidate_3cycle.json'
U=DATA/'candidate_3cycle_LLL_transform.json'
GRAM=DATA/'candidate_3cycle_LLL_gram.json'


def tests():
    selftest()
    proof=verify(CAND,U,GRAM)
    assert proof['basis_transform_det']==1
    assert proof['reduced_gram_determinant']==1
    assert proof['reduced_GS_last_norm_squared']=='1/50'
    G=read_gram(GRAM)
    suff=prefixes(G,depth=1,bound=8)
    assert suff==[[i] for i in range(21)],suff
    killed=Enumerator(G,bound=8,suffix=[20])
    result=killed.run(1000)
    assert result['complete'] and killed.nodes_total==0 and killed.witness is None
    with TemporaryDirectory(prefix='exact_global_test_') as temp:
        td=Path(temp)
        checkpoint=td/'short-run.json'
        e=Enumerator(G,bound=8)
        e.run(1100)
        assert not e.finished and e.witness is None
        a=e.save(checkpoint)
        assert a['certificate_status']=='INCOMPLETE'
        resumed=Enumerator.resume(G,json.loads(checkpoint.read_text()))
        resumed.run(900)
        baseline=Enumerator(G,bound=8)
        baseline.run(2000)
        assert resumed.nodes_total==baseline.nodes_total
        assert resumed.leaves_total==baseline.leaves_total
        def sig(obj):
            return [(f['i'],f['next'],f['active'],f['rn'],f['rd'])
                    for f in obj.frames]
        assert sig(resumed)==sig(baseline)
        assert resumed.x==baseline.x and resumed.shift==baseline.shift
        tampered=json.loads(checkpoint.read_text());tampered['gram_sha256']='0'*64
        try:
            Enumerator.resume(G,tampered)
        except ValueError:pass
        else:raise AssertionError('wrong SHA should fail resume')
        # Exhaust a separate small E8 example and audit 2 suffix reports.
        E=[[2 if i==j else 0 for j in range(8)] for i in range(8)]
        for i,j in [(k,k+1) for k in range(6)]+[(2,7)]:E[i][j]=E[j][i]=-1
        parts=prefixes(E,1,bound=1)
        plist={'status':'COMPLETE_PREFIX_LIST','gram_sha256':digest(E),
               'rank':8,'bound':1,'depth':1,'suffixes':parts}
        (td/'prefixes.json').write_text(json.dumps(plist))
        states=td/'states';states.mkdir()
        for i,suffix in enumerate(parts):
            check=Enumerator(E,bound=1,suffix=suffix)
            assert check.run(100000)['complete']
            assert check.witness is None
            check.save(states/f'{i}.json')
        args=sys.argv;sys.argv=['global_aggregate.py',str(td/'prefixes.json'),
                                str(states),str(td/'coverage.json')]
        try:aggregate_main()
        finally:sys.argv=args
        coverage=json.loads((td/'coverage.json').read_text())
        assert coverage['smallest_nonzero_norm_gt_bound_certified']
        # Restore negative coverage to the intentionally incomplete case.
        not_complete=td/'states'/'0.json'
        d=json.loads(not_complete.read_text());d['finished']=False
        d['certificate_status']='INCOMPLETE';not_complete.write_text(json.dumps(d))
        args=sys.argv;sys.argv=['global_aggregate.py',str(td/'prefixes.json'),
                                str(states),str(td/'coverage2.json')]
        try:
            try:aggregate_main()
            except ValueError:pass
            else:raise AssertionError('bad partition must not validate')
        finally:sys.argv=args
    print('PASS all global 96D proof-method and resumption tests')

if __name__=='__main__':tests()
