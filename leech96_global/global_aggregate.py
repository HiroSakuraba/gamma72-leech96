#!/usr/bin/env python3
"""Check coverage and completion of exact-svp suffix-partition reports.

This verifies that *every* mathematically possible sign-normalized suffix
appears exactly once with matching matrix, bound, rank and complete negative
status. It does not replace re-execution of each branch, nor does it prove
that an externally supplied JSON log is truthful.

Usage:
  python global_aggregate.py Gram.json prefix_list.json results_directory out.json
"""
import json,sys
from pathlib import Path
from exact_svp import read_gram, prefixes, digest

def main():
    if len(sys.argv)!=5:raise SystemExit(__doc__)
    G=read_gram(sys.argv[1])
    src=json.loads(Path(sys.argv[2]).read_text())
    if src['status']!='COMPLETE_PREFIX_LIST':
        raise ValueError('no proof that prefix listing completed')
    if src['gram_sha256']!=digest(G) or src['rank']!=len(G):
        raise ValueError('prefix list has wrong Gram input')
    true_suffixes=prefixes(G,depth=src['depth'],bound=src['bound'])
    if src['suffixes']!=true_suffixes:
        raise ValueError('prefix list omits/adds/misorders mathematically admissible suffixes')
    suffixes={tuple(k) for k in true_suffixes}
    assert len(suffixes)==len(true_suffixes)
    seen={};ignored=0
    for f in Path(sys.argv[3]).rglob('*.json'):
        d=json.loads(f.read_text())
        if d.get('format')!='exact-svp-state':ignored+=1;continue
        if d['gram_sha256']!=src['gram_sha256']:
            raise ValueError(f'wrong matrix in {f}')
        if (d['bound'],d['rank'])!=(src['bound'],src['rank']):
            raise ValueError(f'wrong radius/rank in {f}')
        s=tuple(d['suffix'])
        if len(s)!=src['depth'] or s not in suffixes or s in seen:
            raise ValueError(f'unexpected/duplicate suffix {s} in {f}')
        if not (d['finished'] and not d['stack'] and not d['witness'] and
                d['certificate_status']=='COMPLETE_NO_SHORT_VECTOR'):
            raise ValueError(f'incomplete or positive partition at {f}')
        seen[s]=str(f)
    missing=sorted(suffixes-set(seen))
    result={'status':'COMPLETE_COVERAGE_NO_SHORT_VECTOR' if not missing else 'INCOMPLETE',
            'gram_sha256':src['gram_sha256'],'bound':src['bound'],
            'rank':src['rank'],'suffix_depth':src['depth'],
            'expected_partitions':len(suffixes),'complete_negative_partitions':len(seen),
            'missing_suffixes':[list(x) for x in missing],
            'smallest_nonzero_norm_gt_bound_certified':not missing,
            'warning':'Coverage audit only; individual partition computations must be independently rerunnable and trustworthy.'}
    Path(sys.argv[4]).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if missing:sys.exit(2)

if __name__=='__main__':main()
