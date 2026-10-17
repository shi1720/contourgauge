"""Follow-up ablation: cancellation-preserving bound on all delay instances."""
import json,time,statistics
from pathlib import Path
from fractions import Fraction as Q
from contourgauge import count,verify
from study import ROOT,cases

def main():
    baseline=json.loads((ROOT/'results/study.json').read_text())
    expected={r['name']:r['count'] for r in baseline if r['policy']=='local'}
    rows=[]
    for c in cases():
        if c['family'] not in ['delay','multidelay']:continue
        for policy in ['raw2','global2','local2','jet','oracle9_2']:
            t0=time.perf_counter()
            if policy=='oracle9_2':
                lam=[Q(t['lambda']) for t in c['model']];lo,hi=min(lam),max(lam)
                trials=[];best=None;cert=None
                for k in range(9):
                    s=lo+(hi-lo)*Q(k,8)
                    rr,cc=count(c['model'],c['contour'],s,enclosure_order=2)
                    trials.append(dict(center=str(s),panels=rr['panels']))
                    if best is None or rr['panels']<best['panels']:
                        best,cert=rr,cc
                r=best
            else:
                r,cert=count(c['model'],c['contour'],policy)
                trials=None
            seconds=time.perf_counter()-t0
            assert r['count']==expected[c['name']]
            assert verify(cert,192)['count']==expected[c['name']]
            ts=[]
            for _ in range(5):
                t=time.perf_counter();verify(cert,192);ts.append(time.perf_counter()-t)
            row=dict(name=c['name'],family=c['family'],policy=policy,**r,
                     seconds=seconds,replay_seconds=statistics.median(ts),trials=trials)
            (ROOT/'results/certificates'/f"{c['name']}_{policy}.json").write_text(json.dumps(cert,separators=(',',':')))
            rows.append(row)
            (ROOT/'results/cancellation.json').write_text(json.dumps(rows,indent=2))
            print(c['name'],policy,r['panels'],round(seconds,3),flush=True)

if __name__=='__main__':main()
