"""Offline nine-center fixed-gauge oracle, and certificate replay timing."""
import json,time,statistics
from fractions import Fraction as Q
from pathlib import Path
from contourgauge import count,verify
from study import cases,ROOT

def main():
    rows=[]
    for c in cases():
        if c['family']=='shift':continue
        lam=[Q(t['lambda']) for t in c['model']]
        lo,hi=min(lam),max(lam)
        trials=[];best=None;bestcert=None
        for k in range(9):
            s=lo+(hi-lo)*Q(k,8)
            t=time.perf_counter();r,cert=count(c['model'],c['contour'],s)
            trials.append(dict(center=str(s),panels=r['panels'],seconds=time.perf_counter()-t))
            if best is None or r['panels']<best['panels']:
                best=dict(r,center=str(s));bestcert=cert
        (ROOT/'results/certificates'/f"{c['name']}_oracle9.json").write_text(json.dumps(bestcert,separators=(',',':')))
        rows.append(dict(name=c['name'],family=c['family'],trials=trials,best=best))
        print(c['name'],best['panels'],flush=True)
        (ROOT/'results/global_sweep.json').write_text(json.dumps(rows,indent=2))
    # Exclude common-shift stress for replay comparisons: it is not a fair
    # practical baseline. Median of five independent verifier invocations.
    timing=[]
    for row in json.loads((ROOT/'results/study.json').read_text()):
        if row['family']=='shift':continue
        name,policy=row['name'],row['policy']
        cert=json.loads((ROOT/'results/certificates'/f'{name}_{policy}.json').read_text())
        ts=[]
        for _ in range(5):
            start=time.perf_counter();r=verify(cert,192);ts.append(time.perf_counter()-start)
        timing.append(dict(name=name,policy=policy,panels=r['panels'],seconds=statistics.median(ts),repetitions=5))
    (ROOT/'results/replay_timing.json').write_text(json.dumps(timing,indent=2))

if __name__=='__main__':main()
