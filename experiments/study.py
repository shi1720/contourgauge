"""Deterministic study; writes all inputs, certificates and failures."""
import sys, json, time, platform, random, cmath, math
from pathlib import Path
from fractions import Fraction as Q
from flint import acb, arb, ctx
from contourgauge import count, verify, Unresolved

ROOT=Path(__file__).resolve().parents[1]

def rectangle(x0,x1,y0,y1):
    return [[str(x0),str(y0)],[str(x1),str(y0)],[str(x1),str(y1)],
            [str(x0),str(y1)],[str(x0),str(y0)]]

def term(lam, coeff):
    return {'lambda':str(lam),'poly':[[str(x), '0'] for x in coeff]}

def product_roots(roots, shift=0, rate=1):
    # Exact rational real roots in the w=exp(rate*z) plane.
    c=[Q(1)]
    for w in roots:
        out=[Q(0)]*(len(c)+1)
        for k,v in enumerate(c):
            out[k]-=w*v;out[k+1]+=v
        c=out
    return [term(Q(shift)+rate*k,[v]) for k,v in enumerate(c)]

def exp_expected(roots, rect, rate=Q(1)):
    # Independent ball evaluation of the explicitly known logarithmic lattice.
    ctx.prec=192
    xmin,xmax=map(Q,[rect[0][0],rect[1][0]])
    ymin,ymax=map(Q,[rect[0][1],rect[2][1]])
    n=0
    for w in roots:
        ww=acb(arb(w.numerator)/w.denominator)
        for k in range(-100,101):
            z=(ww.log()+acb(0,2*k*arb.pi()))/(arb(rate.numerator)/rate.denominator)
            inside=(z.real>arb(str(xmin)) and z.real<arb(str(xmax)) and
                    z.imag>arb(str(ymin)) and z.imag<arb(str(ymax)))
            outside=(z.real<arb(str(xmin)) or z.real>arb(str(xmax)) or
                     z.imag<arb(str(ymin)) or z.imag>arb(str(ymax)))
            if not (inside or outside):
                raise ValueError('reference lies on boundary')
            n+=bool(inside)
    # All studied |Im z|<20, rate<=2. Outside |k|<=100 impossible here.
    assert max(abs(float(ymin)),abs(float(ymax)))*float(rate)<100
    return n

def cases():
    out=[]
    # Common-shift stress test; mathematical zeros are unchanged.
    roots=[Q(1,2),Q(2),Q(-1)]
    rect=rectangle(-2,2,-7,7)
    for shift in [0,8,32,128,512]:
        out.append(dict(name=f'shift_{shift}',family='shift',
            model=product_roots(roots,shift),contour=rect,
            expected=exp_expected(roots,rect),shift=shift))
    # Controlled dominance transitions and exponent spans.
    for rate in [Q(1,2),Q(1),Q(2)]:
        for width in [Q(1),Q(3),Q(6)]:
            roots=[Q(1,4),Q(2),Q(-3)]
            rect=rectangle(-width,width,Q(-13,2),Q(13,2))
            out.append(dict(name=f'dominance_{rate}_{width}'.replace('/','-'),
                family='dominance',model=product_roots(roots,0,rate),contour=rect,
                expected=exp_expected(roots,rect,rate)))
    # Close-to-boundary lattice roots at z=0, +-2pi*i,...
    for eps in [Q(1,10),Q(1,100),Q(1,1000),Q(1,10000)]:
        rect=rectangle(-eps,2,-5,5)
        roots=[Q(1),Q(2)]
        out.append(dict(name=f'clearance_{eps.denominator}',family='clearance',
            model=product_roots(roots),contour=rect,expected=exp_expected(roots,rect)))
    rng=random.Random(6402027)
    for i in range(36):
        n=rng.randint(2,5)
        roots=[Q(rng.choice([-1,1])*rng.randint(1,30),rng.randint(2,12)) for _ in range(n)]
        rate=rng.choice([Q(1,2),Q(1),Q(2)])
        rect=rectangle(Q(-7,3),Q(8,3),Q(-17,3),Q(19,3))
        out.append(dict(name=f'lattice_{i:02d}',family='lattice',
            model=product_roots(roots,rng.choice([0,4,16]),rate),contour=rect,
            expected=exp_expected(roots,rect,rate)))
    # Polynomial coefficients: true scalar retarded DDE characteristic.
    for tau in [Q(1,2),Q(1),Q(2),Q(4)]:
        for b in [Q(1,2),Q(2),Q(4)]:
            rect=rectangle(0,3,-20,20)
            out.append(dict(name=f'delay_{tau}_{b}'.replace('/','-'),family='delay',
                model=[term(0,[1,1]),term(-tau,[b])],contour=rect,
                expected=None,tau=str(tau),b=str(b)))
    # Two distinct delays: not reducible to the single-delay crossing formula.
    for i in range(24):
        tau1=Q(rng.randint(2,15),10)
        tau2=Q(rng.randint(16,40),10)
        b1=Q(rng.randint(2,25),10)
        b2=Q(rng.randint(2,25),10)
        out.append(dict(name=f'multidelay_{i:02d}',family='multidelay',
            model=[term(0,[1,1]),term(-tau1,[b1]),term(-tau2,[b2])],
            contour=rectangle(0,5,-12,12),expected=None,
            tau1=str(tau1),tau2=str(tau2),b1=str(b1),b2=str(b2)))
    return out

def quadrature_reference(model,contour):
    # Independent algorithm: integrate f'/f with Arb's adaptive analytic
    # quadrature. No panel disk or exponential gauge is used.
    ctx.prec=192
    def f(z, analytic):
        val=acb(0);der=acb(0)
        for t in model:
            lam=arb(t['lambda'])
            coeff=[acb(arb(x),arb(y)) for x,y in t['poly']]
            p=sum((c*z**k for k,c in enumerate(coeff)),acb(0))
            dp=sum((k*c*z**(k-1) for k,c in enumerate(coeff) if k),acb(0))
            E=(lam*z).exp()
            val+=p*E;der+=(dp+lam*p)*E
        return der/val
    total=acb(0)
    for a,b in zip(contour,contour[1:]):
        aa=acb(arb(a[0]),arb(a[1]));bb=acb(arb(b[0]),arb(b[1]))
        total+=acb.integral(f,aa,bb,rel_tol=arb('1e-35'),abs_tol=arb('1e-35'),
                            eval_limit=200000,depth_limit=40)
    turns=total/acb(0,2*arb.pi())
    n=turns.real.unique_fmpz()
    if n is None or not turns.imag.contains(0):
        raise Unresolved('quadrature reference unresolved: '+str(turns))
    return int(n),str(turns)

def main():
    (ROOT/'results/certificates').mkdir(exist_ok=True)
    allcases=cases(); results=[]
    (ROOT/'results/inputs.json').write_text(json.dumps(allcases,indent=2))
    # Warm imported optimizer before timing.
    count([term(0,[1]),term(1,[1])],rectangle(-1,1,-1,1),'local')
    for case in allcases:
        if case['family'] in ['delay','multidelay']:
            case['expected'],case['reference_ball']=quadrature_reference(case['model'],case['contour'])
        if case['family']=='delay':
            tau,b=float(Q(case['tau'])),float(Q(case['b']))
            exact=0 if b<=1 else 2*sum(tau>(math.acos(-1/b)+2*k*math.pi)/math.sqrt(b*b-1) for k in range(100))
            assert exact==case['expected']
            case['analytic_expected']=exact
        for policy in ['raw','anchor','global','local']:
            row={k:v for k,v in case.items() if k not in ['model','contour']}
            row['policy']=policy
            t0=time.perf_counter()
            try:
                result,cert=count(case['model'],case['contour'],policy,max_attempts=60000)
                row['seconds']=time.perf_counter()-t0
                row.update(result); row['status']='certified'
                assert row['count']==case['expected'],(case['name'],row)
                replay=verify(cert,192)
                assert replay['count']==case['expected']
                row['replay_turns']=replay['turns']
                path=ROOT/'results/certificates'/f"{case['name']}_{policy}.json"
                path.write_text(json.dumps(cert,separators=(',',':')))
            except Unresolved as e:
                row['seconds']=time.perf_counter()-t0
                row['status']='unresolved';row['reason']=str(e)
            results.append(row)
            print(case['name'],policy,row['status'],row.get('count'),row.get('panels'),round(row['seconds'],3),flush=True)
            (ROOT/'results/study.json').write_text(json.dumps(results,indent=2))
    env={'python':sys.version,'platform':platform.platform(),'seed':6402027,
         'precision':128,'replay_precision':192,'safety':'4/5'}
    (ROOT/'results/environment.json').write_text(json.dumps(env,indent=2))

if __name__=='__main__':main()
