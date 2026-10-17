"""Search is heuristic; panel inequalities and integer isolation use Arb.

Inputs are exact rational strings. Polynomial coefficients are in ascending
powers, complex rationals represented by [real, imag]. No floats enter proofs.
"""
from fractions import Fraction as Q
from math import comb, factorial
from flint import arb, acb, ctx


class Unresolved(RuntimeError):
    pass


def real(x):
    q = Q(str(x))
    return arb(q.numerator) / q.denominator


def complex_ball(z):
    return acb(real(z[0]), real(z[1]))


def point(z):
    return tuple(Q(str(x)) for x in z)


def midpoint(a, b):
    return tuple((x+y)/2 for x, y in zip(a, b))


def data(model, a, b):
    a, b = point(a), point(b)
    m = complex_ball(midpoint(a, b))
    d = complex_ball(tuple((y-x)/2 for x,y in zip(a,b)))
    lam = [Q(t['lambda']) for t in model]
    # Midpoint scaling by a nonzero exponential; choose an exponent with largest
    # real growth. This is range reduction, independent of the gauge policy.
    anchor = max(lam, key=lambda l: l*midpoint(a,b)[0])
    coeff = []
    for term, l in zip(model, lam):
        p = [complex_ball(z) for z in term['poly']]
        scale = (real(l-anchor)*m).exp()
        coeff.append([scale * d**k * sum(
            (p[n]*comb(n,k)*m**(n-k) for n in range(k,len(p))), acb(0))
            for k in range(len(p))])
    return lam, coeff, d


def choose_center(model, a, b, policy):
    lam, coeff, d = data(model, a, b)
    if isinstance(policy, (Q, int)):
        return Q(policy)
    if policy == 'raw':
        return Q(0)
    lo, hi = min(lam), max(lam)
    if policy == 'global':
        return (lo+hi)/2
    if policy == 'anchor':
        return lam[0]
    if policy != 'local':
        raise ValueError('unknown policy')
    if lo == hi:
        return lo
    # Optimize log(sum W_j exp(r |lambda_j-s|)). Any candidate is sound:
    # all final comparisons are repeated in ball arithmetic in panel().
    import numpy as np
    from scipy.special import logsumexp
    from scipy.optimize import minimize_scalar
    x = np.array([float(l-lo) for l in lam])
    weights = np.array([float(sum((abs(c) for c in cs), arb(0))) for cs in coeff])
    mask = weights > 0
    if not mask.any():
        return (lo+hi)/2
    logs = np.log(weights[mask]); x=x[mask]
    r = float(abs(d))
    fun = lambda s: float(logsumexp(logs+r*np.abs(x-s)))
    opt = minimize_scalar(fun, bounds=(0,float(hi-lo)), method='bounded',
                          options={'xatol':1e-10, 'maxiter':80})
    candidates = [0., float(hi-lo)/2, float(hi-lo), float(opt.x)] + list(x)
    s = min(candidates,key=fun)
    # Round relative to lo, keeping the common-exponent-shift invariance.
    return lo + Q(round(s*2**30),2**30)


def center_ball(center):
    return complex_ball(center) if isinstance(center,(tuple,list)) else acb(real(center))


def difference_ball(l,center):
    if isinstance(center,(tuple,list)):
        return acb(real(l-Q(center[0])),real(-Q(center[1])))
    return acb(real(l-Q(center)))


def enclosure(model,a,b,center,order=1):
    lam, coeff, d = data(model,a,b)
    F = sum((cs[0] for cs in coeff),acb(0))
    B = arb(0)
    mu = []
    for l, cs in zip(lam,coeff):
        u = difference_ball(l,center)*d
        mu.append(u)
        rad = abs(u)
        B += abs(cs[0])*rad.expm1()
        B += sum((abs(c) for c in cs[1:]),arb(0))*rad.exp()
    if order==2:
        # Preserve cancellation in the first two Taylor coefficients before
        # taking absolute values. The exponential tail is bounded outwards.
        T=[acb(0),acb(0),acb(0)]
        tail=arb(0)
        for cs,u in zip(coeff,mu):
            rad=abs(u)
            for k,c in enumerate(cs):
                if k<=2:
                    for n in range(k,3):
                        T[n]+=c*u**(n-k)/factorial(n-k)
                    tail+=abs(c)*rad.exp()*rad**(3-k)/factorial(3-k)
                else:
                    tail+=abs(c)*rad.exp()
        B2=abs(T[1])+abs(T[2])+tail
        # Selecting either valid enclosure by approximate size is sound.
        if float(B2.upper())<float(B.upper()):B=B2
    elif order!=1:
        raise ValueError('enclosure order must be 1 or 2')
    return F,mu,B,coeff,d


def choose_jet(model,a,b):
    lam,cs,d=data(model,a,b)
    F=sum((c[0] for c in cs),acb(0))
    candidates=[choose_center(model,a,b,'local'),Q(0),lam[0]]
    if not F.contains(0):
        # Subtract the first exponent before division so common shifts do not
        # contaminate the proposal with avoidable cancellation.
        D=sum(((c[1]/d if len(c)>1 else acb(0))+real(l-lam[0])*c[0]
               for l,c in zip(lam,cs)),acb(0))/F
        if D.is_finite():
            candidates.append((lam[0]+Q(round(float(D.real)*2**30),2**30),
                               Q(round(float(D.imag)*2**30),2**30)))
    return min(candidates,key=lambda s:float(enclosure(model,a,b,s,2)[2].upper()))


def panel(model, a, b, center, safety=Q(4,5), order=1):
    """Return certified phase increment and margin, or None if unresolved."""
    F,mu,B,coeff,d=enclosure(model,a,b,center,order)
    if not (B < real(safety)*abs(F)):
        return None
    # G(t)/F lies in |w-1|<1. Principal logs are therefore the same
    # continuous branch, even when the original f turns rapidly.
    ends = []
    for t in [-1,1]:
        G = sum((sum((c*t**k for k,c in enumerate(cs)),acb(0))*(u*t).exp()
                 for cs,u in zip(coeff,mu)),acb(0))
        ratio = G/F
        if not (ratio.real > 0):
            return None
        ends.append(ratio.log().imag)
    delta = 2*(center_ball(center)*d).imag + ends[1]-ends[0]
    return delta, B/abs(F)


def certificate(model, contour, panels, precision=128):
    return {'format':'contourgauge-1','precision':precision,'model':model,
            'contour':[[str(Q(str(x))) for x in z] for z in contour],
            'panels':panels}


def verify(cert, precision=None):
    """Replay a complete, ordered contour partition; trust no saved count."""
    ctx.prec = precision or cert['precision']
    if cert.get('format') != 'contourgauge-1':
        raise ValueError('unknown certificate format')
    order=cert.get('enclosure_order',1)
    if order not in [1,2]:raise ValueError('unknown enclosure order')
    contour = [point(z) for z in cert['contour']]
    if len(contour)<4 or contour[0]!=contour[-1]:
        raise ValueError('contour must close')
    # Experiments use simple CCW rectangles. Restrict the public verifier to
    # exactly that class instead of trusting an unverified geometry claim.
    if len(contour)!=5:
        raise ValueError('only rectangles are supported')
    a,b,c,d,_=contour
    if not (a[1]==b[1] and b[0]==c[0] and c[1]==d[1] and d[0]==a[0]
            and a[0]<b[0] and a[1]<d[1]):
        raise ValueError('expected a simple counterclockwise rectangle')
    groups=[[] for _ in range(4)]
    for p in cert['panels']:
        if not isinstance(p['edge'],int) or not 0<=p['edge']<4:
            raise ValueError('invalid edge')
        groups[p['edge']].append(p)
    total=arb(0); margins=[]
    for e,group in enumerate(groups):
        cursor=Q(0)
        for p in group:
            l,r=Q(p['left']),Q(p['right'])
            if l!=cursor or not l<r<=1:
                raise ValueError('partition gap, overlap, or reversal')
            x=tuple(u+(v-u)*l for u,v in zip(contour[e],contour[e+1]))
            y=tuple(u+(v-u)*r for u,v in zip(contour[e],contour[e+1]))
            s=p['center'] if isinstance(p['center'],list) else Q(p['center'])
            out=panel(cert['model'],x,y,s,order=order)
            if out is None:
                raise Unresolved('panel inequality not established')
            total+=out[0]; margins.append(float(out[1].upper()))
            cursor=r
        if cursor!=1:
            raise ValueError('incomplete edge')
    turns=total/(2*arb.pi())
    n=turns.unique_fmpz()
    if n is None:
        raise Unresolved('integer not isolated')
    return {'count':int(n),'turns':str(turns),'panels':len(margins),
            'max_ratio':max(margins)}


def count(model, contour, policy='local', precision=128, max_depth=28,
          max_attempts=200000, enclosure_order=1):
    """Bisect until every panel satisfies the disk inequality."""
    ctx.prec=precision
    contour=[point(z) for z in contour]
    saved=[]; attempts=0
    if isinstance(policy,str) and policy.endswith('2'):
        policy=policy[:-1];enclosure_order=2
    if policy=='jet':enclosure_order=2
    for e,(start,end) in enumerate(zip(contour,contour[1:])):
        stack=[(Q(0),Q(1),0)]
        while stack:
            l,r,depth=stack.pop(); attempts+=1
            if attempts>max_attempts:
                raise Unresolved('attempt budget exhausted')
            x=tuple(u+(v-u)*l for u,v in zip(start,end))
            y=tuple(u+(v-u)*r for u,v in zip(start,end))
            s=choose_jet(model,x,y) if policy=='jet' else choose_center(model,x,y,policy)
            out=panel(model,x,y,s,order=enclosure_order)
            if out is not None:
                ss=[str(v) for v in s] if isinstance(s,tuple) else str(s)
                saved.append({'edge':e,'left':str(l),'right':str(r),'center':ss})
            else:
                if depth>=max_depth:
                    raise Unresolved('depth exhausted; boundary zero or precision limitation')
                m=(l+r)/2
                stack.extend([(m,r,depth+1),(l,m,depth+1)])
    cert=certificate(model,contour,saved,precision)
    if enclosure_order!=1:cert['enclosure_order']=enclosure_order
    result=verify(cert)
    result['attempts']=attempts
    return result,cert
