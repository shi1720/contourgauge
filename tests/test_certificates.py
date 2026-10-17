import unittest, copy, random
from fractions import Fraction as Q
from contourgauge import count, verify, Unresolved
from contourgauge.core import data, panel, choose_center


def box(x0=-2,x1=2,y0=-2,y1=2):
    return [[str(x0),str(y0)],[str(x1),str(y0)],[str(x1),str(y1)],
            [str(x0),str(y1)],[str(x0),str(y0)]]


class Certificates(unittest.TestCase):
    def test_polynomial_known_roots_and_multiplicity(self):
        # (z-1)^3(z+3), including a multiple interior root.
        m=[{'lambda':'17','poly':[['-3','0'],['8','0'],['-6','0'],['0','0'],['1','0']]}]
        for policy in ['raw','global','local']:
            r,c=count(m,box(),policy)
            self.assertEqual(r['count'],3)
            self.assertEqual(verify(c,192)['count'],3)

    def test_nonzero_exponential_four_panels(self):
        for k in [-10000,0,10000]:
            m=[{'lambda':str(k),'poly':[['1','0']]}]
            r,_=count(m,box(),'local')
            self.assertEqual((r['count'],r['panels']),(0,4))

    def test_common_shift_partition_invariance(self):
        m=[{'lambda':'0','poly':[['-1','0']]},{'lambda':'1','poly':[['1','0']]}]
        C=box(-1,2,-7,7)
        r,c=count(m,C)
        for shift in [-128, Q(10001,7)]:
            mm=copy.deepcopy(m)
            for t in mm:t['lambda']=str(Q(t['lambda'])+shift)
            rr,cc=count(mm,C)
            self.assertEqual(rr['count'],r['count'])
            self.assertEqual([(p['edge'],p['left'],p['right']) for p in c['panels']],
                             [(p['edge'],p['left'],p['right']) for p in cc['panels']])

    def test_mutated_partition_rejected(self):
        m=[{'lambda':'0','poly':[['1','0'],['1','0']]}]
        _,c=count(m,box())
        for mutation in ['gap','overlap','orientation','center','missing']:
            cc=copy.deepcopy(c)
            if mutation=='gap':cc['panels'][0]['left']='1/10'
            elif mutation=='overlap':cc['panels'].insert(0,cc['panels'][0])
            elif mutation=='orientation':cc['contour']=cc['contour'][::-1]
            elif mutation=='center':cc['panels'][0]['center']='100000'
            else:cc['panels'].pop()
            with self.assertRaises((ValueError,Unresolved)):
                verify(cc)

    def test_boundary_zero_unresolved(self):
        m=[{'lambda':'0','poly':[['0','0'],['1','0']]}]
        with self.assertRaises(Unresolved):
            count(m,box(0,2,-1,1),max_depth=10)

    def test_optimizer_majorant_near_grid_minimum(self):
        from flint import arb
        rng=random.Random(640)
        for _ in range(15):
            m=[{'lambda':str(l),'poly':[[str(rng.randint(1,20)),'0']]}
               for l in [-3,0,1,4]]
            a=('0','0');b=('1/10','1/10')
            s=choose_center(m,a,b,'local')
            lam,cs,d=data(m,a,b)
            def bound(v):
                from contourgauge.core import real
                return sum(float(abs(c[0])*(abs(real(l-v)*d)).expm1())
                           for l,c in zip(lam,cs))
            grid=min(bound(Q(-3)+Q(7*k,500)) for k in range(501))
            self.assertLessEqual(bound(s),grid+1e-8)

    def test_cancellation_preserving_certificate(self):
        m=[{'lambda':'0','poly':[['1','0'],['1','0']]},
           {'lambda':'-2','poly':[['4','0']]}]
        r,c=count(m,box(0,3,-20,20),'jet')
        self.assertEqual(r['count'],2)
        self.assertEqual(verify(c,192)['count'],2)
        self.assertEqual(c['enclosure_order'],2)
        c['enclosure_order']=3
        with self.assertRaises(ValueError):verify(c)

    def test_second_order_envelope_contains_sampled_curve(self):
        from flint import acb,arb
        from contourgauge.core import enclosure
        m=[{'lambda':'-2','poly':[['1','2'],['2','-1'],['0','1'],['1','0']]},
           {'lambda':'3','poly':[['-2','1'],['1','1']]}]
        F,mu,B,cs,_=enclosure(m,('1/3','1/4'),('2/3','3/4'),('1/2','1/3'),2)
        for k in range(-19,20):
            t=arb(k)/20
            G=sum((sum((v*t**n for n,v in enumerate(c)),acb(0))*(u*t).exp()
                   for c,u in zip(cs,mu)),acb(0))
            self.assertTrue(abs(G-F)<B)

    def test_jet_shift_invariance(self):
        m=[{'lambda':'0','poly':[['-1','0']]},{'lambda':'1','poly':[['1','0']]}]
        r,c=count(m,box(-1,2,-7,7),'jet')
        mm=copy.deepcopy(m)
        for t in mm:t['lambda']=str(Q(t['lambda'])+100)
        rr,cc=count(mm,box(-1,2,-7,7),'jet')
        self.assertEqual(rr['count'],r['count'])
        self.assertEqual([(p['left'],p['right']) for p in c['panels']],
                         [(p['left'],p['right']) for p in cc['panels']])

if __name__=='__main__':unittest.main()
