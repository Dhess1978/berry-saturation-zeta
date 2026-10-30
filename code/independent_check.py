# Independent check of compute_rigidity.py with different algorithms, from the same text files:
#   primes by trial division + segmented sieve; the prime tail above 1e6 by direct quadrature of int (1 - cos au)/u du (not via Ci);
#   Delta3 by explicit normal equations (numpy.linalg.solve); Delta3 from Berry's Sigma2 by quadrature over r (not term by term);
#   GUE Sigma2 by numerical integration of L - 2 int_0^L (L - u)(sin pi u/pi u)^2 du; Odlyzko phases B ln p mod 2 pi with mpmath (50 digits).
# Sigma2 windows (starts spaced by 0.37) and the Cramer generator are definitions and therefore the same.
# Usage: python3 -I independent_check.py zetazeros100000.txt zeros6.txt zeros3.txt zeros4.txt zeros5.txt rigidity_reference.json
#   prints the largest differences against rigidity_reference.json
import sys, re, json, math, hashlib, numpy as np, mpmath
from scipy.special import sici
from scipy.integrate import quad
F100, Z6, O3, O4, O5, REF = sys.argv[1:7]; mpmath.mp.dps = 50
LG = [0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 80, 100]; LD = [2, 5, 10, 20, 30, 50, 100]; P0 = 10**6; U0 = math.log(P0)
t100 = np.loadtxt(F100, usecols=1); z6 = np.loadtxt(Z6)
Ns = lambda t: t/(2*math.pi)*(np.log(t/(2*math.pi)) - 1) + 7/8
# primes: trial division up to 1000, then a segmented sieve in blocks of 1e5
mala = [p for p in range(2, 1001) if all(p % q for q in range(2, int(p**0.5) + 1))]; PR = list(mala)
for lo in range(1001, P0 + 1, 100000):
    hi = min(lo + 100000, P0 + 1); s = np.ones(hi - lo, bool)
    for p in mala:
        st = max(p*p, ((lo + p - 1)//p)*p)
        if st < hi: s[st - lo::p] = False
    PR += list(np.nonzero(s)[0] + lo)
PR = np.array(PR, float)
def lcg_cramer(seed=20261030):
    s = seed; o = [2.0]
    for n in range(3, P0 + 1):
        s = (1664525*s + 1013904223) % 4294967296
        if (s + 0.5)/4294967296 < 1/math.log(n): o.append(float(n))
    return np.array(o)
CR = lcg_cramer()
def mocn(pr, H):   # (u, w) for p^k <= e^H
    lp = np.log(pr); U, W = [], []
    k = 1
    while True:
        m = k*lp <= H
        if not m.any(): break
        U.append(k*lp[m]); W.append(np.exp(-k*lp[m])/(k*k)); k += 1
    return np.concatenate(U), np.concatenate(W)
def ocas(L, H):
    if H <= U0: return 0.0
    a = 2*math.pi*L/H; v, e = quad(lambda u: (1 - math.cos(a*u))/u, U0, H, limit=2000, epsabs=1e-14, epsrel=1e-13); return v
CACHE = {}
def berry1(L, H, pr=PR, tag='p'):
    k = (tag, H)
    if k not in CACHE: CACHE[k] = mocn(pr, H)
    U, W = CACHE[k]; a = 2*math.pi*L/H; s = float(np.sum(W*(1 - np.cos(a*U)))) + ocas(L, H)
    si = sici(2*math.pi*L)[0]; return s/math.pi**2 + (1 - math.cos(2*math.pi*L))/math.pi**2 + L*(1 - 2/math.pi*si)
def berryB(L, Hs, pr=PR, tag='p'): return float(np.mean([berry1(L, H, pr, tag) for H in Hs]))
def sig2(x, L):
    st = x[0] + 0.37*np.arange(int(np.ceil((x[-1] - L - x[0])/0.37))); st = st[st < x[-1] - L]
    c = np.searchsorted(x, st + L, side='right') - np.searchsorted(x, st, side='right'); return float(np.var(c.astype(float)))
def delta3(x, L, nw=800):
    o = []
    for j in range(nw):
        a = x[0] + (x[-1] - L - x[0])*(j + 0.5)/nw; y = x[(x > a) & (x <= a + L)] - a; y = np.concatenate([[0.0], y, [L]]); n = np.arange(len(y) - 1, dtype=float); dy = np.diff(y)
        I0 = np.sum(n*dy); I1 = np.sum(n*np.diff(y*y)/2); I2 = np.sum(n*n*dy)
        A, B = np.linalg.solve(np.array([[L, L*L/2], [L*L/2, L**3/3]]), np.array([I0, I1])); o.append((I2 - A*I0 - B*I1)/L)
    return float(np.mean(o))
def d3berry(L, H):
    v, e = quad(lambda r: (L**3 - 2*L*L*r + r**3)*berry1(r, H), 0, L, limit=4000, epsabs=1e-12, epsrel=1e-11); return 2/L**4*v
def gueS(L):
    v, e = quad(lambda u: (L - u)*(np.sinc(u))**2, 0, L, limit=4000, epsabs=1e-13); return L - 2*v
def mez(x): d = np.diff(x); d = d/d.mean(); return float(np.mean(d < 0.3)), float(d.std())
def blok(t, x, Hs, LGx, d3=True, cram=True):
    o = {'Hs': [float(h) for h in Hs], 'sig': [sig2(x, L) for L in LGx], 'berry': [berryB(L, Hs) for L in LGx]}
    if cram: o['cram'] = [berryB(L, Hs, CR, 'c') for L in LGx]
    if d3: o['d3'] = [delta3(x, L) for L in LD]
    o['pod03'], o['sdd'] = mez(x); iL = [LGx.index(L) for L in (20, 25, 30, 40, 50)]; o['nas'] = float(np.mean([o['sig'][i] for i in iL])); return o
R = {'nprim': int(len(PR)), 'ncram': int(len(CR)), 'gue': [gueS(L) for L in LG], 'B': {}}
F = json.load(open(REF)); FB = {b['id']: b for b in F['blocks']}; D = {}
def cmp(name, a, b):
    d = float(np.max(np.abs(np.array(a, float) - np.array(b, float)))); D[name] = max(D.get(name, 0.0), d)
for i in (0, 2, 4):
    t = t100[20000*i:20000*(i + 1)]; x = Ns(t); Hs = np.log(t[np.linspace(0, 19999, 48).astype(int)]/(2*math.pi)); o = blok(t, x, Hs, LG); r = FB['lmfdb%d' % i]
    cmp('sigma2', o['sig'], r['sigma2']); cmp('berry', o['berry'], r['berry']); cmp('cramer', o['cram'], r['cramer']); cmp('delta3', o['d3'], r['delta3']); cmp('small_gap_share', [o['pod03']], [r['small_gap_share']]); cmp('heights', o['Hs'], r['heights_H'])
    if i == 0: cmp('delta3_berry (48 heights, quadrature)', [sum(d3berry(L, H) for H in Hs)/len(Hs) for L in (10, 30)], [r['delta3_berry'][2], r['delta3_berry'][4]])
    print('lmfdb', i, flush=True)
for i in (0, 9):
    t = z6[200000*i:200000*(i + 1)]; x = Ns(t); Hs = np.log(t[np.linspace(0, len(t) - 1, 48).astype(int)]/(2*math.pi)); o = blok(t, x, Hs, LG, d3=False); r = FB['zeros6_%d' % i]
    cmp('sigma2', o['sig'], r['sigma2']); cmp('berry', o['berry'], r['berry']); cmp('cramer', o['cram'], r['cramer']); cmp('heights', o['Hs'], r['heights_H']); print('zeros6', i, flush=True)
def nacti(f):
    txt = open(f).read(); B = int(re.search(r'gamma - (\d+)', txt).group(1)); e = int(re.search(r'numbered 10\^(\d+)', txt).group(1))
    d = np.array([float(l) for l in txt.split('\n') if re.fullmatch(r'\s*\d+\.\d+\s*', l)], '<f8'); return B, e, d
for f in (O3, O4, O5):
    B, e, d = nacti(f); H = float(mpmath.log(mpmath.mpf(B)/(2*mpmath.pi))); x = d*H/(2*math.pi) + d*d/(4*math.pi*B); r = FB['odlyzko_1e%d' % e]
    o = blok(None, x, [H], LG[:17]); S = np.arange(len(x)) - x; Sc = S - S.mean()
    pom = [float(2*np.mean(Sc*np.sin(float(mpmath.fmod(mpmath.mpf(B)*mpmath.log(p), 2*mpmath.pi)) + d*math.log(p))))*(-math.pi*math.sqrt(p)) for p in (2, 3, 5, 7, 11, 13)]
    cmp('sigma2', o['sig'], r['sigma2']); cmp('berry', o['berry'], r['berry']); cmp('cramer', o['cram'], r['cramer']); cmp('delta3', o['d3'], r['delta3']); cmp('small_gap_share', [o['pod03']], [r['small_gap_share']])
    cmp('a_p_over_explicit_formula', pom, [r['a_p_over_explicit_formula'][str(p)] for p in (2, 3, 5, 7, 11, 13)])
    if e == 22: cmp('delta3_berry (quadrature, with prime tail)', [d3berry(L, H) for L in (10, 30)], [r['delta3_berry'][2], r['delta3_berry'][4]])
    print('odlyzko', e, flush=True)
cmp('gue_sigma2', R['gue'], F['gue_sigma2'])
x6 = Ns(z6); Hs6 = np.log(z6[np.linspace(0, len(z6) - 1, 96).astype(int)]/(2*math.pi)); cmp('zeros6_whole sigma2(1000)', [sig2(x6, 1000)], [F['zeros6_whole']['sigma2'][3]]); cmp('zeros6_whole berry(1000)', [berryB(1000, Hs6)], [F['zeros6_whole']['berry'][3]])
print('primes up to 1e6:', R['nprim'], '· Cramer primes:', R['ncram'])
for k, v in D.items(): print('max |difference| %-45s %.1e' % (k, v))
print('ALL WITHIN 1e-8' if max(D.values()) < 1e-8 else 'DIFFERENCES ABOVE 1e-8')
