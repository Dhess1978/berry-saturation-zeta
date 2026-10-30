# Number variance and spectral rigidity of the Riemann zeros against Berry's semiclassical formula, Cramer random primes and GUE.
#
# Inputs (text files, not included): the first 100,000 zeros from LMFDB (two columns: index, t), Odlyzko's table of the first 2,001,052 zeros
# (one t per line) and Odlyzko's tables near the zeros number 1e12, 1e21, 1e22 (offsets from a base B, as distributed).
# Output: rigidity_reference.json (all numbers of the article).
#
#   Sigma2(L)   variance of the number of unfolded zeros in a window of length L, window starts spaced by 0.37
#   Delta3(L)   least-squares rigidity, exact for the staircase, mean over 800 uniformly spaced windows
#   Berry       Sigma2_B(L; t) = (1/pi^2) sum_{p^k <= t/2pi} [1 - cos(2 pi L k ln p / H)]/(k^2 p^k) + (1 - cos 2 pi L)/pi^2 + L[1 - (2/pi) Si(2 pi L)],
#               H = ln(t/2pi); primes up to 1e6 explicitly, above by the prime number theorem; mean over 48 heights of a block (1 for Odlyzko)
#   Cramer      the same formula with random primes: n chosen with probability 1/ln n (LCG x -> 1664525 x + 1013904223 mod 2^32, seed 20261030)
#   Delta3 from Berry's Sigma2 exactly term by term: (2/L^4) int_0^L (L^3 - 2L^2 r + r^3)(1 - cos br) dr = 2 A(bL),
#               A(x) = 1/4 - (cos x + 2)/x^2 + 6 sin x/x^3 - 6(1 - cos x)/x^4 (series for x < 2); smooth part by quadrature
#   a(p)        projection of S_n = n - x_n on sin(t_n ln p) (Odlyzko heights: phase B ln p mod 2 pi to 60 digits), ratio to -1/(pi sqrt p)
#
# Usage: python3 -I compute_rigidity.py zetazeros100000.txt zeros6.txt zeros3.txt zeros4.txt zeros5.txt rigidity_reference.json
import sys, re, json, math, time, numpy as np, mpmath
from scipy.special import sici
from scipy.integrate import quad, IntegrationWarning
import warnings; warnings.filterwarnings('ignore', category=IntegrationWarning)   # quad reports round-off at the 1e-13 tolerance
F100, Z6, O3, O4, O5, OUT = sys.argv[1:7]; mpmath.mp.dps = 60
t100 = np.array([float(l.split()[1]) for l in open(F100) if l.strip()]); z6 = np.loadtxt(Z6)
Ns = lambda t: t/(2*np.pi)*(np.log(t/(2*np.pi)) - 1) + 7/8
LG = [0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 80, 100]; LD = [2, 5, 10, 20, 30, 50, 100]; LV = [200, 300, 500, 1000]
P0 = 10**6; U0 = math.log(P0)
def sieve(M):
    s = np.ones(M + 1, bool); s[:2] = False
    for i in range(2, int(M**0.5) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]
PR = sieve(P0).astype(float)
def lcg(seed):
    s = seed & 0xffffffff
    while True:
        s = (1664525*s + 1013904223) & 0xffffffff; yield (s + 0.5)/4294967296
def cramer(M, seed=20261030):
    g = lcg(seed); o = [2.0]
    for n in range(3, M + 1):
        if next(g) < 1/math.log(n): o.append(float(n))
    return np.array(o)
CR = cramer(P0)
def powers(pr, Hmax):
    lp = np.log(pr); u, w = [], []
    for k in range(1, 80):
        m = k*lp <= Hmax
        if not m.any(): break
        u.append(k*lp[m]); w.append(1/(k*k*np.exp(k*lp[m])))
    u = np.concatenate(u); w = np.concatenate(w); o = np.argsort(u, kind='stable'); return u[o], w[o]
UP, WP = powers(PR, 60.0); UC, WC = powers(CR, 60.0)
def berry1(L, H, U=UP, Wt=WP):
    m = U <= H; a = 2*math.pi*L/H; s = float(np.sum(Wt[m]*(1 - np.cos(a*U[m]))))
    if H > U0: s += (math.log(H) - math.log(U0)) - (sici(a*H)[1] - sici(a*U0)[1])
    si, ci = sici(2*math.pi*L); return s/math.pi**2 + (1 - math.cos(2*math.pi*L))/math.pi**2 + L*(1 - 2/math.pi*si)
def gue(L): si, ci = sici(2*math.pi*L); return (math.log(2*math.pi*L) + 0.5772156649015329 + 1 - math.cos(2*math.pi*L) - ci)/math.pi**2 + L*(1 - 2/math.pi*si)
def berryB(L, Hs, U=UP, Wt=WP): return float(np.mean([berry1(L, H, U, Wt) for H in Hs]))
def sig2(x, L):
    st = x[0] + 0.37*np.arange(int(np.ceil((x[-1] - L - x[0])/0.37))); st = st[st < x[-1] - L]
    c = np.searchsorted(x, st + L, 'right') - np.searchsorted(x, st, 'right'); return float(np.mean(c*c.astype(float)) - np.mean(c)**2)
def delta3(x, L, nw=800):
    out = np.empty(nw)
    for j in range(nw):
        a = x[0] + (x[-1] - L - x[0])*(j + 0.5)/nw; j0, j1 = np.searchsorted(x, a, 'right'), np.searchsorted(x, a + L, 'right'); y = np.concatenate([[0.0], x[j0:j1] - a, [L]]); n = np.arange(len(y) - 1, dtype=float)
        d1 = y[1:] - y[:-1]; d2 = (y[1:]**2 - y[:-1]**2)/2; I0 = float(np.sum(n*d1)); I1 = float(np.sum(n*d2)); I2 = float(np.sum(n*n*d1))
        A11, A12, A22 = L**3/3, L**2/2, L; det = A11*A22 - A12*A12; aa = (I1*A22 - I0*A12)/det; bb = (A11*I0 - A12*I1)/det; out[j] = (I2 - aa*I1 - bb*I0)/L
    return float(out.mean())
KS = np.arange(3, 25); KC = np.array([(-1)**(k + 1)*2*(2*k + 1)*(k - 2)/math.factorial(2*k + 2) for k in KS])
def A(x):
    x = np.atleast_1d(np.asarray(x, float)); out = np.empty_like(x); m = x < 2
    xs = x[m]; out[m] = np.sum(KC[None, :]*xs[:, None]**(2*KS[None, :] - 2), axis=1)
    xb = x[~m]; c, s = np.cos(xb), np.sin(xb); out[~m] = 0.25 - (c + 2)/xb**2 + 6*s/xb**3 - 6*(1 - c)/xb**4
    return out
def smooth_d3(L):
    f = lambda r: (L**3 - 2*L*L*r + r**3)*r*(1 - 2/math.pi*sici(2*math.pi*r)[0])
    s, e = quad(f, 0, L, limit=4000, epsabs=1e-13, epsrel=1e-13); return 2/L**4*s
def d3_berry(L, Hs):
    tot = 0.0
    for H in Hs:
        m = UP <= H; s = float(np.sum(WP[m]*A(2*math.pi*UP[m]/H*L)))
        if H > U0: t, e = quad(lambda u: float(A(2*math.pi*u*L/H)[0])/u, U0, H, limit=4000, epsabs=1e-14, epsrel=1e-13); s += t
        tot += 2*s/math.pi**2
    return tot/len(Hs) + 2*float(A(2*math.pi*L)[0])/math.pi**2 + smooth_d3(L)
def block(bid, source, t, x, Hs, LGx=LG):
    sig = [sig2(x, L) for L in LGx]; br = [berryB(L, Hs) for L in LGx]; cr = [berryB(L, Hs, UC, WC) for L in LGx]
    dd = np.diff(x); dd = dd/dd.mean(); m = np.array(sig)
    o = {'id': bid, 'source': source, 'zeros': int(len(x)), 'heights_H': [float(h) for h in Hs], 'sigma2': sig, 'berry': br, 'cramer': cr,
         'small_gap_share': float((dd < 0.3).mean()), 'gap_sd': float(dd.std()), 'delta3': [delta3(x, L) for L in LD], 'delta3_berry': [d3_berry(L, Hs) for L in LD],
         'rms_berry': float(np.sqrt(np.mean(((np.array(br) - m)/m)**2))), 'rms_cramer': float(np.sqrt(np.mean(((np.array(cr) - m)/m)**2)))}
    if t is not None: o['t_min'], o['t_max'] = float(t[0]), float(t[-1])
    return o
R = {'L_sigma2': LG, 'L_delta3': LD, 'L_large': LV, 'prime_cutoff': P0, 'gue_sigma2': [gue(L) for L in LG], 'gue_sigma2_large': [gue(L) for L in LV], 'blocks': []}
t0 = time.time()
for b in range(5):
    t = t100[20000*b:20000*(b + 1)]; Hs = np.log(t[np.linspace(0, len(t) - 1, 48).astype(int)]/(2*np.pi)); R['blocks'].append(block('lmfdb%d' % b, 'LMFDB', t, Ns(t), Hs)); print('lmfdb', b, round(time.time() - t0), flush=True)
for b in range(10):
    t = z6[200000*b:200000*(b + 1)]; Hs = np.log(t[np.linspace(0, len(t) - 1, 48).astype(int)]/(2*np.pi)); R['blocks'].append(block('zeros6_%d' % b, 'Odlyzko zeros6', t, Ns(t), Hs)); print('zeros6', b, round(time.time() - t0), flush=True)
x = Ns(z6); Hs = np.log(z6[np.linspace(0, len(z6) - 1, 96).astype(int)]/(2*np.pi))
R['zeros6_whole'] = {'zeros': int(len(z6)), 'heights_H': [float(h) for h in Hs], 'sigma2': [sig2(x, L) for L in LV], 'berry': [berryB(L, Hs) for L in LV]}
def read_odlyzko(f):
    txt = open(f).read(); B = int(re.search(r'gamma - (\d+)', txt).group(1)); e = int(re.search(r'numbered 10\^(\d+)', txt).group(1))
    d = np.array([float(l) for l in txt.split('\n') if re.fullmatch(r'\s*\d+\.\d+\s*', l)]); return B, e, d
for f in (O3, O4, O5):
    B, e, d = read_odlyzko(f); H = math.log(B/(2*math.pi)); x = d*H/(2*math.pi) + d*d/(4*math.pi*B); S = np.arange(len(x)) - x; Sc = S - S.mean()
    phase = {p: float(mpmath.fmod(mpmath.mpf(B)*mpmath.log(p), 2*mpmath.pi)) for p in (2, 3, 5, 7, 11, 13)}
    a = {p: float(2*np.mean(Sc*np.sin(phase[p] + d*math.log(p)))) for p in phase}
    o = block('odlyzko_1e%d' % e, 'Odlyzko near zero no. 1e%d' % e, None, x, [H], LG[:17])
    o.update({'base_B': str(B), 'zero_number_exponent': e, 'prime_phase': {str(p): v for p, v in phase.items()}, 'a_p': {str(p): v for p, v in a.items()},
              'a_p_over_explicit_formula': {str(p): a[p]*(-math.pi*math.sqrt(p)) for p in a}, 'sd_S': float(Sc.std())})
    R['blocks'].append(o); print('odlyzko', e, round(time.time() - t0), flush=True)
json.dump(R, open(OUT, 'w'), indent=None); print('done', round(time.time() - t0), 's')
