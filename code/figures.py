# Figures of the article (light background, English and Czech), PNG 300 dpi and PDF (vector), from rigidity_reference.json.
# Usage: python3 -I figures.py rigidity_reference.json OUTPUT_FOLDER
import sys, json, math, os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.special import sici
REF, OUT = sys.argv[1:3]; os.makedirs(OUT, exist_ok=True)
J = json.load(open(REF)); LG = np.array(J['L_sigma2']); LD = np.array(J['L_delta3']); G = np.array(J['gue_sigma2'])
def _blk(b): return {'sig': b['sigma2'], 'berry': b['berry'], 'cramer': b['cramer'], 'd3': b['delta3'], 'd3b': b['delta3_berry'], 'pod03': b['small_gap_share'], 'n': b['zeros'], 'H': b['heights_H'],
                     'B': b.get('base_B'), 'pomer': b.get('a_p_over_explicit_formula')}
JB = {b['id']: b for b in J['blocks']}
R = {'lmfdb': [_blk(JB['lmfdb%d' % i]) for i in range(5)], 'z6': [_blk(JB['zeros6_%d' % i]) for i in range(10)], 'odl': {k: _blk(JB['odlyzko_' + k]) for k in ('1e12', '1e21', '1e22')}}
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9, 'legend.fontsize': 8, 'axes.edgecolor': '#555555',
                     'axes.linewidth': 0.8, 'xtick.color': '#333333', 'ytick.color': '#333333', 'axes.grid': True, 'grid.color': '#e3e3e3', 'grid.linewidth': 0.6,
                     'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white', 'legend.frameon': False})
VIR = plt.get_cmap('viridis'); HC = [VIR(x) for x in (0.02, 0.27, 0.50, 0.70, 0.86)]
INK, MUT, GUEC, CRC = '#222222', '#6b6b6b', '#222222', '#b2589a'
def cz(s): return s.replace('.', ',')
def fmt(x, d, en): s = ('%.*f' % (d, x)); return s if en else cz(s)
def hor(e): return '10' + ''.join('⁰¹²³⁴⁵⁶⁷⁸⁹'[int(c)] for c in str(e))
def tlab(t, en):
    e = int(math.floor(math.log10(t))); m = t/10**e; return ('%.1f' % m if en else cz('%.1f' % m)) + '·' + hor(e)
def gue(L): si, ci = sici(2*math.pi*L); return (math.log(2*math.pi*L) + 0.5772156649015329 + 1 - math.cos(2*math.pi*L) - ci)/math.pi**2 + L*(1 - 2/math.pi*si)
def gueD3(L): return (math.log(2*math.pi*L) + 0.5772156649015329 - 1.25)/(2*math.pi**2)
# Berry's saturation level on a grid of heights (single height; primes up to 1e6 explicitly, above by the prime number theorem)
def sito(M):
    s = np.ones(M + 1, bool); s[:2] = False
    for i in range(2, int(M**0.5) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0].astype(float)
PR = sito(10**6); lp = np.log(PR); U, W = [], []
for k in range(1, 80):
    m = k*lp <= 60
    if not m.any(): break
    U.append(k*lp[m]); W.append(np.exp(-k*lp[m])/(k*k))
U = np.concatenate(U); W = np.concatenate(W)
def berry1(L, H):
    m = U <= H; a = 2*math.pi*L/H; s = float(np.sum(W[m]*(1 - np.cos(a*U[m]))))
    u0 = math.log(1e6)
    if H > u0: s += (math.log(H) - math.log(u0)) - (sici(a*H)[1] - sici(a*u0)[1])
    si = sici(2*math.pi*L)[0]; return s/math.pi**2 + (1 - math.cos(2*math.pi*L))/math.pi**2 + L*(1 - 2/math.pi*si)
BL = [R['lmfdb'][i] | {'id': 'l%d' % i, 'zdroj': 'LMFDB'} for i in range(5)] + [R['z6'][i] | {'id': 'z%d' % i, 'zdroj': 'zeros6'} for i in range(10)] + [R['odl'][k] | {'id': 'o' + k[2:], 'zdroj': 'Odlyzko'} for k in ('1e12', '1e21', '1e22')]
for b in BL:
    b['tS'] = math.exp(float(np.mean(b['H'])))*2*math.pi; nL = len(b['sig']); iL = [list(LG[:nL]).index(L) for L in (20, 25, 30, 40, 50)]
    b['nas'] = float(np.mean([b['sig'][i] for i in iL])); b['nasB'] = float(np.mean([b['berry'][i] for i in iL]))
B = {b['id']: b for b in BL}; VYS = ['l2', 'z9', 'o12', 'o21', 'o22']
ZC = {'LMFDB': '#3b6fb6', 'zeros6': '#2a9d6f', 'Odlyzko': '#d07a1a'}
from matplotlib.ticker import FuncFormatter
CZF = FuncFormatter(lambda v, p: ('%g' % round(v, 6)).replace('.', ',').replace('-', '−'))
def carka(fig):
    for ax in fig.axes:
        if ax.get_yscale() == 'linear': ax.yaxis.set_major_formatter(CZF)
        if ax.get_xscale() == 'linear': ax.xaxis.set_major_formatter(CZF)
def uloz(fig, jm):
    if jm.endswith('_cz'): carka(fig)
    fig.savefig(os.path.join(OUT, jm + '.png'), dpi=300, bbox_inches='tight'); fig.savefig(os.path.join(OUT, jm + '.pdf'), bbox_inches='tight'); plt.close(fig)
for en in (False, True):
    sfx = '_en' if en else '_cz'; w = lambda c, e: e if en else c
    # 1) Sigma2(L) at five heights
    fig, ax = plt.subplots(figsize=(6.4, 4.0)); Lf = np.exp(np.linspace(math.log(0.5), math.log(100), 200))
    ax.plot(Lf, [gue(L) for L in Lf], ls='--', c=GUEC, lw=1.2, label=w('GUE (teorie)', 'GUE (theory)'))
    for k, id_ in enumerate(VYS):
        b = B[id_]; n = len(b['sig']); ax.plot(LG[:n], b['berry'], c=HC[k], lw=1.4); ax.plot(LG[:n], b['sig'], 'o', ms=3.6, c=HC[k], mec='white', mew=0.5, label='t ≈ ' + tlab(b['tS'], en))
    ax.set_xscale('log'); ax.set_xticks([0.5, 1, 2, 5, 10, 20, 50, 100]); ax.set_xticklabels([fmt(x, 1, en).rstrip('0').rstrip(',.') for x in [0.5, 1, 2, 5, 10, 20, 50, 100]])
    ax.set_xlabel(w('délka okna L (rozvinuté jednotky)', 'window length L (unfolded units)')); ax.set_ylabel('Σ²(L)'); ax.set_ylim(0, 0.85)
    ax.set_title(w('Number variance nul ζ na pěti výškách: body = nuly, čáry = Berryho vzorec', 'Number variance of the zeta zeros at five heights: dots = zeros, lines = Berry formula'), loc='left')
    ax.legend(loc='upper left', ncol=2); uloz(fig, '01_number_variance_heights' + sfx)
    # 3) saturation level by height
    fig, ax = plt.subplots(figsize=(6.4, 3.8)); ee = np.arange(3, 22.01, 0.25); kr = [np.mean([berry1(L, math.log(10**e/(2*math.pi))) for L in (20, 25, 30, 40, 50)]) for e in ee]
    ax.plot(10**ee, kr, c=INK, lw=1.3, label=w('Berryho vzorec (jedna výška)', 'Berry formula (single height)'))
    nG = float(np.mean([gue(L) for L in (20, 25, 30, 40, 50)])); ax.axhline(nG, ls='--', c=GUEC, lw=1.1, label='GUE')
    for z in ('LMFDB', 'zeros6', 'Odlyzko'):
        bb = [b for b in BL if b['zdroj'] == z]; ax.plot([b['tS'] for b in bb], [b['nas'] for b in bb], 'o', ms=5, c=ZC[z], mec='white', mew=0.6, label=z)
    ax.set_xscale('log'); ax.set_xlabel(w('výška t', 'height t')); ax.set_ylabel(w('průměr Σ²(L) pro L = 20 … 50', 'mean Σ²(L) for L = 20 … 50')); ax.set_ylim(0.25, 0.72)
    ax.set_title(w('Úroveň nasycení roste s výškou (víc prvočísel pod Heisenbergovým časem)', 'The saturation level grows with height (more primes below the Heisenberg time)'), loc='left')
    ax.legend(loc='lower right'); uloz(fig, '03_saturation_level' + sfx)
    # 2) deviation from Berry and from random primes: lowest block and 1e22
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.2), sharey=True)
    for ax, id_, tt in ((axs[0], 'z0', w('zeros6 č. 1–200 000', 'zeros6 no. 1-200 000')), (axs[1], 'o22', w('Odlyzko č. 10²² + 1 …', 'Odlyzko no. 10²² + 1 ...'))):
        b = B[id_]; n = len(b['sig']); s = np.array(b['sig']); ax.axhspan(-1, 1, color='#e9f3ec', lw=0); ax.axhline(0, c=MUT, lw=0.8)
        ax.plot(LG[:n], 100*(s/np.array(b['berry']) - 1), 'o-', ms=3.2, lw=1.1, c=HC[1], label='Berry'); ax.plot(LG[:n], 100*(s/np.array(b['cramer']) - 1), 's--', ms=3.0, lw=1.0, c=CRC, label=w('Cramér (náhodná prvočísla)', 'Cramer (random primes)'))
        ax.set_xscale('log'); ax.set_title(tt + ' (t ≈ ' + tlab(b['tS'], en) + ')', loc='left', fontsize=9); ax.set_xlabel('L')
    axs[0].set_ylabel(w('nuly / vzorec − 1 (%)', 'zeros / formula − 1 (%)')); axs[0].set_ylim(-12, 4); axs[1].legend(loc='lower right')
    for ax in axs: ax.set_xticks([0.5, 1, 2, 5, 10, 20, 50, 100]); ax.set_xticklabels([fmt(x, 1, en).rstrip('0').rstrip(',.') for x in [0.5, 1, 2, 5, 10, 20, 50, 100]])
    fig.suptitle(w('Přesná poloha prvočísel rozhoduje: Berry do 2 %, náhodná prvočísla se stejnou hustotou o 5–10 %', 'The exact position of the primes matters: Berry within 2 %, random primes of the same density miss by 5–10 %'), x=0.01, ha='left', fontsize=10)
    uloz(fig, '02_berry_vs_random_primes' + sfx)
    # 4) Delta3(L)
    fig, ax = plt.subplots(figsize=(6.4, 3.8)); Lf = np.exp(np.linspace(math.log(2), math.log(100), 100)); ax.plot(Lf, [gueD3(L) for L in Lf], ls='--', c=GUEC, lw=1.2, label='GUE (Dyson–Mehta)')
    for k, id_ in enumerate(VYS):
        b = B[id_]; ax.plot(LD, b['d3b'], c=HC[k], lw=1.4); ax.plot(LD, b['d3'], 'o', ms=3.8, c=HC[k], mec='white', mew=0.5, label='t ≈ ' + tlab(b['tS'], en))
    ax.set_xscale('log'); ax.set_xticks(list(LD)); ax.set_xticklabels([str(x) for x in LD]); ax.set_xlabel('L'); ax.set_ylabel('Δ₃(L)')
    ax.set_title(w('Spektrální tuhost: body = nuly, čáry = Δ₃ z Berryho Σ² (Dyson–Mehta)', 'Spectral rigidity: dots = zeros, lines = Δ₃ from the Berry Σ² (Dyson-Mehta)'), loc='left'); ax.legend(loc='upper left'); uloz(fig, '04_delta3' + sfx)
    # 5) small gaps by height
    fig, ax = plt.subplots(figsize=(6.4, 3.4)); ax.axhline(2.7627, ls='--', c=GUEC, lw=1.1, label='GUE 2,76 %' if not en else 'GUE 2.76 %')
    for z in ('LMFDB', 'zeros6', 'Odlyzko'):
        bb = [b for b in BL if b['zdroj'] == z]; p = np.array([b['pod03'] for b in bb]); n = np.array([b['n'] for b in bb]); se = 100*np.sqrt(p*(1 - p)/(n - 1))
        ax.errorbar([b['tS'] for b in bb], 100*p, yerr=se, fmt='o', ms=4.5, c=ZC[z], mec='white', mew=0.6, elinewidth=1, capsize=0, label=z)
    ax.set_xscale('log'); ax.set_xlabel(w('výška t', 'height t')); ax.set_ylabel('P(δ < 0,3) (%)' if not en else 'P(δ < 0.3) (%)')
    ax.set_title(w('Malé mezery se s výškou pomalu blíží GUE (± binomická chyba)', 'Small gaps slowly approach GUE with height (± binomial error)'), loc='left'); ax.legend(loc='lower right'); uloz(fig, '05_small_gaps' + sfx)
    # 6) primes at large heights
    fig, ax = plt.subplots(figsize=(6.4, 3.2)); P = [2, 3, 5, 7, 11, 13]; ax.axhspan(0.99, 1.01, color='#e9f3ec', lw=0); ax.axhline(1, c=MUT, lw=0.8)
    for j, k in enumerate(('1e12', '1e21', '1e22')):
        o = R['odl'][k]; ax.plot(np.arange(6) + (j - 1)*0.16, [o['pomer'][str(p)] for p in P], 'o', ms=5, c=HC[2 + j], mec='white', mew=0.6, label=w('nuly č. ', 'zeros no. ') + hor(k[2:]) + ' (t ≈ ' + tlab(float(o['B']), en) + ')')
    ax.set_xticks(range(6)); ax.set_xticklabels(['p = %d' % p for p in P]); ax.set_ylim(0.96, 1.04); ax.set_ylabel('a(p) / (−1/(π√p))')
    ax.set_title(w('Prvočísla v S(t) na velkých výškách: explicitní formule do 1,2 %', 'Primes in S(t) at large heights: the explicit formula within 1.2 %'), loc='left'); ax.legend(loc='upper right', fontsize=7.5); uloz(fig, '06_primes_large_heights' + sfx)
    # 7) crossover scale from GUE by height
    def Ls(sig, g, Ls_, thr):
        d = 1 - np.array(sig)/np.array(g)
        for i in range(1, len(d)):
            if d[i] > thr >= d[i - 1]: return Ls_[i - 1]*(Ls_[i]/Ls_[i - 1])**((thr - d[i - 1])/(d[i] - d[i - 1]))
        return float('nan')
    fig, ax = plt.subplots(figsize=(6.0, 3.6)); Hs = np.array([float(np.mean(b['H'])) for b in BL])
    for thr, mk, lab in ((0.05, 'o', w('nuly 5 % pod GUE', 'zeros 5 % below GUE')), (0.10, 's', w('nuly 10 % pod GUE', 'zeros 10 % below GUE'))):
        y = [Ls(b['sig'], G[:len(b['sig'])], LG[:len(b['sig'])], thr) for b in BL]; yb = [Ls(b['berry'], G[:len(b['sig'])], LG[:len(b['sig'])], thr) for b in BL]
        ax.plot(Hs, y, mk, ms=4.5, c=HC[1] if thr == 0.05 else HC[3], mec='white', mew=0.5, label=lab); ax.plot(Hs, yb, mk, ms=7, mfc='none', mec=INK, mew=0.7)
    ax.plot([], [], 'o', ms=7, mfc='none', mec=INK, label=w('prázdné značky: Berryho vzorec (stejný práh)', 'open markers: Berry formula (same threshold)'))
    ax.set_xlabel(w('H = ln(t/2π) (Heisenbergův čas)', 'H = ln(t/2π) (Heisenberg time)')); ax.set_ylabel(w('měřítko přechodu L*', 'crossover scale L*'))
    ax.set_title(w('Přechod od GUE je plynulý a jeho měřítko roste s ln(t/2π)', 'The departure from GUE is a smooth crossover whose scale grows with ln(t/2π)'), loc='left'); ax.legend(loc='upper left'); uloz(fig, '07_crossover_scale' + sfx)
print('hotovo', sorted(os.listdir(OUT)))
