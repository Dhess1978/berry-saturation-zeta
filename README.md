# Berry saturation of the Riemann zeros from t ≈ 10⁴ to 10²¹

**David Hess (MAYAN ALFA V1008 LAB v.34.47) · 30.10.2026**

Number variance Σ²(L) and spectral rigidity Δ₃(L) of the non-trivial zeros of the Riemann zeta function in 18 blocks of zeros at heights from t ≈ 10⁴ (LMFDB) through 10⁶ (2,001,052 zeros) to 10²¹ (Odlyzko's tables near the zeros number 10¹², 10²¹, 10²²), compared with Berry's semiclassical formula — which contains no free parameter — with random primes of the same density (Cramér model) and with GUE.

| result | value |
|---|---|
| Berry's formula against the zeros (RMS relative deviation of Σ²(L), L = 0.5 … 100) | 0.26 – 1.21 % |
| the same formula with random primes of the same density | 5.4 – 9.9 % |
| GUE | 15 – 82 % |
| saturation level (mean Σ² over L = 20 … 50) | 0.34 (t ≈ 10⁴) → 0.56 (t ≈ 10²¹) |
| Δ₃ from Berry's Σ² against the measured Δ₃ | within 3.3 % |
| small gaps P(δ < 0.3) | 1.8 % → 2.6 % (GUE 2.76 %) |
| a(p) / (−1/(π√p)) at t ≈ 10²¹, p = 2 … 13 | 1.000 ± 0.012 |

These are known results (Berry 1988); the contribution is their systematic, reproducible verification over 17 orders of magnitude in height. Nothing here bears on the Riemann hypothesis.

## Contents

| folder | what |
|---|---|
| `paper/` | the article: `berry_saturation_zeta_en.pdf`, LaTeX source `main.tex`, Markdown and DOCX versions |
| `figures/` | 7 figures, English (`_en`) and Czech (`_cz`), PNG 300 dpi and PDF (vector) |
| `data/` | `rigidity_reference.json` (all numbers of Sections 4.1–4.7), `sigma2_by_block.csv`, `sieve_reference.json` (Section 4.8) |
| `code/` | `compute_rigidity.py`, `independent_check.py`, `figures.py` |

## Reproducing

The input zeros are not included; download them from the sources cited in the article: the first 100,000 zeros from LMFDB (`zetazeros100000.txt`), Odlyzko's table of the first 2,001,052 zeros (`zeros6`) and Odlyzko's tables near the zeros number 10¹², 10²¹, 10²² (`zeros3`, `zeros4`, `zeros5`).

```
pip install -r requirements.txt
python3 code/compute_rigidity.py zetazeros100000.txt zeros6 zeros3 zeros4 zeros5 rigidity_reference.json      # about 30 s
python3 code/independent_check.py zetazeros100000.txt zeros6 zeros3 zeros4 zeros5 rigidity_reference.json     # all differences below 1e-10
python3 code/figures.py rigidity_reference.json figures/
```

The computations were developed and verified in the protected MAYAN ALFA V1008 laboratory (internal, version 34.47; mayanalfa.com), which is not part of this repository.

## Licence

Text, figures and data: CC BY 4.0 (`LICENSE-CC-BY-4.0.md`). Code: MIT (`LICENSE`).

## Citation

See `CITATION.cff`.
