# ContourGauge

**Local Exponential Gauges for Validated Contour Zero Counts** — Shivam Gupta.

This repository accompanies a computational complex analysis manuscript. It
contains a small validated zero-counting implementation for explicit exponential
polynomials with real rational exponents and complex rational polynomial
coefficients. The public verifier supports simple counterclockwise rectangles.

The method removes a potentially different exponential factor on each contour
panel, proves a zero-excluding disk inequality, and restores the removed phase
analytically. A heuristic scalar optimization proposes gauges. Acceptance,
phase enclosures, partition coverage and integer isolation use Arb ball
arithmetic through python-flint. A failed certificate is **unresolved**, never
interpreted as evidence of a zero.

## Reproduce

Python 3.9 or newer:

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -e . matplotlib==3.9.4
python -m unittest discover -s tests -v
python experiments/study.py
python experiments/global_sweep.py
python experiments/cancellation.py
python experiments/figures.py
python experiments/cancellation_figures.py
```

The main study constructs 360 certificates on 90 deterministic problems and
replays each at 192 bits. References use explicit logarithmic zero lattices or
independent validated integration of `f'/f`. Single-delay examples also agree
with the analytic imaginary-axis crossing count. `global_sweep.py` constructs
the offline best of nine fixed-center comparators and records replay timings.
Timings vary with hardware and load. They are not acceptance conditions.

The separately reported follow-up ablation tests a second-order enclosure on
all 36 delay problems. The `jet` policy proposes a complex logarithmic-derivative
gauge, rounded to exact rational components. Every comparator receives the
same enclosure options. The resulting panel reductions are modest and do not
improve median replay time against the tested best fixed gauges. All 625 saved
main, fixed-sweep and follow-up certificate records pass the final 192-bit
replay audit (`results/final_replay_audit.json`). Nine test methods pass.

Compile the manuscript with `cd paper && tectonic main.tex`. The figures and
numerical table inputs are already included. `paper/result_text.tex` records
the interpretation of the archived results; update it if changing the study.

## Use

```python
from contourgauge import count, verify

model = [
    {'lambda': '0', 'poly': [['1', '0'], ['1', '0']]},
    {'lambda': '-1', 'poly': [['2', '0']]},
]
contour = [['-2','-5'], ['2','-5'], ['2','5'], ['-2','5'], ['-2','-5']]
result, cert = count(model, contour, policy='local')
assert result['count'] == 2
assert verify(cert, precision=192)['count'] == 2
# Optional cancellation-preserving, complex-gauge proposal:
jet_result, jet_cert = count(model, contour, policy='jet')
```

Polynomial coefficients are in ascending powers. Each complex rational is
`[real, imag]`; strings preserve exact input intent. Exponents are rational
strings. Save certificates as JSON, then call `verify` to replay. The claimed
count is recomputed rather than trusted. Verify the model and contour identity
when consuming a certificate supplied by someone else.

## Contents

- `src/contourgauge`: search, analytic panel enclosure, and certificate replay.
- `tests`: mathematical ground-truth, invariance, optimizer and mutation tests.
- `experiments`: deterministic study, fixed-center oracle and plotting scripts.
- `results`: exact inputs, raw results, reference enclosures and panel records.
- `figures`: vector PDF figures and PNG previews.
- `paper`: complete LaTeX source, bibliography and manuscript PDF.
- `RESEARCH_NOTES.md`: literature search, novelty boundary and practical scope.

## What the evidence supports

Local gauges can shorten these disk certificates. Construction can be slower
because the optimizer has a cost. The paper reports this tradeoff directly.
There is no claim of superiority over all argument-principle solvers, no
claim that exponential normalization itself is new, and no guarantee for
physical models beyond their exact encoded mathematical equations.

Search and replay share panel code, and both the verifier and independent
quadrature references share Arb. This is not proof-assistant verification or
independent arithmetic-library validation. Boundary zeros are outside the
certificate's acceptance conditions.

## Licenses

Code: MIT (`LICENSE`). Manuscript, numerical data and figures: CC BY 4.0
(`LICENSE-PAPER`). No third-party source papers are redistributed.
