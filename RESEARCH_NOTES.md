# Research decisions and evidence boundary

Research date: 17 October 2026. Author: Shivam Gupta.

## Target

ICERM Computational Complex Analysis, 4–8 October 2027, Providence.
Official program: https://icerm.brown.edu/program/semester_program_workshop/sp-f27-w2
Application: https://app.icerm.brown.edu/Cube/apply/640

The observed application asks for dates, optional visa information, a funding
choice, a CV PDF under 2.5 MB, a research statement (2,000 characters), and
additional comments (2,000 characters). It does not provide a dedicated paper
upload, paper template, or proceedings track. Participation applications are
open to industry scientists and other researchers. A workshop application is
not a journal or refereed proceedings submission. Presentation selection and
funding remain at the organizers' discretion.

## Prior manuscript overlap

All 78 top-level PDF manuscripts present when the local scan ran were converted
to text and screened for topic overlap (88 PDFs including subdirectories were
present during a later inventory). Earlier work includes applied AI systems,
policy, medicine, historical studies, and exact periodic cellular-automaton
audits. The existing Frobenius Lens manuscript studies finite-field rule
identifiability and Frobenius-time observations, rather than complex contour
zero counts. No prior manuscript found in this review used the present
exponential-panel enclosure or gauge objective. Other tasks may concurrently
add files; the review is a snapshot, not a claim to monitor the whole folder.

## Candidate topics rejected or narrowed

1. AAA approximation plus a Rouché zero-count claim: too close to existing
   rational approximation and certified Nyquist literature to support a broad
   new-method claim without further work.
2. Raw sampled winding plus derivative bounds: certified arc sampling is already
   established, including Apkarian and Noll's sampled Nyquist certificate.
3. Exponential rescaling as a new principle: rejected; analytic normalization and
   rescaling are standard, notably in Evans-function computations.
4. Retained question: can **local** gauges chosen by a particular convex
   enclosure objective compress a replayable disk certificate, while preserving
   the phase by an explicit local correction? This is the manuscript's narrow
   contribution. The underlying complex-analysis ingredients are classical.

## Literature review

Queries included combinations of “argument principle”, “exponential polynomial”,
“normalization”, “scaling”, “preconditioning”, “validated zero counting”,
“winding number”, “local exponential gauge”, “certified Nyquist”, “pole pruning”,
and “Rouché”. Search results were followed to primary publisher pages, author
manuscripts, institutional software pages, and arXiv originals. The review
cannot prove absence of all prior art. It establishes a bounded, attributable
claim rather than a claim that the topic has never been considered.

| Primary source | Relevance and boundary |
| --- | --- |
| Delves & Lyness (1967), DOI 10.1090/S0025-5718-1967-0228165-4 | Classical analytic contour root finding; root counting is not new. |
| Ying & Katz (1988), DOI 10.1007/BF01395882 | Reliable argument-principle counting; certification is not new. |
| Kravanja & Van Barel (2000), DOI 10.1007/BFb0103927 | Established analytic root-finding framework. |
| Petras (2002), DOI 10.1016/S0377-0427(01)00586-6 | Self-validating analytic quadrature. |
| Johansson (2017), DOI 10.1109/TC.2017.2690633 | Ball arithmetic used in the implementation. |
| Johansson (2018), arXiv:1802.07942 | Independent adaptive integration reference. |
| Apkarian & Noll (2018), DOI 10.1002/rnc.4073 | Sampled Nyquist tests with certified between-sample variation; author manuscript section 3 inspected. |
| Barker et al. (2018), DOI 10.1098/rsta.2017.0184 | Evans normalization, analyticity, winding and Rouché bounds; normalization itself is established. |
| Lessard & Mireles James (2020), DOI 10.3934/jcd.2020005 | Validated delay spectra and full unstable-root counts already exist. |
| Chonev, Ouaknine & Worrell (2023), DOI 10.1145/3603543 | Deeper real-zero decision questions; our finite contour claim does not resolve them. |
| Hayes (1950), DOI 10.1112/jlms/s1-25.3.226 | Classical scalar delay characteristic equation. |
| Trefethen & Weideman (2014), DOI 10.1137/130932132 | Analytic quadrature and root-count moments. |
| DDE-BIFTOOL manual, arXiv:1406.7144 | Existing delay stability and bifurcation software. |

## Practical and commercial landscape

MathWorks already documents direct time-delay analysis, frequency margins and
Padé approximations: https://www.mathworks.com/help/control/ug/analyzing-control-systems-with-delays.html
DDE-BIFTOOL already computes stability-determining characteristic roots:
https://twr.cs.kuleuven.be/research/software/delay/ddebiftool.shtml

The potential integration is an auditable, exact-model zero-count record in a
simulation or feedback-design pipeline. Certificate construction may be paid
once, followed by repeated replay during model review or regression checking.
The experiments examine this cost split. They do not estimate market size,
customer willingness to pay, regulatory acceptance, or patentability. The
published mathematical mechanisms are heavily based on classical prior art;
no patentability claim is made.

## Review checklist actually applied

- Every claimed count is computed, with exact rational inputs saved.
- All main certificates replay at higher precision.
- Lattice references retain multiplicities and use ball logarithms.
- Delay references use separate adaptive quadrature of f'/f.
- Single-delay references additionally match analytic crossing thresholds.
- Tail inequalities establish that delay contours cover all unstable roots.
- Corrupted partitions, reversed geometry and invalid centers are rejected.
- A boundary-zero test returns unresolved.
- Construction overhead is reported; panel reductions are not mislabeled as
  universal speedups.
- Fixed-9 is identified as a finite offline oracle, not an optimal global scale.
- No industrial benchmark, customer study, or peer review is claimed.

## Follow-up prompted by the initial study

Small delay-family gains against a sensible fixed scale motivated a second
enclosure. It sums the first two complete Taylor coefficients before taking
absolute values and validates an exponential remainder. A logarithmic-derivative
proposal can cancel the linear coefficient in exact arithmetic; its rounded
rational approximation is verified without assuming exact cancellation.
The follow-up uses all 36 original delay instances, with the improved enclosure
also supplied to every baseline. Jet shortens the records on 26 instances and
ties on ten against the best tested fixed center. Median paired panel ratios
are about 1.06 and 1.07, while median replay ratios are below one. These are
modest proof-length improvements, with no aggregate replay speedup.
The follow-up is explicitly identified as exploratory rather than preregistered.
All 625 archived certificate records replay at 192 bits. Taylor models and
logarithmic derivatives themselves are not claimed as new mathematics.
