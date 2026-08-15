# Reference data

The published values the validation suite compares against, one CSV per table, plus the
reasoning behind every acceptance band. The bands are the part worth reviewing: a residual is
only meaningful against a stated expectation, and an expectation chosen after seeing the
residual is not a test.

Each CSV shares a core schema:

| column | meaning |
| --- | --- |
| `key` | nuclide name, or mass number for a chain yield |
| `value` | the reference number **as its source prints it** |
| `unit` | the unit that number is in |
| `tolerance_rel` | the band, as a fraction |
| `gate` | `gate` asserts in CI; `report` computes and prints without asserting |
| `source` | a key resolved below |
| `note` | why this row deviates, or why it is not gated |

Comparisons are performed in the **source's** units, converting NuSIFT's number rather than the
published one, so nothing in these files is a transcription of a transcription.

`report` rows are not failures being tolerated. They are places where the published quantity and
the computed one are not the same quantity, for a reason named in the row. Widening a band until
such a row passes would state something false; printing it with its cause states something true.

## Sources

**`ninkovic2012`** — Ninković, M. M., & Adrović, F. (2012). *Air Kerma Rate Constants for
Nuclides Important to Gamma Ray Dosimetry and Practical Application.* In F. Adrović (Ed.),
Gamma Radiation (pp. 3–20). InTech. <https://doi.org/10.5772/35029>

Table 1, in µGy·m²/(GBq·h). This recalculation exists because, in its authors' words, published
data are in strong disagreement — which is why it is preferred here over the more commonly
quoted tables. Its Co-60 entry of 309.0 is 13.05 R·cm²/(h·mCi) in the modern roentgen and 13.15
in the pre-1979 one, the latter being the classic 13.2 that older tables give. Same physics,
different decade; see [exposure.md §4](../../docs/exposure.md).

**The cutoff matters more than the choice of table.** The paper counts only photons above
20 keV and excludes bremsstrahlung, stating so in the Table 1 caption and twice in its text.
NuSIFT by default sums the whole spectrum. Compared without that cutoff the two disagree by up
to a factor of ten on the X-ray emitters, none of it physics. The suite applies the same 20 keV
cutoff, which is what the `min_energy_ev` argument on `NuclearData.gamma_constant` is for.

Do not read Table 1's "energy interval from/to" columns as integration bounds — they are the
nuclide's full emitted spectrum, listed for information, and they extend below 20 keV.

**`smith2012`** — Smith, D. S., & Stabin, M. G. (2012). *Exposure rate constants and lead
shielding values for over 1,100 radionuclides.* Health Physics 102(3), 271–291.
<https://doi.org/10.1097/HP.0b013e318235153a>

The second gamma-constant table, in R·cm²/(h·mCi), based on ICRP-107 decay data, with a stated
**15 keV** cutoff — different from Ninković's 20 keV, which is why the cutoff is carried per row
rather than set once. It supplies the nuclides Ninković does not cover.

Two traps in its Table 1, both verified by eye rather than from a text layer, whose exponents
are corrupted by line-wrapping: Rb-86m is printed *above* Rb-86 and the two are easily
transposed, and the Cs-137 entry is by its own footnote the Ba-137m value, not a bare-parent
constant.

That the two tables disagree with each other by about a percent where they overlap — Co-60 is
13.05 R·cm²/(h·mCi) from Ninković against 12.9 here — is the point [exposure.md
§4](../../docs/exposure.md) makes about published constants. NuSIFT lands between them.

**`ame2020`** — Wang, M., Huang, W. J., Kondev, F. G., Audi, G., & Naimi, S. (2021). *The
AME2020 atomic mass evaluation (II). Tables, graphs and references.* Chinese Physics C 45(3),
030003. <https://doi.org/10.1088/1674-1137/abddaf> — data file `mass.mas20`, Atomic Mass Data
Center, <https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt>

Read out of the file rather than a web table. Since the 2019 redefinition of the mole, the molar
mass in g/mol equals the atomic mass in u to within 3×10⁻¹⁰, so the two are compared directly.
One parsing trap worth recording: the integer part of the tabulated mass is `floor(mass in u)`,
not the mass number — Co-60 reads 59.933… because the mass defect pulls it below 60.

**`ensdf2022`** — Evaluated Nuclear Structure Data File, April 2022 snapshot, via the IAEA
Nuclear Data Section Livechart API.
<https://nds.iaea.org/relnsd/vcharthtml/VChartHTML.html>

One evaluation throughout rather than the closest value per nuclide. Mixing compilations would
let a row be quietly rehomed to whichever source happened to agree with the store, which is the
opposite of a test. DDEP (LNHB, <http://www.lnhb.fr/nuclear-data/nuclear-data-table/>) was
cross-checked and agrees on all but three: K-40 by 0.3%, Sr-90 by 0.4%, Ba-133 by 0.1%. Those
disagreements are noted on their rows and are the reason the band is 1% rather than tighter.

**`endf80_cfy`** — Brown, D. A., et al. (2018). *ENDF/B-VIII.0: The 8th Major Release of the
Nuclear Reaction Data Library.* Nuclear Data Sheets 148, 1–142.
<https://doi.org/10.1016/j.nds.2018.02.001>

U-235 neutron-induced fission yields, MAT 9228, MF=8/MT=459, thermal point E = 0.0253 eV, read
out of the tape itself rather than off a web table. Worth knowing: ENDF/B-VIII.0's fission
yields are still the England & Rider (1989) evaluation, unchanged since ENDF/B-VI.

**`peplow2020`** — Peplow, D. E. (2020). *Comparison of Dose Rate Constants.* Health Physics.
<https://doi.org/10.1097/HP.0000000000001136>

Used for the ICRP 116 effective-dose tabulation, and for the definition of a published constant
as a vacuum quantity — a point source in a vacuum, no self-attenuation, no air scatter. That
last point is why the residual on the gamma constants is *not* scatter that an uncollided
calculation omits.

**`glasstone1977`** — Glasstone, S., & Dolan, P. J. (1977). *The Effects of Nuclear Weapons*
(3rd ed.). U.S. Department of Defense and Energy Research and Development Administration.

The t^-1.2 gross fission-product decay rule, and the 1.45×10²³ fissions per kiloton at
180 MeV/fission that `nusift seed-fission` uses.

**`radioactivedecay`** — Malins, A., & Lemoine, T. (2022). *radioactivedecay: A Python package
for radioactive decay calculations.* Journal of Open Source Software 7(71), 3318.
<https://doi.org/10.21105/joss.03318>

ICRP-107 decay data, a matrix-exponential solver, and no shared lineage with NuSIFT or cram.
Version-pinned in `validation/checks.py`; the protocol lives there rather than in a CSV because
it is a procedure, not a published number.

## Why each band is what it is

**Gamma constants — ±8%, and 9 of 35 rows not gated.**
Once the same 20 keV cutoff is applied, most rows land inside 2%. The band is wider than that
because of where the remaining spread comes from: the air table is at its sparsest and most
curved around the Compton minimum near 100 keV, and log-log interpolation of a function that is
turning over there is the least accurate thing the exposure model does. Co-57, whose constant is
carried entirely by 122 and 136 keV photons, is the clearest case at −5.4%.

A row is **not gated** when more than 30% of its above-20-keV constant comes from photons below
100 keV. That criterion is a property of the spectrum, decided before the residual is consulted,
and it is recomputed on every run — the report prints the fraction for every row. Below 10 keV
NuSIFT clamps the air coefficients rather than extrapolating, so those rows are order-of-magnitude
figures by construction and the store census says so too.

**Half-lives — ±1%.**
The staged values clear this by two to three orders of magnitude; ENDF/B-VIII.1 adopts the same
evaluations, so near-exact agreement is expected rather than impressive. The band is set by the
spread *between* compilations, not by the staging error being looked for. What it actually
catches is a unit slip, a misread tape, or a value attached to the wrong isomer.

**Chain yields — ±4%.**
NuSIFT's number is the sum of independent yields over a mass chain at t = 0, which equals the
cumulative yield at the chain's terminus except for delayed-neutron emission — the one process
that moves a nucleus off its chain. That is a real physical offset, not numerical slop: A = 137
sits +2.5% high because I-137 emits a delayed neutron in about 7% of its decays, and A = 85 sits
−2.5% low because the light wing gains from the chain above it. The band accommodates that term.
This is a consistency check on seeding and mass-chain aggregation rather than an independent
measurement, and it is labelled as one.

**ICRP 116 ratios — ±6% on the agreeing rows, ±15% on Am-241.**
This table checks that a caveat is still true. Co-60 and Ba-137m agree with tabulated effective
dose to better than a percent, which is a coincidence of energy near 1 MeV rather than evidence
that NuSIFT's sievert column is effective dose. Am-241 is gated on **diverging** by about a
factor of five. If that row ever started agreeing, the warning in
[exposure.md §6](../../docs/exposure.md#6-units) would have quietly become wrong.

**Way-Wigner — slope within [−1.35, −1.05].**
An empirical fit to gross behaviour, quoted as good to roughly 25% over its validity window, not
an exact exponent. The local slope genuinely moves across the window — measured between −1.07
and −1.23 decade by decade — and the band admits that spread while still rejecting anything
qualitatively wrong, such as a single exponential or a chain that never turns over.

**Cross-code — ±2% single-parent, ±3% mixed.**
Observed worst is 0.82%, so this is roughly a threefold margin. It is not tighter because
ICRP-107 and ENDF/B-VIII.1 genuinely differ: 0.29% on the Cs-137 half-life, and more on some
branchings. It is not looser because an actual error in the chain construction or the solve
would show as tens of percent, not tenths.

One nuclide is held out, with its cause measured rather than assumed: **Xe-131m** disagrees by a
flat −7.8% at every time in every case, which is the ratio of the two evaluations' I-131 →
Xe-131m branchings (0.0108477 against 0.011759 = 0.9225) and not a property of either solver.
Gating it would gate one evaluation against the other. The held-out list is in
`validation/checks.py`; adding to it without first measuring the cause would turn it into a
place for inconvenient failures to go, which is exactly what it must not become.

## Sources considered and rejected

**Unger, L. M., & Trubey, D. K. (1982).** *Specific Gamma-Ray Dose Constants for Nuclides
Important to Dosimetry and Radiological Assessment* (ORNL/RSIC-45/Rev.1).
<https://www.osti.gov/biblio/6246345> — the scanned copy's OCR is unusable for numeric
extraction, rendering 511.0 as "SILO" and destroying exponents. Its quantity is also a
dose-equivalent rate constant rather than air kerma, and its decay data is 1982-vintage. It is
still cited in [exposure.md §5](../../docs/exposure.md) for its argument about folded-in Cs-137
constants, which does not depend on any number being read off the page.

**Smith, D. S., & Stabin, M. G. (2012).** *Exposure rate constants and lead shielding values for
over 1,100 radionuclides.* Health Physics 102(3), 271–291 — the RADAR copy is gone (HTTP 404)
and the host serves a mismatched certificate. Not used rather than cited from memory.

**Atomic weights — ±0.01%.**
By far the tightest band here, and still cleared by three orders of magnitude, because ENDF's
atomic weight ratio and AME2020 are evaluations of the same measurements. It is worth having
anyway: everything expressed per gram passes through this number, and nothing else in the suite
touches it.

## There is no specific-activity table

Specific activity is ln(2)·N_A/(T½·M) — a derivation over two quantities that now each have
their own table here, so a sweep against it would mostly re-test the half-lives twice.

It would also be the least conclusive comparison in the suite, because published
specific-activity tables disagree with each other by more than the staging error being looked
for. Am-241 is tabulated at 3.5 Ci/g by Argonne, 3.43 by several health-physics compilations,
and 3.2 by the DOE New Brunswick Laboratory; the spread is entirely which half-life each
compiler used, 432.2 y against the older ~458 y. Argonne's own fact sheets state they were
computed by scaling from Ra-226 and rounding to two significant figures rather than from the
decay constant at all. Choosing among those would be choosing an answer.

Six keystone nuclides are still checked end to end in
[`tests/validation/test_published_constants.cpp`](../../tests/validation/test_published_constants.cpp)
at 2%, because the quotient is what a user actually reads off a report. That band is set by the
same revision effect: Sr-90 is widely tabulated at 136 Ci/g from the older 29.1 y half-life,
against the 28.79 y the store carries.
