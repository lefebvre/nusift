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
rather than set once. It supplies the nuclides Ninković does not cover. The RADAR web copy of
this table is no longer served (HTTP 404, behind a mismatched certificate), but the Internet
Archive holds the authors' file — capture of 2025-04-19 of
`doseinfo-radar.com/Exposure_Rate_Constants_and_Lead_Shielding_Values 4.pdf` — and every value
here is read from that article's own Table 1, by eye, not from a secondary quotation of it. The
paper states its conventions in one sentence: photons of at least 15 keV and yields of at least
10⁻⁴, bremsstrahlung neglected.

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

**`icrp116`** — ICRP (2010). *Conversion Coefficients for Radiological Protection Quantities for
External Radiation Exposures.* ICRP Publication 116, Ann. ICRP 40(2–5).

Unlike every other source here this one is not a reference to compare against — it is an INPUT.
Table A.1, effective dose per unit fluence for monoenergetic photons in six irradiation
geometries, is transcribed into
[`nusift/exposure/dose_coefficients.cpp`](../../nusift/exposure/dose_coefficients.cpp) and is
what makes a sievert in NuSIFT a sievert. The publication itself is not in the repository; the
table is, printed digit for digit as ICRP prints it, on ICRP's own energy grid, so it can be
checked against a copy of the page line by line.

The transcription is checked without a copy of the page too, and by a route that shares nothing
with it. Publication 116 also gives Table A.2, effective dose per air kerma free-in-air, which it
does not ship: dividing NuSIFT's dose kernel by its air-kerma kernel — one built from this table,
the other from NIST's mass energy-absorption coefficients — has to reproduce A.2, and does, to
better than 1% from 20 keV to 2 MeV. Above that the two drift apart by air's radiative yield,
because ICRP's denominator is kerma and NIST's µ_en/ρ counts only what is absorbed; that
departure is asserted rather than tolerated, in
[`tests/unit/test_point_source.cpp`](../../tests/unit/test_point_source.cpp).

**`icrp119`** — ICRP (2012). *Compendium of Dose Coefficients based on ICRP Publication 60.*
ICRP Publication 119, Ann. ICRP 41(Suppl.).

An INPUT rather than a reference, like `icrp116`: Table A.1's committed effective dose
coefficients for ingested and inhaled particulates by workers are transcribed into
[`data/packs/icrp119-ingestion-worker.csv`](../../data/packs/icrp119-ingestion-worker.csv) and
[`data/packs/icrp119-inhalation-worker-5um.csv`](../../data/packs/icrp119-inhalation-worker-5um.csv).
733 nuclides, every one of which resolves against the ENDF/B-VIII.1 store.

**The extraction could not use the PDF's text layer, and this is worth recording because the
failure is silent.** Every exponent's minus sign is a glyph in a symbol font with no ToUnicode
mapping, so `pdftotext` renders 4.8E-11 as `4.8E11` -- a dose coefficient twenty-two orders of
magnitude too large, in a column where 1e-11 and 1e11 are equally well-formed numbers. The
minus survives in poppler's XML rendering as an empty text fragment in font 9 sitting between
the mantissa and the exponent digits, and the parser reconstructs it from there; 3511 of the
3513 such fragments in the table sit in exactly that position, the other two being a footnote
marker in a different font. Coefficients are then checked against the shape of a dose
coefficient before being accepted, which is what caught the repeated page header leaking into
the column beneath it.

The values are checked against figures published widely enough to be recognised on sight:
Cs-137 ingestion 1.3e-8 Sv/Bq, I-131 ingestion 2.2e-8, Sr-90 ingestion 2.8e-8 at f1 = 0.3,
H-3 (HTO) 1.8e-11, Pu-239 inhalation Type M at 5 um 3.2e-5, Am-241 the same at 2.7e-5.

**`fgr15`** — EPA (2025). *External Exposure to Radionuclides in Air, Water, and Soil.*
Federal Guidance Report No. 15, EPA 402-R-25-001.

An INPUT, like `icrp116` and `icrp119`. Table 4-6's reference person effective dose rate
coefficients for air submersion, adult column, are transcribed into
[`data/packs/fgr15-air-submersion-adult.csv`](../../data/packs/fgr15-air-submersion-adult.csv):
1246 nuclides, every one of which resolves against the ENDF/B-VIII.1 store.

This table extracts cleanly from the PDF's text layer -- the exponents keep their minus signs,
unlike ICRP 119's -- so the checks are on the shape of the result rather than on the glyphs:
every row carries six age columns that parse as numbers, no nuclide appears twice, and the adult
coefficient is the smallest of the six in all 1246 rows, which is both the expected physics
(children take more dose from the same cloud) and a check that the columns were not transposed.
Twenty-four nuclides have a coefficient of exactly zero, which is a real value -- they have no
penetrating emission -- and is carried with that reason in the row rather than as a gap.

FGR-15 states that its coefficients exclude decay products, so the pack declares
`progeny: excluded`, and the publication's own instruction to combine a parent with its progeny
"only after consideration of the equations describing production and decay of daughter
radionuclides over time" is what the decay solve underneath the ranking performs.

**`icrp74`** — ICRP (1996). *Conversion Coefficients for use in Radiological Protection against
External Radiation.* ICRP Publication 74, Ann. ICRP 26(3/4).

An INPUT. Table A.21's conversion coefficients from photon fluence to the ambient dose
equivalent H*(10) are transcribed into
[`data/packs/icrp74-ambient-dose-h10.csv`](../../data/packs/icrp74-ambient-dose-h10.csv), 25
energies from 10 keV to 10 MeV.

**This one was transcribed by eye, and it is the only table here that was.** The copy is a page
scan with no text layer in the body and unusable OCR on the front matter, so there was nothing
to parse; the table was read from the page image. What makes that defensible is that ICRP 74
gives the same quantity twice -- per unit fluence and per unit air kerma -- alongside the air
kerma per unit fluence relating them, so the three columns have to multiply out. They do, to
better than 0.7% at every energy except 1.5 MeV, where the printed H*(10)/Phi of 6.90 pSv cm2
sits 2.0% below the 7.04 its neighbours imply. That row is carried as printed: a transcription
reproduces the page, and the discrepancy is recorded in the pack's own header rather than
quietly corrected. The publication's footnote says the air-kerma column comes from a different
compilation than the protection quantities, which is consistent with some inconsistency between
them.

The independent check is against the other kernel already in the tree: H*(10) is designed to be
a conservative estimate of effective dose for photons below about 10 MeV, and on a caesium field
this pack gives 10.99 Sv/h against ICRP 116's 9.25 Sv/h -- a ratio of 1.19, which is the
published relationship at 662 keV between two curves transcribed from different publications by
different routes.

**`glasstone1977`** — Glasstone, S., & Dolan, P. J. (1977). *The Effects of Nuclear Weapons*
(3rd ed.). U.S. Department of Defense and Energy Research and Development Administration.

The t^-1.2 gross fission-product decay rule, and the 1.45×10²³ fissions per kiloton at
180 MeV/fission that `nusift seed-fission` uses.

**`endf458`** — ENDF-6 section MT458, the fission energy-release partition, as evaluated for
U-235 thermal fission and reproduced in every major library since ENDF/B-VI: delayed betas
6.50 MeV, delayed gammas 6.33 MeV, neutrinos 8.75 MeV, against a total of about 202 MeV.

The delayed beta and gamma terms sum to the recoverable energy released by fission-product
DECAY, which is what `nusift integrate --metric heat` over all time has to reproduce. It is not
staged: NuSIFT fetches only the decay and fission-yield sublibraries, and MT458 lives in the
neutron sublibrary, so this is an outside number rather than one the store could be checked
against itself. That is the point of it.

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
turning over there is the least accurate thing the exposure model does — but its sign is known.
Through the 100–150 keV interval the tabulated μ_en/ρ is convex in log-log, so the chord lies
above the curve and the interpolation can only bias a 122 or 136 keV value *high*, by order a
percent at most. Co-57, carried entirely by those two lines, sits −5.4% against Ninković, and
Tc-99m's 141 keV line lands on the other side of the same source; neither is the interpolation.
The Smith & Stabin row for Co-57, kept as a second cross-source check, settles it: 0.5638
R·cm²/(h·mCi) computed against 0.563 published, +0.1%. The two tables disagree with each other
on Co-57 by five and a half percent, and NuSIFT lands on one of them.

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

**Atomic weights — ±0.01%.**
By far the tightest band here, and still cleared by three orders of magnitude, because ENDF's
atomic weight ratio and AME2020 are evaluations of the same measurements. It is worth having
anyway: everything expressed per gram passes through this number, and nothing else in the suite
touches it.

**Chain yields — ±4%.**
NuSIFT's number is the sum of independent yields over a mass chain at t = 0, which equals the
cumulative yield at the chain's terminus except for delayed-neutron emission — the one process
that moves a nucleus off its chain. That is a real physical offset, not numerical slop: A = 137
sits +2.5% high because I-137 emits a delayed neutron in about 7% of its decays, and A = 85 sits
−2.5% low because the light wing gains from the chain above it. The band accommodates that term.
This is a consistency check on seeding and mass-chain aggregation rather than an independent
measurement, and it is labelled as one.

**ICRP 116 ratios — ±3%.**
This table used to check that a caveat was still true, and now checks an agreement. The sievert
column is computed through ICRP 116's fluence-to-effective-dose kernel in the AP geometry these
coefficients are tabulated in, so it is the same quantity as the reference rather than air kerma
wearing its label, and the expected ratios are 1.00.

The band is the decay data's rather than the kernel's. Co-60 and Ba-137m land within a
thousandth; Am-241 sits 2% high for a measured reason — one of its lines falls below the table's
10 keV floor and takes the floor's coefficient, carrying 1.5% of the total, and excluding it the
ratio is 1.002. The row states that rather than the band absorbing it silently.

The suite also asserts that the air kerma this column no longer reports still diverges from
effective dose by a factor of five for Am-241. Without that, agreement would stop being evidence
the kernel does anything.

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
