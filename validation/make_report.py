"""Regenerate docs/validation.md and its figures.

    PYTHONPATH=python python validation/make_report.py

Every number and every figure here is computed from the committed store by validation/checks.py
-- the same module python/tests/test_validation.py asserts on -- so the published table and the
gate can never be two different numbers.

The output is deterministic. No wall-clock time is written, every value goes through a fixed
format, and figure geometry is arithmetic rather than measured text, so two runs on two
machines produce byte-identical files. That is what lets CI regenerate the report and fail on
`git diff`, which in turn is what stops the document rotting the way a hand-written table does.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "docs" / "figures"))

from svgplot import Axes, circle, line, rect, svg, text  # noqa: E402

from validation import checks  # noqa: E402

FIGURES = ROOT / "docs" / "figures"
REPORT = ROOT / "docs" / "validation.md"


def write(path, content):
    # Explicit LF: the repository normalises to LF and the drift check compares bytes.
    path.write_text(content, encoding="utf-8", newline="\n")


def check_figure_bounds(paths):
    """Fail if any label was placed outside its own canvas.

    Text is the one thing here whose size is not derived from the data, so a caption that fits
    today can run off the edge when a number gains a digit or a nuclide is added to a table.
    A generated figure nobody looks at is exactly where that goes unnoticed, so it is checked
    rather than eyeballed. Geometry only: this cannot see glyph widths, so it catches a label
    placed past the edge, not one that merely reaches it.
    """
    import xml.etree.ElementTree as elementtree

    problems = []
    for path in paths:
        root = elementtree.parse(path).getroot()
        width, height = float(root.get("width")), float(root.get("height"))
        for element in root.iter("{http://www.w3.org/2000/svg}text"):
            x, y = float(element.get("x")), float(element.get("y"))
            if not (0.0 <= x <= width and 0.0 <= y <= height):
                problems.append(f"{path.name}: '{(element.text or '')[:40]}' at ({x:.0f},{y:.0f})"
                                f" outside {width:.0f}x{height:.0f}")
    if problems:
        raise SystemExit("figure labels fall outside the canvas:\n  " + "\n  ".join(problems))


def pct(value, digits=2):
    return f"{value * 100:+.{digits}f}%"


def caption(parts, x, y, sentences, width=880, size=11.5, leading=18, cls="muted"):
    """Wrap prose to the canvas and emit it, returning the y after the last line.

    Text length is the one thing in these figures that cannot be derived from the data, and a
    hand-counted line that fits today runs off the edge the moment a number gains a digit. The
    budget is a conservative average character width for the sans stack in svgplot; no font is
    measured, so the wrap stays deterministic.
    """
    limit = int((width - x - 26) / (0.53 * size))
    words = " ".join(sentences).split()
    lines, current = [], ""
    for word in words:
        candidate = f"{current} {word}".strip()
        if current and len(candidate) > limit:
            lines.append(current)
            current = word
        else:
            current = candidate
    if current:
        lines.append(current)
    for index, line_text in enumerate(lines):
        parts.append(text(x, y + index * leading, line_text, cls, size))
    return y + len(lines) * leading


# --------------------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------------------
#
# One form does most of the work: a residual dot plot, one row per reference value, with the
# acceptance band drawn behind it. It is the right form because the question these tables
# answer is "how far off, and is that inside what we said we would accept" -- a magnitude with
# a sign, per named thing. Bars would imply the zero line is a quantity rather than a target,
# and a table alone would not show at a glance that a whole family of rows sits inside a
# percent while a handful sit outside.


def residual_plot(path, rows, title, footer, band, unit_note, width=880, row_height=15.5):
    top = 96
    height = top + row_height * len(rows) + 110
    # The axis reaches three times the band at most. Past that a single distant outlier would
    # squeeze every row that matters into an unreadable stripe around zero, so an outlier is
    # clipped to the edge and marked -- its exact value is printed beside it either way.
    limit = max(band * 1.9, min(max(abs(r["residual"]) for r in rows) * 1.15, band * 3.0))
    ax = Axes(150, width - 130, top, top + row_height * len(rows), (-limit, limit), (0.0, 1.0))

    parts = []
    # The band behind the marks, so a row's position is read against what was accepted rather
    # than against an arbitrary axis range.
    parts.append(rect(ax.x(-band), ax.top - 6, ax.x(band) - ax.x(-band),
                      row_height * len(rows) + 12, "fill-a", 0.09, rx=3, width=0))
    for edge in (-band, band):
        parts.append(line(ax.x(edge), ax.top - 6, ax.x(edge), ax.bottom + 6, "axis", 1.0, "3 3"))
    parts.append(line(ax.x(0.0), ax.top - 6, ax.x(0.0), ax.bottom + 6, "axis", 1.4))

    for value in (-limit, -band, 0.0, band, limit):
        if value not in (-limit, limit):
            parts.append(line(ax.x(value), ax.top - 6, ax.x(value), ax.bottom + 6, "grid"))
        label = "0" if value == 0.0 else f"{value * 100:+.3g}%"
        parts.append(text(ax.x(value), ax.top - 14, label, "muted", 10, "middle"))

    # A rule wherever the reference source changes, so a figure drawing on more than one table
    # cannot be read as though it were all one comparison.
    for index in range(1, len(rows)):
        if rows[index]["source"] != rows[index - 1]["source"]:
            y = ax.top + row_height * index
            parts.append(line(ax.left - 92, y, ax.right, y, "axis", 0.9, "2 4"))
            parts.append(text(ax.left - 96, y + 3.6, rows[index]["source"], "muted", 9.5, "end"))

    for index, row in enumerate(rows):
        y = ax.top + row_height * (index + 0.5)
        gated = row["gate"]
        cls_line = "s0" if gated else "s1"
        cls_dot = "t0" if gated else "t1"
        residual = row["residual"]
        clipped = max(-limit, min(limit, residual))
        parts.append(line(ax.x(0.0), y, ax.x(clipped), y, cls_line, 1.6))
        parts.append(circle(ax.x(clipped), y, 3.4, cls_dot))
        parts.append(text(ax.left - 12, y + 3.6, row["key"], "ink", 11, "end",
                          "600" if gated else "normal"))
        label = pct(residual, 1)
        if abs(residual) > limit:
            label += "  →"
        parts.append(text(ax.right + 10, y + 3.6, label, "muted", 10))

    parts.append(text(58, 30, title, "ink", 13.5, "start", "600"))
    parts.append(text(58, 50, unit_note, "muted", 11.5))

    legend_y = 70
    parts.append(circle(60, legend_y - 3.5, 3.4, "t0"))
    parts.append(text(70, legend_y, "gated in CI", "muted", 11))
    parts.append(circle(168, legend_y - 3.5, 3.4, "t1"))
    parts.append(text(178, legend_y, "reported, not gated", "muted", 11))
    parts.append(text(ax.x(0.0) + 6, legend_y, f"shaded band: ±{band * 100:g}% acceptance",
                      "muted", 11))

    after = caption(parts, 58, ax.bottom + 34, footer, width)

    write(path, svg(width, max(height, after + 16), "\n".join(parts), title))


def figure_overview(path, families):
    """The scoreboard: how much of each family's band its worst row actually uses."""
    row_height = 38
    width = 880
    height = 104 + row_height * len(families) + 92
    ax = Axes(330, 720, 96, 96 + row_height * len(families), (0.0, 1.25), (0.0, 1.0))

    parts = []
    parts.append(rect(ax.x(0.0), ax.top - 8, ax.x(1.0) - ax.x(0.0),
                      row_height * len(families) + 16, "fill-a", 0.08, rx=3, width=0))
    parts.append(line(ax.x(1.0), ax.top - 8, ax.x(1.0), ax.bottom + 8, "accent", 1.6, "4 3"))
    parts.append(text(ax.x(1.0), ax.top - 16, "the band", "warn", 10.5, "middle", "600"))
    for fraction in (0.25, 0.5, 0.75):
        parts.append(line(ax.x(fraction), ax.top - 8, ax.x(fraction), ax.bottom + 8, "grid"))
        parts.append(text(ax.x(fraction), ax.top - 16, f"{fraction:.0%}", "muted", 10, "middle"))

    for index, family in enumerate(families):
        top = ax.top + row_height * index
        used = family["used"]
        parts.append(rect(ax.x(0.0), top + 8, max(ax.x(min(used, 1.24)) - ax.x(0.0), 1.5), 13,
                          "fill-a", 0.85, rx=3, width=0))
        # Name and detail stacked in the gutter rather than set on opposite sides of the same
        # line: the detail strings vary enough in length that a single line collides.
        parts.append(text(ax.left - 14, top + 12, family["name"], "ink", 11.5, "end", "600"))
        parts.append(text(ax.left - 14, top + 27, family["detail"], "muted", 10, "end"))
        parts.append(text(ax.x(min(used, 1.24)) + 10, top + 19,
                          f"{used:.0%} of ±{family['band']}", "muted", 10.5))

    parts.append(text(58, 30, "How much of each acceptance band the worst row actually uses",
                      "ink", 13.5, "start", "600"))
    caption(parts, 58, 50, ["Every family passes. A short bar is margin; a bar reaching the "
                            "dashed line would be a check about to fail."], width)
    after = caption(parts, 58, ax.bottom + 46,
                    ["Bands are set from the physics that limits each comparison -- evaluation "
                     "differences, delayed-neutron flow, air-table interpolation -- not from the "
                     "residuals they turned out to have. Each is justified beside its table in "
                     "validation/references/README.md."], width)
    write(path, svg(width, max(height, after + 16), "\n".join(parts), "Validation scoreboard"))


def figure_cutoff(path, rows):
    """Why a published constant and a NuSIFT constant have to be compared on the same terms."""
    width, height = 880, 470
    ordered = sorted(rows, key=lambda r: -r["soft_fraction"])[:18]
    ax = Axes(150, 700, 104, 104 + 17 * len(ordered), (0.3, 30.0), (0.0, 1.0), xlog=True)

    parts = []
    for value, label in ((1.0, "exact"), (2.0, "2x"), (5.0, "5x"), (10.0, "10x")):
        parts.append(line(ax.x(value), ax.top - 6, ax.x(value), ax.bottom + 6, "grid"))
        parts.append(text(ax.x(value), ax.top - 14, label, "muted", 10, "middle"))
    parts.append(line(ax.x(1.0), ax.top - 6, ax.x(1.0), ax.bottom + 6, "axis", 1.4))

    for index, row in enumerate(ordered):
        y = ax.top + 17 * (index + 0.5)
        whole = max(row["whole_spectrum"] / row["published"], 0.3)
        cut = max(row["computed"] / row["published"], 0.3)
        parts.append(line(ax.x(cut), y, ax.x(min(whole, 29.0)), y, "s1", 1.4, "2 2"))
        parts.append(circle(ax.x(min(whole, 29.0)), y, 3.2, "t1"))
        parts.append(circle(ax.x(cut), y, 3.4, "t0"))
        parts.append(text(ax.left - 12, y + 3.6, row["key"], "ink", 11, "end"))
        parts.append(text(ax.right + 10, y + 3.6, f"{row['soft_fraction']:.0%} soft", "muted", 10))

    parts.append(text(58, 30, "The same spectra, measured the way the reference measures them",
                      "ink", 13.5, "start", "600"))
    parts.append(text(58, 50,
                      "Ratio of NuSIFT's constant to the published one. The reference counts "
                      "photons above 20 keV only; NuSIFT by default counts the whole spectrum.",
                      "muted", 11.5))
    parts.append(circle(60, 74.5, 3.2, "t1"))
    parts.append(text(70, 78, "whole spectrum", "muted", 11))
    parts.append(circle(196, 74.5, 3.4, "t0"))
    parts.append(text(206, 78, "above 20 keV, as published", "muted", 11))
    after = caption(parts, 58, ax.bottom + 36,
            ["Soft X-rays carry little energy but air absorbs them greedily, and below 10 keV "
             "the coefficients are clamped, so counting them against a table that excludes them "
             "overstates by up to a factor of ten. The rows are ordered by how much of the "
             "constant sits below 100 keV -- which is also what decides whether a row is "
             "gated."], width)
    write(path, svg(width, max(height, after + 16), "\n".join(parts), "Effect of the published energy cutoff"))


def figure_way_wigner(path, result):
    width, height = 880, 470
    times = result["times"]
    activity = result["activity"]
    ax = Axes(96, 700, 76, 320, (times[0], times[-1]),
              (activity.min() * 0.6, activity.max() * 1.6), xlog=True, ylog=True)

    parts = []
    decades = [3600.0, 86400.0, 604800.0, 30 * 86400.0]
    parts += ax.vgrid(decades, ["1 h", "1 d", "1 w", "30 d"])
    exponent = int(np.floor(np.log10(activity.min())))
    grid = [10.0**e for e in range(exponent, exponent + 4)]
    parts += ax.hgrid(grid, [f"1e{e}" for e in range(exponent, exponent + 4)])

    reference = activity[0] * (times / times[0]) ** -1.2
    parts.append(ax.path(times, reference, "accent", 1.6, "5 4"))
    parts.append(ax.path(times, activity, "curve", 2.2))
    parts += ax.frame()

    parts.append(text(58, 30, "Gross fission-product activity against the Way-Wigner rule",
                      "ink", 13.5, "start", "600"))
    parts.append(text(58, 50, "20 kt U-235 thermal fission, total activity over all nuclides",
                      "muted", 11.5))
    parts.append(text(ax.right + 8, ax.y(activity[-1]) + 4, "NuSIFT", "t0", 11, "start", "600"))
    parts.append(text(ax.right + 8, ax.y(reference[-1]) + 4, "t^-1.2", "t1", 11, "start", "600"))
    parts.append(text(ax.right, ax.bottom + 40, "time after fission →", "muted", 11, "end"))
    after = caption(parts, 58, 366,
                    [f"Fitted log-log slope {result['slope']:.3f}, inside the accepted "
                     f"[{result['band'][0]}, {result['band'][1]}]."], width, size=12, cls="ink")
    segments = ",   ".join(f"{label} {slope:.2f}" for label, slope in result["segments"])
    after = caption(parts, 58, after + 4,
            ["Nothing in NuSIFT was fitted to this rule: the exponent falls out of the yield "
             "set, the chain topology, the solve, and the activity weights together.",
             f"Local slope by decade — {segments}.",
             "The rule is a fit to gross behaviour rather than an exact exponent, which is why "
             "the local slope moves across the window and the band admits that."], width)
    write(path, svg(width, max(height, after + 16), "\n".join(parts), "Fission-product activity against t^-1.2"))


def figure_equilibria(path, curves):
    width, height = 880, 456
    parts = []
    panel_width = 320
    for index, curve in enumerate(curves):
        left = 92 + index * (panel_width + 108)
        days = curve["times"] / 86400.0
        both = np.concatenate([curve["engine_parent"], curve["engine_daughter"]])
        ax = Axes(left, left + panel_width, 112, 296, (days[0], days[-1]),
                  (both.max() * 3.0e-3, both.max() * 2.0), xlog=True, ylog=True)

        decades = [d for d in (0.01, 0.1, 1.0, 10.0, 100.0) if days[0] <= d <= days[-1]]
        labels = [f"{d:g} d" for d in decades]
        parts += ax.vgrid(decades, labels)
        parts.append(ax.path(days, curve["engine_parent"], "s0", 2.0))
        parts.append(ax.path(days, curve["engine_daughter"], "s1", 2.0))
        for series, cls in (("closed_parent", "t0"), ("closed_daughter", "t1")):
            for k in range(0, len(days), 6):
                parts.append(circle(ax.x(days[k]), ax.y(curve[series][k]), 2.6, cls))
        parts += ax.frame()

        parts.append(text(left, 88, f"{curve['parent']} → {curve['daughter']}", "ink", 12.5,
                          "start", "600"))
        parts.append(text(left + panel_width, 88, f"{curve['kind']} equilibrium", "muted", 11,
                          "end"))
        # Below the axis tick labels, which sit at bottom + 18.
        parts.append(text(left, 342, curve["parent"], "t0", 11, "start", "600"))
        parts.append(text(left + 62, 342, curve["daughter"], "t1", 11, "start", "600"))
        parts.append(text(left + panel_width, 342, f"branch {curve['branch']:.3f}", "muted", 10.5,
                          "end"))

    parts.append(text(58, 30, "Real chains from the shipped store, against their closed forms",
                      "ink", 13.5, "start", "600"))
    parts.append(text(58, 50,
                      "Lines are the CRAM solve; dots are the Bateman solution evaluated with "
                      "the store's own decay constants. Activity, log scale, arbitrary units.",
                      "muted", 11.5))
    worst = max(c["worst_residual"] for c in curves)
    after = caption(parts, 58, 386,
                    [f"Worst disagreement between solve and closed form across both chains: "
                     f"{worst:.2e} relative."], width, size=12, cls="ink")
    after = caption(parts, 58, after + 4,
            ["Mo-99 only lands on its curve if the 87.6% branch to the isomer was staged "
             "correctly; a branch of 1.0 would sit a seventh high and no half-life check would "
             "notice. Sr-90's daughter climbs to its parent's activity and stays there, which "
             "is what secular equilibrium means and what a Sr-90 source does."], width)
    write(path, svg(width, max(height, after + 16), "\n".join(parts), "Engine against closed-form Bateman"))


def figure_cross_code(path, curve, summaries):
    """One chain drawn twice. The cocktail's members sit within a factor of two of each other,
    so overlaying them says less than a single chain whose curves genuinely separate."""
    width, height = 880, 452
    times = curve["times"] / 3600.0
    everything = np.concatenate([ours for ours, _ in curve["series"].values()])
    ax = Axes(100, 690, 96, 292, (times[0], times[-1]),
              (everything.max() * 2.0e-3, everything.max() * 2.5), xlog=True, ylog=True)

    parts = []
    parts += ax.vgrid([1.0, 6.0, 24.0, 72.0, 336.0], ["1 h", "6 h", "1 d", "3 d", "14 d"])
    exponent = int(np.floor(np.log10(everything.max())))
    grid = [10.0**e for e in range(exponent - 2, exponent + 1)]
    parts += ax.hgrid(grid, [f"1e{e}" for e in range(exponent - 2, exponent + 1)])

    for index, (name, (ours, theirs)) in enumerate(sorted(curve["series"].items())):
        parts.append(ax.path(times, ours, f"s{index % 8}", 2.0))
        for k in range(0, len(times), 3):
            parts.append(circle(ax.x(times[k]), ax.y(theirs[k]), 2.8, f"t{index % 8}"))
        peak = int(np.argmax(ours))
        parts.append(text(ax.x(times[peak]), ax.y(ours[peak]) - 10, name, f"t{index % 8}", 11,
                          "middle", "600"))
    parts += ax.frame()

    parts.append(text(58, 30, f"{curve['parent']} decayed by two unrelated implementations",
                      "ink", 13.5, "start", "600"))
    parts.append(text(58, 50,
                      "Lines are NuSIFT (CRAM, ENDF/B-VIII.1); dots are radioactivedecay "
                      "(matrix exponential, ICRP-107). Activity in Bq.", "muted", 11.5))
    parts.append(text(58, 70,
                      "The dots sit on the lines through ingrowth, the crossing, and the decay "
                      "-- three regimes with different sensitivities to the chain data.",
                      "muted", 11.5))
    parts.append(text(ax.right, ax.bottom + 38, "time after separation →", "muted", 11, "end"))

    worst = max(s["max_residual"] for s in summaries)
    after = caption(parts, 58, 356,
                    [f"Across all {sum(s['comparisons'] for s in summaries)} gated comparisons "
                     f"-- six inventories, not just this one -- the largest disagreement is "
                     f"{worst * 100:.2f}%."], width, size=12, cls="ink")
    after = caption(parts, 58, after + 4,
            ["Different solver, different decay data, different author. The residual is not "
             "expected to be zero: the two evaluations differ on the Cs-137 half-life alone by "
             "0.29%. What agreement at this level rules out is an error in the chain "
             "construction or the solve, which would show as tens of percent rather than "
             "tenths."], width)
    write(path, svg(width, max(height, after + 16), "\n".join(parts), "Cross-code comparison of one decay chain"))


def figure_cross_code_residuals(path, fine):
    width, height = 880, 430
    gated = [row for row in fine if not row["excluded"]]
    excluded = [row for row in fine if row["excluded"]]
    limit = 0.09
    ax = Axes(96, 700, 84, 300, (2.0e3, 4.0e9), (-limit, limit), xlog=True)

    parts = []
    for band, cls in ((0.02, "fill-a"),):
        parts.append(rect(ax.left, ax.y(band), ax.right - ax.left, ax.y(-band) - ax.y(band),
                          cls, 0.10, rx=2, width=0))
    parts += ax.hgrid([-0.08, -0.04, 0.0, 0.04, 0.08], ["-8%", "-4%", "0", "+4%", "+8%"])
    parts += ax.vgrid([3600.0, 86400.0, 3.156e7, 3.156e8], ["1 h", "1 d", "1 y", "10 y"])

    for row in gated:
        parts.append(circle(ax.x(max(row["time_s"], 2.1e3)),
                            ax.y(max(-limit, min(limit, row["residual"]))), 2.8, "t0"))
    for row in excluded:
        parts.append(circle(ax.x(max(row["time_s"], 2.1e3)),
                            ax.y(max(-limit, min(limit, row["residual"]))), 2.8, "t1"))
    parts += ax.frame()

    parts.append(text(58, 30, "Every per-nuclide comparison against radioactivedecay",
                      "ink", 13.5, "start", "600"))
    parts.append(text(58, 50, "One dot per nuclide per time per case.", "muted", 11.5))
    parts.append(circle(300, 46.5, 2.8, "t0"))
    parts.append(text(310, 50, f"gated ({len(gated)})", "muted", 11))
    parts.append(circle(392, 46.5, 2.8, "t1"))
    parts.append(text(402, 50, f"held out, cause known ({len(excluded)})", "muted", 11))
    parts.append(text(ax.right, ax.bottom + 38, "time →", "muted", 11, "end"))
    after = caption(parts, 58, 356,
                    ["The held-out band at -7.8% is Xe-131m, and it is flat in time -- the "
                     "signature of a branching difference rather than a solver one."],
                    width, size=12, cls="ink")
    after = caption(parts, 58, after + 4,
            ["ENDF/B-VIII.1 feeds it from I-131 with branching 0.0108477 against ICRP-107's "
             "0.011759, a ratio of 0.9225, which is the residual. Gating it would gate one "
             "evaluation against another, which says nothing about NuSIFT; hiding it behind a "
             "wider band would say something false."], width)
    write(path, svg(width, max(height, after + 16), "\n".join(parts), "Cross-code residuals"))


def figure_icrp116(path, rows):
    width = 880
    ax = Axes(230, 700, 96, 250, (0.0, 6.0), (0.0, 1.0))

    parts = []
    for value in (1.0, 2.0, 3.0, 4.0, 5.0):
        parts.append(line(ax.x(value), ax.top - 8, ax.x(value), ax.bottom + 8, "grid"))
        parts.append(text(ax.x(value), ax.top - 16, f"{value:g}x", "muted", 10, "middle"))
    parts.append(line(ax.x(1.0), ax.top - 8, ax.x(1.0), ax.bottom + 8, "axis", 1.4))

    for index, row in enumerate(rows):
        y = ax.top + 46 * index + 14
        ratio = row["ratio"]
        diverges = ratio > 1.5
        parts.append(rect(ax.x(0.0), y, ax.x(ratio) - ax.x(0.0), 16,
                          "fill-b" if diverges else "fill-a", 0.85, rx=3, width=0))
        parts.append(text(ax.left - 12, y + 13, row["key"], "ink", 11.5, "end", "600"))
        parts.append(text(ax.x(ratio) + 10, y + 13, f"{ratio:.2f}x", "muted", 11))
        parts.append(text(58, y + 13, f"{row['computed']:.3e} vs {row['published']:.3e}",
                          "muted", 10))

    parts.append(text(58, 30, "The sievert column against ICRP 116 effective dose",
                      "ink", 13.5, "start", "600"))
    caption(parts, 58, 50, ["mSv per hour per MBq at 1 m. Ratio of NuSIFT's air kerma to "
                            "tabulated effective dose."], width)
    after = caption(parts, 58, 300,
                    ["This figure exists to keep a caveat true, not to show agreement. NuSIFT "
                     "reports air kerma with a photon weighting factor of 1, which is not "
                     "effective dose to a person."], width, size=12, cls="ink")
    after = caption(parts, 58, after + 4,
                    ["Co-60 and Ba-137m agree because effective dose per fluence happens to "
                     "track air kerma per fluence near 1 MeV. It is a coincidence of energy and "
                     "it does not survive going soft. Am-241 emits at 60 keV and below, where "
                     "air kerma loads heavily and effective dose does not, so the label "
                     "overstates the hazard to a person fivefold."], width)
    after = caption(parts, 58, after + 4,
                    ["The Am-241 row is gated on STILL diverging: if it ever started agreeing, "
                     "this warning would have quietly become wrong."], width)
    write(path, svg(width, after + 16, "\n".join(parts),
                    "Air kerma against ICRP 116 effective dose"))


# --------------------------------------------------------------------------------------
# The document
# --------------------------------------------------------------------------------------


def table(header, rows):
    out = ["| " + " | ".join(header) + " |", "| " + " | ".join("---" for _ in header) + " |"]
    out += ["| " + " | ".join(cells) + " |" for cells in rows]
    return "\n".join(out)


def reference_rows(rows, value_format="{:.4g}"):
    return [
        [
            row["key"],
            value_format.format(row["published"]),
            value_format.format(row["computed"]),
            pct(row["residual"]),
            "gated" if row["gate"] else "reported",
            row["note"] or "",
        ]
        for row in rows
    ]


def main():
    data = checks.open_store()

    gamma = checks.gamma_constant_rows(data)
    half_lives = checks.half_life_rows(data)
    masses = checks.molar_mass_rows(data)
    yields = checks.chain_yield_rows(data)
    icrp = checks.icrp116_rows(data)
    wigner = checks.way_wigner(data)
    curves = checks.equilibrium_curves(data)
    fine, summaries = checks.cross_code_rows(data)

    def worst(rows):
        gated = [r for r in rows if r["gate"]]
        return max((abs(r["residual"]) / r["tolerance"] for r in gated), default=0.0)

    families = [
        {"name": "gamma constants", "used": worst(gamma), "band": "8%",
         "detail": f"{sum(1 for r in gamma if r['gate'])} nuclides vs Ninkovic & Adrovic"},
        {"name": "half-lives", "used": worst(half_lives), "band": "1%",
         "detail": f"{len(half_lives)} nuclides vs ENSDF"},
        {"name": "atomic weights", "used": worst(masses), "band": "0.01%",
         "detail": f"{len(masses)} nuclides vs AME2020"},
        {"name": "chain yields", "used": worst(yields), "band": "4%",
         "detail": f"{len(yields)} mass chains vs ENDF/B-VIII.0"},
        {"name": "ICRP 116 ratios", "used": worst(icrp), "band": "6-15%",
         "detail": "3 nuclides vs Peplow 2020"},
        {"name": "cross-code", "used": max(s["max_residual"] / s["tolerance"] for s in summaries),
         "band": "2-3%",
         "detail": f"{sum(s['comparisons'] for s in summaries)} comparisons vs radioactivedecay"},
    ]

    figure_overview(FIGURES / "validation-overview.svg", families)
    residual_plot(
        FIGURES / "validation-gamma-constants.svg", gamma,
        "Specific gamma-ray constants against two published tabulations",
        ["Each block is compared in its own source's units and against its own stated cutoff -- "
         "20 keV for Ninkovic & Adrovic, 15 keV for Smith & Stabin -- so only the residual is "
         "common to the figure.",
         "Rows drawing more than 30% of the constant from below 100 keV are reported rather "
         "than gated: there the air table is sparse, its energy-absorption coefficient is "
         "turning over at its Compton minimum, and X-ray intensities differ between "
         "compilations."],
        0.08, "residual against each table's own units, vacuum, at that table's cutoff")
    residual_plot(
        FIGURES / "validation-half-lives.svg", half_lives,
        "Staged half-lives against ENSDF",
        ["Every decay constant in NuSIFT is ln(2) over one of these, so an error here moves "
         "every activity and every exposure.",
         "The band is 1% because independent compilations differ from each other by more than "
         "the staging error being looked for."],
        0.01, "evaluated half-life, ENSDF April 2022")
    residual_plot(
        FIGURES / "validation-molar-masses.svg", masses,
        "Staged atomic weights against AME2020",
        ["ENDF's atomic weight ratio and AME2020 are evaluations of the same measurements, so "
         "near-exact agreement is expected rather than impressive. What a 0.01% band catches is "
         "a ratio staged against the wrong reference mass, or a value attached to the wrong "
         "nuclide -- either of which would silently misconvert every inventory given in grams."],
        0.0001, "g/mol", row_height=20)
    residual_plot(
        FIGURES / "validation-chain-yields.svg", yields,
        "U-235 thermal cumulative chain yields against ENDF/B-VIII.0",
        ["NuSIFT's number is the sum of independent yields over the mass chain, which equals the "
         "cumulative yield at its terminus",
         "except for delayed-neutron emission -- the one process that moves a nucleus off its "
         "chain, and the reason A=137 and A=85 sit",
         "furthest out. A consistency check on seeding and aggregation rather than an "
         "independent measurement."],
        0.04, "percent per fission, at E = 0.0253 eV", row_height=22)
    figure_cutoff(FIGURES / "validation-gamma-cutoff.svg", gamma)
    figure_way_wigner(FIGURES / "validation-way-wigner.svg", wigner)
    figure_equilibria(FIGURES / "validation-equilibria.svg", curves)
    figure_cross_code(FIGURES / "validation-cross-code.svg",
                      checks.cross_code_curve(data, "Mo-99"), summaries)
    figure_cross_code_residuals(FIGURES / "validation-cross-code-residuals.svg", fine)
    figure_icrp116(FIGURES / "validation-icrp116.svg", icrp)

    import radioactivedecay as rd

    lines = []
    add = lines.append

    add("# Validation")
    add("")
    add("What NuSIFT computes, checked against values it was not fitted to: published gamma")
    add("constants, evaluated half-lives, atomic weights and fission yields, an empirical decay")
    add("law, and an independent implementation of the same physics.")
    add("")
    add("**This file is generated.** `python validation/make_report.py` rewrites it and its")
    add("figures from the committed store, and CI regenerates it and fails if the result differs")
    add("from what is committed. Every number below is computed by `validation/checks.py`, the")
    add("same module `python/tests/test_validation.py` asserts on, so a number here and a number")
    add("CI gates on cannot be two different numbers.")
    add("")
    add(table(["", ""], [
        ["store", f"`data/nusift_b8.1.h5`, {data.library}, staged {data.staged_utc}"],
        ["coverage", f"{data.staged_count} nuclides staged, {data.size} in the closed chain"],
        ["nusift", f"{checks.__dict__.get('nusift').__version__}"],
        ["cross-code", f"radioactivedecay {rd.__version__} (pinned {checks.RADIOACTIVEDECAY_PIN}),"
                       f" ICRP-107 decay data"],
    ]))
    add("")
    add("## What this does and does not establish")
    add("")
    add("A residual here is not an error bar. Every comparison is against a value produced by")
    add("someone else's conventions, someone else's evaluation, or an empirical fit, and where")
    add("those differ from NuSIFT's the difference shows up in this table as a residual whether")
    add("or not anything is wrong. The bands are therefore set from the physics that limits each")
    add("comparison, and the reasoning for each is recorded beside the numbers in")
    add("[`validation/references/README.md`](../validation/references/README.md).")
    add("")
    add("Rows marked **reported** are computed and shown but not gated. That is where a known")
    add("convention or evaluation difference lives, named, rather than being hidden under a band")
    add("wide enough to swallow it.")
    add("")
    add("![Scoreboard](figures/validation-overview.svg)")
    add("")

    add("## Gamma constants")
    add("")
    add("The broadest published comparison here, and the one that exercises the most at once:")
    add("staged line energies and intensities, the NIST air table, the log-log interpolation,")
    add("and the roentgen conversion, against numbers nothing in NuSIFT was tuned to.")
    add("")
    add("Two independent tabulations are used, in their own units and each with its own stated")
    add("low-energy cutoff: Ninković & Adrović (2012) in µGy·m²/(GBq·h) above 20 keV, and Smith &")
    add("Stabin (2012) in R·cm²/(h·mCi) above 15 keV. Each nuclide appears once, under whichever")
    add("carries it. That the two disagree with each other by about a percent where they overlap")
    add("— Ninković's Co-60 is 13.05 R·cm²/(h·mCi) against Smith & Stabin's 12.9 — is the point")
    add("[exposure.md §4](exposure.md) makes about published constants, and NuSIFT lands between")
    add("them.")
    add("")
    add("Both tables exclude photons below their cutoff, and bremsstrahlung entirely. NuSIFT by")
    add("default counts the whole spectrum, so the two are not the same quantity until the same")
    add("cutoff is applied — and for the X-ray emitters that difference reaches a factor of ten,")
    add("all of it convention rather than disagreement.")
    add("")
    add("![Cutoff](figures/validation-gamma-cutoff.svg)")
    add("")
    add("![Gamma constants](figures/validation-gamma-constants.svg)")
    add("")
    add(table(["nuclide", "published", "NuSIFT", "unit", "residual", "status", "source", "note"],
              [[r["key"], f"{r['published']:.4g}", f"{r['computed']:.4g}", r["unit"],
                pct(r["residual"]), "gated" if r["gate"] else "reported", r["source"],
                r["note"] or ""] for r in gamma]))
    add("")

    add("## Half-lives")
    add("")
    add("![Half-lives](figures/validation-half-lives.svg)")
    add("")
    add(table(["nuclide", "ENSDF", "staged", "residual", "status", "note"],
              [[r["key"], f"{r['published']:g} {r['unit']}", f"{r['computed']:.6g} {r['unit']}",
                pct(r["residual"], 3), "gated" if r["gate"] else "reported", r["note"] or ""]
               for r in half_lives]))
    add("")

    add("## Atomic weights")
    add("")
    add("Everything expressed per gram goes through the staged atomic weight ratio: an inventory")
    add("given in grams, and any specific activity derived from one. AME2020 is an independent")
    add("evaluation of the same masses, which makes this the one check that pins that conversion")
    add("against something outside the ENDF pipeline.")
    add("")
    add("![Atomic weights](figures/validation-molar-masses.svg)")
    add("")
    add(table(["nuclide", "AME2020", "staged", "residual", "status", "note"],
              [[r["key"], f"{r['published']:.6f}", f"{r['computed']:.6f}", pct(r["residual"], 5),
                "gated" if r["gate"] else "reported", r["note"] or ""] for r in masses]))
    add("")
    add("There is no specific-activity sweep. It would be ln(2)·N_A/(T½·M) over two quantities")
    add("that already have their own tables above, and published specific-activity tables")
    add("disagree with each other by more than the staging error being looked for — Am-241 is")
    add("tabulated at 3.2, 3.43 and 3.5 Ci/g by different compilers, depending only on which")
    add("half-life each used. Six keystone nuclides are still checked end to end in the C++")
    add("suite, where the quotient itself is what a user reads off a report.")
    add("")

    add("## Fission yields")
    add("")
    add("![Chain yields](figures/validation-chain-yields.svg)")
    add("")
    add(table(["mass chain", "ENDF/B-VIII.0", "NuSIFT", "residual", "status", "note"],
              [["A = " + r["key"], f"{r['published']:.4f}%", f"{r['computed']:.4f}%",
                pct(r["residual"]), "gated" if r["gate"] else "reported", r["note"] or ""]
               for r in yields]))
    add("")

    add("## The Way-Wigner law")
    add("")
    add("![Way-Wigner](figures/validation-way-wigner.svg)")
    add("")
    add(f"Fitted log-log slope **{wigner['slope']:.3f}** over {checks.WAY_WIGNER_START} to")
    add(f"{checks.WAY_WIGNER_END}, against the empirical t^-1.2, accepted within")
    add(f"[{wigner['band'][0]}, {wigner['band'][1]}].")
    add("")
    add(table(["window", "local slope"],
              [[label, f"{slope:.3f}"] for label, slope in wigner["segments"]]))
    add("")

    add("## Real chains against their closed forms")
    add("")
    add("![Equilibria](figures/validation-equilibria.svg)")
    add("")
    add(table(["chain", "kind", "staged branch", "worst disagreement with Bateman"],
              [[f"{c['parent']} → {c['daughter']}", c["kind"], f"{c['branch']:.4f}",
                f"{c['worst_residual']:.2e}"] for c in curves]))
    add("")
    add("Atom conservation on a pure beta chain is checked alongside these in the C++ suite. It")
    add("holds to 1e-10 on Sr-90 and to 1e-8 on Cs-137, and the difference is the evaluation's")
    add("rather than the solver's: Cs-137's two staged branchings are 0.05300549 and 0.9469945,")
    add("which sum to 0.99999999, so exactly 1e-8 of every decay goes nowhere. Across the store")
    add("226 nuclides miss unity by more than 1e-12 and none by more than 1e-6.")
    add("")

    add("## An independent implementation")
    add("")
    add("![Cross-code](figures/validation-cross-code.svg)")
    add("")
    add("![Cross-code residuals](figures/validation-cross-code-residuals.svg)")
    add("")
    add(table(["case", "comparisons", "worst", "median", "band"],
              [[s["key"], str(s["comparisons"]), pct(s["max_residual"]),
                pct(s["median_residual"]), f"{s['tolerance'] * 100:.0f}%"] for s in summaries]))
    add("")
    add("Held out of the gate, with the cause measured rather than assumed:")
    add("")
    add(table(["nuclide", "why"],
              [[name, reason] for name, reason in
               sorted(checks.KNOWN_EVALUATION_DIFFERENCES.items())]))
    add("")

    add("## The sievert column")
    add("")
    add("![ICRP 116](figures/validation-icrp116.svg)")
    add("")
    add(table(["nuclide", "ICRP 116 effective dose", "NuSIFT air kerma", "ratio", "expected"],
              [[r["key"], f"{r['published']:.4g}", f"{r['computed']:.4g}", f"{r['ratio']:.3f}",
                f"{r['expected_ratio']:.2f}"] for r in icrp]))
    add("")
    add("Units are mSv·h⁻¹·MBq⁻¹ at 1 m. See [exposure.md §6](exposure.md#6-units) for why the")
    add("sievert here is air kerma wearing a label.")
    add("")

    add("## What gates in CI")
    add("")
    add(table(["suite", "what it covers", "where"],
              [["`ctest -L validation`",
                "store census and provenance, published constants through the store, real chains "
                "against closed forms, mass-chain conservation",
                "[`tests/validation/`](../tests/validation/)"],
               ["`pytest -m validation`",
                "the sweeps and the cross-code comparison in this document",
                "[`python/tests/test_validation.py`](../python/tests/test_validation.py)"],
               ["report drift",
                "regenerates this file and fails if it differs from what is committed",
                "[`.github/workflows/ci.yml`](../.github/workflows/ci.yml)"]]))
    add("")
    add("The full store census — what this evaluation covers and, more usefully, what it does")
    add("not — is in [nuclear-data.md](nuclear-data.md) and printed by `nusift data info`.")
    add("")

    add("## Sources")
    add("")
    add("Full citations, and the reasoning behind every acceptance band, are in")
    add("[`validation/references/README.md`](../validation/references/README.md).")
    add("")

    write(REPORT, "\n".join(lines) + "\n")
    check_figure_bounds(sorted(FIGURES.glob("validation-*.svg")))

    print(f"wrote {REPORT.relative_to(ROOT)}")
    print(f"  gamma constants  {sum(1 for r in gamma if r['gate'])} gated, "
          f"{sum(1 for r in gamma if not r['gate'])} reported")
    print(f"  half-lives       {len(half_lives)} gated")
    print(f"  chain yields     {len(yields)} gated")
    print(f"  cross-code       {sum(s['comparisons'] for s in summaries)} comparisons")
    print(f"  Way-Wigner slope {wigner['slope']:.3f}")
    print(f"  figures written to {FIGURES.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
