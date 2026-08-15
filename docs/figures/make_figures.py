"""Regenerate the data-driven figures in docs/figures.

    python docs/figures/make_figures.py        # numpy only

Nothing here is drawn freehand from memory. Each figure is built from one of:

  * the NIST air-coefficient table, parsed straight out of
    nusift/exposure/air_coefficients.cpp, so a figure cannot drift from the values the
    code interpolates;
  * verbatim CLI output committed under data/ (see data/README.md for the commands);
  * float64 arithmetic on a closed-form exponential, where the point is the conditioning
    of a subtraction rather than any particular nuclide.

The structural diagrams -- the ones with no measurements in them -- live in the docs as
mermaid or as hand-authored SVG, not here.
"""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import numpy as np

OUT = Path(__file__).parent
DATA = OUT / "data"
ROOT = OUT.parent.parent

# Cs-137: the nuclide the README's example inventory leads with.
HALF_LIFE_Y = 30.08
YEAR_S = 365.25 * 86400.0
LAMBDA = math.log(2.0) / (HALF_LIFE_Y * YEAR_S)
N0 = 1.0e20
EPS = np.finfo(float).eps
AIR_DENSITY = 1.205  # kg/m^3, the default geometry


# --------------------------------------------------------------------------------------
# SVG scaffolding. Every figure paints an explicit background so it stays legible whether
# or not the viewer honours prefers-color-scheme; the media query then adapts it to a dark
# page rather than leaving a light card sitting in one.
# --------------------------------------------------------------------------------------

STYLE = """
  .bg     { fill: #fbfbfa; stroke: #d8d6d1; }
  .ink    { fill: #1c1b19; }
  .muted  { fill: #6b6862; }
  .axis   { stroke: #b4b1ab; fill: none; }
  .grid   { stroke: #e6e4df; fill: none; }
  .curve  { stroke: #2f6f8f; fill: none; }
  .fill-a { fill: #2f6f8f; }
  .fill-b { fill: #b4551f; }
  .accent { stroke: #b4551f; fill: none; }
  .warn   { fill: #b4551f; }
  .band   { fill: #b4551f; }
  .box    { fill: none; stroke: #b4b1ab; }
  .s0 { stroke: #2f6f8f; fill: none; } .t0 { fill: #2f6f8f; }
  .s1 { stroke: #b4551f; fill: none; } .t1 { fill: #b4551f; }
  .s2 { stroke: #4d7c4d; fill: none; } .t2 { fill: #4d7c4d; }
  .s3 { stroke: #7a5aa0; fill: none; } .t3 { fill: #7a5aa0; }
  .s4 { stroke: #a03a5a; fill: none; } .t4 { fill: #a03a5a; }
  .s5 { stroke: #7a6a3a; fill: none; } .t5 { fill: #7a6a3a; }
  .s6 { stroke: #3a7a7a; fill: none; } .t6 { fill: #3a7a7a; }
  .s7 { stroke: #8a5a3a; fill: none; } .t7 { fill: #8a5a3a; }
  text    { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif; }
  @media (prefers-color-scheme: dark) {
    .bg     { fill: #12161c; stroke: #2b3239; }
    .ink    { fill: #e6edf3; }
    .muted  { fill: #9aa4ae; }
    .axis   { stroke: #4a545e; }
    .grid   { stroke: #232a31; }
    .curve  { stroke: #6cb6d9; }
    .fill-a { fill: #6cb6d9; }
    .fill-b { fill: #e08a4c; }
    .accent { stroke: #e08a4c; }
    .warn   { fill: #e08a4c; }
    .band   { fill: #e08a4c; }
    .box    { stroke: #4a545e; }
    .s0 { stroke: #6cb6d9; } .t0 { fill: #6cb6d9; }
    .s1 { stroke: #e08a4c; } .t1 { fill: #e08a4c; }
    .s2 { stroke: #7fbf7f; } .t2 { fill: #7fbf7f; }
    .s3 { stroke: #b49ae0; } .t3 { fill: #b49ae0; }
    .s4 { stroke: #e07f9f; } .t4 { fill: #e07f9f; }
    .s5 { stroke: #c9b56a; } .t5 { fill: #c9b56a; }
    .s6 { stroke: #6ac9c9; } .t6 { fill: #6ac9c9; }
    .s7 { stroke: #d19a72; } .t7 { fill: #d19a72; }
  }
"""


def svg(width, height, body, title):
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'width="{width}" height="{height}" role="img" aria-label="{title}">\n'
        f"<style>{STYLE}</style>\n"
        f'<rect class="bg" x="0.5" y="0.5" width="{width - 1}" height="{height - 1}" rx="6"/>\n'
        f"{body}\n</svg>\n"
    )


def esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def text(x, y, s, cls="ink", size=12, anchor="start", weight="normal"):
    return (
        f'<text class="{cls}" x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
        f'text-anchor="{anchor}" font-weight="{weight}">{esc(s)}</text>'
    )


def polyline(pts, cls, width=1.8, dash=None):
    if len(pts) < 2:
        return ""
    d = " ".join(f"{x:.2f},{y:.2f}" for x, y in pts)
    stroke = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<polyline class="{cls}" points="{d}" stroke-width="{width}" stroke-linejoin="round"{stroke}/>'


def rect(x, y, w, h, cls, opacity=1.0, rx=0, width=1.2):
    return (
        f'<rect class="{cls}" x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
        f'rx="{rx}" fill-opacity="{opacity}" stroke-width="{width}"/>'
    )


def line(x1, y1, x2, y2, cls, width=1.2, dash=None):
    stroke = f' stroke-dasharray="{dash}"' if dash else ""
    return (
        f'<line class="{cls}" x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
        f'stroke-width="{width}"{stroke}/>'
    )


class Axes:
    """Minimal linear/log axis mapper with a frame, grid, and ticks."""

    def __init__(self, left, right, top, bottom, xlim, ylim, xlog=False, ylog=False):
        self.left, self.right, self.top, self.bottom = left, right, top, bottom
        self.xlog, self.ylog = xlog, ylog
        self.x0, self.x1 = (math.log10(v) for v in xlim) if xlog else xlim
        self.y0, self.y1 = (math.log10(v) for v in ylim) if ylog else ylim

    def x(self, v):
        v = math.log10(v) if self.xlog else v
        return self.left + (self.right - self.left) * (v - self.x0) / (self.x1 - self.x0)

    def y(self, v):
        if self.ylog:
            v = math.log10(max(v, 10.0**self.y0))
        return self.bottom - (self.bottom - self.top) * (v - self.y0) / (self.y1 - self.y0)

    def frame(self):
        return [
            line(self.left, self.bottom, self.right, self.bottom, "axis", 1.4),
            line(self.left, self.top, self.left, self.bottom, "axis", 1.4),
        ]

    def hgrid(self, values, labels=None, size=10.5):
        out = []
        for i, v in enumerate(values):
            y = self.y(v)
            out.append(line(self.left, y, self.right, y, "grid"))
            if labels:
                out.append(text(self.left - 9, y + 3.5, labels[i], "muted", size, "end"))
        return out

    def vgrid(self, values, labels=None, size=10.5):
        out = []
        for i, v in enumerate(values):
            x = self.x(v)
            out.append(line(x, self.top, x, self.bottom, "grid"))
            if labels:
                out.append(text(x, self.bottom + 18, labels[i], "muted", size, "middle"))
        return out

    def path(self, xs, ys, cls, width=1.8, dash=None):
        return polyline([(self.x(a), self.y(b)) for a, b in zip(xs, ys)], cls, width, dash)


# --------------------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------------------

def air_table():
    """The NIST dry-air table, parsed out of the source that uses it."""
    source = (ROOT / "nusift" / "exposure" / "air_coefficients.cpp").read_text(encoding="utf-8")
    body = source[source.index("kAir[] = {") : source.index("};", source.index("kAir[] = {"))]
    rows = re.findall(r"\{\s*([0-9.e+-]+),\s*([0-9.e+-]+),\s*([0-9.e+-]+)\s*\}", body)
    table = np.array([[float(a) for a in row] for row in rows])
    if len(table) < 20:
        raise SystemExit(f"parsed only {len(table)} air-coefficient rows; the table moved")
    return table[:, 0], table[:, 1], table[:, 2]


AIR_E, AIR_MU, AIR_MUEN = air_table()


def air_interp(energy_ev, values):
    """Log-log interpolation, clamped -- the same rule air_coefficients.cpp applies."""
    e = np.clip(np.asarray(energy_ev, dtype=float), AIR_E[0], AIR_E[-1])
    return np.exp(np.interp(np.log(e), np.log(AIR_E), np.log(values)))


def read_csv(name):
    with (DATA / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def cumulative(t):
    return N0 / LAMBDA * (1.0 - np.exp(-LAMBDA * t))


def stable_interval(t1, dt):
    return N0 / LAMBDA * np.exp(-LAMBDA * t1) * -np.expm1(-LAMBDA * dt)


# --------------------------------------------------------------------------------------
# Figure -- the interval as a difference of two areas measured from zero.
# --------------------------------------------------------------------------------------

def figure_area(path):
    w, h = 760, 330
    ax = Axes(66, 700, 40, 250, (0.0, 60.0 * YEAR_S), (0.0, N0))
    t1, t2 = 12.0 * YEAR_S, 30.0 * YEAR_S
    ts = np.linspace(0.0, 60.0 * YEAR_S, 400)

    parts = []
    for frac in (0.25, 0.5, 0.75, 1.0):
        y = ax.y(frac * N0)
        parts.append(line(ax.left, y, ax.right, y, "grid"))

    def area(a, b, cls, opacity):
        seg = np.concatenate(([a], ts[(ts >= a) & (ts <= b)], [b]))
        pts = [(ax.x(a), ax.bottom)]
        pts += [(ax.x(t), ax.y(N0 * math.exp(-LAMBDA * t))) for t in seg]
        pts += [(ax.x(b), ax.bottom)]
        d = " ".join(f"{x:.2f},{y:.2f}" for x, y in pts)
        return f'<polygon class="{cls}" points="{d}" fill-opacity="{opacity}"/>'

    parts.append(area(0.0, t1, "fill-a", 0.16))
    parts.append(area(t1, t2, "fill-b", 0.42))
    parts.append(ax.path(ts, N0 * np.exp(-LAMBDA * ts), "curve", 2.0))

    for t, label in ((t1, "t₁"), (t2, "t₂")):
        parts.append(
            line(ax.x(t), ax.y(N0 * math.exp(-LAMBDA * t)), ax.x(t), ax.bottom, "axis", 1.2, "3 3")
        )
        parts.append(text(ax.x(t), ax.bottom + 18, label, "ink", 13, "middle", "600"))

    parts += ax.frame()
    parts.append(text(ax.left - 10, ax.top + 4, "n₀", "muted", 11, "end"))
    parts.append(text(ax.left - 10, ax.bottom + 4, "0", "muted", 11, "end"))
    parts.append(text(ax.right, ax.bottom + 34, "cooling time →", "muted", 11, "end"))
    parts.append(text(ax.left, 22, "atoms of one nuclide, n(τ)", "muted", 12))

    mid_a = ax.x(t1 * 0.42)
    parts.append(text(mid_a, ax.bottom - 26, "G(t₁)", "ink", 13, "middle", "600"))
    parts.append(text(mid_a, ax.bottom - 10, "cumulative from 0", "muted", 10, "middle"))
    mid_b = (ax.x(t1) + ax.x(t2)) / 2
    parts.append(text(mid_b, ax.bottom - 20, "the answer", "ink", 12, "middle", "600"))
    parts.append(text(mid_b, ax.bottom - 6, "atom·seconds", "muted", 10, "middle"))

    parts.append(text(ax.left, 288, "∫", "ink", 15, "start", "600"))
    parts.append(text(ax.left + 12, 288, "n(τ) dτ  over [t₁, t₂]  =  G(t₂) − G(t₁)", "ink", 13))
    parts.append(
        text(
            ax.left,
            310,
            "Both areas come from one solve each. No grid is laid inside the window, so no half-life can fall between its points.",
            "muted",
            11.5,
        )
    )
    path.write_text(svg(w, h, "\n".join(parts), "The interval integral as a difference of areas"), encoding="utf-8")


# --------------------------------------------------------------------------------------
# Figure -- what the subtraction costs, measured.
# --------------------------------------------------------------------------------------

def figure_cancellation(path):
    w, h = 760, 430
    ax = Axes(78, 690, 46, 320, (1e-3, 1e9), (1e-17, 1e-2), xlog=True, ylog=True)

    t1 = 30.0 * YEAR_S
    dts = np.logspace(-3, 9.0, 400)
    exact = stable_interval(t1, dts)
    naive = cumulative(t1 + dts) - cumulative(t1)
    rel = np.abs(naive - exact) / exact
    retained = (cumulative(t1 + dts) - cumulative(t1)) / cumulative(t1 + dts)

    parts = []
    decades = [10.0**e for e in range(-17, -1, 3)]
    parts += ax.hgrid(decades, [f"1e{e}" for e in range(-17, -1, 3)])
    parts += ax.vgrid([10.0**e for e in range(-3, 10, 3)])

    idx = int(np.argmax(retained >= 1e-8))
    dt_guard = dts[idx]
    xg = ax.x(dt_guard)
    parts.append(rect(ax.left, ax.top, xg - ax.left, ax.bottom - ax.top, "band", 0.09, width=0))
    parts.append(line(xg, ax.top, xg, ax.bottom, "accent", 1.6, "5 4"))

    parts.append(ax.path(dts, EPS / retained, "accent", 1.6, "4 3"))
    parts.append(ax.path(dts, rel, "curve", 1.6))
    parts += ax.frame()

    for dt, label in ((1e-3, "1 ms"), (1.0, "1 s"), (3600.0, "1 h"), (86400.0, "1 d"),
                      (2.63e6, "1 mo"), (3.156e7, "1 y"), (3.156e8, "10 y")):
        x = ax.x(dt)
        parts.append(line(x, ax.bottom, x, ax.bottom + 5, "axis", 1.2))
        parts.append(text(x, ax.bottom + 19, label, "muted", 10.5, "middle"))

    parts.append(text(ax.left, 24, "relative error of G(t₂) − G(t₁), window starting at t₁ = 30 y", "ink", 13, "start", "600"))
    parts.append(text(ax.right, ax.bottom + 40, "window width t₂ − t₁ →", "muted", 11, "end"))
    parts.append(text(xg - 10, ax.top + 18, "guard fires", "warn", 11.5, "end", "600"))
    parts.append(text(xg - 10, ax.top + 33, "re-solve from t₁", "muted", 10.5, "end"))
    parts.append(text(xg + 10, ax.top + 18, "plain subtraction kept", "muted", 11, "start"))

    lx = ax.x(10.0**4.2)
    parts.append(line(lx, ax.top + 76, lx + 26, ax.top + 76, "curve", 1.8))
    parts.append(text(lx + 34, ax.top + 80, "measured in float64", "muted", 11))
    parts.append(line(lx, ax.top + 96, lx + 26, ax.top + 96, "accent", 1.8, "4 3"))
    parts.append(text(lx + 34, ax.top + 100, "ε ÷ retained fraction", "muted", 11))

    parts.append(text(ax.left, 372, "The floor the guard buys: retained ≥ 1e-8 bounds this curve at about 1e-8.", "ink", 12))
    parts.append(text(ax.left, 392, "A one-second window here retains 7.3e-10 of G(t₂) and loses roughly half the mantissa; the guard turns it into a third solve.", "muted", 11.5))
    parts.append(text(ax.left, 410, "Closed-form G in float64 — a CRAM solve adds its own error, so this is a floor on the loss, not an estimate of it.", "muted", 11.5))

    path.write_text(svg(w, h, "\n".join(parts), "Relative error of the naive difference against window width"), encoding="utf-8")
    print(f"  guard threshold at dt = {dt_guard:.3g} s")


# --------------------------------------------------------------------------------------
# Figure -- the air coefficients, and why they cannot be a single number.
# --------------------------------------------------------------------------------------

def figure_air_coefficients(path):
    w, h = 760, 418
    ax = Axes(80, 610, 46, 300, (1e4, 1e7), (1e-3, 1.0), xlog=True, ylog=True)

    grid = np.logspace(4, 7, 300)
    parts = []
    parts += ax.hgrid([1e-3, 1e-2, 1e-1, 1.0], ["0.001", "0.01", "0.1", "1"])
    parts += ax.vgrid([1e4, 1e5, 1e6, 1e7], ["10 keV", "100 keV", "1 MeV", "10 MeV"])

    parts.append(ax.path(grid, air_interp(grid, AIR_MU), "s0", 2.0))
    parts.append(ax.path(grid, air_interp(grid, AIR_MUEN), "s1", 2.0))
    for e, mu, muen in zip(AIR_E, AIR_MU, AIR_MUEN):
        parts.append(f'<circle class="t0" cx="{ax.x(e):.1f}" cy="{ax.y(mu):.1f}" r="2.1"/>')
        parts.append(f'<circle class="t1" cx="{ax.x(e):.1f}" cy="{ax.y(muen):.1f}" r="2.1"/>')
    parts += ax.frame()

    parts.append(text(ax.left, 24, "dry-air mass coefficients, NIST tabulation", "ink", 13, "start", "600"))
    parts.append(text(ax.left - 62, ax.top - 12, "m²/kg", "muted", 11))
    parts.append(text(ax.right, ax.bottom + 36, "photon energy →", "muted", 11, "end"))

    parts.append(text(ax.x(2.3e4), ax.y(0.30), "μ/ρ", "t0", 13, "start", "600"))
    parts.append(text(ax.x(2.3e4), ax.y(0.30) + 14, "attenuation along the path", "muted", 10))
    parts.append(text(ax.x(1.1e5), ax.y(0.0055), "μ_en/ρ", "t1", 13, "start", "600"))
    parts.append(text(ax.x(1.1e5), ax.y(0.0055) + 14, "energy absorbed at the point", "muted", 10))

    # Where real decay photons actually sit.
    for e, label in ((6.617e5, "Ba-137m 662 keV"), (1.3325e6, "Co-60 1333 keV")):
        parts.append(line(ax.x(e), ax.top, ax.x(e), ax.bottom, "axis", 1.0, "2 4"))
    parts.append(text(ax.x(6.617e5), ax.top + 14, "Ba-137m", "muted", 10, "middle"))
    parts.append(text(ax.x(6.617e5), ax.top + 26, "662 keV", "muted", 10, "middle"))
    parts.append(text(ax.x(1.3325e6), ax.top + 14, "Co-60", "muted", 10, "middle"))
    parts.append(text(ax.x(1.3325e6), ax.top + 26, "1333 keV", "muted", 10, "middle"))

    span = air_interp(1e4, AIR_MU) / air_interp(1e7, AIR_MU)
    parts.append(text(ax.left, 358, f"μ/ρ spans a factor of {span:.0f} across the tabulated range, so it cannot come out of the sum over lines.", "ink", 12))
    parts.append(text(ax.left, 378, "Points are NIST's own energy grid; the curve is the log-log interpolation the code applies between them.", "muted", 11.5))
    parts.append(text(ax.left, 396, "Outside 10 keV – 10 MeV both coefficients are clamped to the end value, which is reported by `nusift data info`.", "muted", 11.5))

    path.write_text(svg(w, h, "\n".join(parts), "Air mass attenuation and mass energy-absorption coefficients against photon energy"), encoding="utf-8")
    print(f"  mu/rho span 10 keV to 10 MeV: {span:.1f}x")


def figure_air_attenuation(path):
    w, h = 760, 418
    ax = Axes(80, 600, 46, 300, (0.0, 20.0), (0.3, 1.0))

    parts = []
    parts += ax.hgrid([0.4, 0.6, 0.8, 1.0], ["40%", "60%", "80%", "100%"])
    parts += ax.vgrid([0, 5, 10, 15, 20], ["0", "5 m", "10 m", "15 m", "20 m"])

    d = np.linspace(0.0, 20.0, 200)
    lines = [
        (3.0e4, "30 keV", "s4"),
        (8.0e4, "80 keV", "s3"),
        (6.617e5, "662 keV  (Ba-137m)", "s0"),
        (1.3325e6, "1333 keV  (Co-60)", "s2"),
    ]
    transmissions = {}
    for energy, label, cls in lines:
        mu = float(air_interp(energy, AIR_MU)) * AIR_DENSITY
        t = np.exp(-mu * d)
        transmissions[label] = t[-1]
        parts.append(ax.path(d, t, cls, 2.0))
        parts.append(text(ax.right + 8, ax.y(t[-1]) + 4, label, cls.replace("s", "t"), 11, "start", "600"))
    parts += ax.frame()

    parts.append(text(ax.left, 24, "fraction of emitted photons surviving the air path", "ink", 13, "start", "600"))
    parts.append(text(ax.right - 60, ax.bottom + 36, "distance from the source →", "muted", 11, "end"))

    lo = transmissions["30 keV"]
    hi = transmissions["1333 keV  (Co-60)"]
    parts.append(text(ax.left, 358, f"At 20 m a 30 keV line has lost {(1 - lo) * 100:.0f}% of its photons and a 1333 keV line {(1 - hi) * 100:.0f}%.", "ink", 12))
    parts.append(text(ax.left, 378, "The curves are not parallel, so their ratio depends on distance — which is exactly why one per-nuclide", "muted", 11.5))
    parts.append(text(ax.left, 396, "exposure constant cannot be right at more than one distance, and why the store keeps whole spectra.", "muted", 11.5))

    path.write_text(svg(w, h, "\n".join(parts), "Air transmission against distance for four photon energies"), encoding="utf-8")
    print(f"  transmission at 20 m: 30 keV {lo:.3f}, 1333 keV {hi:.3f}")


# --------------------------------------------------------------------------------------
# Figure -- who leads, and when that changes. Real ENDF/B-VIII.1 output.
# --------------------------------------------------------------------------------------

def figure_dominance(path):
    rows = read_csv("fission-exposure-timeline.csv")
    windows = read_csv("fission-dominance-windows.csv")

    times = sorted({float(r["time_s"]) for r in rows})
    series = {}
    for r in rows:
        series.setdefault(r["contributor"], {})[float(r["time_s"])] = float(r["fraction"])

    leaders = [win["label"] for win in windows]
    # Draw the contributors that actually lead, plus the next most prominent, so the plot
    # shows the competition rather than only the winners.
    extra = sorted(series, key=lambda k: -max(series[k].values()))
    ordered = leaders + [k for k in extra if k not in leaders][:2]

    w, h = 900, 470
    strip_top, strip_h = 46, 34
    ax = Axes(74, 792, 112, 336, (times[0], times[-1]), (0.0, 1.0), xlog=True)

    parts = []
    parts += ax.hgrid([0.0, 0.25, 0.5, 0.75, 1.0], ["0", "25%", "50%", "75%", "100%"])
    decades = [60, 600, 3600, 86400, 864000, 3.156e7, 3.156e8, 3.156e9]
    labels = ["1 m", "10 m", "1 h", "1 d", "10 d", "1 y", "10 y", "100 y"]
    keep = [(v, l) for v, l in zip(decades, labels) if times[0] <= v <= times[-1]]
    parts += ax.vgrid([v for v, _ in keep], [l for _, l in keep])

    # The leader strip: one band per dominance window, labels staggered so the narrow ones fit.
    for i, win in enumerate(windows):
        x0, x1 = ax.x(float(win["start_s"])), ax.x(float(win["end_s"]))
        cls = f"t{i % 8}"
        parts.append(rect(x0, strip_top, x1 - x0, strip_h, cls, 0.22, rx=2, width=0))
        parts.append(line(x0, strip_top, x0, ax.bottom, "axis", 0.9, "2 4"))
        row = i % 2
        parts.append(
            text((x0 + x1) / 2, strip_top + 14 + row * 13, win["label"], cls, 10.5, "middle", "600")
        )
    parts.append(rect(ax.x(times[0]), strip_top, ax.x(times[-1]) - ax.x(times[0]), strip_h, "box", 0, rx=2))
    parts.append(text(ax.left, strip_top - 8, "leads the exposure", "muted", 11))

    for i, name in enumerate(ordered):
        pts = sorted(series[name].items())
        cls = f"s{i % 8}"
        parts.append(ax.path([t for t, _ in pts], [f for _, f in pts], cls, 1.9))
        best_t, best_f = max(pts, key=lambda p: p[1])
        parts.append(text(ax.x(best_t), ax.y(best_f) - 7, name, cls.replace("s", "t"), 10, "middle", "600"))
    parts += ax.frame()

    parts.append(text(ax.left, 24, "share of the total exposure rate — 20 kt U-235 thermal fission, ENDF/B-VIII.1", "ink", 13, "start", "600"))
    parts.append(text(ax.right, ax.bottom + 40, "time after fission →", "muted", 11, "end"))
    parts.append(text(ax.left, 396, f"{len(windows)} leadership changes over {len(times)} log-spaced samples; the boundaries are interpolated crossings, not sample times.", "ink", 12))
    parts.append(text(ax.left, 416, "A share this concentrated is a late-time property: at one minute the leader holds 11%, at a century it holds essentially all of it.", "muted", 11.5))
    parts.append(text(ax.left, 434, "Curves are the per-time top 8, so a contributor's line stops where it drops out of that set.", "muted", 11.5))

    path.write_text(svg(w, h, "\n".join(parts), "Exposure share against time for the dominant fission products, with dominance windows"), encoding="utf-8")
    print(f"  dominance windows: {len(windows)} over {len(times)} samples")


# --------------------------------------------------------------------------------------
# Figure -- the same source, the same instant, two metrics.
# --------------------------------------------------------------------------------------

def figure_metrics(path):
    activity = read_csv("fission-1d-activity.csv")
    exposure = read_csv("fission-1d-exposure.csv")

    w, h = 760, 420
    top, row_h = 96, 34
    left_x, right_x = 250, 510

    parts = []
    parts.append(text(380, 26, "20 kt U-235 fission, 1 day after — the same atoms, ranked twice", "ink", 13, "middle", "600"))
    parts.append(text(left_x, 66, "by activity", "ink", 12.5, "middle", "600"))
    parts.append(text(left_x, 82, "Bq", "muted", 10.5, "middle"))
    parts.append(text(right_x, 66, "by exposure at 1 m", "ink", 12.5, "middle", "600"))
    parts.append(text(right_x, 82, "R/h", "muted", 10.5, "middle"))

    colour = {}
    for i, r in enumerate(activity):
        colour[r["contributor"]] = i % 8
    for r in exposure:
        colour.setdefault(r["contributor"], len(colour) % 8)

    def place(rows, x, anchor, dx):
        out = {}
        for i, r in enumerate(rows):
            y = top + i * row_h
            name = r["contributor"]
            cls = f"t{colour[name] % 8}"
            out[name] = y
            out_label = f"{name}   {float(r['fraction']) * 100:.1f}%"
            out.setdefault("_parts", []).append(text(x + dx, y + 4, out_label, cls, 12, anchor, "600"))
        return out

    la = place(activity, left_x, "end", 46)
    le = place(exposure, right_x, "start", -46)
    parts += la.pop("_parts") + le.pop("_parts")

    for name, y0 in la.items():
        if name in le:
            y1 = le[name]
            cls = f"s{colour[name] % 8}"
            x0, x1 = left_x + 58, right_x - 58
            parts.append(
                f'<path class="{cls}" d="M {x0} {y0} C {(x0 + x1) / 2} {y0}, {(x0 + x1) / 2} {y1}, {x1} {y1}" '
                f'stroke-width="1.6" fill="none" stroke-opacity="0.75"/>'
            )
        else:
            parts.append(text(left_x + 66, y0 + 4, "—", "muted", 12))

    for name, y1 in le.items():
        if name not in la:
            parts.append(text(right_x - 66, y1 + 4, "—", "muted", 12, "end"))

    note = 380
    only_exposure = [n for n in le if n not in la]
    only_activity = [n for n in la if n not in le]
    parts.append(text(60, note - 20, f"{len(only_activity)} of the top 8 by activity do not appear in the top 8 by exposure, and {len(only_exposure)} appear only there.", "ink", 12))
    parts.append(text(60, note, "A pure beta emitter can lead a decay count and contribute no exposure at all; a nuclide can lead exposure through photons its", "muted", 11.5))
    parts.append(text(60, note + 18, "daughter emits. Ranking by the wrong one is not a rounding error — it names a different nuclide.", "muted", 11.5))

    path.write_text(svg(w, h, "\n".join(parts), "The same fission source at one day ranked by activity and by exposure"), encoding="utf-8")
    print(f"  activity-only {sorted(only_activity)}, exposure-only {sorted(only_exposure)}")


if __name__ == "__main__":
    print("figures:")
    figure_area(OUT / "interval-area.svg")
    figure_cancellation(OUT / "cancellation-error.svg")
    figure_air_coefficients(OUT / "air-coefficients.svg")
    figure_air_attenuation(OUT / "air-attenuation.svg")
    figure_dominance(OUT / "dominance-timeline.svg")
    figure_metrics(OUT / "metric-divergence.svg")
    print(f"wrote to {OUT}")
