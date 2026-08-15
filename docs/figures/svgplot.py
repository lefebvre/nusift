"""SVG scaffolding shared by the figure scripts.

Both docs/figures/make_figures.py and validation/make_report.py draw with this, so a
methodology figure and a validation figure cannot end up in different visual languages.

Everything here is deterministic: no font metrics are measured, no timestamps are written,
and every coordinate is emitted through a fixed format. Two runs on two machines produce
byte-identical SVGs, which is what lets the generated figures be committed and diff-checked.
"""

from __future__ import annotations

import math

# --------------------------------------------------------------------------------------
# Every figure paints an explicit background so it stays legible whether or not the viewer
# honours prefers-color-scheme; the media query then adapts it to a dark page rather than
# leaving a light card sitting in one.
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


def esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def svg(width, height, body, title):
    # The accessible label is escaped like any other text. An unescaped ampersand here -- an
    # author's name, most likely -- makes the whole file fail to parse as XML, and a browser
    # renders an error page rather than the figure.
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'width="{width}" height="{height}" role="img" aria-label="{esc(title)}">\n'
        f"<style>{STYLE}</style>\n"
        f'<rect class="bg" x="0.5" y="0.5" width="{width - 1}" height="{height - 1}" rx="6"/>\n'
        f"{body}\n</svg>\n"
    )


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


def circle(x, y, r, cls):
    return f'<circle class="{cls}" cx="{x:.1f}" cy="{y:.1f}" r="{r}"/>'


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
