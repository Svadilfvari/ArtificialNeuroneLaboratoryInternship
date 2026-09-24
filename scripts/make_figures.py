#!/usr/bin/env python3
"""Regenerate every figure in docs/figures from the data in data/.

    pip install -r scripts/requirements.txt
    python scripts/make_figures.py

Each figure is written twice (``*-light`` and ``*-dark``) so the README can
follow GitHub's theme. The script also prints the metrics quoted in the README.

The LIF spike generator below is a line-by-line port of matlab/LIF.m, and the
polynomial fits use the same scaling as matlab/main.m (I_ex in nA, f_spike in
kHz, P_rms in nW), so the numbers match the MATLAB run.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "docs" / "figures"

ORDER = 15  # polynomial order used in matlab/main.m

THEMES = {
    "light": dict(
        surface="#ffffff", ink="#1f2328", muted="#59636e", grid="#e6e8eb", line="#d0d7de",
        model="#2a78d6", cadence="#eb6834", alt="#1baf7a", tint="#eaf2fc", tint2="#fdeee8",
        myelin="#c9d1d9",
    ),
    "dark": dict(
        surface="#0d1117", ink="#e6edf3", muted="#9198a1", grid="#21262d", line="#30363d",
        model="#3987e5", cadence="#d95926", alt="#199e70", tint="#11263f", tint2="#2d1a12",
        myelin="#3d444d",
    ),
}


# --------------------------------------------------------------------- data
def load_operating_points():
    """Cadence sweep of the FS eNeuron: I_ex (A), P_rms (W), f_spike (Hz)."""
    d = np.genfromtxt(DATA / "TB_Ferreira2020_FS_MLneuron_PLS_FoM.csv", delimiter=",", skip_header=2)
    return d[:, 0], d[:, 1], d[:, 2]


def load_transient():
    """Cadence transient at I_ex = 30 pA: time (ms), Vout (mV)."""
    d = np.genfromtxt(DATA / "Vout_LIF_30p.csv", delimiter=",", skip_header=1)
    return d[:, 0] * 1e3, d[:, 1] * 1e3


def polyfit(x, y, n):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # high orders are ill-conditioned, as in MATLAB
        return np.polyfit(x, y, n)


def lif(dc, fe, fsp, vpp):
    """Port of matlab/LIF.m -> (time in ms, membrane potential in mV)."""
    te = 1 / fe
    period = 1 / (3.40 * fsp) * 1e4
    t = np.arange(0, 1 + te / 2, te)
    ton = dc * period
    a1, a2 = 0.01, 100.0
    lam1 = np.log(1 + vpp / a1) / ton
    lam2 = np.log(1 - vpp / a2) / (ton - period)
    b1, b2 = -a1, -a2
    v = np.zeros_like(t)
    x, y, z, a = 0.0, ton, period, 1.0
    for i, ti in enumerate(t):
        if ti >= y and v[i] == 0:
            if a == 1:
                a = ti
            x, y, z = x + a, y + a, z + a
        if x <= ti <= y:
            v[i] = a1 * np.exp(lam1 * (ti - x)) + b1
        elif y < ti <= z:
            v[i] = a2 * np.exp(-lam2 * (ti - x)) + b2
    return t, v


def spike_times(t, v, level):
    idx = np.where((v[:-1] < level) & (v[1:] >= level))[0]
    return t[idx]


# ------------------------------------------------------------------ styling
def style(th):
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 12,
        "figure.facecolor": th["surface"],
        "axes.facecolor": th["surface"],
        "savefig.facecolor": th["surface"],
        "axes.edgecolor": th["line"],
        "axes.labelcolor": th["muted"],
        "axes.titlecolor": th["ink"],
        "axes.titlesize": 13.5,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.titlepad": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": th["grid"],
        "grid.linewidth": 0.8,
        "xtick.color": th["muted"],
        "ytick.color": th["muted"],
        "xtick.labelcolor": th["muted"],
        "ytick.labelcolor": th["muted"],
        "legend.frameon": False,
        "legend.labelcolor": th["ink"],
        "lines.solid_capstyle": "round",
        "xtick.major.pad": 6,
        "ytick.major.pad": 6,
        "lines.solid_joinstyle": "round",
    })


def save(fig, name, theme):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}-{theme}.png", dpi=180, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)


# ------------------------------------------------------------------ figures
def fig_fits(th, theme, iex, prms, fsp):
    x = iex * 1e9
    series = [
        ("Spike frequency  f$_{spike}$", "kHz", fsp * 1e-3),
        ("Power consumption  P$_{rms}$", "nW", prms * 1e9),
    ]
    xs = np.linspace(x.min(), x.max(), 800)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, (title, unit, y) in zip(axes, series):
        y2 = np.polyval(polyfit(x, y, 2), xs)
        y15 = np.polyval(polyfit(x, y, ORDER), xs)
        ax.scatter(x, y, s=22, color=th["cadence"], linewidth=0, zorder=2, label="Cadence operating points")
        ax.plot(xs, y2, color=th["alt"], lw=2, zorder=3, label="MATLAB fit, order 2")
        ax.plot(xs, y15, color=th["model"], lw=2, zorder=4, label=f"MATLAB fit, order {ORDER}")
        ax.set_title(title)
        ax.set_xlabel("Excitation current  I$_{ex}$ (nA)")
        ax.set_ylabel(unit)
        ax.set_xlim(0, 20.5)
        ax.set_xticks([0, 5, 10, 15, 20])
        ax.set_ylim(bottom=0)
        ax.tick_params(length=0)
    axes[0].legend(loc="lower right", handlelength=1.6)
    fig.tight_layout(w_pad=3)
    save(fig, "fits", theme)


def fig_error(th, theme, iex, prms, fsp):
    x = iex * 1e9
    orders = np.arange(1, 21)
    series = [
        ("f$_{spike}$ fit error", fsp * 1e-3),
        ("P$_{rms}$ fit error", prms * 1e9),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
    for ax, (title, y) in zip(axes, series):
        span = y.max() - y.min()
        err = [np.sqrt(np.mean((np.polyval(polyfit(x, y, n), x) - y) ** 2)) / span * 100 for n in orders]
        ax.plot(orders, err, color=th["model"], lw=2, zorder=2)
        for n, col in ((2, th["alt"]), (ORDER, th["model"])):
            e = err[n - 1]
            ax.scatter([n], [e], s=70, color=col, edgecolor=th["surface"], linewidth=2, zorder=3)
            ax.annotate(f"order {n}: {e:.2g} %", (n, e), xytext=(8, 8), textcoords="offset points",
                        color=th["ink"], fontsize=11.5)
        ax.set_yscale("log")
        ax.set_title(title)
        ax.set_xlabel("Polynomial order N")
        ax.set_xticks([1, 5, 10, 15, 20])
        ax.tick_params(length=0)
    axes[0].set_ylabel("Normalised RMSE (% of range)")
    axes[0].set_yticks([0.5, 1, 2, 5, 10, 20])
    axes[0].yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    axes[0].yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    fig.tight_layout(w_pad=3)
    save(fig, "error_vs_order", theme)


def fig_spike_model(th, theme, fit_f):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw=dict(width_ratios=[1, 1.25]))

    # left: anatomy of one spike (I_ex = 30 pA), finer time step so the peak reaches V_pp
    t, v = lif(0.5, 20000, fit_f(0.03), 90)
    m = t <= 0.13
    ax1.plot(t[m], v[m], color=th["model"], lw=2)
    ton = t[np.argmax(v[t < 0.07])]
    vmax = 90.0
    ax1.axhline(vmax, color=th["muted"], lw=0.8, ls=(0, (3, 3)), zorder=0)
    ax1.text(0.001, vmax + 2.5, "V$_{pp}$ = 90 mV", color=th["muted"], va="bottom", ha="left", fontsize=11)
    ax1.annotate("", xy=(0, -9), xytext=(ton, -9),
                 arrowprops=dict(arrowstyle="<->", color=th["muted"], lw=1, shrinkA=0, shrinkB=0))
    ax1.text(ton / 2, -15, "T$_{on}$", color=th["muted"], ha="center", va="top")
    ax1.annotate("integrate:\n$V = A_1 e^{\\lambda_1 t} + B_1$", xy=(0.036, 12), xytext=(0.003, 50),
                 color=th["ink"], fontsize=11,
                 arrowprops=dict(arrowstyle="-", color=th["muted"], lw=0.8))
    ax1.annotate("fire", xy=(ton, vmax), xytext=(ton + 0.008, vmax - 6), color=th["ink"], fontsize=11,
                 va="top")
    ax1.annotate("reset", xy=(ton + 0.0008, 3), xytext=(ton + 0.012, 22), color=th["ink"], fontsize=11,
                 arrowprops=dict(arrowstyle="-", color=th["muted"], lw=0.8))
    ax1.set_xlim(-0.004, 0.14)
    ax1.set_ylim(-26, 100)
    ax1.set_yticks([0, 25, 50, 75, 100])
    ax1.set_title("One spike of the MATLAB model")
    ax1.set_xlabel("Time (ms)")
    ax1.set_ylabel("Membrane potential (mV)")
    ax1.tick_params(length=0)

    # right: rate coding, stronger input -> more spikes
    currents = [0.03, 0.05, 0.1]  # nA
    offset = 110
    for k, i_na in enumerate(currents):
        t, v = lif(0.5, 2000, fit_f(i_na), 90)
        m = t <= 0.5
        base = (len(currents) - 1 - k) * offset
        ax2.plot(t[m], v[m] + base, color=th["model"], lw=1.4)
        n = len(spike_times(t[m], v[m], 45))
        ax2.text(0.505, base + 40, f"{i_na * 1e3:.0f} pA\n{n} spikes", color=th["ink"], fontsize=11,
                 va="center", ha="left")
    ax2.set_yticks([])
    ax2.spines["left"].set_visible(False)
    ax2.grid(axis="y", visible=False)
    ax2.set_xlim(0, 0.5)
    ax2.set_title("Frequency coding: more current, more spikes")
    ax2.set_xlabel("Time (ms)")
    ax2.tick_params(length=0)
    fig.tight_layout(w_pad=3)
    save(fig, "spike_model", theme)


def fig_validation(th, theme, fit_f, stats):
    tc, vc = load_transient()
    t, v = lif(0.5, 2000, fit_f(0.03), 90)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 5.6), gridspec_kw=dict(height_ratios=[1, 1]))
    for ax, lo, hi in ((ax1, 0, 1), (ax2, 0.30, 0.42)):
        mc = (tc >= lo) & (tc <= hi)
        mm = (t >= lo) & (t <= hi)
        ax.plot(tc[mc], vc[mc], color=th["cadence"], lw=1.6, label="Cadence (transistor level)")
        ax.plot(t[mm], v[mm], color=th["model"], lw=1.6, label="MATLAB LIF model")
        ax.set_xlim(lo, hi)
        ax.set_ylim(0, 100)
        ax.set_ylabel("V (mV)")
        ax.tick_params(length=0)
    ax1.set_title("Spike train at I$_{ex}$ = 30 pA")
    ax1.set_title(f"{stats['n_model']} vs {stats['n_cad']} spikes  ·  period {stats['p_err']:+.1f} %  ·  "
                  f"peak {(stats['pk_model'] / stats['pk_cad'] - 1) * 100:+.1f} %",
                  loc="right", fontsize=11.5, fontweight="normal", color=th["muted"])
    ax2.set_title("Zoom 0.30 – 0.42 ms")
    ax2.set_xlabel("Time (ms)")
    ax2.legend(loc="upper left", bbox_to_anchor=(0.17, 1.0), handlelength=1.6)
    fig.tight_layout(h_pad=2)
    save(fig, "validation_30pA", theme)


# ------------------------------------------------------------ SVG diagrams
def sub(base, s):
    """SVG text with a subscript, e.g. sub('I', 'ex')."""
    return f'{base}<tspan font-size="75%" dy="3">{s}</tspan><tspan dy="-3">\u200b</tspan>'


IEX, PRMS, FSPK = sub("I", "ex"), sub("P", "rms"), sub("f", "spike")


def svg_pipeline(th):
    W, H = 1000, 300
    cards = [
        ("1", "Cadence", ["Transistor-level simulation", "of the FS eNeuron circuit"], "chip"),
        ("2", "Dataset", ["309 operating points", f"{IEX} \u2192 {PRMS}, {FSPK}"], "table"),
        ("3", "Polynomial fit", ["polyfit, order N = 15", f"f({IEX}) and P({IEX})"], "curve"),
        ("4", "LIF spike generator", ["exponential charge + reset", "→ spike train V(t)"], "spikes"),
    ]
    cw, ch, gap, x0, y0 = 214, 150, 30, 26, 78
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
        f'font-family="-apple-system, Segoe UI, Helvetica, Arial, sans-serif">',
        f'<rect width="{W}" height="{H}" rx="12" fill="{th["surface"]}"/>',
        '<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
        f'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{th["muted"]}"/></marker></defs>',
    ]
    # brackets: slow part / fast part
    def bracket(xa, xb, label, sub, color):
        y = 44
        out.append(f'<path d="M{xa},{y + 12} V{y} H{xb} V{y + 12}" fill="none" stroke="{color}" stroke-width="1.5"/>')
        out.append(f'<text x="{(xa + xb) / 2}" y="{y - 10}" text-anchor="middle" font-size="14" '
                   f'font-weight="600" fill="{th["ink"]}">{label} '
                   f'<tspan font-weight="400" fill="{th["muted"]}">{sub}</tspan></text>')

    xs = [x0 + k * (cw + gap) for k in range(4)]
    bracket(xs[0], xs[1] + cw, "Cadence", "— run once, ≈ 1 day per neuron", th["cadence"])
    bracket(xs[2], xs[3] + cw, "MATLAB", f"\u2014 any {IEX}, &lt; 1 s", th["model"])

    for (num, title, lines, icon), x in zip(cards, xs):
        stage_col = th["cadence"] if num in "12" else th["model"]
        tint = th["tint2"] if num in "12" else th["tint"]
        out.append(f'<rect x="{x}" y="{y0}" width="{cw}" height="{ch}" rx="10" fill="{tint}" '
                   f'stroke="{th["line"]}"/>')
        out.append(f'<circle cx="{x + 22}" cy="{y0 + 24}" r="12" fill="{stage_col}"/>')
        out.append(f'<text x="{x + 22}" y="{y0 + 29}" text-anchor="middle" font-size="13" '
                   f'font-weight="700" fill="#ffffff">{num}</text>')
        out.append(f'<text x="{x + 42}" y="{y0 + 29}" font-size="15" font-weight="700" '
                   f'fill="{th["ink"]}">{title}</text>')
        for k, line in enumerate(lines):
            out.append(f'<text x="{x + 16}" y="{y0 + 118 + k * 18}" font-size="12.5" '
                       f'fill="{th["muted"]}">{line}</text>')
        # icon area: x+16 .. x+cw-16, y0+44 .. y0+96
        ix, iy, iw, ih = x + 16, y0 + 46, cw - 32, 48
        if icon == "chip":
            cx = ix + iw / 2
            out.append(f'<rect x="{cx - 26}" y="{iy + 6}" width="52" height="36" rx="4" fill="none" '
                       f'stroke="{stage_col}" stroke-width="2"/>')
            for k in range(4):
                px = cx - 18 + k * 12
                out.append(f'<line x1="{px}" y1="{iy}" x2="{px}" y2="{iy + 6}" stroke="{stage_col}" stroke-width="2"/>')
                out.append(f'<line x1="{px}" y1="{iy + 42}" x2="{px}" y2="{iy + 48}" stroke="{stage_col}" stroke-width="2"/>')
            out.append(f'<text x="{cx}" y="{iy + 29}" text-anchor="middle" font-size="11" '
                       f'font-weight="600" fill="{stage_col}">SPICE</text>')
        elif icon == "table":
            cx = ix + iw / 2
            tw, rh = 120, 10
            for r in range(5):
                fill = stage_col if r == 0 else "none"
                op = "0.85" if r == 0 else "1"
                out.append(f'<rect x="{cx - tw / 2}" y="{iy + r * rh - 1}" width="{tw}" height="{rh}" '
                           f'fill="{fill}" fill-opacity="{op}" stroke="{th["line"]}"/>')
            for c in (1, 2):
                xx = cx - tw / 2 + c * tw / 3
                out.append(f'<line x1="{xx}" y1="{iy - 1}" x2="{xx}" y2="{iy + 5 * rh - 1}" stroke="{th["line"]}"/>')
        elif icon == "curve":
            pts = []
            for k in range(41):
                u = k / 40
                pts.append((ix + 30 + u * (iw - 60), iy + ih - 4 - (ih - 10) * (1 - np.exp(-4 * u)) / (1 - np.exp(-4))))
            d = "M" + " L".join(f"{a:.1f},{b:.1f}" for a, b in pts)
            out.append(f'<path d="{d}" fill="none" stroke="{stage_col}" stroke-width="2.5"/>')
            for k in (0, 4, 9, 15, 23, 32, 40):
                a, b = pts[k]
                out.append(f'<circle cx="{a:.1f}" cy="{b + (3 if k % 2 else -3):.1f}" r="3" fill="{th["cadence"]}"/>')
        elif icon == "spikes":
            d = []
            n, px, base, top = 5, ix + 14, iy + ih - 4, iy + 4
            per = (iw - 28) / n
            for k in range(n):
                sx = px + k * per
                seg = [f"{sx + u / 10 * per:.1f},{base - (base - top) * (np.exp(3 * u / 10) - 1) / (np.exp(3) - 1):.1f}"
                       for u in range(11)]
                d.append(("M" if k == 0 else "L") + " L".join(seg) + f" L{sx + per:.1f},{base}")
            out.append(f'<path d="{" ".join(d)}" fill="none" stroke="{stage_col}" stroke-width="2"/>')
    # arrows between cards
    for k in range(3):
        xa, xb = xs[k] + cw + 3, xs[k + 1] - 3
        out.append(f'<line x1="{xa}" y1="{y0 + ch / 2}" x2="{xb}" y2="{y0 + ch / 2}" stroke="{th["muted"]}" '
                   f'stroke-width="1.8" marker-end="url(#ah)"/>')
    # validation loop
    yb = y0 + ch + 34
    xa, xb = xs[3] + cw / 2, xs[0] + cw / 2
    out.append(f'<path d="M{xa},{y0 + ch + 3} V{yb} H{xb} V{y0 + ch + 5}" fill="none" stroke="{th["muted"]}" '
               f'stroke-width="1.5" stroke-dasharray="5 4" marker-end="url(#ah)"/>')
    out.append(f'<rect x="{(xa + xb) / 2 - 215}" y="{yb - 13}" width="430" height="26" rx="13" fill="{th["surface"]}"/>')
    out.append(f'<text x="{(xa + xb) / 2}" y="{yb + 5}" text-anchor="middle" font-size="13" fill="{th["ink"]}">'
               f'Validation: MATLAB spike train vs Cadence transient at {IEX} = 30 pA</text>')
    out.append("</svg>")
    return "\n".join(out)


def svg_neuron_circuit(th):
    W, H = 1000, 400
    ink, muted, line, blue, orange = th["ink"], th["muted"], th["line"], th["model"], th["cadence"]
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
        f'font-family="-apple-system, Segoe UI, Helvetica, Arial, sans-serif">',
        f'<rect width="{W}" height="{H}" rx="12" fill="{th["surface"]}"/>',
        '<defs><marker id="a2" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
        f'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{muted}"/></marker></defs>',
        f'<text x="30" y="40" font-size="16" font-weight="700" fill="{ink}">Biological neuron</text>',
        f'<text x="560" y="40" font-size="16" font-weight="700" fill="{ink}">LIF equivalent circuit</text>',
        f'<line x1="520" y1="60" x2="520" y2="320" stroke="{line}"/>',
    ]

    def badge(x, y, n, col=None):
        col = col or ink
        out.append(f'<circle cx="{x}" cy="{y}" r="11" fill="{col}"/>')
        out.append(f'<text x="{x}" y="{y + 4.5}" text-anchor="middle" font-size="12.5" font-weight="700" '
                   f'fill="{th["surface"]}">{n}</text>')

    def label(x, y, text, anchor="start", size=12.5, col=None):
        out.append(f'<text x="{x}" y="{y}" font-size="{size}" fill="{col or muted}" text-anchor="{anchor}">{text}</text>')

    # ---- neuron (left): dendrites on the left, soma, hillock, myelinated axon, terminals
    sx, sy = 190, 190
    dend = [
        "M158,176 C130,158 108,128 80,108", "M110,136 C100,116 104,96 96,80", "M94,121 C70,122 58,130 46,140",
        "M156,190 C126,190 104,196 70,192",
        "M158,206 C128,222 106,248 78,268", "M106,244 C90,245 70,238 54,232", "M92,256 C90,278 94,294 88,308",
    ]
    for d in dend:
        out.append(f'<path d="{d}" fill="none" stroke="{blue}" stroke-width="4" stroke-linecap="round"/>')
    out.append(f'<path d="M{sx - 36},{sy} C{sx - 38},{sy - 34} {sx + 10},{sy - 44} {sx + 30},{sy - 20} '
               f'C{sx + 42},{sy - 8} {sx + 42},{sy + 12} {sx + 30},{sy + 22} '
               f'C{sx + 8},{sy + 44} {sx - 34},{sy + 34} {sx - 36},{sy} Z" '
               f'fill="{th["tint"]}" stroke="{blue}" stroke-width="3"/>')
    out.append(f'<circle cx="{sx - 4}" cy="{sy}" r="11" fill="{blue}" fill-opacity="0.35"/>')
    for cx in (sx - 18, sx + 4):  # ion channels (pores) on the lower membrane
        cy = sy + 33
        out.append(f'<g transform="translate({cx},{cy})">'
                   f'<rect x="-6" y="-7" width="4.5" height="14" rx="2" fill="{orange}"/>'
                   f'<rect x="1.5" y="-7" width="4.5" height="14" rx="2" fill="{orange}"/></g>')
    ax0 = sx + 38
    out.append(f'<path d="M{sx + 30},{sy - 14} L{ax0 + 24},{sy - 4} L{ax0 + 24},{sy + 4} L{sx + 30},{sy + 16} Z" '
               f'fill="{th["tint"]}" stroke="{blue}" stroke-width="3" stroke-linejoin="round"/>')
    out.append(f'<line x1="{ax0 + 24}" y1="{sy}" x2="425" y2="{sy}" stroke="{blue}" stroke-width="4"/>')
    for k in range(5):
        out.append(f'<rect x="{ax0 + 36 + k * 34}" y="{sy - 9}" width="26" height="18" rx="9" fill="{th["myelin"]}"/>')
    for d in (f"M425,{sy} C440,{sy - 5} 450,{sy - 25} 465,{sy - 35}", f"M425,{sy} L470,{sy}",
              f"M425,{sy} C440,{sy + 5} 450,{sy + 25} 465,{sy + 35}"):
        out.append(f'<path d="{d}" fill="none" stroke="{blue}" stroke-width="3.5" stroke-linecap="round"/>')
    for bx, by in ((467, sy - 36), (473, sy), (467, sy + 36)):
        out.append(f'<circle cx="{bx}" cy="{by}" r="5" fill="{blue}"/>')
    # action potential travelling along the axon
    b0, yb = 330, sy - 28
    sp = (f"M{b0},{yb} L{b0 + 18},{yb} C{b0 + 26},{yb} {b0 + 28},{yb - 42} {b0 + 34},{yb - 42} "
          f"C{b0 + 40},{yb - 42} {b0 + 40},{yb} {b0 + 50},{yb} L{b0 + 66},{yb}")
    out.append(f'<path d="{sp}" fill="none" stroke="{orange}" stroke-width="2.5"/>')
    out.append(f'<line x1="{b0 + 72}" y1="{yb}" x2="{b0 + 100}" y2="{yb}" stroke="{orange}" stroke-width="2" marker-end="url(#a2)"/>')
    label(b0, yb - 54, "action potential")
    # input
    out.append(f'<line x1="22" y1="96" x2="66" y2="104" stroke="{muted}" stroke-width="1.5" marker-end="url(#a2)"/>')
    label(22, 88, "input")
    # badges + labels
    badge(110, 70, 1)
    label(126, 75, "dendrites")
    badge(sx - 4, sy - 66, 2)
    out.append(f'<line x1="{sx - 4}" y1="{sy - 55}" x2="{sx - 4}" y2="{sy - 42}" stroke="{muted}" stroke-width="1.2"/>')
    label(sx + 12, sy - 61, "cell membrane")
    badge(sx - 7, sy + 68, 3, orange)
    out.append(f'<line x1="{sx - 7}" y1="{sy + 57}" x2="{sx - 7}" y2="{sy + 42}" stroke="{muted}" stroke-width="1.2"/>')
    label(sx + 9, sy + 73, "ion channels")
    badge(ax0 + 14, sy - 30, 4)
    label(ax0 + 30, sy - 25, "axon hillock")
    label(330, sy + 32, "axon")

    # ---- circuit (right)
    top, bot = 115, 275
    xl, xr = 575, 800
    out.append(f'<path d="M{xl},{top} H{xr} M{xl},{bot} H{xr}" stroke="{ink}" stroke-width="2" fill="none"/>')
    gx = 690
    out.append(f'<line x1="{gx}" y1="{bot}" x2="{gx}" y2="{bot + 14}" stroke="{ink}" stroke-width="2"/>')
    for k, w in enumerate((22, 14, 6)):
        out.append(f'<line x1="{gx - w}" y1="{bot + 14 + k * 5}" x2="{gx + w}" y2="{bot + 14 + k * 5}" stroke="{ink}" stroke-width="2"/>')
    ym = (top + bot) / 2
    # 1: current source
    cx = xl
    out.append(f'<path d="M{cx},{top} V{ym - 22} M{cx},{ym + 22} V{bot}" stroke="{ink}" stroke-width="2" fill="none"/>')
    out.append(f'<circle cx="{cx}" cy="{ym}" r="22" fill="{th["surface"]}" stroke="{ink}" stroke-width="2"/>')
    out.append(f'<line x1="{cx}" y1="{ym + 12}" x2="{cx}" y2="{ym - 12}" stroke="{blue}" stroke-width="2.5" marker-end="url(#a2)"/>')
    out.append(f'<text x="{cx - 30}" y="{ym + 5}" font-size="14" text-anchor="end" fill="{ink}" font-style="italic">'
               f'I<tspan font-size="10" dy="4">ex</tspan></text>')
    badge(cx - 34, ym - 34, 1)
    # 3: leak resistor
    rx = 645
    zz = " ".join(f"L{rx + (8 if k % 2 == 0 else -8)},{ym - 25 + k * 10}" for k in range(6))
    out.append(f'<path d="M{rx},{top} V{ym - 30} {zz} L{rx},{ym + 35} V{bot}" fill="none" stroke="{ink}" '
               f'stroke-width="2" stroke-linejoin="round"/>')
    out.append(f'<text x="{rx + 16}" y="{ym + 5}" font-size="14" fill="{ink}" font-style="italic">R</text>')
    badge(rx + 24, ym - 40, 3, orange)
    # 2: membrane capacitor
    cpx = 725
    out.append(f'<path d="M{cpx},{top} V{ym - 7} M{cpx - 20},{ym - 7} H{cpx + 20} M{cpx - 20},{ym + 7} H{cpx + 20} '
               f'M{cpx},{ym + 7} V{bot}" fill="none" stroke="{ink}" stroke-width="2"/>')
    out.append(f'<text x="{cpx + 26}" y="{ym + 5}" font-size="14" fill="{ink}" font-style="italic">C</text>')
    badge(cpx + 32, ym - 30, 2)
    out.append(f'<text x="{cpx - 28}" y="{top - 10}" font-size="14" fill="{ink}" font-style="italic">V(t)</text>')
    # reset switch
    swx = xr
    out.append(f'<path d="M{swx},{top} V{ym - 20} M{swx},{bot} V{ym + 22}" fill="none" stroke="{ink}" stroke-width="2"/>')
    out.append(f'<circle cx="{swx}" cy="{ym - 20}" r="3" fill="{ink}"/><circle cx="{swx}" cy="{ym + 22}" r="3" fill="{ink}"/>')
    out.append(f'<line x1="{swx}" y1="{ym + 22}" x2="{swx - 18}" y2="{ym - 14}" stroke="{ink}" stroke-width="2"/>')
    label(swx + 8, ym + 5, "reset")
    # 4: threshold comparator
    bx, by, bw, bh = 850, 93, 118, 44
    out.append(f'<line x1="{xr}" y1="{top}" x2="{bx}" y2="{top}" stroke="{ink}" stroke-width="2"/>')
    out.append(f'<rect x="{bx}" y="{by}" width="{bw}" height="{bh}" rx="8" fill="{th["tint"]}" stroke="{blue}" stroke-width="2"/>')
    out.append(f'<text x="{bx + bw / 2}" y="{by + 27}" text-anchor="middle" font-size="14" fill="{ink}">'
               f'<tspan font-style="italic">V</tspan> ≥ <tspan font-style="italic">V</tspan>'
               f'<tspan font-size="10" dy="4">th</tspan><tspan dy="-4"> ?</tspan></text>')
    badge(bx + bw - 2, by - 2, 4)
    ox, oy = bx + 20, 190
    out.append(f'<line x1="{bx + bw / 2}" y1="{by + bh}" x2="{bx + bw / 2}" y2="{oy - 30}" stroke="{ink}" '
               f'stroke-width="2" marker-end="url(#a2)"/>')
    spk = (f"M{ox},{oy + 22} L{ox + 22},{oy + 22} C{ox + 32},{oy + 22} {ox + 36},{oy - 18} {ox + 40},{oy - 18} "
           f"L{ox + 40},{oy + 22} L{ox + 78},{oy + 22}")
    out.append(f'<path d="{spk}" fill="none" stroke="{orange}" stroke-width="2.5"/>')
    label(ox + 39, oy + 44, "spike out", anchor="middle")
    out.append(f'<path d="M{bx + 8},{by + bh} V{250} H{swx + 10}" fill="none" stroke="{muted}" stroke-width="1.5" '
               f'stroke-dasharray="4 4" marker-end="url(#a2)"/>')
    label((bx + swx) / 2 + 14, 244, "fire", anchor="middle", size=11.5)

    # ---- legend row: badge -> biological part = circuit element
    items = [
        ("1", ink, "input current → source I<tspan font-size='9' dy='3'>ex</tspan>"),
        ("2", ink, "cell membrane → capacitor C"),
        ("3", orange, "ion channels → leak resistor R"),
        ("4", ink, "axon hillock → threshold, fire, reset"),
    ]
    out.append(f'<line x1="30" y1="{H - 62}" x2="{W - 30}" y2="{H - 62}" stroke="{line}"/>')
    colw = (W - 60) / 4
    for k, (n, col, text) in enumerate(items):
        x = 30 + k * colw + 12
        badge(x, H - 32, n, col)
        label(x + 18, H - 27.5, text, col=ink)
    out.append("</svg>")
    return "\n".join(out)


# --------------------------------------------------------------------- main
def main():
    iex, prms, fsp = load_operating_points()
    x = iex * 1e9
    p_f = polyfit(x, fsp * 1e-3, ORDER)
    p_p = polyfit(x, prms * 1e9, ORDER)

    def fit_f(i_na):  # Hz, like polyval(P, Iex)*1e3 in main.m
        return np.polyval(p_f, i_na) * 1e3

    # ---- metrics quoted in the README
    f_hat = np.polyval(p_f, x)
    p_hat = np.polyval(p_p, x)
    f, p = fsp * 1e-3, prms * 1e9
    nrmse = lambda y, yh: np.sqrt(np.mean((yh - y) ** 2)) / (y.max() - y.min()) * 100  # noqa: E731
    medrel = lambda y, yh: np.median(np.abs((yh - y) / y)) * 100  # noqa: E731
    tc, vc = load_transient()
    t, v = lif(0.5, 2000, fit_f(0.03), 90)
    sm, sc = spike_times(t, v, 45), spike_times(tc, vc, vc.max() / 2)
    stats = dict(
        n_model=len(sm), n_cad=len(sc),
        p_model=np.diff(sm).mean() * 1e3, p_cad=np.diff(sc).mean() * 1e3,
        pk_model=v.max(), pk_cad=vc.max(),
    )
    stats["p_err"] = (stats["p_model"] / stats["p_cad"] - 1) * 100
    print(f"operating points      : {len(x)}  ({x.min() * 1e3:.1f} pA - {x.max():.1f} nA)")
    print(f"order {ORDER} f_spike fit : NRMSE {nrmse(f, f_hat):.2f} %, median |rel err| {medrel(f, f_hat):.2f} %")
    print(f"order {ORDER} P_rms fit   : NRMSE {nrmse(p, p_hat):.2f} %, median |rel err| {medrel(p, p_hat):.2f} %")
    print(f"order 2 f_spike / P_rms NRMSE: {nrmse(f, np.polyval(polyfit(x, f, 2), x)):.1f} % / "
          f"{nrmse(p, np.polyval(polyfit(x, p, 2), x)):.1f} %")
    print("validation at 30 pA   : " + ", ".join(f"{k}={val:.3g}" for k, val in stats.items()))

    for theme, th in THEMES.items():
        style(th)
        fig_fits(th, theme, iex, prms, fsp)
        fig_error(th, theme, iex, prms, fsp)
        fig_spike_model(th, theme, fit_f)
        fig_validation(th, theme, fit_f, stats)
        (OUT / f"pipeline-{theme}.svg").write_text(svg_pipeline(th), encoding="utf-8")
        (OUT / f"neuron_to_circuit-{theme}.svg").write_text(svg_neuron_circuit(th), encoding="utf-8")
    print(f"figures written to {OUT.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()
