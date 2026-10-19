"""Shared design system for the v0.5 main figures (JACC: BTS).
Load in every figure script with exec(open(path).read()); then call apply_jbts_style().
Source: figure_storyboard_v05.json -> design_system. Do not redefine these per figure.
"""
import matplotlib as mpl
from matplotlib import font_manager

MM = 1 / 25.4
WIDTH_MM = 183
WIDTH_IN = WIDTH_MM * MM

SIZES = dict(panel_letter=14, panel_title=10, axis_label=10, tick=9,
             in_panel_stats=10, legend=9, min_gene_label=8)

COHORT = {
    "SCP3342": "#8E3B7E",      # HFpEF snRNA-seq (circle markers)
    "Hahn": "#CDA9F2",         # Hahn biopsy, same HFpEF program (square, bulk)
    "GSE290367": "#3A9A8C",    # DCM atlas (triangle family)
    "MAGNet": "#9A7034",       # DCM bulk (square)
    "GTEx": "#3D5A80",         # ventilator-managed deaths
}
GSE_STUDY = {"Chaffin": "#1F6F66", "Koenig": "#3A9A8C", "Reichart": "#8FD1C6"}
STUDY_MARKER = {"Chaffin": "o", "Koenig": "s", "Reichart": "^"}
COHORT_MARKER = {"SCP3342": "o", "Hahn": "s", "GSE290367": "^", "MAGNet": "s", "GTEx": "D"}
CONTROL = "#A6A6A6"            # donor / non-failing controls, every cohort; GTEx other deaths

GC_CLASS = {"induced": "#E0823D", "repressed": "#2D6A9F", "insensitive": "#BDBDBD",
            "ambiguous": "#FFFFFF", "not_measured": "#EDEDED"}
AMBIGUOUS_EDGE = "#7F7F7F"
DIVERGING = ["#2D6A9F", "#F7F7F7", "#E0823D"]   # GC classes / in vitro logFC heatmaps
DISEASE_CMAP = "RdBu_r"                          # disease-effect heatmaps; never red-green
TEXT = "#222222"
GRID = "#E6E6E6"
POSTHOC_TAG = dict(color="#7F7F7F", style="italic")

# Encodings (storyboard): filled = base model; open = + GC score; half-filled = + GC score + macrophage content.
# Fig 4 shapes: triangle = +circadian, square = +cAMP, diamond = both.
# Score versions: circle = 3-gene GC score; small grey diamond = 8-gene GC response score.
# Significance: print exact values; never encode by fill; no stars.
# Test badges: black = passed; outlined x = failed; half = not met; dashed = not evaluated.


def diverging_cmap(name="gc_div"):
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(name, DIVERGING)


def apply_jbts_style():
    fams = {f.name for f in font_manager.fontManager.ttflist}
    fam = "Arial" if "Arial" in fams else ("Liberation Sans" if "Liberation Sans" in fams else "DejaVu Sans")
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": [fam, "Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": SIZES["tick"], "axes.titlesize": SIZES["panel_title"], "axes.labelsize": SIZES["axis_label"],
        "xtick.labelsize": SIZES["tick"], "ytick.labelsize": SIZES["tick"], "legend.fontsize": SIZES["legend"],
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out", "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "legend.frameon": False,
        "axes.edgecolor": TEXT, "axes.labelcolor": TEXT, "xtick.color": TEXT, "ytick.color": TEXT, "text.color": TEXT,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.dpi": 600, "figure.dpi": 150, "lines.linewidth": 1.0, "patch.linewidth": 0.6,
    })
    return fam


def panel_label(ax, letter, x=-0.12, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=SIZES["panel_letter"], fontweight="bold",
            va="bottom", ha="left", color=TEXT)


def save_figure(fig, stem):
    """Writes stem.pdf (vector), stem.tif (600 dpi RGB, LZW) and stem.png (300 dpi preview)."""
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{stem}.tif", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
