"""Supplementary Figure S1: colorimetric LAMP (proof of principle with laboratory strains), uidA | papGII, 0-35 minutes.

Each photo (two rows of eight tubes: row E = uidA primers, row P = papGII primers; 8 samples, see STRAINS) is cut into the two rows. The tube rows are
located automatically (centres of the coloured liquid), each row is cut to the same window (8 tubes + margin) and scaled to the same size; the cut-outs
are stored in data/colorimetric/ (JPEG) so that the figure can be rebuilt without the original photos.

Usage: python figS1_colorimetric.py [--photos DIR]   (DIR with C_<min>_min.HEIC; only needed to create the cut-outs; macOS: HEIC converted with sips)
Output: data/colorimetric/crop_<min>_<E|P>.jpg and results/figures/figS1.{pdf,png,svg}
"""
import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy import ndimage as ndi

from lamp_common import DATA, RESULTS

CROPS = DATA / "colorimetric"
OUT = RESULTS / "figures"
TIMES = [0, 10, 15, 20, 25, 30, 35]
CROP_W, CROP_H = 1400, 500                      # px of each cut-out (aspect of 8 tubes + margin)
STRAINS = [("$\\it{E. coli}$ 700223-20 ($\\it{papGII}$)"), ("$\\it{E. coli}$ ATCC-41-3 ($\\it{papGII}$)"), ("$\\it{E. coli}$ AVS0171 ($\\it{papGII}$)"),
           ("$\\it{E. coli}$ EC6804 ($\\it{papGIII}$)"), ("$\\it{E. coli}$ EC6792 ($\\it{papGIII}$)"), ("$\\it{E. coli}$ 700636-20 (no $\\it{papG}$)"),
           ("$\\it{K. pneumoniae}$ (no $\\it{papG}$)"), "NTC"]
mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
                     "pdf.fonttype": 42, "svg.fonttype": "none", "svg.hashsalt": "lamp-figS1"})
LETTER = {"fontsize": 16, "fontweight": "bold", "family": "DejaVu Sans"}


def tube_rows(img):
    """Centres (x, y) of the liquid of the 8 tubes in the upper (E) and lower (P) row of a photo."""
    hsv = np.asarray(img.convert("HSV")).astype(float) / 255
    for smin in (0.2, 0.12, 0.09):                                     # pale (early) tubes need a lower saturation threshold
        m = (hsv[..., 1] > smin) & (hsv[..., 2] > 0.4)
        m = ndi.binary_closing(ndi.binary_opening(m, structure=np.ones((7, 7))), structure=np.ones((15, 15)))
        lab, _ = ndi.label(m)
        pts = []
        for i, sl in enumerate(ndi.find_objects(lab), 1):
            a = int((lab[sl] == i).sum())
            if 4000 < a < 150000:
                cy, cx = ndi.center_of_mass(lab == i); pts.append((cx, cy))
        pts = np.array(pts)
        if len(pts) >= 12:
            break
    y_split = np.median(pts[:, 1])
    rows = []
    for sel in (pts[:, 1] < y_split, pts[:, 1] >= y_split):
        r = pts[sel][np.argsort(pts[sel][:, 0])]
        p0 = np.median(np.diff(r[:, 0]))
        idx = np.round((r[:, 0] - r[0, 0]) / p0).astype(int)
        pitch, x0 = np.polyfit(idx, r[:, 0], 1)                         # x = x0 + i * pitch
        y = np.polyval(np.polyfit(idx, r[:, 1], 1), [0, 7]).mean()
        rows.append((x0, pitch, y))
    return rows


def make_crops(photos):
    CROPS.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as d:
        for t in TIMES:
            src = Path(photos) / f"C_{t}_min.HEIC"
            jpg = Path(d) / f"C_{t}.jpg"
            subprocess.run(["sips", "-s", "format", "jpeg", "-s", "formatOptions", "95", str(src), "--out", str(jpg)], capture_output=True, check=True)
            img = Image.open(jpg).convert("RGB")
            for name, (x0, pitch, y) in zip("EP", tube_rows(img)):
                box = (x0 - 0.95 * pitch, y - 2.55 * pitch, x0 + 7 * pitch + 0.95 * pitch, y + 0.55 * pitch)
                box = (box[0], box[1], box[2], box[1] + (box[2] - box[0]) * CROP_H / CROP_W)
                img.crop(tuple(int(round(v)) for v in box)).resize((CROP_W, CROP_H), Image.LANCZOS).save(CROPS / f"crop_{t}_{name}.jpg", quality=92)
            print("cropped", t)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--photos"); a = ap.parse_args()
    if a.photos: make_crops(a.photos)
    OUT.mkdir(parents=True, exist_ok=True)
    mm = 1 / 25.4
    W, left, col_w, gap = 192.0, 22.0, 78.0, 10.0
    row_h = col_w * CROP_H / CROP_W
    top, head, leg = 10.0, 12.0, 34.0
    H = top + head + len(TIMES) * (row_h + 2.0) + leg
    fig = plt.figure(figsize=(W * mm, H * mm))
    for c, (name, gene, letter) in enumerate((("E", "$\\it{uidA}$ ($\\it{E. coli}$)", "(A)"), ("P", "$\\it{papGII}$", "(B)"))):
        x0 = left + c * (col_w + gap)
        fig.text((x0 - (left - 1 if c == 0 else 4)) / W, 1 - 1.5 / H, letter, va="top", **LETTER)
        fig.text((x0 + col_w / 2) / W, 1 - 6 / H, gene, ha="center", va="center", fontsize=11)
        for r, t in enumerate(TIMES):
            y = top + head + r * (row_h + 2.0)
            ax = fig.add_axes([x0 / W, 1 - (y + row_h) / H, col_w / W, row_h / H])
            ax.imshow(Image.open(CROPS / f"crop_{t}_{name}.jpg")); ax.axis("off")
            if c == 0:
                fig.text((left - 2) / W, 1 - (y + row_h / 2) / H, f"{t} min", ha="right", va="center", fontsize=10, fontweight="bold")
        # tube numbers above the first row (tube centres from the 0-min cut-out geometry: 0.95 pitch margin on each side)
        for k in range(8):
            xx = x0 + col_w * (0.95 + k) / (7 + 1.9)
            fig.text(xx / W, 1 - (top + head - 1.5) / H, str(k + 1), ha="center", va="bottom", fontsize=9)
    # strain legend
    ly = top + head + len(TIMES) * (row_h + 2.0) + 4
    for k, s in enumerate(STRAINS):
        col, row = divmod(k, 4)
        fig.text((left + col * 80) / W, 1 - (ly + row * 6.0) / H, f"{k + 1}.  {s}", ha="left", va="top", fontsize=9.5)
    for ext in ("pdf", "png", "svg"):
        meta = {"pdf": {"CreationDate": None, "Creator": "figS1_colorimetric.py"}, "svg": {"Date": None}}.get(ext)
        fig.savefig(OUT / f"figS1.{ext}", dpi=600 if ext == "png" else None, metadata=meta)


if __name__ == "__main__":
    main()
