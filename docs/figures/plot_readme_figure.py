#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
"""Figure of the README: six frames of the three kinks computed by run_graph.m on the graph.

Run sgg_readme_data.m first (MATLAB, or GNU Octave), then, from this folder,
    python3 plot_readme_figure.py
which writes sgg-kinks.png. The style is the stamped copy manuscript_style.py beside it.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import TwoSlopeNorm
from PIL import Image
from scipy.io import loadmat

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import manuscript_style as hs  # noqa: E402

TIMES = (0.0, 7.0, 10.0, 16.0, 25.0, 33.0)    # instants shown, among the stored frames
LOW, HIGH = -4.2 * np.pi, 2.1 * np.pi         # ends of the colour scale, white at u = 0

data = loadmat(HERE / "sgg_frames.mat")
frames = data["frames"]                       # frames[i, j, k]: node i of edge j at time tkeep[k]
stored = data["tkeep"].ravel()
X = data["X"]                                 # X[i, :, j]: position of node i of edge j in the plane
V = data["V"]                                 # V[:, v]: position of vertex v
chosen = [int(np.argmin(np.abs(stored - t))) for t in TIMES]

hs.configure(family="libertine")
plt.rcParams["savefig.dpi"] = 265
plt.rcParams["figure.facecolor"] = hs.PAPER
plt.rcParams["savefig.facecolor"] = hs.PAPER

norm = TwoSlopeNorm(vcenter=0.0, vmin=LOW, vmax=HIGH)
cmap = matplotlib.colormaps[hs.SIGNED_COLORMAP]
fig, axes = plt.subplots(2, 3, figsize=(hs.WIDTH, 4.35))
fig.get_layout_engine().set(wspace=0.04, hspace=0.06)
for n, (axis, k) in enumerate(zip(axes.flat, chosen)):
    for j in range(X.shape[2]):
        points = X[:, :, j]
        segments = np.stack([points[:-1], points[1:]], axis=1)
        values = 0.5 * (frames[:-1, j, k] + frames[1:, j, k])
        axis.add_collection(LineCollection(segments, colors=hs.SLATE, linewidths=3.7,
                                           capstyle="round", zorder=1))
        edge = LineCollection(segments, cmap=cmap, norm=norm, linewidths=2.9,
                              capstyle="round", zorder=2)
        edge.set_array(values)
        axis.add_collection(edge)
    axis.scatter(V[0], V[1], s=14, color=hs.INK, zorder=3)
    axis.set_xlim(-12.6, 12.6)
    axis.set_ylim(-12.6, 12.6)
    axis.set_aspect("equal")
    axis.axis("off")
    hs.panel(axis, rf"({'abcdef'[n]}) $t = {stored[k]:g}$")
    if n == 0:
        offsets = {0: (0.0, -1.9), 1: (0.0, 1.55), 2: (1.6, 1.15), 3: (-1.6, 1.15)}
        for v, (dx, dy) in offsets.items():
            axis.text(V[0, v] + dx, V[1, v] + dy, rf"$v_{v + 1}$", ha="center", va="center",
                      fontsize=plt.rcParams["legend.fontsize"], color=hs.INK)

bar = hs.colorbar(fig, edge, list(axes.flat), r"solution $u$", shrink=0.62)
bar.set_ticks([-4 * np.pi, -3 * np.pi, -2 * np.pi, -np.pi, 0.0, np.pi, 2 * np.pi],
              labels=[r"$-4\pi$", r"$-3\pi$", r"$-2\pi$", r"$-\pi$", r"$0$", r"$\pi$", r"$2\pi$"])
out = HERE / "sgg-kinks.png"
hs.save(fig, out)
# The background is opaque: without its alpha channel, and compressed a second time without
# loss, the file is smaller.
Image.open(out).convert("RGB").save(out, optimize=True)
print("written", out, "frames at t =", [round(float(stored[k]), 2) for k in chosen])
