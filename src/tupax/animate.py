'''
Gray-Scott animation from *.npz.
Author: Tsiry Avisoa Randrianasolo
'''
# ---------------------
# Standard library imports
# ---------------------
import os
import shutil

# import scipy.io

# ---------------------
# Matplotlib setup
# ---------------------
import matplotlib

matplotlib.use("Agg")  # Use a non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

# ---------------------
# Other imports
# ---------------------
import numpy as np
import argparse
import shutil


# ---------------------
# Main simulation function
# ---------------------
def make_video(
    output_animations_dir: str = "output_animations",
    npz_path: str = "output_solution/reference_solution.npz",
    out_path: str = "gray_scott.mp4",
    fps: int = 15,
    dpi: int = 200,
    cmap: str = "turbo",
    use_fixed_limits: bool = True,
    vmin_fixed: float = 0.0,
    vmax_fixed: float = 1.0,
    which: str = "both",
):
    """Create a video (MP4 or GIF) from snapshots.npz produced by the simulator.

    Parameters
    ----------
    npz_path : str
        Path to NPZ file containing arrays t, u, v and metadata.
    out_path : str
        Output filename (.mp4 recommended; .gif if you prefer).
    which : {'both', 'u', 'v'}
        Whether to animate both fields side-by-side, or only one.
    fps : int
        Frames per second in the output video.
    dpi : int
        Dots per inch (resolution) of the output video.
    cmap : str
        Colormap to use (see matplotlib colormaps).
    use_fixed_limits : bool
        Whether to use fixed color limits (vmin_fixed, vmax_fixed) or dynamic limits.
    vmin_fixed : float
        Fixed minimum value for color limits (only if use_fixed_limits is True).
    vmax_fixed : float
        Fixed maximum value for color limits (only if use_fixed_limits is True).
    """

    print("Running...", end="\r")

    # ---------------------
    # Output directories
    # ---------------------
    os.makedirs(output_animations_dir, exist_ok=True)

    # ---------------------
    # Load data
    # ---------------------
    data = np.load(npz_path)
    t = data["t"]
    u = data["u"]
    v = data["v"]
    NT = u.shape[0]

    # ---------------------
    # Set up figure and animation
    # ---------------------
    if use_fixed_limits:  # use fixed vmin, vmax
        vmin, vmax = vmin_fixed, vmax_fixed
    else:  # use dynamic vmin, vmax
        vmin = float(np.nanmin([u.min(), v.min()]))
        vmax = float(np.nanmax([u.max(), v.max()]))

    if which == "both":  # animate both fields side-by-side

        fig, axes = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)

        ims = [
            axes[0].imshow(u[0], origin="lower", cmap=cmap, vmin=vmin, vmax=vmax),
            axes[1].imshow(v[0], origin="lower", cmap=cmap, vmin=vmin, vmax=vmax),
        ]

        axes[0].set_title(f"u at t = {t[0]:.2f}")
        axes[1].set_title(f"v at t = {t[0]:.2f}")

        for ax in axes:
            ax.axis("off")

        cbar = fig.colorbar(
            ims[0], ax=axes[0], orientation="vertical", fraction=0.03, pad=0.04
        )
        cbar.set_label("Concentration")

        cbar = fig.colorbar(
            ims[1], ax=axes[1], orientation="vertical", fraction=0.03, pad=0.04
        )
        cbar.set_label("Concentration")
    else:  # animate only one field (u or v)
        field = u if which == "u" else v
        fig, ax = plt.subplots(1, 1, figsize=(7, 5), constrained_layout=True)
        ims = [ax.imshow(field[0], origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)]
        ax.set_title(f"{which} ({t[0]:.2f} tu)")
        ax.axis("off")
        axes = [ax]

        cbar = fig.colorbar(
            ims[0], ax=axes[0], orientation="vertical", fraction=0.03, pad=0.04
        )
        cbar.set_label("Concentration")

    def update(i):  # update function for animation
        if which == "both":
            ims[0].set_data(u[i])
            ims[1].set_data(v[i])
            ims[0].set_clim(vmin, vmax)
            ims[1].set_clim(vmin, vmax)
            axes[0].set_title(f"Species u ({t[i]:.2f} tu)")
            axes[1].set_title(f"Species v ({t[i]:.2f} tu)")
            return ims
        else:
            ims[0].set_data((u if which == "u" else v)[i])
            ims[0].set_clim(vmin, vmax)
            axes[0].set_title(f"Species {which} ({t[i]:.2f} tu)")
            return ims

    # Create animation object
    anim = FuncAnimation(fig, update, frames=NT, interval=1000 / fps, blit=False)

    # Save to file (MP4)
    ffmpeg_bin = shutil.which("ffmpeg")

    if ffmpeg_bin:
        plt.rcParams['animation.ffmpeg_path'] = ffmpeg_bin
    else:
        plt.rcParams['animation.ffmpeg_path'] = r"C:\ffmpeg\bin\ffmpeg.exe"

    writer = FFMpegWriter(fps=fps, bitrate=2000)
    fname = os.path.join(
        output_animations_dir,
        out_path if out_path.endswith(".mp4") else "gray_scott.mp4",
    )
    # Check for name conflicts and add _00, _01, ... if needed
    if os.path.exists(fname):
        root, ext = os.path.splitext(fname)
        i = 0
        while True:
            candidate = f"{root}_{i:02d}{ext}"
            if not os.path.exists(candidate):
                fname = candidate
                break
            i += 1
    anim.save(fname, writer=writer, dpi=dpi)

    print("\nDone.")


def main():  # command-line interface

    p = argparse.ArgumentParser(description="Create a video from reference_solution.npz")
    p.add_argument("--output-animations-dir", default="output_animations")
    p.add_argument("--npz-path", nargs="?", default="output_solution/reference_solution.npz")
    p.add_argument("--out", default="gray_scott.mp4")
    p.add_argument("--fps", type=int, default=15)
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--cmap", default="jet")
    p.add_argument("--use-fixed-limits", action="store_true")
    p.add_argument("--vmin", type=float, default=0.0)
    p.add_argument("--vmax", type=float, default=1.0)
    p.add_argument("--which", choices=["both", "u", "v"], default="both")
    args = p.parse_args()

    # Call make_video function
    make_video(
        npz_path=args.npz_path,
        out_path=args.out,
        fps=args.fps,
        dpi=args.dpi,
        cmap=args.cmap,
        use_fixed_limits=args.use_fixed_limits,
        vmin_fixed=args.vmin,
        vmax_fixed=args.vmax,
        which=args.which,
    )

    

# ---------------------
# Entry point for script
# ---------------------
if __name__ == "__main__":

    main()
