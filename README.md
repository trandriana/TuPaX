# Gray–Scott IMEX–FiPy

This code implements a semi-implicit (IMEX) finite-volume scheme for the Gray–Scott reaction–diffusion model using FiPy.

---

## Installation

### Using conda (recommended)
```bash
conda create -n grayscott python=3.11 -y
conda activate grayscott
conda install -c conda-forge fipy matplotlib scipy ffmpeg -y
pip install -e .

## Quick start

Run a short simulation and save snapshots:
```bash
gray-scott-compute --nx 128 --T 500 --steps 4000 --n-image-snapshots 5


src/gray_scott/
  compute.py      # runs the simulation (CLI: gray-scott-compute)
  animate.py      # builds MP4/GIF from snapshots (CLI: gray-scott-animate)
  __init__.py     # exposes API functions

examples/
  run_default.sh  # simple demo to compute
  make_video.sh   # simple demo to animate

output_images/      # PNG frames (created at runtime)
output_solution/    # NPZ/MAT snapshots + optional MP4 (created at runtime)

---


