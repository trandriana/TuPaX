# Gray–Scott IMEX–FiPy

This code implements a semi-implicit (IMEX) finite-volume scheme for the Gray–Scott reaction–diffusion model using FiPy.

---

## Project structure
gray_scott_to_ZENODO/
├── src/
│   └── gray_scott/
│       ├── __init__.py
│       ├── compute.py  # runs the simulation (CLI: gray-scott-compute)
│       ├── animate.py  # builds MP4 from snapshots (CLI: gray-scott-animate)
│       └── (maybe utils.py)
│
├── examples/
│   ├── run_default.sh
│   └── make_video.sh
│
├── output_images/        # PNG frames created when you run compute
├── output_solution/      # NPZ/MAT snapshots created when you run compute
│
├── README.md
├── LICENSE
├── pyproject.toml
├── CITATION.cff
└── .gitignore

## Installation

### Using conda (recommended)

```bash
conda create -n grayscott python=3.11 -y
conda activate grayscott
conda install -c conda-forge fipy matplotlib scipy ffmpeg -y
pip install -e .
```

## Quick start

### Run a short simulation and save snapshots:
```bash
gray-scott-compute [options]
```
Common options:

--nx INT grid resolution per side (default: 256)

--L FLOAT domain size (default: 1.0)

--du FLOAT diffusion for u (default: 1.6e-5)

--dv FLOAT diffusion for v (default: 0.8e-5)

--F FLOAT feed rate (default: 0.025)

--k FLOAT kill rate (default: 0.06)

--T FLOAT final time (default: 2000.0)

--steps INT number of time steps (default: 16000)

--n-image-snapshots INT how many PNG frames to save (default: 5)

--output-images-dir PATH (default: output_images)

--output-solution-dir PATH (default: output_solution)


### Make a video from NPZ snapshots
```
gray-scott-animate [options]
```
Common options:

--which {both,u,v} animate both fields or just one

--use-fixed-limits keep a constant color scale (recommended)

--vmin, --vmax color scale (e.g., 0.0 and 1.0)

--cmap jet colormap (default: jet)





