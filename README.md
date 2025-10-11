[![DOI](https://zenodo.org/badge/1074261727.svg)](https://doi.org/10.5281/zenodo.17328875)

# Gray-Scott IMEX-FiPy

This code implements a semi-implicit (IMEX) finite-volume scheme for the Gray–Scott reaction–diffusion model using FiPy.

We solve the coupled reaction-diffusion equations

$$\partial_t u = d_u \Delta u - u v^2 + F (1 - u),$$ $$\partial_t v = d_v \Delta v + u v^2 - (F + k) v,$$

with Neumann boundary conditions (the defaul on FiPy), where:

| Symbol | Meaning |
|:-------:|:--------|
| $u, v$ | species concentrations |
| $d_u, d_v$ | diffusion coefficients |
| $F, k$ | feed and kill rates |


---

## Project structure
```
gray_scott_IMEX/
├── src/gray_scott/
│   ├── __init__.py
│   ├── compute.py        # CLI: gray-scott-compute (run simulation)
│   └── animate.py        # CLI: gray-scott-animate (build MP4 from NPZ)
│
├── examples/
│   ├── run_default.sh
│   └── make_video.sh
│
├── output_images/        # PNG frames (created at runtime)
├── output_solution/      # NPZ/MAT snapshots, MP4 (created at runtime)
│
├── README.md
├── LICENSE
├── pyproject.toml
├── CITATION.cff
└── .gitignore
```


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
gray-scott-compute [OPTIONS]
```
| Option                       | Meaning                       | Default           |
| ---------------------------- | ----------------------------- | ----------------- |
| `--nx INT`                   | Grid resolution per side      | `256`             |
| `--L FLOAT`                  | Domain size                   | `1.0`             |
| `--du FLOAT`                 | Diffusion coefficient for *u* | `1.6e-5`          |
| `--dv FLOAT`                 | Diffusion coefficient for *v* | `0.8e-5`          |
| `--F FLOAT`                  | Feed rate                     | `0.025`           |
| `--k FLOAT`                  | Kill rate                     | `0.06`            |
| `--T FLOAT`                  | Final time                    | `2000.0`          |
| `--steps INT`                | Number of time steps          | `16000`           |
| `--n-image-snapshots INT`    | How many PNGs to save         | `5`               |
| `--output-images-dir PATH`   | Directory for PNGs            | `output_images`   |
| `--output-solution-dir PATH` | Directory for NPZ/MAT         | `output_solution` |

Example
```
gray-scott-compute --nx 128 --T 2000 --steps 4000 --n-image-snapshots 5
```

### Make a video from NPZ snapshots
```
gray-scott-animate <path/to/snapshots.npz> [OPTIONS]
```
| Option               | Meaning                                 | Default                          |
| -------------------- | --------------------------------------- | -------------------------------- |
| `--which {both,u,v}` | Animate both fields or only one         | `both`                           |
| `--use-fixed-limits` | Keep constant color scale across frames | off                              |
| `--vmin FLOAT`       | Color min                               | `0.0`                            |
| `--vmax FLOAT`       | Color max                               | `1.0`                            |
| `--cmap NAME`        | Matplotlib colormap                     | `turbo`                          |
| `--fps INT`          | Frames per second                       | `15`                             |
| `--dpi INT`          | Render DPI                              | `200`                            |
| `--out PATH`         | Output file (`.mp4`)          | `output_solution/gray_scott.mp4` |

Example:
```
gray-scott-animate output_solution/snapshots.npz \
  --out output_solution/gray_scott.mp4 --fps 15 \
  --which both --use-fixed-limits --vmin 0.0 --vmax 1.0
```
