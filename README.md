<img src="https://github.com/user-attachments/assets/007e1a3c-e4e8-4b56-8442-7302773b8c9c"
           alt="TuPaX logo" width="256" height="256">     
[![DOI](https://zenodo.org/badge/1074261727.svg)](https://doi.org/10.5281/zenodo.17328875)
[![arXiv](https://img.shields.io/badge/arXiv-2508.18910-b31b1b.svg)](https://doi.org/10.48550/arXiv.2508.18910)
[![arXiv](https://img.shields.io/badge/arXiv-2510.03972-b31b1b.svg)](https://doi.org/10.48550/arXiv.2510.03972)
[![Powered by FiPy](https://img.shields.io/badge/Powered%20by-FiPy-1f425f.svg)](https://github.com/usnistgov/fipy)

# Turing Pattern eXperiments

This code implements a **discrete data assimilation algorithm** for the Gray--Scott reaction-diffusion model, using a semi-implicit (IMEX) finite-volume scheme built with FiPy. The continuous form of the data assimilation algorithm reads

$$
\partial_t \tilde{u} = d_u \Delta \tilde{u} - \tilde{u} \tilde{v}^2 + F(1-\tilde{u}) + \mu_u (\mathcal{I}_H\tilde{u} - \mathcal{I}_H u ), \qquad
\partial_t \tilde{v} = d_v \Delta \tilde{v} + \tilde{u} \tilde{v}^2 - (F+k)\tilde{v} + \mu_v (\mathcal{I}_H\tilde{v} - \mathcal{I}_H v ),
$$

with Neumann boundary conditions and a finite volume interpolant $\mathcal{I}_H$ that characterizes the coarse observation we have on the Truth $(u,v)$.
<div align="center">
           
| Symbol | Meaning |
|:------:|:--------|
| $u, v$ | species concentrations (Truth) |
| $\tilde{u}, \tilde{v}$ | reconstructed concentrations  |
| $d_u, d_v$ | diffusion coefficients |
| $F, k$ | feed and kill rates |
| $H$ | observation resolution |
| $(\mu_u, \mu_v)$ | nudging gain |

</div >

The framework supports a **multigrid (multiresolution) approach**, where the observations are defined on a **coarse grid** (low resolution) and the model state is reconstructed on a **fine grid** (high resolution). This enables the recovery of fine-scale Gray--Scott patterns from sparse or low-resolution observations.

The underlying IMEX solver can also be used on its own to simulate the Gray--Scott system by setting $\mu_u = \mu_v = 0$.

The data assimilation module updates the model state at discrete, possibly sparse, times using coarse or noisy observations, while the IMEX solver advances the forecast between updates. It can therefore both **generate synthetic data** and **reconstruct the full state** from partial measurements.




---

## Project structure
```
TuPaX/
├── tests
│   └── test_import.py
├── src
│   ├── tupax
│   │   ├── animate.py
│   └── imex
│       ├── reconstruct.py
│       ├── compute.py
├── pyproject.toml
├── output_solution
├── output_images
├── output_animations
├── examples
│   ├── run_default.sh
│   └── make_video.sh
├── README.md
├── MANIFEST.in
├── LICENSE
└── CITATION.cff
```


## Installation

```bash
conda create -n env-tupax
conda activate env-tupax
conda install -c conda-forge fipy matplotlib scipy ffmpeg
pip install -e .
```

## Quick start

### Generate synthetic data:
```bash
imex-compute [OPTIONS]
```
| Option                       | Meaning                       | Default           |
| ---------------------------- | ----------------------------- | ----------------- |
| `--nx INT`                   | Grid resolution per side      | `240`             |
| `--L FLOAT`                  | Domain size                   | `1.0`             |
| `--du FLOAT`                 | Diffusion coefficient for *u* | `1.6e-5`          |
| `--dv FLOAT`                 | Diffusion coefficient for *v* | `0.8e-5`          |
| `--F FLOAT`                  | Feed rate                     | `0.025`           |
| `--k FLOAT`                  | Kill rate                     | `0.06`            |
| `--T FLOAT`                  | Final time                    | `4000.0`          |
| `--steps INT`                | Number of time steps          | `8000`            |
| `--n-solution-snapshots`     | Number of snapshots           | `8000`            |
| `--n-image-snapshots INT`    | How many PNGs to save         | `5`               |
| `--output-images-dir PATH`   | Directory for PNGs            | `output_images`   |
| `--output-solution-dir PATH` | Directory for NPZ             | `output_solution` |

**Example:**
```
imex-compute
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
| `--out PATH`         | Output file (`.mp4`)                    | `output_solution/gray_scott.mp4` |

Example:
```
gray-scott-animate output_solution/snapshots.npz \
  --out output_solution/gray_scott.mp4 --fps 15 \
  --which both --use-fixed-limits --vmin 0.0 --vmax 1.0
```




### Reconstruction of the Gray--Scott dynamics

<img width="300" height="300" alt="image" src="https://github.com/user-attachments/assets/b53e5e51-d948-4bc4-9501-4af511a1b1d5"  />
<img width="300" height="300" alt="image" src="https://github.com/user-attachments/assets/a45224b9-8a9d-43f6-a9a1-fbb9af2930b5" />
<img width="300" height="300" alt="image" src="https://github.com/user-attachments/assets/c580caf9-e038-4389-811a-e2434c7d0877" />






