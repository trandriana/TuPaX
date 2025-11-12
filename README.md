<div align="center">
           <img src="https://github.com/user-attachments/assets/007e1a3c-e4e8-4b56-8442-7302773b8c9c"
           alt="TuPaX logo" width="256" height="256">   

</div >

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
├── tests
│   └── test_import.py
├── src
│   ├── tupax
│   │   └── animate.py
│   └── imex
│       ├── reconstruct.py
│       └── compute.py
├── output_solution
├── output_images
├── output_animations
├── pyproject.toml
├── README.md
├── LICENSE
└── CITATION.cff
```


## Installation

```bash
conda create -n grayscott python=3.11 -y
conda activate grayscott
conda install -c conda-forge fipy matplotlib scipy ffmpeg -y
pip install -e .
```

## Quick start

### Generate synthetic data with `imex-compute`:
```bash
imex-compute [OPTIONS]
```
  
           
| Option                       | Meaning                       | Default           |
| ---------------------------- | ----------------------------- | ----------------- |
| `--nx INT`                   | Grid resolution per side      | `240`             |
| `--L FLOAT`                  | Domain size                   | `1.0`             |
| `--du FLOAT`                 | Diffusion coefficient for *u* | `1.6e-5`          |
| `--dv FLOAT`                 | Diffusion coefficient for *v* | `0.8e-5`          |
| `--F FLOAT`                  | Feed rate                     | `0.037`           |
| `--k FLOAT`                  | Kill rate                     | `0.060`            |
| `--T FLOAT`                  | Final time                    | `4000.0`          |
| `--steps INT`                | Number of time steps          | `8000`            |
| `--n-solution-snapshots`     | Number of snapshots           | `8000`            |
| `--n-image-snapshots INT`    | How many PNGs to save         | `4`               |
| `--output-images-dir PATH`   | Directory for PNGs            | `output_images`   |
| `--output-solution-dir PATH` | Directory for NPZ             | `output_solution` |




**Example:**
We generate synthetic data, referred to as the truth, by simulating the system using the default imex-compute parameters.

```
imex-compute
```
This process generates visual outputs for quick verification and saves the numerical results as `reference_solution.npz` in the folder `output_solution/`.


### Reconstruction of the Gray--Scott system with `imex–reconstruct`:
```bash
imex-reconstruct [OPTIONS]
```

 
           
| Option                       | Meaning                                                        | Default                                  |
| ---------------------------- | -------------------------------------------------------------- | ---------------------------------------- |
| `--npz-path PATH`            | Path to the NPZ file containing the reference (truth) solution | `output_solution/reference_solution.npz` |
| `--output-solution-dir PATH` | Directory to save reconstructed solutions                      | `output_solution`                        |
| `--output-images-dir PATH`   | Directory to save reconstruction images and error plots        | `output_images`                          |
| `--obs-nx INT`               | Number of coarse grid points for observations                  | `24`                                     |
| `--mu-u FLOAT`               | Nudging coefficient for species *u*                            | `0.0`                                    |
| `--mu-v FLOAT`               | Nudging coefficient for species *v*                            | `1.0`                                    |
| `--n-solution-snapshots INT` | Number of solution snapshots to save                           | `80`                                     |
| `--n-image-snapshots INT`    | Number of images to save during the reconstruction             | `5`                                      |
| `--reconstruct-time FLOAT`   | Time (in tu) to start the reconstruction                       | `0.0`                                    |
| `--reconstruct-freq INT`     | Frequency (in steps) of reconstruction updates                 | `0`                                      |     
         


**Example:**
We reconstruct the Gray--Scott system using the default parameters:
```
imex-reconstruct
```
This process loads `output_solution/reference_solution.npz` (the synthetic truth generated by `imex-compute`),
runs the reconstruction with default nudging coefficients, and saves the numerical results as `observed_solution.npz` and `reconstructed_solution.npz` in the folder `output_solution/`.

**Coarse observation:**

![gray_scott_04](https://github.com/user-attachments/assets/5f3626b0-0187-4b5f-bfe2-cdd72aa51114)


**Fine scale reconstructed solution:**

![gray_scott_03](https://github.com/user-attachments/assets/ae3a04aa-ad2e-43b5-adcf-bcbcb2ab4b71)



The above command produces also error plots between the reconstructed solution and the reference solution:

<img width="300" height="300" alt="image" src="https://github.com/user-attachments/assets/b53e5e51-d948-4bc4-9501-4af511a1b1d5" >
 

### Create an animation of the Gray–Scott system with `imex-animate`
```
imex-animate [OPTIONS]
```
| Option                         | Meaning                                                       | Default                                  |
| :----------------------------- | :------------------------------------------------------------ | :--------------------------------------- |
| `--output-animations-dir PATH` | Directory to save the generated video                         | `output_animations`                      |
| `--npz-path PATH`              | Path to the `.npz` file containing the simulated data         | `output_solution/reference_solution.npz` |
| `--out PATH`                   | Output filename (MP4 or GIF)                                  | `gray_scott.mp4`                         |
| `--fps INT`                    | Frames per second of the video                                | `15`                                     |
| `--dpi INT`                    | Resolution of the output video (dots per inch)                | `200`                                    |
| `--cmap STR`                   | Colormap to use (any valid Matplotlib colormap)               | `jet`                                    |
| `--use-fixed-limits`           | Use fixed color scale limits (`--vmin`, `--vmax`)             | *(disabled by default)*                  |
| `--vmin FLOAT`                 | Fixed minimum color scale value                               | `0.0`                                    |
| `--vmax FLOAT`                 | Fixed maximum color scale value                               | `1.0`                                    |
| `--which {both,u,v}`           | Choose whether to animate both fields (*u*, *v*), or only one | `both`                                   |

**Example**. To generate an animation of both species from the default reconstructed solution:
```
imex-animate --npz-path output_solution/reconstructed_solution.npz
```
or observed solution:
```
imex-animate --npz-path output_solution/reconstructed_solution.npz
```
This process produces an MP4 animation named `gray_scott.mp4` or `gray_scott_01.mp4`  in the folder `output_animations/`.






