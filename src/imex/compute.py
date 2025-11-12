'''
Gray-Scott FiPy simulation with a semi-implicit (IMEX) scheme.
Author: Tsiry Avisoa Randrianasolo
'''

# ---------------------
# Standard library imports
# ---------------------
import os
import scipy.io

# ---------------------
# Matplotlib setup
# ---------------------
import matplotlib

matplotlib.use('Agg')  # Use a non-interactive backend
import matplotlib.pyplot as plt

# ---------------------
# FiPy imports
# ---------------------
from fipy import CellVariable, Grid2D, TransientTerm, DiffusionTerm 
from fipy.solvers.scipy import LinearPCGSolver  # PCG = Preconditioned Conjugate Gradient

# ---------------------
# Other imports
# ---------------------
import numpy as np
import argparse


# ===============================================================
# Main simulation function
# ===============================================================
def run_simulation(
    output_images_dir: str = 'output_images',
    output_solution_dir: str = 'output_solution',
    L: float = 1.0,
    nx: int = 240,
    du: float = 1.6e-5,
    dv: float = 0.8e-5,
    F: float = 0.037,
    k: float = 0.060,
    T: float = 4000.0,
    steps: int = 8000,
    n_solution_snapshots: int = 8000,
    n_image_snapshots: int = 5,
):
    '''
    Run the Gray-Scott simulation with IMEX scheme.

    1. Set up the domain and mesh.
    2. Define variables and initial conditions.
    3. Define the equations using FiPy.
    4. Time-stepping loop to solve the equations.
    5. Save snapshots and images at specified intervals.

     Parameters:
     ----------
    - output_images_dir: Directory to save output images.
    - output_solution_dir: Directory to save solution snapshots.
    - L: Length of the domain (assumed square).
    - nx: Number of grid points in each direction.
    - du: Diffusion coefficient for species u.
    - dv: Diffusion coefficient for species v.
    - F: Feed rate.
    - k: Kill rate.
    - T: Total time to simulate.
    - steps: Number of time steps.
    - n_solution_snapshots: Number of solution snapshots to save.
    - n_image_snapshots: Number of image snapshots to save.

    '''

    # ---------------------
    # Domain and mesh
    # ---------------------
    dx = L / nx
    mesh = Grid2D(dx=dx, dy=dx, nx=nx, ny=nx)
    x, y = mesh.cellCenters[0], mesh.cellCenters[1]

    # ---------------------
    # Temporal discretization
    # ---------------------
    t0 = 0.0
    t = np.linspace(t0, T, steps + 1)
    dt = t[1] - t[0]

    # ---------------------
    # Variables
    # ---------------------
    u = CellVariable(name='u', mesh=mesh)
    v = CellVariable(name='v', mesh=mesh)
    u_old = CellVariable(mesh=mesh)
    v_old = CellVariable(mesh=mesh)

    # ---------------------
    # Initial conditions
    # ---------------------
    u.setValue(1.0)
    v.setValue(0.0)
    perturb = (abs(x - L / 2) <= 0.1) & (abs(y - L / 2) <= 0.10)
    u.setValue(0.50, where=perturb)
    v.setValue(0.25, where=perturb)

    # ---------------------
    # Define equations
    # ---------------------
    u_old.setValue(u.value.copy())
    v_old.setValue(v.value.copy())

    eq_u = TransientTerm(var=u) == DiffusionTerm(coeff=du, var=u) + (
        -u_old * v_old * v_old + F * (1.0 - u_old)
    )
    eq_v = TransientTerm(var=v) == DiffusionTerm(coeff=dv, var=v) + (
        u_old * v_old * v_old - (F + k) * v_old
    )

    # ---------------------
    # Output directories
    # ---------------------
    os.makedirs(output_images_dir, exist_ok=True)
    os.makedirs(output_solution_dir, exist_ok=True)

    # ---------------------
    # Function to save images
    # ---------------------
    def save_images(u_var, v_var, step_idx):
        plt.figure(figsize=(16, 5))
        plt.subplot(1, 2, 1)
        # plt.imshow(u_var.value.reshape((nx, nx)), vmin=0.0, vmax=1.0, cmap='turbo', origin='lower')
        plt.pcolormesh(
            u_var.value.reshape((nx, nx)), cmap='jet', vmin=0, vmax=1, shading='auto'
        )
        plt.axis('equal')
        plt.axis('off')
        plt.colorbar(label='Concentration u', fraction=0.03, pad=0.04)
        plt.title(f'Species u - {t[step_idx]:7.2f} tu')

        plt.subplot(1, 2, 2)
        # plt.imshow(v_var.value.reshape((nx, nx)), vmin=0.0, vmax=1.0, cmap='turbo', origin='lower')
        plt.pcolormesh(
            v_var.value.reshape((nx, nx)), cmap='jet', vmin=0, vmax=1, shading='auto'
        )
        plt.axis('equal')
        plt.axis('off')
        plt.colorbar(label='Concentration v', fraction=0.03, pad=0.04)
        plt.title(f'Species v - {t[step_idx]:7.2f} tu')

        # plt.tight_layout()  # avoids overlapping titles/colorbars
        fname = os.path.join(
            output_images_dir, f'GS_concentration_{int(t[step_idx]):05d}.png'
        )
        plt.savefig(fname, dpi=300, bbox_inches='tight', transparent=False)
        plt.close()  # Close the figure to free memory

    # ---------------------
    # Time sampling for saving
    # ---------------------
    
    times_solution = np.linspace(t0, T, n_solution_snapshots+1)
    times_images = np.linspace(t0, T, n_image_snapshots+1)

    # initial save
    save_images(u, v, 0)
    snap_u = [u.value.copy().reshape((nx, nx), order='F')]
    snap_v = [v.value.copy().reshape((nx, nx), order='F')]

    # solver
    solver = LinearPCGSolver(tolerance=3e-16, iterations=1000)

    print('Running...')

    for step in range(1, steps + 1):

        u_old.setValue(u.value.copy())
        v_old.setValue(v.value.copy())

        eq_u.solve(var=u, dt=dt, solver=solver)
        eq_v.solve(var=v, dt=dt, solver=solver)

        if t[step] in times_images:
            save_images(u, v, step)

        if t[step] in times_solution:

            snap_u.append(u.value.copy().reshape((nx, nx), order='F'))
            snap_v.append(v.value.copy().reshape((nx, nx), order='F'))

            print(f'Step {step:5d}/{steps}, t = {t[step]:07.2f}')

    print('Saving...')
    # ---------------------
    # Save snapshots + metadata
    # ---------------------
    metadata = dict(  # to be saved alongside the snapshots
        L=L,
        nx=nx,
        du=du,
        dv=dv,
        F=F,
        k=k,
        T=T,
        steps=steps,
        n_solution_snapshots=n_solution_snapshots,
        n_image_snapshots=n_image_snapshots,
    )

    # Save to NPZ for Python
    np.savez_compressed(
        os.path.join(output_solution_dir, 'reference_solution.npz'),
        t=np.array(times_solution),
        u=np.array(snap_u),
        v=np.array(snap_v),
        **metadata,  # injects the metadata as arrays
    )

    # Save to MAT for MATLAB
    # scipy.io.savemat(
    #     os.path.join(output_solution_dir, 'snapshots.mat'),
    #     {
    #         't': np.array(times_solution),
    #         'u': np.array(snap_u),
    #         'v': np.array(snap_v),
    #         **metadata,  # adds same metadata in MATLAB file
    #     },
    #     do_compression=True,
    # )

    print('Done.')


# ===============================================================
# Command-line interface
# ===============================================================
def main():
    p = argparse.ArgumentParser(
        description='Gray-Scott FiPy simulation with IMEX scheme'
    )
    p.add_argument('--output-images-dir', default='output_images')
    p.add_argument('--output-solution-dir', default='output_solution')
    p.add_argument('--L', type=float, default=1.0)
    p.add_argument('--nx', type=int, default=240)
    p.add_argument('--du', type=float, default=1.6e-5)
    p.add_argument('--dv', type=float, default=0.8e-5)
    p.add_argument('--F', type=float, default=0.037)
    p.add_argument('--k', type=float, default=0.060)
    p.add_argument('--T', type=float, default=4000.0)
    p.add_argument('--steps', type=int, default=8000)
    p.add_argument('--n-solution-snapshots', type=int, default=8000)
    p.add_argument('--n-image-snapshots', type=int, default=4)
    args = p.parse_args()

    run_simulation(
        output_images_dir=args.output_images_dir,
        output_solution_dir=args.output_solution_dir,
        L=args.L,
        nx=args.nx,
        du=args.du,
        dv=args.dv,
        F=args.F,
        k=args.k,
        T=args.T,
        steps=args.steps,
        n_solution_snapshots=args.n_solution_snapshots,
        n_image_snapshots=args.n_image_snapshots,
    )


# ===============================================================
# Entry point
# ===============================================================
if __name__ == '__main__':

    main()