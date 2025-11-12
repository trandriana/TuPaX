'''
Reconstruction of the Gray-Scott system with multiscale nudging 
Underlying scheme: semi-implicit (IMEX / FiPy)
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
import scipy.sparse
from scipy.spatial import cKDTree

# ===============================================================
# Helper functions for R, P matrices 
# ===============================================================

def bilinear_weights(fine_centers, pt, nx):
    '''Return indices and bilinear interpolation weights for a point.'''
    x, y = pt
    X = fine_centers[:, 0].reshape((nx, nx), order='F')
    Y = fine_centers[:, 1].reshape((nx, nx), order='F')

    i = np.searchsorted(X[:, 0], x) - 1
    j = np.searchsorted(Y[0, :], y) - 1
    i = np.clip(i, 0, nx - 2)
    j = np.clip(j, 0, nx - 2)

    x1, x2 = X[i, j], X[i + 1, j]
    y1, y2 = Y[i, j], Y[i, j + 1]
    dx, dy = x2 - x1, y2 - y1

    if dx == 0 or dy == 0:
        return [], []

    w11 = (x2 - x) * (y2 - y) / (dx * dy)
    w21 = (x - x1) * (y2 - y) / (dx * dy)
    w12 = (x2 - x) * (y - y1) / (dx * dy)
    w22 = (x - x1) * (y - y1) / (dx * dy)

    idx11 = i + j * nx
    idx21 = (i + 1) + j * nx
    idx12 = i + (j + 1) * nx
    idx22 = (i + 1) + (j + 1) * nx
    return [idx11, idx21, idx12, idx22], [w11, w21, w12, w22]

def build_restriction_prolongation_matrix(mesh_fine, mesh_coarse, nx):
    '''Build multiscale restriction (R) and prolongation (P) operators.'''
    fine_centers = mesh_fine.cellCenters.value.T
    coarse_centers = mesh_coarse.cellCenters.value.T

    H = float(mesh_coarse.cellVolumes[0] ** 0.5)
    q_offsets = np.array(
        [[-H / 2, -H / 2],
         [ H / 2, -H / 2],
         [-H / 2,  H / 2],
         [ H / 2,  H / 2]]
    )

    rows, cols, data = [], [], []
    for j, c in enumerate(coarse_centers):
        for pt_offset in q_offsets:
            q_point = c + pt_offset
            fine_ids, weights = bilinear_weights(fine_centers, q_point, nx)
            for fi, w in zip(fine_ids, weights):
                rows.append(j)
                cols.append(fi)
                data.append(0.25 * w)
    R = scipy.sparse.coo_matrix((data, (rows, cols)),
                                shape=(len(coarse_centers), len(fine_centers))).tocsr()

    # simple nearest neighbor prolongation
    tree = cKDTree(coarse_centers)
    _, nearest = tree.query(fine_centers, p=np.inf)
    P = scipy.sparse.coo_matrix(
        (np.ones_like(nearest), (np.arange(len(fine_centers)), nearest)),
        shape=(len(fine_centers), len(coarse_centers))
    ).tocsr()
    return R, P

# ===============================================================
# Main simulation function with assimilation
# ===============================================================
def run_simulation(
    npz_path='output_solution/reference_solution.npz',
    output_solution_dir='output_solution',
    output_images_dir='output_images',
    obs_nx=24,
    mu_u=0.0,
    mu_v=1.0,
    n_solution_snapshots=80,
    n_image_snapshots=5,
    reconstruct_time=0.0,
    reconstruct_freq=0,
):
     
    '''
    Reconstruct u,v fields from sparse observations (multiscale nudging).

    1. Load reference solution data.
    2. Set up fine and coarse meshes.
    3. Build restriction (R) and prolongation (P) matrices.
    4. Define variables and initial conditions.
    5. Define local reaction terms with nudging.
    6. Time-stepping loop to solve the equations with nudging.
    7. Save snapshots and images at specified intervals.

    Parameters:
    ----------
    - npz_path: Path to NPZ file containing reference solution.
    - output_solution_dir: Directory to save output solution snapshots.
    - output_images_dir: Directory to save output images.
    - obs_nx: Number of coarse grid points in each direction for observations.
    - mu_u: Nudging coefficient for species u.
    - mu_v: Nudging coefficient for species v.
    - n_solution_snapshots: Number of solution snapshots to save.
    - n_image_snapshots: Number of image snapshots to save.
    - reconstruct_time: Time to start reconstruction (in time units).
    - reconstruct_freq: Frequency of reconstruction updates (in time steps).
    '''

    # ---------------- Load reference data ----------------
    ref = np.load(npz_path)
    L = float(ref['L'])
    nx = int(ref['nx'])
    du, dv, F, k = map(float, (ref['du'], ref['dv'], ref['F'], ref['k']))
    T = float(ref['T'])
    steps = int(ref['steps'])
    t_ref = ref['t']
    u_ref, v_ref = ref['u'], ref['v']

    # ---------------- Mesh setup ----------------
    dx = L / nx
    mesh = Grid2D(dx=dx, dy=dx, nx=nx, ny=nx)
    x, y = mesh.cellCenters
    dt = T / steps

    obs_mesh = Grid2D(dx=L / obs_nx, dy=L / obs_nx, nx=obs_nx, ny=obs_nx)

    # ---------------- Build R, P matrices ----------------
    # R = restriction, P = prolongation
    R, P = build_restriction_prolongation_matrix(mesh, obs_mesh, nx)

    # ---------------- Variables ----------------
    tu = CellVariable(name='tu', mesh=mesh, hasOld=True)
    tv = CellVariable(name='tv', mesh=mesh, hasOld=True)
    tu.setValue(1.0)
    tv.setValue(0.0)
    perturb = (abs(x - L / 2 - 0.25) < 0.1) & (abs(y - L / 2 - 0.25) < 0.1)
    tu.setValue(0.6, where=perturb)
    tv.setValue(0.15, where=perturb)

    # ---------------- Local reaction terms ----------------
    def Rtu(tu, tv, Ih_u, Ih_tu, mu): return -tu*tv**2 + F*(1.0 - tu) + mu*(Ih_u - Ih_tu)
    def Rtv(tu, tv, Ih_v, Ih_tv, mu): return  tu*tv**2 - (F + k)*tv + mu*(Ih_v - Ih_tv)
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
        plt.title('Species u at t = %7.2f' % t_ref[step_idx])

        plt.subplot(1, 2, 2)
        # plt.imshow(v_var.value.reshape((nx, nx)), vmin=0.0, vmax=1.0, cmap='turbo', origin='lower')
        plt.pcolormesh(
            v_var.value.reshape((nx, nx)), cmap='jet', vmin=0, vmax=1, shading='auto'
        )
        plt.axis('equal')
        plt.axis('off')
        plt.colorbar(label='Concentration v', fraction=0.03, pad=0.04)
        plt.title('Species v at t = %7.2f' % t_ref[step_idx])

        # plt.tight_layout()  # avoids overlapping titles/colorbars
        fname = os.path.join(
            output_images_dir, f'GS_concentration_{int(t_ref[step_idx]):05d}.png'
        )
        plt.savefig(fname, dpi=300, bbox_inches='tight', transparent=False)
        plt.close()  # Close the figure to free memory

    # ---------------- Initialize arrays ----------------
    # A small tolerance solver is required otherwise the round-off errors destabilize the synchronization
    solver = LinearPCGSolver(tolerance=3e-16, iterations=1000)
    snap_tu, snap_tv, snap_Ih_u, snap_Ih_v = [], [], [], []

    error_u, error_v = [], []
    u, v = u_ref[0,:,:].flatten(order='F'), v_ref[0,:,:].flatten(order='F')
    e0u = np.linalg.norm(tu.value - u)
    e0v = np.linalg.norm(tv.value - v)
    error_u.append(e0u)
    error_v.append(e0v)

    # ---------------------
    # Time sampling for saving
    # ---------------------
    times_solution = np.linspace(t_ref[0], t_ref[-1], n_solution_snapshots+1)
    times_images = np.linspace(t_ref[0], t_ref[-1], n_image_snapshots+1)
    activate_u = 0.0
    activate_v = 0.0
    
    # initial 
    # Restrict and prolong
    Ih_u = P @ (R @ u)
    Ih_v = P @ (R @ v)
    Ih_tu = P @ (R @ tu.value)
    Ih_tv = P @ (R @ tv.value)

    t_rec = [t_ref[0]]
    snap_tu.append(tu.value.reshape((nx, nx), order='F'))
    snap_tv.append(tv.value.reshape((nx, nx), order='F'))
    snap_Ih_u.append(Ih_u.reshape((nx, nx), order='F'))
    snap_Ih_v.append(Ih_v.reshape((nx, nx), order='F'))
    # ---------------- Time loop ----------------
    print(f'Running reconstruction ({steps} steps)...')
    for step in range(1, steps+1):
        activate_v = float((t_ref[step] >= reconstruct_time) and (step % 2000 >= reconstruct_freq))

        # The reconstruction starts only at reconstruct_time (by default 0.0)
        # Source terms
        Stu = CellVariable(mesh=tu.mesh, value=Rtu(tu.value, tv.value, Ih_u, Ih_tu, activate_u*mu_u))
        Stv = CellVariable(mesh=tv.mesh, value=Rtv(tu.value, tv.value, Ih_v, Ih_tv, activate_v*mu_v))

        # IMEX step
        eq_tu = TransientTerm(var=tu) == DiffusionTerm(coeff=du, var=tu) + Stu
        eq_tv = TransientTerm(var=tv) == DiffusionTerm(coeff=dv, var=tv) + Stv
        eq_tu.solve(var=tu, dt=dt, solver=solver)
        eq_tv.solve(var=tv, dt=dt, solver=solver)

        tu.updateOld()
        tv.updateOld()

        u, v = u_ref[step,:,:].flatten(order='F'), v_ref[step,:,:].flatten(order='F')

        # Restrict and prolong
        Ih_u = P @ (R @ u)
        Ih_v = P @ (R @ v)
        Ih_tu = P @ (R @ tu.value)
        Ih_tv = P @ (R @ tv.value)

        # Errors
        e_u = np.linalg.norm(u - tu.value) / e0u
        e_v = np.linalg.norm(v - tv.value) / e0v
        error_u.append(e_u)
        error_v.append(e_v)

        if t_ref[step] in times_images:
            save_images(tu, tv, step)

        # Save occasional snapshots
        if t_ref[step] in times_solution:
            print(f'Step {step:5d}/{steps}, Error u={e_u:.2e}, v={e_v:.2e}')

            t_rec.append(t_ref[step])
            snap_tu.append(tu.value.copy().reshape((nx, nx), order='F'))
            snap_tv.append(tv.value.copy().reshape((nx, nx), order='F'))
            snap_Ih_u.append(Ih_u.reshape((nx, nx), order='F'))
            snap_Ih_v.append(Ih_v.reshape((nx, nx), order='F'))

    print('Saving...')
    # ---------------- Save results ----------------
    os.makedirs(output_solution_dir, exist_ok=True)

    np.savez_compressed(
        os.path.join(output_solution_dir, 'reconstructed_solution.npz'),
        t=t_rec,
        u=np.array(snap_tu),
        v=np.array(snap_tv),
        error_u=np.array(error_u),
        error_v=np.array(error_v),
        obs_nx=obs_nx, mu_u=mu_u, mu_v=mu_v,
    )

    np.savez_compressed(
        os.path.join(output_solution_dir, 'observed_solution.npz'),
        t=t_rec,
        u=np.array(snap_Ih_u),
        v=np.array(snap_Ih_v),
        obs_nx=obs_nx
    )

    # --- Plot error evolution ---
    # import matplotlib

    # matplotlib.use('Agg')
    plt.figure()
    plt.plot(t_ref, error_u, label='u L2-error')
    plt.plot(t_ref, error_v, label='v L2-error')
    plt.xlabel('Time unit (tu)')
    plt.ylabel('L2-Error')
    plt.yscale('log')
    plt.legend()
    plt.title('Error Over Time (Log Scale)')
    plt.grid(True, which='both', ls='--')
     # plt.tight_layout()  # avoids overlapping titles/colorbars
    fname = os.path.join(output_images_dir, f'error_plot.png')

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

    plt.savefig(fname, dpi=300, bbox_inches='tight', transparent=False)
    plt.close()  # Close the figure to free memory

    # scipy.io.savemat(
    #     os.path.join(output_solution_dir, 'reconstructed_solution.mat'),
    #     {
    #         't': t_ref,
    #         'tu': np.array(snap_tu),
    #         'tv': np.array(snap_tv),
    #         'Ih_u': np.array(snap_Ih_u),
    #         'Ih_v': np.array(snap_Ih_v),
    #         'error_u': np.array(error_u),
    #         'error_v': np.array(error_v),
    #         'obs_nx': obs_nx,
    #         'mu_u': mu_u,
    #         'mu_v': mu_v,
    #     },
    #     do_compression=True,
    # )

    print('Reconstruction done.')


# ===============================================================
# Command-line interface
# ===============================================================
def main():
    p = argparse.ArgumentParser(description='Gray-Scott DA (multiscale nudging)')
    p.add_argument('--npz-path', default='output_solution/reference_solution.npz')
    p.add_argument('--output-solution-dir', default='output_solution')
    p.add_argument('--output-images-dir', default='output_images')
    p.add_argument('--obs-nx', type=int, default=24)
    p.add_argument('--mu-u', type=float, default=0.0)
    p.add_argument('--mu-v', type=float, default=1.0)
    p.add_argument('--n-solution-snapshots', type=int, default=80)
    p.add_argument('--n-image-snapshots', type=int, default=5)
    p.add_argument('--reconstruct-time', type=float, default=0.0)
    p.add_argument('--reconstruct-freq', type=int, default=0)
    args = p.parse_args()

    run_simulation(
        npz_path=args.npz_path,
        output_solution_dir=args.output_solution_dir,
        output_images_dir=args.output_images_dir,
        obs_nx=args.obs_nx,
        mu_u=args.mu_u,
        mu_v=args.mu_v,
        n_solution_snapshots=args.n_solution_snapshots,
        n_image_snapshots=args.n_image_snapshots,
        reconstruct_time=args.reconstruct_time,
        reconstruct_freq=args.reconstruct_freq,
    )

# ===============================================================
# Entry point
# ===============================================================
if __name__ == '__main__':

    main()