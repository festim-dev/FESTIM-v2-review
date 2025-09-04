import festim as F
import numpy as np
import dolfinx
import ufl
from mpi4py import MPI


def festim_sim_v2_disc_nietsche(n):
    # Create the mesh
    fenics_mesh = dolfinx.mesh.create_unit_square(MPI.COMM_WORLD, n, n)

    # Create the FESTIM model
    my_model = F.HydrogenTransportProblemDiscontinuous()

    my_model.mesh = F.Mesh(fenics_mesh)

    D_left, D_right = 2, 5  # diffusion coeffs
    S_left = 3
    S_right = 6

    mat_left = F.Material(D_0=D_left, E_D=0, K_S_0=S_left, E_K_S=0)
    mat_right = F.Material(D_0=D_right, E_D=0, K_S_0=S_right, E_K_S=0)
    left_volume = F.VolumeSubdomain(
        id=1, material=mat_left, locator=lambda x: x[0] < 0.5 + 1e-10
    )
    right_volume = F.VolumeSubdomain(
        id=2, material=mat_right, locator=lambda x: x[0] >= 0.5 - 1e-10
    )

    left_surface = F.SurfaceSubdomain(id=1, locator=lambda x: np.isclose(x[0], 0.0))
    right_surface = F.SurfaceSubdomain(id=2, locator=lambda x: np.isclose(x[0], 1.0))
    top_left_surface = F.SurfaceSubdomain(
        id=3, locator=lambda x: np.logical_and(np.isclose(x[1], 1.0), x[0] <= 0.5)
    )
    bottom_left_surface = F.SurfaceSubdomain(
        id=4, locator=lambda x: np.logical_and(np.isclose(x[1], 0.0), x[0] <= 0.5)
    )
    top_right_surface = F.SurfaceSubdomain(
        id=5, locator=lambda x: np.logical_and(np.isclose(x[1], 1.0), x[0] > 0.5)
    )
    bottom_right_surface = F.SurfaceSubdomain(
        id=6, locator=lambda x: np.logical_and(np.isclose(x[1], 0.0), x[0] > 0.5)
    )

    my_model.subdomains = [
        left_volume,
        right_volume,
        left_surface,
        right_surface,
        top_left_surface,
        bottom_left_surface,
        top_right_surface,
        bottom_right_surface,
    ]

    my_model.surface_to_volume = {
        left_surface: left_volume,
        right_surface: right_volume,
        top_left_surface: left_volume,
        bottom_left_surface: left_volume,
        top_right_surface: right_volume,
        bottom_right_surface: right_volume,
    }

    my_model.interfaces = [F.Interface(id=7, subdomains=[left_volume, right_volume])]
    my_model.method_interface = "nietsche"

    H = F.Species("H", subdomains=my_model.volume_subdomains)
    my_model.species = [H]

    def exact_solution_left(mod):
        return lambda x: (
            1 + mod.sin(2 * mod.pi * (x[0] + 0.25)) + mod.cos(2 * mod.pi * x[1])
        )

    def exact_solution_right(mod):
        return lambda x: S_right / S_left * exact_solution_left(mod)(x)

    exact_solution_left_ufl = exact_solution_left(ufl)
    exact_solution_right_ufl = exact_solution_right(ufl)

    # source terms
    def f_left(x):
        return -ufl.div(D_left * ufl.grad(exact_solution_left_ufl(x)))

    def f_right(x):
        return -ufl.div(D_right * ufl.grad(exact_solution_right_ufl(x)))

    my_model.sources = [
        F.ParticleSource(f_left, volume=left_volume, species=H),
        F.ParticleSource(f_right, volume=right_volume, species=H),
    ]

    # boundary conditions
    my_model.boundary_conditions = [
        F.FixedConcentrationBC(
            subdomain=surface, value=exact_solution_left_ufl, species=H
        )
        for surface in [left_surface, top_left_surface, bottom_left_surface]
    ] + [
        F.FixedConcentrationBC(
            subdomain=surface, value=exact_solution_right_ufl, species=H
        )
        for surface in [right_surface, top_right_surface, bottom_right_surface]
    ]

    my_model.temperature = 500.0

    my_model.settings = F.Settings(atol=1e-10, rtol=1e-10, transient=False)

    my_model.exports = [
        F.VTXSpeciesExport(
            filename="results/disc_nietsche_l.bp", field=H, subdomain=left_volume
        ),
        F.VTXSpeciesExport(
            filename="results/disc_nietsche_r.bp", field=H, subdomain=right_volume
        ),
    ]

    my_model.initialise()
    my_model.run()
