from foam2dolfinx import OpenFOAMReader
import festim as F
from mpi4py import MPI
from dolfinx.mesh import create_rectangle
import numpy as np
from pathlib import Path
from dolfinx import fem
import zipfile


def read_velocity_field():
    zip_path = Path("data/cavity.zip")
    extract_path = Path("data/cavity")

    # if the extract path does not exist, create it
    if not extract_path.exists():
        # Unzip the file
        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            zip_ref.extractall(extract_path)

        # make a cavity.foam file
        open("data/cavity/cavity/cavity.foam", "w").close()

    # Read the velocity field from the OpenFOAM reader
    ldc_reader = OpenFOAMReader(filename="data/cavity/cavity/cavity.foam", cell_type=12)
    velocity_field = ldc_reader.create_dolfinx_function(t=2.5, name="U")

    return velocity_field


def build_festim_model(
    velocity_field: fem.Function | None = None,
    results_filename="results/mobile_concentration.bp",
) -> F.HydrogenTransportProblem:
    """_summary_

    Args:
        mesh: mesh to use for the model
        velocity_field: velocity field. Defaults to None.
        results_filename: Filename of the exported field. Defaults to "results/example.bp".

    Returns:
        F.HydrogenTransportProblem: FESTIM model object
    """
    my_model = F.HydrogenTransportProblem()

    fenics_mesh = create_rectangle(MPI.COMM_WORLD, [[0, 0], [0.1, 0.1]], [50, 50])

    my_mesh = F.Mesh(fenics_mesh)
    my_model.mesh = my_mesh

    # define species
    H = F.Species("H", mobile=True)
    my_model.species = [H]

    # define subdomains
    my_mat = F.Material(D_0=5e-04, E_D=0)

    fluid = F.VolumeSubdomain(id=1, material=my_mat)
    left = F.SurfaceSubdomain(id=2, locator=lambda x: np.isclose(x[0], 0.0))
    right = F.SurfaceSubdomain(id=3, locator=lambda x: np.isclose(x[0], 0.1))
    bottom = F.SurfaceSubdomain(id=4, locator=lambda x: np.isclose(x[1], 0.0))
    top = F.SurfaceSubdomain(id=5, locator=lambda x: np.isclose(x[1], 0.1))
    my_model.subdomains = [fluid, left, right, bottom, top]

    # define boundary conditions
    my_model.boundary_conditions = [
        F.FixedConcentrationBC(subdomain=left, value=0, species=H),
        F.FixedConcentrationBC(subdomain=right, value=0, species=H),
        F.FixedConcentrationBC(subdomain=bottom, value=0, species=H),
        F.FixedConcentrationBC(subdomain=top, value=0, species=H),
    ]

    # define temperature
    my_model.temperature = 500

    # define sources
    my_model.sources = [
        F.ParticleSource(volume=fluid, species=H, value=10),
    ]

    if velocity_field:
        # define advection term
        my_model.advection_terms = [
            F.AdvectionTerm(
                velocity=lambda t: velocity_field(t), species=H, subdomain=fluid
            )
        ]

    # define settings
    my_model.settings = F.Settings(atol=1e-10, rtol=1e-10, transient=False)

    # define exports
    my_model.exports = [F.VTXSpeciesExport(filename=f"{results_filename}", field=H)]

    return my_model


def test_case(velocity_field: fem.Function):
    # case without field
    model_1 = build_festim_model(
        velocity_field=None, results_filename="results/mobile_concentration_standard.bp"
    )
    model_1.initialise()
    print("Running model without advection...")
    model_1.run()

    # case with velocity field
    model_2 = build_festim_model(
        velocity_field=velocity_field,
        results_filename="results/mobile_concentration_advection.bp",
    )
    model_2.initialise()
    print("Running model with advection...")
    model_2.run()


if __name__ == "__main__":
    """
    Run the test case.
    """
    # Generate the mesh every time
    velocity_field = read_velocity_field()

    test_case(velocity_field=velocity_field)
