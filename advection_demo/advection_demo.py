from mpi4py import MPI
import festim as F
from basix.ufl import element
from dolfinx import fem
import dolfinx.mesh
import numpy as np
from dolfinx.io import VTXWriter
import ufl


def build_velocity_field_function(
    mesh: dolfinx.mesh.Mesh, magnitude: float = 1.0, export_field: bool = False
) -> fem.Function:
    """
    Build a velocity field function for advection.
    This function is a placeholder and should be replaced with actual velocity logic.
    """
    el = element("Lagrange", mesh.topology.cell_name(), 2, shape=(mesh.geometry.dim,))
    V = fem.functionspace(mesh, el)
    velocity = fem.Function(V)

    def constant_velocity(x):
        return (
            np.full_like(x[0], float(magnitude)),
            np.full_like(x[1], 0.0),
        )

    velocity.interpolate(constant_velocity)

    if export_field:
        writer = VTXWriter(MPI.COMM_WORLD, "results/velocity_field.bp", velocity, "BP5")
        writer.write(t=0)

    return velocity


def build_festim_model(
    D_value: float,
    results_filename="results/example.bp",
    with_velocity_field=True,
) -> F.HydrogenTransportProblem:
    """
    Build FESTIM model using the provided mesh and tags.
    This function is a placeholder for the FESTIM model setup.
    """

    my_model = F.HydrogenTransportProblem()

    # define mesh and meshtags
    rectangle_width = 20.0
    rectangle_height = 10.0
    mesh = dolfinx.mesh.create_rectangle(
        MPI.COMM_WORLD,
        [[0.0, 0.0], [rectangle_width, rectangle_height]],
        n=[int(rectangle_width * 20), int(rectangle_height * 20)],
    )
    my_model.mesh = F.Mesh(mesh=mesh)

    # define subdomains
    my_mat = F.Material(D_0=D_value, E_D=0)
    vol = F.VolumeSubdomain(id=1, material=my_mat)
    left = F.SurfaceSubdomain(id=1, locator=lambda x: np.isclose(x[0], 0.0))
    my_model.subdomains = [vol, left]

    # define species
    H = F.Species(name="H", mobile=True)
    my_model.species = [H]

    # define initial condition
    def initial_condition(x, a=5):
        """Initial condition: Gaussian centered in the left part of the rectangle"""
        pos_x = rectangle_width / 5
        pos_y = rectangle_height / 2
        return ufl.exp(-a * ((x[0] - pos_x) ** 2 + (x[1] - pos_y) ** 2))

    my_model.initial_conditions = [
        F.InitialConcentration(value=initial_condition, species=H, volume=vol),
    ]

    # define temperature
    my_model.temperature = 500

    if with_velocity_field:
        velocity_field = build_velocity_field_function(
            mesh=my_model.mesh.mesh,
            magnitude=1,
            export_field=True,
        )
        # define advection terms
        my_model.advection_terms = [
            F.AdvectionTerm(species=H, velocity=velocity_field, subdomain=vol),
        ]

    # define exports
    my_model.exports = [
        F.VTXSpeciesExport(
            filename=f"{results_filename}.bp",
            field=H,
            subdomain=vol,
            times=[5, 10],
        ),
    ]

    dt = F.Stepsize(initial_value=0.02)
    my_model.settings = F.Settings(
        atol=1e-10,
        rtol=1e-10,
        transient=True,
        final_time=10,
        stepsize=dt,
    )

    return my_model


def test_case():
    """
    Test case to demonstrate the mesh generation and model setup.
    This function is a placeholder for actual test logic.
    """
    print("Running test case...")
    for D, filename in zip(
        [1.5, 1e-2],
        ["results/diff_case", "results/advec_case"],
    ):
        model = build_festim_model(
            D_value=D,
            results_filename=filename,
        )

        print(f"Running model with D={D} and results saved to {filename}")

        model.initialise()
        model.run()


def export_initial_condition():
    """
    Export initial condition for the model.
    This function is a placeholder for actual export logic.
    """
    my_model = build_festim_model(D_value=1)

    my_model.initialise()

    print("Exporting initial condition...")

    writer = VTXWriter(
        MPI.COMM_WORLD, "results/initial_condition.bp", my_model.u_n, "BP5"
    )
    writer.write(t=0)


if __name__ == "__main__":
    # model = build_festim_model(
    #     D_value=1e-02,
    #     results_filename="results/test",
    #     with_velocity_field=True,
    # )

    # model.initialise()

    # model.run()

    export_initial_condition()

    test_case()
