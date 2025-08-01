"""
GMSH mesh generation script for a rectangle with embedded circular domain.
This script creates a mesh with:
- Rectangle domain: 1x4 dimensions (background domain, tag=1)
- Circular domain: centered at [0.5, 0.5], diameter=0.4 (special domain, tag=2)
- Converts GMSH mesh to DOLFINx compatible format
"""

import gmsh
from mpi4py import MPI
from dolfinx.io import gmshio
import festim as F
from basix.ufl import element
from dolfinx import fem
import dolfinx.mesh
import numpy as np
from dolfinx.io import VTXWriter


def create_mesh_with_circle() -> tuple[
    dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags
]:
    """
    Create a GMSH mesh with rectangle containing a circular domain.

    Returns:
        mesh: DOLFINx mesh object
        cell_tags: Volume tags (1=rectangle, 2=circle)
        facet_tags: Boundary tags
    """

    # Initialize GMSH
    gmsh.initialize()
    gmsh.model.add("rectangle_with_circle")

    # Define geometry parameters
    rect_width = 20
    rect_height = 10
    circle_center_x = 4
    circle_center_y = 5
    circle_radius = 0.1  # diameter 0.4 / 2

    # Mesh characteristic length (adjust for finer/coarser mesh)
    mesh_size = 0.025

    # Create rectangle points
    p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, mesh_size)
    p2 = gmsh.model.geo.addPoint(rect_width, 0.0, 0.0, mesh_size)
    p3 = gmsh.model.geo.addPoint(rect_width, rect_height, 0.0, mesh_size)
    p4 = gmsh.model.geo.addPoint(0.0, rect_height, 0.0, mesh_size)

    # Create rectangle lines
    l1 = gmsh.model.geo.addLine(p1, p2)  # bottom
    l2 = gmsh.model.geo.addLine(p2, p3)  # right
    l3 = gmsh.model.geo.addLine(p3, p4)  # top
    l4 = gmsh.model.geo.addLine(p4, p1)  # left

    # Create rectangle curve loop and surface
    rect_loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])

    # Create circle
    circle_center = gmsh.model.geo.addPoint(
        circle_center_x, circle_center_y, 0.0, mesh_size / 2
    )

    # Create circle points for better control
    p_circle_right = gmsh.model.geo.addPoint(
        circle_center_x + circle_radius, circle_center_y, 0.0, mesh_size / 2
    )
    p_circle_top = gmsh.model.geo.addPoint(
        circle_center_x, circle_center_y + circle_radius, 0.0, mesh_size / 2
    )
    p_circle_left = gmsh.model.geo.addPoint(
        circle_center_x - circle_radius, circle_center_y, 0.0, mesh_size / 2
    )
    p_circle_bottom = gmsh.model.geo.addPoint(
        circle_center_x, circle_center_y - circle_radius, 0.0, mesh_size / 2
    )

    # Create circle arcs
    arc1 = gmsh.model.geo.addCircleArc(p_circle_right, circle_center, p_circle_top)
    arc2 = gmsh.model.geo.addCircleArc(p_circle_top, circle_center, p_circle_left)
    arc3 = gmsh.model.geo.addCircleArc(p_circle_left, circle_center, p_circle_bottom)
    arc4 = gmsh.model.geo.addCircleArc(p_circle_bottom, circle_center, p_circle_right)

    # Create circle curve loop
    circle_loop = gmsh.model.geo.addCurveLoop([arc1, arc2, arc3, arc4])

    # Create surfaces
    circle_surface = gmsh.model.geo.addPlaneSurface([circle_loop])
    rect_surface = gmsh.model.geo.addPlaneSurface(
        [rect_loop, circle_loop]
    )  # rectangle minus circle

    # Synchronize the model
    gmsh.model.geo.synchronize()

    # Add physical groups for domain identification
    # Physical surfaces (volumes in 2D)
    gmsh.model.addPhysicalGroup(2, [rect_surface], tag=1, name="rectangle")
    gmsh.model.addPhysicalGroup(2, [circle_surface], tag=2, name="circle")

    # Physical curves (boundaries)
    gmsh.model.addPhysicalGroup(1, [l1], tag=1, name="bottom")
    gmsh.model.addPhysicalGroup(1, [l2], tag=2, name="right")
    gmsh.model.addPhysicalGroup(1, [l3], tag=3, name="top")
    gmsh.model.addPhysicalGroup(1, [l4], tag=4, name="left")
    gmsh.model.addPhysicalGroup(
        1, [arc1, arc2, arc3, arc4], tag=5, name="circle_boundary"
    )

    # Generate 2D mesh
    gmsh.model.mesh.generate(2)

    # Optional: optimize the mesh
    gmsh.model.mesh.optimize("Netgen")

    # Convert to DOLFINx
    mesh, cell_tags, facet_tags = gmshio.model_to_mesh(
        gmsh.model, MPI.COMM_WORLD, 0, gdim=2
    )

    # Clean up GMSH
    gmsh.finalize()

    return mesh, cell_tags, facet_tags


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
            np.full_like(x[0], 0.0),
        )

    velocity.interpolate(constant_velocity)

    if export_field:
        writer = VTXWriter(MPI.COMM_WORLD, "results/velocity_field.bp", velocity, "BP5")
        writer.write(t=0)

    return velocity


def build_festim_model(
    mesh: dolfinx.mesh.Mesh,
    cell_tags: dolfinx.mesh.MeshTags,
    facet_tags: dolfinx.mesh.MeshTags,
    velocity_field: fem.Function,
    D_value: float,
    results_filename="results/example.bp",
) -> F.HydrogenTransportProblem:
    """
    Build FESTIM model using the provided mesh and tags.
    This function is a placeholder for the FESTIM model setup.
    """

    my_model = F.HydrogenTransportProblem()

    # define mesh and meshtags
    my_model.mesh = F.Mesh(mesh=mesh)
    my_model.facet_meshtags = facet_tags
    my_model.volume_meshtags = cell_tags

    # define subdomains
    my_mat = F.Material(D_0=D_value, E_D=0)
    vol_rec = F.VolumeSubdomain(id=1, material=my_mat)
    circle = F.VolumeSubdomain(id=2, material=my_mat)
    my_model.subdomains = [vol_rec, circle]

    # define species
    H = F.Species(name="H", mobile=True, subdomains=[vol_rec, circle])
    my_model.species = [H]

    # define initial conditions
    my_model.initial_conditions = [
        F.InitialConcentration(value=1, species=H, volume=circle),
    ]

    # define temperature
    my_model.temperature = 500

    if velocity_field is not None:
        # define advection terms
        my_model.advection_terms = [
            F.AdvectionTerm(species=H, velocity=velocity_field, subdomain=vol_rec),
            F.AdvectionTerm(species=H, velocity=velocity_field, subdomain=circle),
        ]

    # define exports
    my_model.exports = [
        F.VTXSpeciesExport(
            filename=f"{results_filename}.bp",
            field=H,
            subdomain=vol_rec,
            times=[5.0, 10.0],
        ),
    ]

    # define settings
    dt = F.Stepsize(initial_value=0.1)
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
            mesh=mesh,
            cell_tags=cell_tags,
            facet_tags=facet_tags,
            velocity_field=velocity_field,
            D_value=D,
            results_filename=filename,
        )

        print(f"Running model with D={D} and results saved to {filename}")

        model.initialise()
        model.run()


def export_initial_condition(
    mesh: dolfinx.mesh.Mesh,
    cell_tags: dolfinx.mesh.MeshTags,
    facet_tags: dolfinx.mesh.MeshTags,
):
    """
    Export initial condition for the model.
    This function is a placeholder for actual export logic.
    """
    my_model = build_festim_model(
        mesh=mesh,
        cell_tags=cell_tags,
        facet_tags=facet_tags,
        velocity_field=None,
        D_value=1e-2,
    )

    my_model.initialise()

    print("Exporting initial condition...")

    writer = VTXWriter(
        MPI.COMM_WORLD, "results/initial_condition.bp", my_model.u_n, "BP5"
    )
    writer.write(t=0)


if __name__ == "__main__":
    print("Creating GMSH mesh with rectangle and circular domain...")

    # Generate the mesh every time
    mesh, cell_tags, facet_tags = create_mesh_with_circle()

    velocity_field = build_velocity_field_function(
        mesh, magnitude=1.0, export_field=True
    )

    # standard
    # model = build_festim_model(
    #     mesh,
    #     cell_tags,
    #     facet_tags,
    #     velocity_field=velocity_field,
    #     D_value=1e-02,
    #     results_filename="results/advec_case.bp",
    # )

    # model.initialise()
    # model.run()

    export_initial_condition(mesh, cell_tags, facet_tags)

    test_case()
