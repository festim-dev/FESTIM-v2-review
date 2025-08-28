import openmc
import openmc_data_downloader as odd
from openmc2dolfinx import StructuredGridReader
import festim as F
import numpy as np
import dolfinx
from mpi4py import MPI

dim = 60


def run_openmc_example():
    # Materials

    lithium = openmc.Material(name="lithium")
    lithium.set_density("g/cc", 0.534)
    lithium.add_element("Li", 1.0)

    mats = openmc.Materials([lithium])

    odd.download_cross_section_data(
        mats,
        libraries=["FENDL-3.1d"],
        set_OPENMC_CROSS_SECTIONS=True,
        particles=["neutron"],
    )

    # Geometry

    cube_surface = openmc.model.RectangularParallelepiped(
        -dim, dim, -dim, dim, -dim, dim
    )
    region = -cube_surface
    cell = openmc.Cell(region=region, fill=lithium)

    vacuum_surf = openmc.Sphere(r=dim * 2, boundary_type="vacuum")
    vacuum_region = +cube_surface & -vacuum_surf
    vacuum = openmc.Cell(region=vacuum_region, fill=None)

    geometry = openmc.Geometry([cell, vacuum])

    # Tallies
    tally = openmc.Tally(name="tritium_production")
    tally.scores = ["(n,Xt)"]
    mesh = openmc.RegularMesh()
    mesh.dimension = [30, 30, 15]
    mesh.lower_left = [-dim, -dim, -dim]
    mesh.upper_right = [dim, dim, dim]
    tally.filters = [openmc.MeshFilter(mesh)]

    # cell_filter = openmc.CellFilter(cell)
    # tally.filters = [cell_filter]

    tallies = openmc.Tallies([tally])

    # Settings

    source = openmc.IndependentSource()
    source_pos_z = dim + 10
    source.space = openmc.stats.Point((0, 0, source_pos_z))
    source.angle = openmc.stats.Isotropic()
    source.energy = openmc.stats.Discrete([14.1e6], [1.0])

    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.source = source
    settings.batches = 10
    settings.particles = 10000

    my_model = openmc.Model(
        geometry=geometry, settings=settings, materials=mats, tallies=tallies
    )

    my_model.run(apply_tally_results=True)

    mesh.write_data_to_vtk(
        "tritium_production_mesh.vtk", {"tritium_production": tally.mean}
    )


def festim_model(use_openmc_mesh=True):
    """
    Args:
        use_openmc_mesh (bool, optional): if True, FESTIM will use the OpenMC mesh, otherwise
            we make our own.
    """
    # read openmc vtk
    reader = StructuredGridReader("tritium_production_mesh.vtk")
    reader.create_dolfinx_mesh()

    my_model = F.HydrogenTransportProblem()

    tritium = F.Species("T")

    if use_openmc_mesh:
        my_model.mesh = F.Mesh(reader.dolfinx_mesh)
    else:
        refined_mesh = dolfinx.mesh.create_box(
            MPI.COMM_WORLD,
            points=[(-dim, -dim, -dim), (dim, dim, dim)],
            n=(30, 30, 70),
            cell_type=dolfinx.mesh.CellType.tetrahedron,
        )
        my_model.mesh = F.Mesh(refined_mesh)

    boundary = F.SurfaceSubdomain(id=1)

    my_mat = F.Material(D_0=1, E_D=0)

    volume = F.VolumeSubdomain(id=1, material=my_mat)

    my_model.subdomains = [boundary, volume]

    my_model.boundary_conditions = [
        F.FixedConcentrationBC(subdomain=boundary, value=0.0, species=tritium)
    ]

    if use_openmc_mesh:
        source_term = reader.create_dolfinx_function("tritium_production")
    else:
        original_source_term = reader.create_dolfinx_function("tritium_production")
        V = dolfinx.fem.functionspace(refined_mesh, ("DG", 0))
        source_term = dolfinx.fem.Function(V)
        F.helpers.nmm_interpolate(source_term, original_source_term)
        # export to VTX
        writer = dolfinx.io.VTXWriter(
            MPI.COMM_WORLD, "results/source.bp", [source_term], "BP5"
        )
        writer.write(0.0)

    my_model.sources = [
        F.ParticleSource(
            value=source_term,
            volume=volume,
            species=tritium,
        )
    ]

    my_model.species = [tritium]

    my_model.temperature = 300

    my_model.settings = F.Settings(
        atol=1e-10,
        rtol=1e-10,
        transient=False,
    )

    my_model.exports = [
        F.VTXSpeciesExport(filename="results/tritium_conc.bp", field=tritium)
    ]

    my_model.initialise()
    my_model.run()


if __name__ == "__main__":
    run_openmc_example()
    festim_model(use_openmc_mesh=False)
