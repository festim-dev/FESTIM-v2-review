import openmc
import openmc_data_downloader as odd
import matplotlib.pyplot as plt

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

dim = 60
cube_surface = openmc.model.RectangularParallelepiped(-dim, dim, -dim, dim, -dim, dim)
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
