from mpi4py import MPI
from petsc4py import PETSc

import basix
import numpy as np
import tqdm.autonotebook
from dolfinx import fem, mesh
from dolfinx.fem.petsc import (
    NonlinearProblem,
)
from dolfinx.nls.petsc import NewtonSolver
import ufl
from dolfinx.io import VTXWriter


# define geometry
indices = np.linspace(0, 3e-04, num=1000)
domain = ufl.Mesh(basix.ufl.element("Lagrange", "interval", 1, shape=(1,)))
mesh_points = np.reshape(indices, (len(indices), 1))
indexes = np.arange(mesh_points.shape[0])
cells = np.stack((indexes[:-1], indexes[1:]), axis=-1)
my_mesh = mesh.create_mesh(MPI.COMM_WORLD, cells, mesh_points, domain)
vdim = my_mesh.topology.dim
fdim = vdim - 1

# define function space and measures
element_CG = basix.ufl.element(
    basix.ElementFamily.P,
    my_mesh.basix_cell(),
    1,
    basix.LagrangeVariant.equispaced,
)
element_DG = basix.ufl.element(
    "DG",
    my_mesh.basix_cell(),
    1,
    basix.LagrangeVariant.equispaced,
)
mixed_element = basix.ufl.mixed_element([element_CG, element_DG])
V = fem.functionspace(my_mesh, mixed_element)

num_cells = my_mesh.topology.index_map(vdim).size_local
cells = np.arange(num_cells, dtype=np.int32)
markers = np.full(num_cells, 1, dtype=np.int32)
mesh_entities = mesh.locate_entities(
    my_mesh,
    vdim,
    lambda x: np.logical_and(x[0] >= indices[0], x[0] <= indices[1]),
)
markers[mesh_entities] = 1
mesh_tags_volumes = mesh.meshtags(my_mesh, vdim, cells, markers)
num_facets = my_mesh.topology.index_map(fdim).size_local
mesh_facet_indices = np.arange(num_facets, dtype=np.int32)
tags_facets = np.full(num_facets, 0, dtype=np.int32)
indices_left = mesh.locate_entities_boundary(my_mesh, 0, lambda x: np.isclose(x[0], 0))
tags_facets[indices_left] = 1
indices_right = mesh.locate_entities_boundary(
    my_mesh, 0, lambda x: np.isclose(x[0], indices[-1])
)
tags_facets[indices_right] = 2
mesh_tags_facets = mesh.meshtags(my_mesh, fdim, mesh_facet_indices, tags_facets)
dx = ufl.Measure("dx", domain=my_mesh, subdomain_data=mesh_tags_volumes)

# define solutions
u = fem.Function(V)
u_n = fem.Function(V)
v = ufl.TestFunction(V)
mobile, trapped = ufl.split(u)
mobile_n, trapped_n = ufl.split(u_n)
mobile_test, trapped_test = ufl.TestFunctions(V)
mobile_sub_function_space = V.sub(0)
trapped_sub_function_space = V.sub(1)
mobile_sub_function = u.sub(0)
trapped_sub_function = u.sub(1)

# define boundary conditions
left_facets = mesh_tags_facets.find(1)
left_dofs = fem.locate_dofs_topological(mobile_sub_function_space, fdim, left_facets)
bc_left = fem.dirichletbc(
    fem.Constant(my_mesh, PETSc.ScalarType(1e12)),
    left_dofs,
    mobile_sub_function_space,
)
right_facets = mesh_tags_facets.find(2)
right_dofs = fem.locate_dofs_topological(mobile_sub_function_space, fdim, right_facets)
bc_right = fem.dirichletbc(
    fem.Constant(my_mesh, PETSc.ScalarType(0)),
    right_dofs,
    mobile_sub_function_space,
)
bcs = [bc_left, bc_right]

# define variational formulation
temperature = 500
k_B = 8.6173303e-5
D = 1.90e-7 * ufl.exp(-0.2 / k_B / temperature)
dt = 1 / 10
final_time = 20
k = 3.8e-17 * ufl.exp(-0.2 / k_B / temperature)
p = 1e13 * ufl.exp(-1.2 / k_B / temperature)
n_traps = 1e19
F = ufl.dot(D * ufl.grad(mobile), ufl.grad(mobile_test)) * dx(1)
F += ((mobile - mobile_n) / dt) * mobile_test * dx(1)
F += ((trapped - trapped_n) / dt) * trapped_test * dx(1)
F += (k * mobile * (n_traps - trapped) - p * trapped) * mobile_test * dx(1)
F += -(k * mobile * (n_traps - trapped) - p * trapped) * trapped_test * dx(1)

# build solver
problem = NonlinearProblem(F, u, bcs=bcs)
solver = NewtonSolver(MPI.COMM_WORLD, problem)
solver.rtol = 1e-10
solver.atol = 1e-10
solver.max_it = 30
ksp = solver.krylov_solver
ksp.setType("preonly")
ksp.getPC().setType("lu")
ksp.getPC().setFactorSolverType("mumps")
ksp.setErrorIfNotConverged(True)
writer_mobile = VTXWriter(
    MPI.COMM_WORLD,
    "results/fenicsx_mobile.bp",
    [mobile_sub_function],
    "BP5",
)
writer_trapped = VTXWriter(
    MPI.COMM_WORLD,
    "results/fenicsx_trapped.bp",
    [trapped_sub_function],
    "BP5",
)
t = 0
progress = tqdm.autonotebook.tqdm(
    desc="Solving H transport problem", total=final_time, unit_scale=True
)
while t < final_time:
    progress.update(float(dt))
    t += float(dt)

    solver.solve(u)

    writer_mobile.write(float(t))
    writer_trapped.write(float(t))

    u_n.x.array[:] = u.x.array[:]

progress.refresh()
progress.close()

writer_mobile.close()
writer_trapped.close()
