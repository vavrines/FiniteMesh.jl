using Test
import MeshArt

#MeshArt.read_mesh("t1.msh")

using PyCall
meshio = pyimport("meshio")
m0 = meshio.read("t1.msh")
cells = m0.cells
points = m0.points