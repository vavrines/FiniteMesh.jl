using Test
import MeshArt

#MeshArt.read_mesh("t1.msh")

using PyCall
meshio = pyimport("meshio")
m0 = meshio.read("t1.msh")
cells = m0.cells[end][1]
#points = m0.points
#=
py"""
import meshio

def read_mesh(file):
    m0 = meshio.read(file)

    points = m0.points
    cells = m0.cells
    for cell in cells:
        if not(cell[0] in ["line", "vertex"]):
            _cells = cell[1] + 1

    return _cells, points
"""
cells, nodes = py"read_mesh"("t1.msh")=#