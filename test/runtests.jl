using Test
import MeshArt

#MeshArt.read_mesh("t1.msh")

#using PyCall
#meshio = pyimport("meshio")
#meshio.read("t1.msh")
#cells = m0.cells[end][1]
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


cd(@__DIR__)
cellid, nodes = MeshArt.read_mesh("t1.msh")
edgeNodes, edgeCells, cellNeighbors = MeshArt.mesh_connectivity_2D(cellid)
cellType = MeshArt.mesh_cell_type(cellNeighbors)
cellArea = MeshArt.mesh_area_2D(nodes, cellid)
cellCenter = MeshArt.mesh_center_2D(nodes, cellid)
edgeCenter = MeshArt.mesh_edge_center(nodes, edgeNodes)
cellEdges = MeshArt.mesh_cell_edge(cellid, edgeCells)

cellid, nodes = MeshArt.read_mesh("square.msh")
edgeNodes, edgeCells, cellNeighbors = MeshArt.mesh_connectivity_2D(cellid)
cellArea = MeshArt.mesh_area_2D(nodes, cellid)