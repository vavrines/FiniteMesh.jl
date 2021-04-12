using Test
import MeshArt

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
