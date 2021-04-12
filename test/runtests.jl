using Test
import Mesher

cd(@__DIR__)
cellid, nodes = Mesher.read_mesh("t1.msh")
edgeNodes, edgeCells, cellNeighbors = Mesher.mesh_connectivity_2D(cellid)
cellType = Mesher.mesh_cell_type(cellNeighbors)
cellArea = Mesher.mesh_area_2D(nodes, cellid)
cellCenter = Mesher.mesh_center_2D(nodes, cellid)
edgeCenter = Mesher.mesh_edge_center(nodes, edgeNodes)
cellEdges = Mesher.mesh_cell_edge(cellid, edgeCells)

cellid, nodes = Mesher.read_mesh("square.msh")
edgeNodes, edgeCells, cellNeighbors = Mesher.mesh_connectivity_2D(cellid)
cellArea = Mesher.mesh_area_2D(nodes, cellid)
