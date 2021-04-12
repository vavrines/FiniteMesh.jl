using Test, FiniteMesh

cd(@__DIR__)
cellid, nodes = read_mesh("t1.msh")
edgeNodes, edgeCells, cellNeighbors = mesh_connectivity_2D(cellid)
cellType = mesh_cell_type(cellNeighbors)
cellArea = mesh_area_2D(nodes, cellid)
cellCenter = mesh_center_2D(nodes, cellid)
edgeCenter = mesh_edge_center(nodes, edgeNodes)
cellEdges = mesh_cell_edge(cellid, edgeCells)

cellid, nodes = read_mesh("square.msh")
edgeNodes, edgeCells, cellNeighbors = mesh_connectivity_2D(cellid)
cellArea = mesh_area_2D(nodes, cellid)
