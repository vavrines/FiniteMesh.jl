using Test, FiniteMesh

cd(@__DIR__)
cells, nodes = read_mesh("t1.msh")
cellid = extract_cell(cells)
edgeNodes, edgeCells, cellNeighbors = mesh_connectivity_2D(cellid)
cellType = mesh_cell_type(cellNeighbors)
cellArea = mesh_area_2D(nodes, cellid)
cellCenter = mesh_center_2D(nodes, cellid)
edgeCenter = mesh_face_center(nodes, edgeNodes)
cellEdges = mesh_cell_face(cellid, edgeCells)
edgeType = mesh_face_type(edgeCells, cellType)

cells, nodes = read_mesh("square.msh")
cellid = extract_cell(cells)
edgeNodes, edgeCells, cellNeighbors = mesh_connectivity_2D(cellid)
cellArea = mesh_area_2D(nodes, cellid)
