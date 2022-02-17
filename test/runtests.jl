using Test, FiniteMesh

cd(@__DIR__)
mesh = Mesh("../mesh/Gmsh/t1.msh")

cells, nodes = read_mesh("../mesh/Gmsh/square.msh")
cellid = extract_cell(cells)
edgeNodes, edgeCells = mesh_face_connectivity_2D(cellid)
cellArea = mesh_cell_area_2D(nodes, cellid)
mesh_cell_neighbor_2D(cellid, edgeNodes, edgeCells)

cells, nodes = read_mesh("../mesh/SU2/naca0012.su2")
add_group!(cells, "../mesh/SU2/naca0012.su2")

rg = [0.1, 0.3, 0.5, 0.7, 0.9]
x = hcat(rg, rg, rg)
rg = [0.1 0.3 0.5]
y = vcat(rg, rg, rg, rg, rg)

cells, points = unstructured_grid(cat(x, y, dims = 3), [0.0, 1.0], [0.0, 0.6])

_points = zeros(4, 2)
_points[2, :] = [1.0, 0.0]
_points[3, :] = [1.0, 1.0]
_points[4, :] = [0.0, 1.0]

cells = triangulate(_points)

FiniteMesh.unit_normal(rand(3), rand(3), rand(3))

cellNeighborsReg = regularize_cell_neighbor(mesh.cellid, mesh.cellNeighbors)
FiniteMesh.check_cell_neighbor(mesh.cellid, cellNeighborsReg)

cellFacesReg = regularize_cell_face(mesh.cellid, mesh.cellFaces, mesh.facePoints)
FiniteMesh.check_cell_face(mesh.cellid, cellFacesReg, mesh.facePoints)
