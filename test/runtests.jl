using Test, FiniteMesh

cd(@__DIR__)
cells, nodes = read_mesh("../mesh/Gmsh/t1.msh")
mesh_connectivity_2D(cells, nodes)

cells, nodes = read_mesh("../mesh/Gmsh/square.msh")
cellid = extract_cell(cells)
edgeNodes, edgeCells = mesh_face_connectivity_2D(cellid)
cellArea = mesh_area_2D(nodes, cellid)

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
