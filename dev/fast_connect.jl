using FiniteMesh

cd(@__DIR__)
cells, points = read_mesh("../mesh/SU2/checkerboard.su2")

cellid = extract_cell(cells)
facePoints, faceCells = mesh_face_connectivity_2D(cellid)
