import FiniteMesh

cells, points = FiniteMesh.read_mesh("../mesh/SU2/naca0012.su2")

FiniteMesh.add_group!(cells, "../mesh/SU2/naca0012.su2")

cells, points = FiniteMesh.read_mesh("../test/square.msh")
