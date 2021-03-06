module FiniteMesh

using Base.Threads: @threads
using DocStringExtensions
using ProgressMeter
using PyCall
using LinearAlgebra

export FM
export read_mesh
export Mesh, Cells, extract_cell
export mesh_connectivity_2D,
       mesh_face_connectivity_2D,
       mesh_face_center,
       mesh_face_type,
       mesh_face_area_2D,
       mesh_cell_neighbor_2D,
       mesh_cell_type,
       mesh_cell_center,
       mesh_cell_area_2D,
       mesh_cell_face,
       mesh_cell_normals_2D
export add_group!, su2_group!
export unstructured_index, unstructured_grid, triangulate
export unit_normal
export regularize_cell_neighbor, regularize_cell_face

const FM = FiniteMesh

include("tools.jl")
include("struct.jl")
include("connectivity.jl")
include("group.jl")
include("transform.jl")
include("regularize.jl")
include("check.jl")

"""
$(SIGNATURES)

Read mesh file

* @return cells: node ids inside cells
* @return points: are saved with 3D coordinates (z=0 for 2D case)
"""
function read_mesh(file::T) where {T<:AbstractString}
    py"""
    import meshio
    def read(file):
        m0 = meshio.read(file)
        points = m0.points
        cells = m0.cells
        keys = []
        vals = []
        for cell in cells:
            keys.append(cell.type)
            vals.append(cell.data)
        return points, keys, vals
    """

    points, keys, vals = py"read"(file)
    for val in vals
        val .+= 1 # python index is zero-based
    end

    cells = Cells(keys, vals)

    return cells, points
end

end # module
