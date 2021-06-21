"""
Struct of mesh information

"""
struct Mesh{
    A,
    B<:AbstractMatrix{<:AbstractFloat},
    C<:AbstractMatrix{<:Integer},
    D<:AbstractVector{<:Integer},
    E<:AbstractVector{<:AbstractFloat},
    F<:AbstractArray{<:AbstractFloat,3},
}
    cells::A # all information: cell, line, vertex
    points::B # locations of vertices
    
    cellid::C # node indices of elements
    cellType::D # inner/boundary cell
    cellNeighbors::C # neighboring cells id
    cellFaces::C # cell edges id
    cellCenter::B # cell center location
    cellArea::E # cell size
    cellNormals::F # unit normal vectors of cell

    facePoints::C # point ids affiliated to face
    faceCells::C # ids of neighbor cells
    faceCenter::B # center location of face
    faceType::D
    faceArea::E # face area
end

function Mesh(file::T) where T<:AbstractString
    cells, points = read_mesh(file)
    p = mesh_connectivity_2D(cells, points)

    return Mesh((cells, points)..., p...)
end


"""
    struct Cells{T1,T2}
        type::T1
        index::T2
    end

Cell connectivity information
"""
struct Cells{T1,T2}
    type::T1
    index::T2
end


"""
    extract_cell(cells::Cells)

Extract cell id list from struct
"""
function extract_cell(cells::Cells)
    for i in eachindex(cells.type)
        if (cells.type[i] in ["triangle", "quad"])
            return cells.index[i]
        end
    end
end
