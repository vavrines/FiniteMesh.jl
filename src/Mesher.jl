module Mesher

using PyCall

export read_mesh,
       mesh_connectivity_2D,
       mesh_cell_type,
       mesh_center_2D,
       mesh_area_2D,
       mesh_edge_center,
       mesh_cell_edge

"""
    read_mesh(file::T) where {T<:AbstractString}

Read mesh file

* @return cells: node ids inside cells
* @return points: are saved with 3D coordinates (z=0 for 2D case)

"""
function read_mesh(file::T) where {T<:AbstractString}
    meshio = pyimport("meshio")

    py"""
    import meshio

    def read(file):
        m0 = meshio.read(file)

        points = m0.points
        celllist = m0.cells
        for cell in celllist:
            if not(cell[0] in ["line", "vertex"]):
                cells = cell[1] + 1

        return cells, points
    """

    m0 = meshio.read(file)

    return py"read"(file)
end


"""
    mesh_connectivity_2D(_cells::AbstractArray{<:Integer,2})

Compute connectivity of 2D unstructured mesh
"""
function mesh_connectivity_2D(cells::T) where {T<:AbstractArray{<:Integer,2}}
    nNodesPerCell = size(cells, 2)
    nCells = size(cells, 1)
    nEdgesMax = nNodesPerCell * nCells

    tmpEdgeNodes = -ones(Int, nEdgesMax, 2)
    tmpEdgeCells = -ones(Int, nEdgesMax, 2)

    counter = 0
    for i = 1:nCells, k = 1:nNodesPerCell
        isNewEdge = true
        for j = 1:counter
            if tmpEdgeNodes[j, :] == [cells[i, k], cells[i, k%nNodesPerCell+1]] ||
               tmpEdgeNodes[j, :] == [cells[i, k%nNodesPerCell+1], cells[i, k]]
                isNewEdge = false
                tmpEdgeCells[j, 2] = i
            end
        end
        if isNewEdge
            counter += 1
            tmpEdgeNodes[counter, 1] = cells[i, k]
            tmpEdgeNodes[counter, 2] = cells[i, k%nNodesPerCell+1]
            tmpEdgeCells[counter, 1] = i
        end
    end

    nEdges = counter
    edgeNodes = tmpEdgeNodes[1:nEdges, :]
    edgeCells = tmpEdgeCells[1:nEdges, :]

    cellNeighbors = -ones(Int, nCells, nNodesPerCell)
    for i = 1:nCells, k = 1:nNodesPerCell, j = 1:nEdges
        if edgeNodes[j, 1] == cells[i, k] &&
           edgeNodes[j, 2] == cells[i, k%nNodesPerCell+1] ||
           edgeNodes[j, 1] == cells[i, k%nNodesPerCell+1] && edgeNodes[j, 2] == cells[i, k]
            if edgeCells[j, 1] != i && edgeCells[j, 2] == i
                cellNeighbors[i, k] = edgeCells[j, 1]
            elseif edgeCells[j, 1] == i && edgeCells[j, 2] != i
                cellNeighbors[i, k] = edgeCells[j, 2]
            else
                throw("wrong info in neighboring cells of edge")
            end
        end
    end

    return edgeNodes, edgeCells, cellNeighbors
end


"""
    mesh_cell_type(cellNeighbors::T) where {T<:AbstractArray{<:Integer,2}}

Compute types of elements
- 0: inner
- 1: boundary
"""
function mesh_cell_type(cellNeighbors::T) where {T<:AbstractArray{<:Integer,2}}
    cellid = zeros(eltype(cellNeighbors), size(cellNeighbors, 1))
    for i in axes(cellNeighbors, 1)
        if -1 in cellNeighbors[i, :]
            cellid[i] = 1
        end
    end

    return cellid
end


"""
    mesh_area_2D(nodes::AbstractArray{<:AbstractFloat,2}, cells::AbstractArray{<:Int,2})

Compute areas of 2D elements
"""
function mesh_area_2D(
    nodes::X,
    cells::Y,
) where {X<:AbstractArray{<:AbstractFloat,2},Y<:AbstractArray{<:Integer,2}}

    ΔS = zeros(size(cells, 1))

    if size(cells, 2) == 3 # triangular mesh
        for i in eachindex(ΔS)
            ΔS[i] = abs(
                (
                    nodes[cells[i, 1], 1] *
                    (nodes[cells[i, 2], 2] - nodes[cells[i, 3], 2]) +
                    nodes[cells[i, 2], 1] *
                    (nodes[cells[i, 3], 2] - nodes[cells[i, 1], 2]) +
                    nodes[cells[i, 3], 1] * (nodes[cells[i, 1], 2] - nodes[cells[i, 2], 2])
                ) / 2,
            )
        end
    elseif size(cells, 2) == 4 # quadrilateral mesh
        for i in eachindex(ΔS)
            d1 = [
                nodes[cells[i, 1], 1] - nodes[cells[i, 2], 1],
                nodes[cells[i, 1], 2] - nodes[cells[i, 2], 2],
            ]
            d2 = [
                nodes[cells[i, 2], 1] - nodes[cells[i, 3], 1],
                nodes[cells[i, 2], 2] - nodes[cells[i, 3], 2],
            ]
            d3 = [
                nodes[cells[i, 3], 1] - nodes[cells[i, 4], 1],
                nodes[cells[i, 3], 2] - nodes[cells[i, 4], 2],
            ]
            d4 = [
                nodes[cells[i, 4], 1] - nodes[cells[i, 1], 1],
                nodes[cells[i, 4], 2] - nodes[cells[i, 1], 2],
            ]

            a = sqrt(d1[1]^2 + d1[2]^2)
            b = sqrt(d2[1]^2 + d2[2]^2)
            c = sqrt(d3[1]^2 + d3[2]^2)
            d = sqrt(d4[1]^2 + d4[2]^2)
            T = 0.5 * (a + b + c + d)

            alpha = acos((d4[1] * d1[1] + d4[2] * d1[2]) / (a * d))
            beta = acos((d2[1] * d3[1] + d2[2] * d3[2]) / (b * c))

            ΔS[i] = sqrt(
                (T - a) * (T - b) * (T - c) * (T - d) -
                a * b * c * d * cos(0.5 * (alpha + beta)) * cos(0.5 * (alpha + beta)),
            )
        end
    end

    return ΔS

end


"""
    mesh_center_2D(nodes::AbstractArray{<:AbstractFloat,2}, cells::AbstractArray{<:Integer,2})

Compute central points of 2D elements
"""
function mesh_center_2D(
    nodes::X,
    cells::Y,
) where {X<:AbstractArray{<:AbstractFloat,2},Y<:AbstractArray{<:Integer,2}}

    cellMidPoints = zeros(size(cells, 1), 2)
    for i in axes(cellMidPoints, 1) # nCells
        for j in axes(cells, 2) # nNodesPerCell
            cellMidPoints[i, :] .+= nodes[cells[i, j], 1:2]
        end
    end
    cellMidPoints ./= size(cells, 2)

    return cellMidPoints

end


"""
    mesh_edge_center(nodes::AbstractArray{<:AbstractFloat,2}, edgeNodes::AbstractArray{<:Integer,2})

Compute central points of cell edges
"""
function mesh_edge_center(
    nodes::X,
    edgeNodes::Y,
) where {X<:AbstractArray{<:AbstractFloat,2},Y<:AbstractArray{<:Integer,2}}
    
    edgeCenter = zeros(size(edgeNodes, 1), size(nodes, 2))
    for i in axes(edgeCenter, 1)
        id1 = edgeNodes[i, 1]
        id2 = edgeNodes[i, 2]
        @. edgeCenter[i, :] = 0.5 * (nodes[id1, :] + nodes[id2, :])
    end

    return edgeCenter

end


"""
    mesh_edge_center(nodes::AbstractArray{<:AbstractFloat,2}, edgeNodes::AbstractArray{<:Integer,2})

Compute central points of cell edges
"""
function mesh_cell_edge(cells::X, edgeCells::Y) where {X<:AbstractArray{<:Integer,2},Y<:AbstractArray{<:Integer,2}}
    ncell = size(cells, 1)
    vv = [Int[] for i in 1:ncell]
    for i in axes(edgeCells, 1)
        if edgeCells[i, 1] != -1
            push!(vv[edgeCells[i, 1]], i)
        end
        if edgeCells[i, 2] != -1
            push!(vv[edgeCells[i, 2]], i)
        end
    end

    cellEdges = zero(cells)
    for i in axes(cellEdges, 1)
        cellEdges[i, :] .= vv[i]
    end

    return cellEdges
end

end # module
