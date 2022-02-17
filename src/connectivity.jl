"""
$(SIGNATURES)

Compute connectivity of 2D unstructured mesh
"""
function mesh_connectivity_2D(cells, points::AbstractMatrix{T}) where {T<:Real}
    cellid = extract_cell(cells)

    facePoints, faceCells = mesh_face_connectivity_2D(cellid)
    cellNeighbors = mesh_cell_neighbor_2D(cellid, facePoints, faceCells)
    cellFaces = mesh_cell_face(cellid, faceCells)
    cellType = mesh_cell_type(cellNeighbors)
    cellArea = mesh_cell_area_2D(points, cellid)
    cellCenter = mesh_cell_center(points, cellid)
    cellNormals = mesh_cell_normals_2D(points, cellid, cellCenter)
    faceCenter = mesh_face_center(points, facePoints)
    faceType = mesh_face_type(faceCells, cellType)
    faceArea = mesh_face_area_2D(points, facePoints)

    return cellid,
    cellType,
    cellNeighbors,
    cellFaces,
    cellCenter,
    cellArea,
    cellNormals,
    facePoints,
    faceCells,
    faceCenter,
    faceType,
    faceArea
end

# ------------------------------------------------------------
# Face-indexed connectivity
# ------------------------------------------------------------

"""
$(SIGNATURES)

Compute affliated points and neighbor cells of faces
"""
function mesh_face_connectivity_2D(cells::AbstractMatrix{T}) where {T<:Integer}
    nCells = size(cells, 1)
    nNodesPerCell = size(cells, 2)
    nEdgesMax = nCells * nNodesPerCell

    tmpEdgeNodes = -ones(Int, nEdgesMax, 2)
    tmpEdgeCells = -ones(Int, nEdgesMax, 2)

    counter = 0
    @inbounds for i = 1:nCells
        for k = 1:nNodesPerCell
            endpoints = [cells[i, k], cells[i, k%nNodesPerCell+1]]
            isNewEdge = true
            for j = 1:counter
                if tmpEdgeNodes[j, 1] ∈ endpoints && tmpEdgeNodes[j, 2] ∈ endpoints
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
    end

    nEdges = counter
    edgeNodes = tmpEdgeNodes[1:nEdges, :]
    edgeCells = tmpEdgeCells[1:nEdges, :]

    return edgeNodes, edgeCells
end


"""
$(SIGNATURES)

Compute central points of cell faces
"""
function mesh_face_center(
    nodes::AbstractMatrix{T1},
    edgeNodes::AbstractMatrix{T2},
) where {T1<:Real,T2<:Integer}
    edgeCenter = zeros(size(edgeNodes, 1), size(nodes, 2))
    @inbounds for i in axes(edgeCenter, 1)
        pids = edgeNodes[i, :]

        tmp = zeros(size(nodes, 2))
        for i in eachindex(pids)
            tmp .+= nodes[pids[i], :]
        end
        tmp ./= length(pids)

        @. edgeCenter[i, :] = tmp
    end

    return edgeCenter
end

"""
$(SIGNATURES)

Compute type of faces

- 0: inner
- 1: boundary
"""
function mesh_face_type(
    edgeCells::AbstractMatrix{T},
    cellType::AbstractVector{T},
) where {T<:Integer}
    edgeType = zeros(eltype(edgeCells), size(edgeCells, 1))

    @inbounds for i in axes(edgeCells, 1)
        i1 = edgeCells[i, 1]
        i2 = edgeCells[i, 2]
        if cellType[i1] == 0 && cellType[i2] == 0
            edgeType[i] = 0
        else
            edgeType[i] = 1
        end
    end

    return edgeType
end


"""
$(SIGNATURES)

Compute area of faces
"""
function mesh_face_area_2D(points, facePoints)
    faceArea = zeros(size(facePoints, 1))

    @inbounds for i in eachindex(faceArea)
        faceArea[i] = norm(points[facePoints[i, 1], 1:2] .- points[facePoints[i, 2], 1:2])
    end

    return faceArea
end

# ------------------------------------------------------------
# Cell-indexed connectivity
# ------------------------------------------------------------

"""
$(SIGNATURES)

Compute neighbor cells of cells
"""
function mesh_cell_neighbor_2D(cells::AbstractMatrix{T}) where {T<:Integer}
    nCells = size(cells, 1)
    nNodesPerCell = size(cells, 2)

    prog = Progress(nCells, 5, "Computing mesh connectivity...", 50)
    cellNeighbors = -ones(Int, nCells, nNodesPerCell)
    @inbounds @threads for i = 1:nCells
        for k = 1:nNodesPerCell
            for j = 1:nCells
                if j == i
                    continue
                end
                
                if length(intersect(cells[j, :], [cells[i, k], cells[i, k%nNodesPerCell+1]])) == 2
                    cellNeighbors[i, k] = j
                    break
                end
            end
        end
        next!(prog)
    end

    #=cellNeighbors = -ones(Int, nCells, nNodesPerCell)
    @showprogress for i = 1:nCells, k = 1:nNodesPerCell
        for j = 1:nEdges
            if length(intersect(edgeNodes[j, :], [cells[i, k], cells[i, k%nNodesPerCell+1]])) ==
            2 # shared face
                idx = findall(x -> x != i, edgeCells[j, :]) |> first # idx takes value 1 or 2
                cellNeighbors[i, k] = edgeCells[j, idx]
                break
            end
        end
    end=#

    # cellNeighbors are already regularized, i.e.,
    # face 1 -> neighbor cell 1
    # face 2 -> neighbor cell 2
    # face 3 -> neighbor cell 3
    # so we can get rid of this part
    #=
    cellNeighbors_regular = -ones(Int, nCells, nNodesPerCell)
    for i = 1:nCells
        cids = cellNeighbors[i, :]

        fpids = [cells[i, j:j+1] for j = 1:nNodesPerCell-1]
        push!(fpids, [cells[i, nNodesPerCell], cells[i, 1]])

        for j = 1:nNodesPerCell
            fpid = fpids[j]

            for k in eachindex(cids)
                if cids[k] != -1 && length(intersect(fpid, cells[cids[k], :])) == 2
                    cellNeighbors_regular[i, j] = cids[k]
                end
            end
        end
    end
    =#

    return cellNeighbors
end

function mesh_cell_neighbor_2D(
    cells::AbstractMatrix{T},
    edgeNodes::AbstractMatrix{T},
    edgeCells::AbstractMatrix{T},
) where {T<:Integer}
    nCells = size(cells, 1)
    nNodesPerCell = size(cells, 2)
    nEdges = size(edgeNodes, 1)

    prog = Progress(nCells, 5, "Computing mesh connectivity...", 50)
    cellNeighbors = -ones(Int, nCells, nNodesPerCell)
    @inbounds for i = 1:nCells
        for k = 1:nNodesPerCell
            for j = 1:nEdges
                if length(intersect(edgeNodes[j, :], [cells[i, k], cells[i, k%nNodesPerCell+1]])) ==
                2 # shared face
                    idx = findall(x -> x != i, edgeCells[j, :]) |> first # idx takes value 1 or 2
                    cellNeighbors[i, k] = edgeCells[j, idx]
                    break
                end
            end
        end
        next!(prog)
    end

    # cellNeighbors are already regularized, i.e.,
    # face 1 -> neighbor cell 1
    # face 2 -> neighbor cell 2
    # face 3 -> neighbor cell 3
    # so we can get rid of this part
    #=
    cellNeighbors_regular = -ones(Int, nCells, nNodesPerCell)
    for i = 1:nCells
        cids = cellNeighbors[i, :]

        fpids = [cells[i, j:j+1] for j = 1:nNodesPerCell-1]
        push!(fpids, [cells[i, nNodesPerCell], cells[i, 1]])

        for j = 1:nNodesPerCell
            fpid = fpids[j]

            for k in eachindex(cids)
                if cids[k] != -1 && length(intersect(fpid, cells[cids[k], :])) == 2
                    cellNeighbors_regular[i, j] = cids[k]
                end
            end
        end
    end
    =#

    return cellNeighbors
end


"""
$(SIGNATURES)

Compute types of elements

- 0: inner
- 1: boundary
"""
function mesh_cell_type(cellNeighbors::AbstractMatrix{T}) where {T<:Integer}
    cellid = zeros(eltype(cellNeighbors), size(cellNeighbors, 1))
    @inbounds for i in axes(cellNeighbors, 1)
        if -1 in cellNeighbors[i, :]
            cellid[i] = 1
        end
    end

    return cellid
end


"""
$(SIGNATURES)

Compute areas of 2D elements
"""
function mesh_cell_area_2D(
    nodes::AbstractMatrix{T1},
    cells::AbstractMatrix{T2},
) where {T1<:Real,T2<:Integer}

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
$(SIGNATURES)

Compute central points of elements
"""
function mesh_cell_center(
    nodes::AbstractMatrix{T1},
    cells::AbstractMatrix{T2},
) where {T1<:AbstractFloat,T2<:Integer}

    cellMidPoints = zeros(size(cells, 1), size(nodes, 2))
    @inbounds for i in axes(cellMidPoints, 1) # nCells
        for j in axes(cells, 2) # nNodesPerCell
            cellMidPoints[i, :] .+= nodes[cells[i, j], :]
        end
    end
    cellMidPoints ./= size(cells, 2)

    return cellMidPoints

end


"""
$(SIGNATURES)

Compute surrounding faces of cell
"""
function mesh_cell_face(
    cells::AbstractMatrix{T},
    edgeCells::AbstractMatrix{T},
) where {T<:Integer}

    ncell = size(cells, 1)
    vv = [Int[] for i = 1:ncell]
    @inbounds for i in axes(edgeCells, 1)
        if edgeCells[i, 1] != -1
            push!(vv[edgeCells[i, 1]], i)
        end
        if edgeCells[i, 2] != -1
            push!(vv[edgeCells[i, 2]], i)
        end
    end

    cellEdges = zero(cells)
    @inbounds for i in axes(cellEdges, 1)
        cellEdges[i, :] .= vv[i]
    end

    return cellEdges

end


"""
$(SIGNATURES)

Compute unit normal vectors of cells
"""
function mesh_cell_normals_2D(points, cells, cellCenter)
    ncell = size(cells, 1)
    np = size(cells, 2)
    cell_normal = zeros(ncell, np, 2)

    @inbounds for i = 1:ncell
        pids = [cells[i, j:j+1] for j = 1:np-1]
        push!(pids, [cells[i, np], cells[i, 1]])

        for j = 1:np
            cell_normal[i, j, :] .=
                unit_normal(points[pids[j][1], :], points[pids[j][2], :])
        end

        p = [(points[pids[j][1], :] .+ points[pids[j][2], :]) ./ 2 for j = 1:np]

        for j = 1:np
            if dot(cell_normal[i, j, :], p[j][1:2] - cellCenter[i, 1:2]) < 0
                cell_normal[i, j, :] .= -cell_normal[i, j, :]
            end
        end
    end

    return cell_normal
end
