"""
$(SIGNATURES)

check if cellNeighbors are regularized, e.g., in a triangle,
    
- face 1 -> neighbor cell 1
- face 2 -> neighbor cell 2
- face 3 -> neighbor cell 3
"""
function check_cell_neighbor(cells, cellNeighbors)
    nCells = size(cells, 1)
    nNodesPerCell = size(cells, 2)

    for i in 1:nCells
        cids = cellNeighbors[i, :]

        fpids = [cells[i, j:j+1] for j in 1:nNodesPerCell-1]
        push!(fpids, [cells[i, nNodesPerCell], cells[i, 1]])

        for j in 1:nNodesPerCell
            nid = cids[j]
            if nid != -1
                @assert length(intersect(fpids[j], cells[nid, :])) == 2 "neighbor is not regularized at ($i, $j)"
            end
        end
    end

    return true
end

"""
$(SIGNATURES)

check if cellFaces are regularized, e.g., in a triangle,
    
- face 1 -> points 1 & 2
- face 2 -> points 2 & 3
- face 3 -> points 3 & 1
"""
function check_cell_face(cells, cellFaces, facePoints)
    nCells = size(cells, 1)
    nNodesPerCell = size(cells, 2)

    for i in 1:nCells
        fpids = [cells[i, j:j+1] for j in 1:nNodesPerCell-1]
        push!(fpids, [cells[i, nNodesPerCell], cells[i, 1]])

        for j in 1:nNodesPerCell
            @assert sort(facePoints[cellFaces[i, j], :]) == sort(fpids[j]) "face is not regularized at ($i, $j)"
        end
    end

    return true
end
