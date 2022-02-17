"""
$(SIGNATURES)

regularize neighboring cell IDs, e.g., in a triangle,
    
- face 1 -> neighbor cell 1
- face 2 -> neighbor cell 2
- face 3 -> neighbor cell 3
"""
function regularize_cell_neighbor(cells, cellNeighbors)
    nCells = size(cells, 1)
    nNodesPerCell = size(cells, 2)

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

    return cellNeighbors_regular
end


"""
$(SIGNATURES)

regularize affiliated face IDs, e.g., in a triangle,
    
- face 1 -> points 1 & 2
- face 2 -> points 2 & 3
- face 3 -> points 3 & 1
"""
function regularize_cell_face(cells, cellFaces, facePoints)
    nCells = size(cells, 1)
    nNodesPerCell = size(cells, 2)

    cellFacesReg = zeros(Int, nCells, nNodesPerCell)
    for i = 1:nCells
        fids = cellFaces[i, :]

        pids = [cells[i, j:j+1] for j = 1:nNodesPerCell-1]
        push!(pids, [cells[i, nNodesPerCell], cells[i, 1]])

        for j = 1:nNodesPerCell
            pid = pids[j]

            for fid in fids
                if sort(pid) == sort(facePoints[fid, :])
                    cellFacesReg[i, j] = fid
                end
            end
        end
    end

    return cellFacesReg
end
