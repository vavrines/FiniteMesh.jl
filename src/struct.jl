"""
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
