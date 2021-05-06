function add_group!(cells::Cells, file::AbstractString)
    func = file[end-2:end] * "_group!"
    eval(Symbol(func))(cells, file)
end

function su2_group!(cells::Cells, file::AbstractString)
    lines = readlines(file)
    ids = Int[]
    for i in 1:length(lines)-1
        if lines[i][1:4] * lines[i+1][1:4] == "MARKMARK"
            push!(ids, i)
        end
    end

    tags = [[lines[id]] for id in ids]
    for idx in eachindex(tags)
        for i = ids[idx]+2:length(lines)
            lines[i][1:4] == "MARK" && break
            push!(tags[idx], lines[i])
        end
    end

    for idx in eachindex(tags)
        key = split(tags[idx][1], "= ")[2] |> String
        
        value = zeros(Int, length(tags[idx])-1, length(split(tags[idx][2], "\t"))-1)
        for i in axes(value, 1)
            value[i, :] = parse.(Int, split(tags[idx][i+1], "\t")[2:end]) .+ 1
        end

        push!(cells.type, key)
        push!(cells.index, value)
    end
end
