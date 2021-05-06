import FiniteMesh

lines = readlines("../mesh/SU2/naca0012.su2")
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

split(tags[1][2], "\t")









for i = 1:3
    @show i
    if i==2
        break
    end
end


for line in lines
    if line[1:10] == "MARKER_TAG"
    end

end

parse(Float64, lines[3][1] |> string)


parse(Float64, "n")




split(lines[3], "\t")