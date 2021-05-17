unstructured_index(i::T, j::T, nx::T) where {T<:Integer} = (j-1) * nx + i

function unstructured_grid(coords::T, corners::T) where {T<:AbstractArray{<:Real,3}}
    nx = size(coords, 1)
    ny = size(coords, 2)
    ne = nx * ny
    np = (nx+1) * (ny+1)

    points = zeros(np, 2)
    for i = 1:nx+1, j = 1:ny+1
        idx = unstructured_index(i, j, nx+1)
        points[idx, :] .= corners[i, j, :]
    end

    cells = zeros(Int, ne, 4)
    for i = 1:nx, j = 1:ny
        idx = unstructured_index(i, j, nx)
        cells[idx, 1] = unstructured_index(i, j, nx+1)
        cells[idx, 2] = unstructured_index(i+1, j, nx+1)
        cells[idx, 3] = unstructured_index(i+1, j+1, nx+1)
        cells[idx, 4] = unstructured_index(i, j+1, nx+1)
    end

    return cells, points
end

function unstructured_grid(x::T, y::T, x0, x1, y0, y1) where {T<:AbstractMatrix}
    Δx = zero(x[:, 1])
    Δx[1] = 2.0 * (x[1, 1] - x0)
    for i = 2:length(Δx)-1
        Δx[i] = 2.0 * (x[i, 1] - x[i-1, 1] - Δx[i-1])
    end
    Δx[end] = 2.0 * (x1 - x[end, 1])

    Δy = zero(y[1, :])
    Δy[1] = 2.0 * (y[1, 1] - y0)
    for i = 2:length(Δy)-1
        Δy[i] = 2.0 * (y[1, i] - y[1, i-1] - Δy[i-1])
    end
    Δy[end] = 2.0 * (y1 - y[1, end])

    points = zeros(length(Δx)+1, length(Δy)+1, 2)
    for i in axes(points, 1), j in axes(points, 2)
        if i <= length(Δx) && j <= length(Δy)
            points[i, j, 1] = x[i, j] - Δx[i]/2
            points[i, j, 2] = y[i, j] - Δy[j]/2
        elseif i == length(Δx)+1 && j <= length(Δy)
            points[i, j, 1] = x[i-1, j] + Δx[i-1]/2
            points[i, j, 2] = y[i-1, j] - Δy[j]/2
        elseif i <= length(Δx) && j == length(Δy)+1
            points[i, j, 1] = x[i, j-1] - Δx[i]/2
            points[i, j, 2] = y[i, j-1] + Δy[j-1]/2
        elseif i == length(Δx)+1 && j == length(Δy)+1
            points[i, j, 1] = x[i-1, j-1] + Δx[i-1]/2
            points[i, j, 2] = y[i-1, j-1] + Δy[j-1]/2
        end
    end

    unstructured_grid(cat(x, y, dims=3), points)
end

function unstructured_grid(coords::T, xrange, yrange) where {T<:AbstractArray{<:Real,3}}
    xcoords = coords[:, :, 1]
    ycoords = coords[:, :, 2]
    x0, x1 = xrange
    y0, y1 = yrange

    unstructured_grid(xcoords, ycoords, x0, x1, y0, y1)
end
