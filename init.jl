# Same as FLatticeExact, but get_odometers returns matrices instead of dictionaries.
if !isdefined(@__MODULE__, :FLatticeExact)
    include(joinpath(@__DIR__,"src","FLatticeExact.jl"))
end

getFracs(word) = FLatticeExact.getFracs(word)

function _odometer_matrix(tile)
    xmin,xmax = extrema(real,keys(tile))
    ymin,ymax = extrema(imag,keys(tile))
    result = Matrix{Union{BigInt,Float64}}(undef,Int(xmax-xmin+3),Int(ymax-ymin+3))
    fill!(result,NaN)
    for (z,h) in tile
        result[Int(real(z)-xmin+2),Int(imag(z)-ymin+2)] = h
    end
    result
end

"""Return the standard and alternate tiles as matrices. Points outside the tile are NaN,
and there is a border of NaN around each matrix."""
function get_odometers(word)
    standards,alternates = FLatticeExact.get_odometers(word)
    map(_odometer_matrix,standards),map(_odometer_matrix,alternates)
end
