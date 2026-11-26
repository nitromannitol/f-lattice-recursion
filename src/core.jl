# A tile is a dictionary from lattice points to integer heights.
# Points outside the tile are not in the dictionary.
const Site = Complex{BigInt}
const Tile = Dict{Site,BigInt}
pt(x::Integer, y::Integer) = complex(BigInt(x), BigInt(y))
_exactint(x::Integer) = BigInt(x)
function _exactint(x::Rational)
    denominator(x) == 1 || throw(ArgumentError("nonintegral lattice component: $x"))
    return BigInt(numerator(x))
end
_site(z::Number) = pt(_exactint(real(z)), _exactint(imag(z)))
function _require_nonempty(t::Tile, which::AbstractString)
    isempty(t) && throw(ArgumentError("empty tile: $which"))
    return t
end
# ll is the leftmost point of the bottom row, ur is the rightmost point of the top row.
function ll(t::Tile)
    _require_nonempty(t, "LL")
    return argmin(z -> (imag(z),real(z)), keys(t))
end
function ur(t::Tile)
    _require_nonempty(t, "UR")
    return argmax(z -> (imag(z),real(z)), keys(t))
end
"""Shift a tile by v and add the linear function with the given slope, evaluated at the new position."""
function place(t::Tile, v, slope)
    _require_nonempty(t,"place")
    w = _site(v); a = _site(slope)
    out = Tile(); sizehint!(out,length(t))
    if iszero(a)
        for (z,h) in t
            out[z+w] = h
        end
    else
        ax,ay = real(a),imag(a)
        for (z,h) in t
            u=z+w
            out[u] = h + ax*real(u) + ay*imag(u)
        end
    end
    return out
end
"""Shift a tile so that its smallest x and smallest y coordinates are 0."""
function trim(t::Tile)
    _require_nonempty(t,"trim")
    return place(t,-pt(minimum(real,keys(t)), minimum(imag,keys(t))),0)
end
"""Shift a tile so that the leftmost point of its bottom row is at 0."""
llcenter(t::Tile) = place(t,-ll(t),0)
function maxzero(t::Tile)
    _require_nonempty(t,"maxzero")
    m=maximum(values(t))
    return Tile(z => h-m for (z,h) in t)
end
"""Shift the incoming tile up or down so that it touches the old tile from above, then take the
minimum of the two. The tiles must overlap."""
function touch(old::Tile, incoming::Tile; context="")
    _require_nonempty(old,"old in $context"); _require_nonempty(incoming,"incoming in $context")
    delta = nothing
    small = length(old) < length(incoming) ? old : incoming
    for z in keys(small)
        if haskey(old,z) && haskey(incoming,z)
            v=incoming[z]-old[z]
            if isnothing(delta) || v<delta
                delta=v
            end
        end
    end
    isnothing(delta) && throw(ArgumentError("tiles do not overlap: $context"))
    # After the shift the old tile is the minimum on the overlap, so only new points are added.
    out=copy(old)
    for (z,h) in incoming
        if !haskey(out,z)
            out[z]=h-delta
        end
    end
    return out
end
"""Write a tile as lines of x, y and height, with the leftmost point of the bottom row at (0,0)."""
function write_tile(io::IO,t::Tile)
    a=ll(t)
    println(io,"# x\ty\theight")
    for z in sort!(collect(keys(t));by=z->(real(z),imag(z)))
        u=z-a
        println(io,real(u),'\t',imag(u),'\t',t[z])
    end
    return nothing
end
