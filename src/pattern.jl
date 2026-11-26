"""The periodic extension of a tile: one period of heights, the two period constants and the Laplacian."""
struct PeriodicPattern
    u::Site
    v::Site
    gradient::Site
    area::BigInt
    cu::BigInt
    cv::BigInt
    heights::Tile
    labels::Dict{Site,Int8}
    horizontal_parity::Int
end

_cross(z,w) = real(z)*imag(w)-imag(z)*real(w)
_dot(z,w) = real(z)*real(w)+imag(z)*imag(w)

# z = r + i*u + j*v, with r in the half-open fundamental parallelogram.
function _period_coordinates(z,u,v,area)
    i = fld(_cross(z,v),area)
    j = fld(_cross(u,z),area)
    (z-i*u-j*v, i, j)
end

function _period_term(r,j,gradient,v)
    j*_dot(gradient,r) + div(j*(j-1),2)*_dot(gradient,v)
end

function _exact_character(e1,e2)
    a,b,c = e1
    d,e,f = e2
    den = a*e-b*d
    x,rx = divrem(c*e-b*f,den)
    y,ry = divrem(a*f-c*d,den)
    iszero(rx) && iszero(ry) || throw(ArgumentError("period constants are not integers"))
    x,y
end

"""
    periodic_pattern(tile, fraction; horizontal_parity=nothing)

Extend a tile periodically and compute its Laplacian. `fraction` is the slope of
the tile. Throws an error if the tile does not cover a full period or if its
heights are not periodic. If both checkerboards give a valid Laplacian and the
results differ, pass horizontal_parity = 0 or 1 to pick one.
"""
function periodic_pattern(tile::Tile, fraction::Rational; horizontal_parity=nothing)
    0 < fraction < 1 || throw(ArgumentError("fraction must be between 0 and 1"))
    n,d = BigInt(numerator(fraction)),BigInt(denominator(fraction))
    u = (isodd(n+d) ? 2 : 1)*pt(d,n)
    v,gradient = pt(n-d,n+d),pt(n,-d)
    area = _cross(u,v)
    firsts = Dict{Site,NTuple{3,BigInt}}()
    e1 = constants = nothing
    for (z,h) in tile
        r,i,j = _period_coordinates(z,u,v,area)
        h0 = h-_period_term(r,j,gradient,v)
        if !haskey(firsts,r)
            firsts[r] = (i,j,h0)
        elseif isnothing(constants)
            i0,j0,g = firsts[r]
            eq = (i-i0,j-j0,h0-g)
            if isnothing(e1)
                e1 = eq
            elseif e1[1]*eq[2] != e1[2]*eq[1]
                constants = _exact_character(e1,eq)
            end
        end
    end
    length(firsts)==area || throw(ArgumentError("tile covers $(length(firsts)) of the $area points in a period"))
    isnothing(constants) && throw(ArgumentError("could not find both period constants"))
    cu,cv = constants
    heights = Tile(r=>g-i*cu-j*cv for (r,(i,j,g)) in firsts)
    function at(z)
        r,i,j = _period_coordinates(z,u,v,area)
        heights[r]+i*cu+j*cv+_period_term(r,j,gradient,v)
    end
    all(at(z)==h for (z,h) in tile) || throw(ArgumentError("heights are not periodic"))
    horizontal_parity in (nothing,0,1) || throw(ArgumentError("horizontal_parity must be 0 or 1"))
    phases = isnothing(horizontal_parity) ? (0,1) : (Int(horizontal_parity),)
    result = nothing
    for phase in phases
        labels = Dict{Site,Int8}()
        for (r,h) in heights
            step = iseven(real(r)+imag(r)+phase) ? 1 : im
            lap = at(r-step)+at(r+step)-2h
            lap in (-1,0) || break
            labels[r] = Int8(lap)
        end
        length(labels)==area || continue
        if !isnothing(result)
            result.labels==labels && continue
            throw(ArgumentError("both checkerboards work and give different answers; pass horizontal_parity"))
        end
        result = PeriodicPattern(u,v,gradient,area,cu,cv,heights,labels,phase)
    end
    isnothing(result) && throw(ArgumentError("Laplacian is not -1 or 0 everywhere"))
    result
end

"""Height of the periodic extension at the point z."""
function periodic_height(p::PeriodicPattern,z::Number)
    r,i,j = _period_coordinates(_site(z),p.u,p.v,p.area)
    p.heights[r]+i*p.cu+j*p.cv+_period_term(r,j,p.gradient,p.v)
end

"""Laplacian (-1 or 0) of the periodic extension at the point z."""
function pattern_value(p::PeriodicPattern,z::Number)
    r,_,_ = _period_coordinates(_site(z),p.u,p.v,p.area)
    p.labels[r]
end
