# Base case tiles: the empty word and words of all 2s or all 3s.

"One of the base case tiles, with size parameter d."
function _base_shape(kind::Symbol, d::Integer)
    d = big(d)
    staircase = kind in (:stair_even, :stair_odd, :stair_even_alt, :stair_odd_alt)
    ladder = kind in (:ladder_odd, :ladder_even, :ladder_odd_alt, :ladder_even_alt)
    ((staircase && isodd(d) && d >= 3) || (ladder && iseven(d) && d >= 4)) ||
        throw(ArgumentError("invalid shape or size"))
    T = Tile()
    if kind == :stair_even
        for x in (2-d):d, y in big(0):(d+1)
            z = -div(y*(y-1),2) + min(x-1,0)
            if 1 <= x+y <= d+2
                T[pt(x,y)] = z
            end
            if (x == 0 && y == 0) || (x == 2 && y == d+1)
                T[pt(x,y)] = z+1
            end
        end
    elseif kind == :stair_odd
        for x in big(0):(d+1), y in big(1):(2d-1)
            z = -fld((x-y)^2,4) + min(d-y,0)
            if -1 <= y-x <= d
                T[pt(x,y)] = z
            end
            if (x == 0 && y == d+1) || (x == d+1 && y == d-1)
                T[pt(x,y)] = z+1
            end
        end
    elseif kind == :stair_even_alt
        for x in (4-2d):d, y in big(-1):(d+1)
            z = -div(y*(y+1),2) + min(0,2-x) + min(0,d-2+x) + max(x+y-2,0)
            if 2-d <= x+y <= d+2 && 0 <= y <= d
                T[pt(x,y)] = z
            end
            if (y == -1 && 1-d <= x <= 0) || (y == d+1 && 4-d <= x <= 3)
                T[pt(x,y)] = z
            end
            if (x+y == 2 && 0 < y < d) ||
               (x == 1-d && (y == 0 || y == -1)) ||
               (x == 3 && (y == d || y == d+1))
                T[pt(x,y)] = z+1
            end
        end
    elseif kind == :stair_odd_alt
        for x in big(-1):(d+1), y in big(1):(3d-3)
            z = -fld((x-y)^2,4) + min(d-y-1,0) + min(2d-1-y,0) + max(y-x-d+1,0)
            if 0 <= x <= d && -1 <= y-x <= 2d-1
                T[pt(x,y)] = z
            end
            if y-x == d-1 && 0 < x < d
                T[pt(x,y)] = z+1
            end
            if (x == -1 && d+1 <= y <= 2d) || (x == d+1 && d-2 <= y <= 2d-3)
                T[pt(x,y)] = z
            end
            if (x == 0 && y == 2d) || (x == d && y == d-2) ||
               (x == -1 && y == 2d) || (x == d+1 && y == d-2)
                T[pt(x,y)] = z+1
            end
        end
    elseif kind == :ladder_odd
        for x in (2-d):(2d), y in big(0):(d+2)
            z = -div(y*(y-1),2) + min(0,x-1) + max(y+x-d-2,0) + min(d+1-x,0)
            if 1 <= x+y <= 2d+3 && !((y == 0 && x >= d) || (y == d+2 && x <= 2))
                T[pt(x,y)] = z
            end
            if (x == 0 && y == 0) || (x == d+2 && y == d+2) ||
               (x+y == d+2 && 1 < y < d+1)
                T[pt(x,y)] = z+1
            end
        end
    elseif kind == :ladder_even
        for x in big(-1):(d+1), y in big(1):(3d-1)
            z = -fld((x-y)^2,4) + min(0,d-y) + max(0,y-x-d) + min(0,2d-y)
            if -1 <= y-x <= 2d+1 && !((x == -1 && y <= d+1) || (x == d+1 && y >= 2d-1))
                T[pt(x,y)] = z
            end
            if (x == -1 && y == 2d+1) || (x == d+1 && y == d-1) ||
               (y-x == d && 0 < x < d)
                T[pt(x,y)] = z+1
            end
        end
    elseif kind == :ladder_odd_alt
        for x in (1-d):(2d-2), y in big(-1):(d+2)
            z1 = -div(y*(y+1),2)
            z2 = min(x,0)
            z3 = max(x+y-d,0)
            z = z1+z2+z3+min(d-x-1,0)
            if -1 <= x+y <= 2d+1 && !((-1 <= y <= 0 && x >= d) || (d+1 <= y <= d+2 && x <= -1))
                T[pt(x,y)] = z
                if x+y == d && 0 < x < d-1
                    # no min(d-x-1,0) term on this diagonal
                    T[pt(x,y)] = z1+z2+z3+1
                end
            end
            if (x == -1 && y == -1) || (x == d && y == d+2)
                T[pt(x,y)] = z+1
            end
        end
    elseif kind == :ladder_even_alt
        for x in big(-2):(d+1), y in big(1):(3d-2)
            z = -fld((x-y)^2,4) + min(0,d-y) + min(2d-1-y,0) + max(y-x-d,0)
            if -1 <= y-x <= 2d+1 && !((x <= -1 && y <= d-1) || (x >= d && y >= 2d))
                T[pt(x,y)] = z
                if y-x == d && 0 < x < d-1
                    T[pt(x,y)] = z+1
                end
            end
            if (x == -2 && y == 2d) || (x == d+1 && y == d-1)
                T[pt(x,y)] = z+1
            end
        end
    else
        throw(ArgumentError("unknown shape $kind"))
    end
    return trim(T)
end

"The small tile for the fraction 0 or 1, placed at the leftmost point of the bottom row of T."
function _base_endpoint(T::Tile, endpoint::Symbol)
    isempty(T) && throw(ArgumentError("empty tile"))
    z0 = ll(T)
    B = Tile()
    if endpoint == :zero
        for x in 0:1, y in 0:2
            B[z0+pt(x,y)] = big(y == 2 ? -1 : 0)
        end
    elseif endpoint == :one
        for x in 0:1, y in 0:1
            B[z0+pt(x,y)] = big(x == 0 && y == 1 ? -1 : 0)
        end
    else
        throw(ArgumentError("endpoint must be :zero or :one"))
    end
    return trim(B)
end

"""
    base_quad(word)

Return the standard and alternate tiles for the empty word or a word of all 2s
or all 3s, in the order odd, even, odd parent, even parent.
"""
function base_quad(word::AbstractString)
    all(c -> c == '2' || c == '3', word) || throw(ArgumentError("word must be empty, all 2s or all 3s"))
    k = big(length(word))
    if isempty(word)
        p = _base_shape(:stair_odd,3)
        q = _base_shape(:stair_even,3)
        p0,q0 = _base_endpoint(p,:zero),_base_endpoint(q,:one)
        return ((p,q,p0,q0),(_base_shape(:stair_odd_alt,3),_base_shape(:stair_even_alt,3),copy(p0),copy(q0)))
    elseif all(==('2'),word)
        p = _base_shape(:stair_odd,2k+3)
        q = _base_shape(:ladder_even,2k+2)
        p0,q0 = _base_shape(:stair_odd,2k+1),_base_endpoint(q,:one)
        return ((p,q,p0,q0),(_base_shape(:stair_odd_alt,2k+3),_base_shape(:ladder_even_alt,2k+2),_base_shape(:stair_odd_alt,2k+1),copy(q0)))
    elseif all(==('3'),word)
        p = _base_shape(:ladder_odd,2k+2)
        q = _base_shape(:stair_even,2k+3)
        p0,q0 = _base_endpoint(p,:zero),_base_shape(:stair_even,2k+1)
        return ((p,q,p0,q0),(_base_shape(:ladder_odd_alt,2k+2),_base_shape(:stair_even_alt,2k+3),copy(p0),_base_shape(:stair_even_alt,2k+1)))
    else
        throw(ArgumentError("word must be empty, all 2s or all 3s"))
    end
end
