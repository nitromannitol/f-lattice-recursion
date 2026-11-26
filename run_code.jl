# julia --project=. run_code.jl 3312 patterns.svg
# Draws the Laplacian of the two standard and two alternate tiles of a word.
include(joinpath(@__DIR__,"src","FLatticeExact.jl"))
using .FLatticeExact

function draw_panel(out,tile,fraction,family,role,x,y)
    p = periodic_pattern(tile,fraction)
    xmin,xmax = extrema(real,keys(tile))
    ymin,ymax = extrema(imag,keys(tile))
    width,height = xmax-xmin+1,ymax-ymin+1
    println(out,"<text x='$x' y='$(y+20)'>$family $role · $(numerator(fraction))/$(denominator(fraction))</text>")
    println(out,"<svg x='$x' y='$(y+32)' width='314' height='264' viewBox='0 0 $width $height' preserveAspectRatio='xMidYMid meet' shape-rendering='crispEdges'>")
    for (value,color) in ((0,"#c5c5c5"),(-1,"#171717"))
        print(out,"<path fill='$color' d='")
        for z in keys(tile)
            if pattern_value(p,z)==value
                print(out,"M",real(z)-xmin," ",ymax-imag(z),"h1v1h-1z")
            end
        end
        println(out,"'/>")
    end
    println(out,"</svg>")
end

function render(word,path)
    standards,alternates = get_odometers(word)
    fractions = getFracs(word)
    svg = IOBuffer()
    println(svg,"<svg xmlns='http://www.w3.org/2000/svg' width='704' height='686' viewBox='0 0 704 686' font-family='sans-serif' font-size='17'>")
    println(svg,"<title>F-lattice patterns for $(isempty(word) ? "the empty word" : word)</title><rect width='704' height='686' fill='white'/>")
    for (f,(family,tiles)) in enumerate((("Standard",standards),("Alternate",alternates)))
        for (i,role) in enumerate(("odd","even"))
            draw_panel(svg,tiles[i],fractions[i],family,role,24+340*(i-1),12+312*(f-1))
        end
    end
    println(svg,"<text x='24' y='661' font-size='14'>Black = −1     ·     Gray = 0     ·     White = outside the tile</text></svg>")
    write(path,take!(svg))
    println("Wrote ",abspath(path))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)<=2 || error("usage: julia run_code.jl [word] [patterns.svg]")
    word = isempty(ARGS) ? "3312" : ARGS[1]
    path = length(ARGS)==2 ? ARGS[2] : "patterns.svg"
    render(word,path)
elseif isinteractive() && isempty(PROGRAM_FILE)
    render("3312","patterns.svg")
end
