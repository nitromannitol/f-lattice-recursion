@testset "Periodic extension" begin
    for (word,counts) in (("",nothing),("3312",nothing),
            ("1"*repeat("2",37),(42125,82011,82028,82028)),
            ("1"*repeat("2",44),(58701,114757,114774,114774)))
        standards,alternates = get_odometers(word)
        fractions = getFracs(word)
        for (f,tiles) in enumerate((standards,alternates)), j in 1:2
            tile = tiles[j]
            if !isnothing(counts)
                @test length(tile)==counts[2*(f-1)+j]
            end
            field = periodic_pattern(tile,fractions[j])
            @test length(field.labels)==field.area
            @test all(periodic_height(field,z)==h for (z,h) in tile)
            z = -31+23im
            @test periodic_height(field,z+field.u)-periodic_height(field,z)==field.cu
            @test periodic_height(field,z+field.v)-periodic_height(field,z)==
                real(field.gradient)*real(z)+imag(field.gradient)*imag(z)+field.cv
            @test pattern_value(field,z)==pattern_value(field,z+3*field.u-2*field.v)
            step = iseven(real(z)+imag(z)+field.horizontal_parity) ? 1 : im
            @test pattern_value(field,z)==periodic_height(field,z+step)+periodic_height(field,z-step)-2*periodic_height(field,z)
        end
    end
end

module MatrixInterface
using Test
include("../init.jl")
@testset "Matrix interface" begin
    m,a = get_odometers("23")
    em,ea = FLatticeExact.get_odometers("23")
    @test getFracs("23")==FLatticeExact.getFracs("23")
    for (matrices,tiles) in ((m,em),(a,ea)), j in 1:4
        arr,t = matrices[j],tiles[j]
        xmin,ymin = minimum(real,keys(t)),minimum(imag,keys(t))
        @test all(arr[Int(real(z)-xmin+2),Int(imag(z)-ymin+2)]==h for (z,h) in t)
        @test count(isfinite,arr)==length(t)
        @test all(isnan,arr[1,:]) && all(isnan,arr[end,:]) && all(isnan,arr[:,1]) && all(isnan,arr[:,end])
    end
end
end
