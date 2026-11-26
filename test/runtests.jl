using Test
using FLatticeExact
using FLatticeExact: _recbirth, rational_approximation, pt, place, touch, llcenter, ur
const X = FLatticeExact
@testset "Basic operations" begin
    v=big(2)^125; h=big(2)^140+1
    T=Tile(pt(0,0)=>h,pt(1,0)=>h+1,pt(-2,1)=>-h)
    U=place(T,pt(-v,v+1),pt(v,2v))
    @test length(T)==length(U)==3
    for (z,hz) in T
        k=z+pt(-v,v+1)
        @test U[k]==hz+v*real(k)+2v*imag(k)
    end
    @test X.ll(T)==pt(0,0) && ur(T)==pt(-2,1)
    @test X.ll(X.trim(T))==pt(2,0)
    # The old tile keeps its values where the two tiles overlap.
    B=Tile(pt(0,0)=>h+7, pt(1,0)=>h+9, pt(2,0)=>h+10)
    R=touch(T,B;context="unit example")
    @test R[pt(1,0)]==h+1 && R[pt(2,0)]==h+3
    @test length(R)==4 && length(T)==3
    @test_throws ArgumentError touch(T,place(T,pt(v,0),0);context="disjoint")
    @test_throws ArgumentError touch(Tile(),T)
    @test_throws ArgumentError place(T,1//2,0)
    @test_throws ArgumentError get_odometers("12𝟛")
    @test denominator(getFracs(repeat("1",100))[1])>typemax(Int64)
end
@testset "Saved tiles" begin
    dir=joinpath(@__DIR__,"fixtures")
    fixtures=sort(filter(f->endswith(f,".tsv"),readdir(dir)))
    @test length(fixtures)==31
    cache=Dict{String,Any}()
    for file in fixtures
        w,kind,i=split(first(splitext(file)),'_')
        w=(w=="empty" ? "" : w)
        m,a=get!(()->get_odometers(w),cache,w)
        actual=(kind=="std" ? m : a)[parse(Int,i)]
        expected=Tile()
        for row in eachline(joinpath(dir,file))
            startswith(row,"#") && continue
            isempty(strip(row)) && continue
            x,y,h=parse.(BigInt,split(row))
            expected[pt(x,y)]=h
        end
        @test actual==expected
    end
end
include("recursion_test.jl")
include("pattern_test.jl")
