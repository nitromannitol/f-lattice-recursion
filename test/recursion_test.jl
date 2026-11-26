using Test
@testset "Tile sizes" begin
    expected=Dict(""=>(21,21,30,30),"1"=>(93,93,110,110),
      "2"=>(45,57,72,58),"3"=>(57,45,58,72),
      "11"=>(499,499,586,586),"12"=>(221,291,308,308),
      "13"=>(291,221,308,308),"112"=>(1233,1633,1720,1720),
      "113"=>(1633,1233,1720,1720))
    for (w,ns) in expected
        m,a=get_odometers(w)
        @test (length(m[1]),length(m[2]),length(a[1]),length(a[2]))==ns
        @test all(T->all(v->v isa BigInt,values(T)),(m...,a...))
    end
    @test_throws ArgumentError get_odometers("10")
    @test_throws ArgumentError get_odometers("abc")
    # _recbirth and rational_approximation should agree.
    words=[""]
    for depth in 1:4
        layer=[w*string(c) for w in words if length(w)==depth-1 for c in ('1','2','3')]
        append!(words,layer)
    end
    for w in words, odd in (true,false)
        f=getFracs(w)[odd ? 3 : 4]
        if 0<f<1
            @test _recbirth(f,w,odd)==rational_approximation(f)
        end
    end
end
