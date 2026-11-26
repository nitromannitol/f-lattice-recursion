# julia --project=. examples/export.jl 23 alt odd > tile.tsv
using FLatticeExact
length(ARGS)==3 || error("usage: julia --project=. examples/export.jl WORD std|alt odd|even")
word,kind,which=ARGS
kind in ("std","alt") || error("use std or alt")
which in ("odd","even") || error("use odd or even")
m,a=get_odometers(word)
write_tile(stdout,(kind=="std" ? m : a)[which=="odd" ? 1 : 2])
