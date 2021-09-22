## run this first 
using Plots

using FileIO,Colors, IndirectArrays #to produce pictures


include("tile_words.jl")
include("build_odometers_init.jl")
include("rational_recursion.jl")
include("tile_decomp/tile_offsets.jl")
include("tile_decomp/annotate_alt_tile.jl")
include("tile_decomp/annotate_std_tile.jl")
include("tile_decomp/annotate_corrected_tile.jl")
include("tile_decomp/plot_alt_tile.jl")
include("base_case.jl")
include("odometer_recursion.jl")
include("gen_figs_init.jl")