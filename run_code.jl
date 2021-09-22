## example useage of the code
include("init.jl")

## do not be alarmed, the code takes a few seconds to run 

str = "3312"; #change to any string in the alphabet 1,2,3

mats,alt_mats = get_odometers(str)
p1,q1,p0,q0 = getFracs(str); 
isTile = false; #change to true to tile
isDecomp = true; #change to false to not plot the tile decomposition

plt1 = full_plot_std(mats[1], p1, isTile = isTile, isDecomp = isDecomp); plot!(title = "std -$(p1)"); 
plt2 = full_plot_std(mats[2], q1, isTile = isTile, isDecomp = isDecomp); plot!(title = "std -$(q1)");
plt3 = full_plot_alt(alt_mats[1], p1, isTile = isTile, isDecomp = isDecomp); plot!(title = "alt -$(p1)");
plt4 = full_plot_alt(alt_mats[2], q1, isTile = isTile, isDecomp = isDecomp); plot!(title = "alt -$(q1)");
display(plot(plt1,plt2,plt3,plt4))
