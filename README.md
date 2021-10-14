# f-lattice-recursion


This repository contains code to produce tile odometers and tiles following the algorithm in 
https://nitromannitol.github.io/manuscripts/farey_recursion.pdf


Instructions

1. install Julia 1.5.0 and install Plots with GrPlots backend, FileIO,Colors, IndirectArrays
	
	https://julialang.org/downloads/

2. 	start Julia and run 
  
		>include("PATHTOFILE/init.jl")
		


3.  compute standard and alternate tile odometers by running, for example

		>str = "31"
		>mats,alt_mats = get_odometers(str)
		>p1,q1,p0,q0 = getFracs(str); 
		>full_plot_std(mats[1], p1) 

--- 
	- see run_code.jl for more sample code


