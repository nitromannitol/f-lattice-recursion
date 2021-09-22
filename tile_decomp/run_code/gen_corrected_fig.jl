SS = "111111"
textsize = 8; 
for str in [SS*"23",SS*"3", SS*"32", SS*"33"], str in ["1"*str, str], isFlip in [false, true]
	p0,q0,p1,q1 = getFracs(str); 
	W1 = count_full(str, '1')%2
	if(W1 == 0 && !isFlip)
		ss = "00"
	elseif(W1 == 1 && !isFlip)
		ss = "10"
	elseif(W1 == 0 && isFlip)
		ss = "01"
	elseif(W1 == 1 && isFlip)
		ss = "11"
	end
	for p in [p0,q0]
		plot(aspect_ratio = true, axis = :off, ticks = :none, legend=false)
		plt2 = plot_and_annotate_corrected_tile!(p, isFlip = isFlip, textsize = textsize);
		plt1 = plot_and_annotate_std_tile(p, textsize = textsize,isAnnotate=true,isLong=true);
		ss = ss*"_"*str[end];
		println(ss)
		folder_dir = "even_tile"; if(isOdd(p)) folder_dir = "odd_tile"; end
		savefig(plot(plt2, legend=false,dpi=512), "annotated_tiles/corrected_tiles/$(folder_dir)/$(ss)")
		savefig(plot(plt1, legend=false,dpi=512), "annotated_tiles/corrected_tiles/$(folder_dir)/$(ss[1])")
	end
end