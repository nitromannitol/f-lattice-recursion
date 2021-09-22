
SS = "23311111"
width = 1; 
isAnnotate=true;
isRecurse = false;
isDouble = false;
num = 1; 
for str in [SS*"1", SS*"2", SS*"3", SS*"11",SS*"12",SS*"13"]
	p0,q0,p1,q1 = getFracs(str); 
	W1 = count_full(str, '1')%2
	if(W1 == 0) ss = "e" else ss = "o"; end
	ss = ss*str[end];
	println(ss)

	textsize = 10; 


	#dont draw the doubled tile, but label them individually
	isDouble = false; 
	pltp_std = plot_and_annotate_std_tile(p0, width = width, isDouble = isDouble, textsize = textsize,isAnnotate=isAnnotate)
	pltq_std = plot_and_annotate_std_tile(q0, width = width, isDouble = isDouble, textsize = textsize,isAnnotate=isAnnotate)

	pltp_alt = plot_and_annotate_alt_tile(p0, textsize = textsize, width = width,isAnnotate=isAnnotate);
	pltq_alt = plot_and_annotate_alt_tile(q0, textsize = textsize, width = width,isAnnotate=isAnnotate);


	savefig(plot(pltp_std, pltq_std, legend=false,dpi=512), "annotated_tiles/std_tile/$(ss)")
	savefig(plot(pltp_alt, pltq_alt, legend=false,dpi=512), "annotated_tiles/alt_tile/$(ss)")


	## save the single versions
	savefig(plot(pltp_std, legend=false,dpi=512), "annotated_tiles/single/std_tile/$(ss)_p")
	savefig(plot(pltq_std, legend=false,dpi=512), "annotated_tiles/single/std_tile/$(ss)_q")

	savefig(plot(pltp_alt, legend=false,dpi=512), "annotated_tiles/single/alt_tile/$(ss)_p")
	savefig(plot(pltq_alt, legend=false,dpi=512), "annotated_tiles/single/alt_tile/$(ss)_q")



	#double decomposition
	isDoubleDecomp=true
	textsize=8;

	pltp_alt = double_decomp_plot_and_annotate_alt_tile(p0, textsize=textsize, isDoubleDecomp=isDoubleDecomp,width=width, isAnnotate=isAnnotate);
	pltq_alt = double_decomp_plot_and_annotate_alt_tile(q0, textsize=textsize, isDoubleDecomp=isDoubleDecomp,width=width, isAnnotate=isAnnotate);

	pltp_std = double_decomp_plot_and_annotate_std_tile(p0, isAnnotate=isAnnotate,isRecurse = isRecurse, width = width, isDouble = isDouble, num = num, textsize=textsize)
	pltq_std = double_decomp_plot_and_annotate_std_tile(q0, isAnnotate=isAnnotate,isRecurse = isRecurse, width = width, isDouble = isDouble, num = num, textsize=textsize)


	savefig(plot(pltp_std, pltq_std, legend=false,dpi=512), "annotated_tiles/std_tile/double_decomp/$(ss)")
	savefig(plot(pltp_alt, pltq_alt, legend=false,dpi=512), "annotated_tiles/alt_tile/double_decomp/$(ss)")


	## save the single versions
	savefig(plot(pltp_std, legend=false,dpi=512), "annotated_tiles/single/std_tile/double/$(ss)_p")
	savefig(plot(pltq_std, legend=false,dpi=512), "annotated_tiles/single/std_tile/double/$(ss)_q")

	savefig(plot(pltp_alt, legend=false,dpi=512), "annotated_tiles/single/alt_tile/double/$(ss)_p")
	savefig(plot(pltq_alt, legend=false,dpi=512), "annotated_tiles/single/alt_tile/double/$(ss)_q")
end