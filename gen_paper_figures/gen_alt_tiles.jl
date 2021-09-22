## generate tile figures for paper 
## only modified for even-first, even tiles

function double_decomp_plot_and_annotate_alt_tile!(p::Rational{Int64}; annotate_center=true,width = 1, offSet = 0, textsize = 10, num = 1, isAnnotate = true, isRecurse = false, isDoubleDecomp=false, isPlotParent=false)
	if(p == 0//1 || p == 1//1) return plot_and_annotate_single_alt!(p, offSet=offSet, textsize=textsize, width=width,isAnnotate=isAnnotate, num = num); end
	str = rational_approximation(p,9999)[1]; 
	if(all2(str) || all3(str)) 
		return plot_and_annotate_single_alt!(p, offSet=offSet, textsize=textsize, width=width,isAnnotate=isAnnotate, num = num)
	end
	p0,q0,p1,q1 = getFracs(str); 
	W1 = count_full(str, '1')%2


	vp1,vp2 = vf(p1); vq1,vq2 = vf(q1)
	tp1, tq1 = getTile(p1), getTile(q1)

	tp1,tq1 = tp1.+offSet, tq1 .+offSet



	ano,anx = "+","-";
	#label each part of a doubled tile
	ano1,ano2 = ano*",1",ano*",2";
	anx1,anx2 = anx*",1",anx*",2";

	if(isAnnotate == false)
		ano=anx=ano1=ano2=anx1=anx2=false;
		annotate_center = false; 
	end

	O1,O2,O3,O4,Od = getAltTileOffSets(p, str);


	isLCorrect = false;
	s = str[end]; 
	if(s== '2' && isOdd(p) || s == '3' && isEven(p))  isLCorrect = true; end


	if(isLCorrect)

		O_ = O1+offSet;
		plot_and_annotate_corrected_tile!(q1, offSet = O_, isFlip = false, isLong = true, annotate_center=annotate_center,width = width,  textsize = textsize, num = num+1, isAnnotate = isAnnotate, isDoubleDecomp=isDoubleDecomp, isPlotParent=isPlotParent)

		O_ = O2+Od+offSet;
		plot_and_annotate_corrected_tile!(q1, offSet = O_, isFlip = true, isLong = true, annotate_center=annotate_center,width = width,  textsize = textsize, num = num+1, isAnnotate = isAnnotate, isDoubleDecomp=isDoubleDecomp, isPlotParent=isPlotParent)


		if(bool4)
			plot_and_annotate(tq1, q1, O = O1+(Od), label = q1, textsize = textsize,num = num, ann = ano2)  
			plot_and_annotate(tq1, q1, O = O2, label = q1, textsize = textsize,num = num, ann = anx1)  
		else
			plot_and_annotate_std_tile!(q1, offSet = offSet+O1+Od, textsize = textsize,num = num+1, width=width)  
			plot_and_annotate_std_tile!(q1, offSet = offSet+O2, textsize = textsize,num = num+1, width=width)  
		end

	else
		if(bool1)

		
			plot_and_annotate_std_tile!(q1, offSet = offSet+O1, textsize = textsize,num = num+1, width=width)  
			plot_and_annotate_std_tile!(q1, offSet = offSet+O2+Od, textsize = textsize,num = num+1, width=width)  

			if(bool5)
				plot_and_annotate(tq1, q1, O = O1+(Od), label = q1, textsize = textsize,num = num, ann = ano2)  
				plot_and_annotate(tq1, q1, O = O2, label = q1, textsize = textsize,num = num, ann = anx1)  


			else
				plot_and_annotate_std_tile!(q1, offSet = offSet+O1+Od, textsize = textsize,num = num+1, width=width)  
				plot_and_annotate_std_tile!(q1, offSet = offSet+O2, textsize = textsize,num = num+1, width=width)  

			end

		else
			plot_and_annotate(tq1, q1, O = O1, label = q1, textsize = textsize,num = num, ann = ano1)  
			plot_and_annotate(tq1, q1, O = O1+(Od), label = q1, textsize = textsize,num = num, ann = ano2)  
			plot_and_annotate(tq1, q1, O = O2, label = q1, textsize = textsize,num = num, ann = anx1)  
			plot_and_annotate(tq1, q1, O = O2+(Od), label = q1, textsize = textsize,num = num, ann = anx2)  

		end


	end

	if(bool3)
		plot_and_annotate(tp1, p1, O = O3, label = p1, textsize = textsize,num = num, ann = ano)  
		plot_and_annotate(tp1, p1, O = O4, textsize = textsize,num = num, ann = anx)
	else
		plot_and_annotate_std_tile!(p1, offSet = offSet+O3, textsize = textsize,num = num+1, width=width)  
		plot_and_annotate_std_tile!(p1, offSet = offSet+O4, textsize = textsize,num = num+1, width=width)  
	end

	#double decomposition
	if(isLCorrect) return plot!(); end #double overlap case

	##otherwise plot alterante
	if(bool2)
		plot_and_annotate_alt_tile!(p1, offSet=offSet+get_even_offset(p1,q1,str), textsize=textsize, width=width,isAnnotate=annotate_center, num = num+1)
	else
		plot_and_annotate_single_alt!(p1, offSet=offSet+get_even_offset(p1,q1,str), textsize=textsize, width=width,isAnnotate=annotate_center, num = num)
	end

	return plot!()


end
bool1 = true; 
bool2 = true; 
bool3 = true; 
bool4  =false
bool5 = true

str = "11111111"; 
tt = 5


p1,q1,p0,q0 = getFracs(str); 
plt1 = double_decomp_plot_and_annotate_alt_tile(q1, isAnnotate=true,textsize = tt); plot!(legend=false);

p1,q1,p0,q0 = getFracs(str*"2"); 
plt2 = double_decomp_plot_and_annotate_alt_tile(q1, isAnnotate=true,textsize = tt); plot!(legend=false);

p1,q1,p0,q0 = getFracs(str*"3"); 
plt3_ = double_decomp_plot_and_annotate_alt_tile(q1, isAnnotate=true,textsize = tt); plot!(legend=false);
plt3 = plot_and_annotate_alt_tile(q1, isAnnotate=true,textsize = tt); plot!(legend=false);


p1,q1,p0,q0 = getFracs(str*"23"); 
plt4_ = double_decomp_plot_and_annotate_alt_tile(q1, isAnnotate=true,textsize = tt); plot!(legend=false);
plt4 = plot_and_annotate_alt_tile(q1, isAnnotate=true,textsize = tt); plot!(legend=false);



k = 4
p1,q1,p0,q0 = getFracs(str*"3"^k); 
plt5_ = double_decomp_plot_and_annotate_alt_tile(q1, isAnnotate=true,textsize = tt); plot!(legend=false);
plt5 = plot_and_annotate_alt_tile(q1, isAnnotate=true,textsize = tt); plot!(legend=false);


p1,q1,p0,q0 = getFracs(str*"23"*"3"^(k-1)); 
plt6_ = double_decomp_plot_and_annotate_alt_tile(q1, isAnnotate=true,textsize = tt); plot!(legend=false);
plt6 = plot_and_annotate_alt_tile(q1, isAnnotate=true,textsize = tt); plot!(legend=false);

plot(plt5_,plt3_)
 # savefig(plot(plt1,dpi=512), "s_1.png")
 # savefig(plot(plt2,dpi=512), "s_2.png")
 # savefig(plot(plt3,dpi=512), "s_13.png")
 # savefig(plot(plt4,dpi=512), "s_23.png")
 # savefig(plot(plt3_,dpi=512), "s_13_.png")
 # savefig(plot(plt4_,dpi=512), "s_23_.png")
