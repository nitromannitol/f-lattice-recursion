## generate tile figures for paper 
## only modified for even-first, odd tiles
function double_decomp_plot_and_annotate_std_tile!(p::Rational{Int64}; width = 1, offSet = 0, isDouble = false, textsize = 10, num = 1, isAnnotate = true, isRecurse = false)
	if(p == 0//1 || p == 1//1) return plot!(); end
	str = rational_approximation(p,9999)[1]; 
	t0 = getTile(p); 

	if(all2(str) || all3(str)) 
		return plot_and_annotate(t0, p, O = offSet, label = p, textsize = textsize,num = num, width=width)  
	end

#hi 
	p0,q0,p1,q1 = getFracs(str); 
	W1 = count_full(str, '1')%2;


	vp1,vp2 = vf(p1); vq1,vq2 = vf(q1)

	tp1, tq1 = getTile(p1), getTile(q1)

	t0 = t0 .+ offSet
	tp1,tq1 = tp1.+offSet, tq1 .+offSet


	##double the tiles 
	tp1d = [tp1; tp1 .+ vp2]
	tq1d = [tq1; tq1 .+ vq1]


	ano,anx = "+","-";
	#label each part of a doubled tile
	ano1,ano2 = ano*",1",ano*",2";
	anx1,anx2 = anx*",1",anx*",2";

	if(isAnnotate == false)
		ano=anx=ano1=ano2=anx1=anx2=false;
	end

	O1,O2,O3,O4 = getTileOffSets(p, str);


	#determine whether or not there is an L correction
	isLCorrect = false;
	s = str[end]; 
	if(isEven(p) && s == '2' || isOdd(p) && s == '3') isLCorrect = true; end

	## set some auxillary terms
	annotate_center = isAnnotate;
	isDoubleDecomp=isDouble
	isPlotParent=false



	if(isLCorrect)
		if(W1 == 0)
			#L corrections
			O_ = O1 + offSet+vq1; 
			O_ = O_+ (vp2+ (2vp1-vq1))
			plot_and_annotate_corrected_tile!(q1, offSet = O_, isFlip = true,annotate_center=annotate_center,width = width,  textsize = textsize, num = num+1, isAnnotate = isAnnotate, isDoubleDecomp=isDoubleDecomp, isPlotParent=isPlotParent)
			O_ = O2 + offSet; 
			O_ = O_- (vp2+ (2vp1-vq1))
			plot_and_annotate_corrected_tile!(q1, offSet = O_, isFlip = false,annotate_center=annotate_center,width = width,  textsize = textsize, num = num+1, isAnnotate = isAnnotate, isDoubleDecomp=isDoubleDecomp, isPlotParent=isPlotParent)

			#standard tiles
			plot_and_annotate_std_tile!(q1, offSet = offSet+O1, textsize = textsize,num = num+1, width=width)  
			plot_and_annotate_std_tile!(q1, offSet = offSet+O2+vq1, textsize = textsize,num = num+1, width=width)  
		end
	else
		#plt1 
		if(bool1)
			plot_and_annotate(tq1, q1, O = O1, label = q1, textsize = textsize,num = num, ann = ano1, width=width)  
			plot_and_annotate(tq1, q1, O = O2+vq1, textsize = textsize,num = num, ann = anx2, width=width)  

			plot_and_annotate_std_tile!(q1, offSet = offSet+O1+vq1, textsize = textsize,num = num+1, width=width)  
			plot_and_annotate_std_tile!(q1, offSet = offSet+O2, textsize = textsize,num = num+1, width=width)  
		else
			plot_and_annotate(tq1, q1, O = O1, label = q1, textsize = textsize,num = num, ann = ano1, width=width)  
			plot_and_annotate(tq1, q1, O = O1+vq1, textsize = textsize,num = num, ann = ano2, width=width)  
			plot_and_annotate(tq1, q1, O = O2, label = q1, textsize = textsize,num = num, ann = anx1, width=width)  
			plot_and_annotate(tq1, q1, O = O2+vq1, textsize = textsize,num = num, ann = anx2, width=width)  

		end
	end

	if(!bool1)
		plot_and_annotate_std_tile!(p1, offSet = offSet+O3, textsize = textsize,num = num+1, width=width)  
		plot_and_annotate_std_tile!(p1, offSet = offSet+O4, textsize = textsize,num = num+1, width=width)  
	else
		plot_and_annotate(tp1, p1, O = O3, label = p1, textsize = textsize,num = num, ann = ano, width=width)  
		plot_and_annotate(tp1, p1, O = O4, textsize = textsize,num = num, ann = anx, width=width)  

	end


	return plot!()


end



## 4 cases
str = "11111111"; 

tt = 5

global bool1 = true ;
#toggle to get correct double decomp for plt1 -true or plt2 -false must be true for plt3,plt4








p1,q1,p0,q0 = getFracs(str); 
plt1 = double_decomp_plot_and_annotate_std_tile(p1, isAnnotate=true,textsize = tt); plot!(legend=false);


p1,q1,p0,q0 = getFracs(str*"2"); 
plt2 = double_decomp_plot_and_annotate_std_tile(p1, isAnnotate=true,textsize = tt); plot!(legend=false);

p1,q1,p0,q0 = getFracs(str*"3"); 
plt3_ = double_decomp_plot_and_annotate_std_tile(p1, isAnnotate=true,textsize = tt); plot!(legend=false);
plt3 = plot_and_annotate_std_tile(p1, isAnnotate=true,textsize = tt); plot!(legend=false);

p1,q1,p0,q0 = getFracs(str*"23"); 
plt4_ = double_decomp_plot_and_annotate_std_tile(p1, isAnnotate=true,textsize = tt); plot!(legend=false);
plt4 = plot_and_annotate_std_tile(p1, isAnnotate=true,textsize = tt); plot!(legend=false);






#savefig(plot(plt1,dpi=512), "s_1.png")
#savefig(plot(plt2,dpi=512), "s_2.png")
# savefig(plot(plt3,dpi=512), "s_13.png")
# savefig(plot(plt4,dpi=512), "s_23.png")
# savefig(plot(plt3_,dpi=512), "s_13_.png")
# savefig(plot(plt4_,dpi=512), "s_23_.png")




k = 4
p1,q1,p0,q0 = getFracs(str*"3"^k); 
plt5_ = double_decomp_plot_and_annotate_std_tile(p1, isAnnotate=true,textsize = tt); plot!(legend=false);
plt5 = plot_and_annotate_std_tile(p1, isAnnotate=true,textsize = tt); plot!(legend=false);


p1,q1,p0,q0 = getFracs(str*"23"*"3"^(k-1)); 
plt6_ = double_decomp_plot_and_annotate_std_tile(p1, isAnnotate=true,textsize = tt); plot!(legend=false);
plt6 = plot_and_annotate_std_tile(p1, isAnnotate=true,textsize = tt); plot!(legend=false);

plot(plt5_,plt3_)