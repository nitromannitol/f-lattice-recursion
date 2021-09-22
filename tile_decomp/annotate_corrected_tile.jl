#uses functions in entire folder
function plot_and_annotate_corrected_tile_odd!(p::Rational{Int64}; isLong = false, isFlip = false, annotate_center=true,width = 1, offSet = 0, textsize = 10, num = 1, isAnnotate = true, isDoubleDecomp=false, isPlotParent=false)
	if(p == 0//1 || p == 1//1) return plot!(); end

	str = rational_approximation(p,9999)[1]; # be careful
	W1 = count_full(str, '1')%2

	O1,O2,O3,O4 = getTileOffSets(p, str);
	p3,q3,p4,q4 = getFracs(str); vp1,vp2 = vf(p4); vq1,vq2 = vf(q4)


	##initial tile offset 
	if(W1==0)
		if(isFlip)
			O = O2+vq1
		else
			O = O1; 
		end
	else
		if(isFlip)
			O = O2
		else
			O = O1+vq1; 
		end
	end



	k = count_end(str, '2')
	str = str[1:end-k]
	p0,q0,p1,q1 = getFracs(str)

	## use the first two
	tp1, tq1 = getTile(p1), getTile(q1)
	tp1,tq1 = tp1.+offSet, tq1 .+offSet
	vp1,vp2 = vf(p1); vq1,vq2 = vf(q1)


	##number times alternate repeats
	num_iters = 2*(k+1); 
	#if(length(str) > 0 && str[end] == '3' && !isLong) num_iters-=1; end
	if(!isLong) num_iters-=1; end #test

	K1,K2,K3,K4,A1,A2,A3,A4 = get_odd_l_offsets(str, isFlip)

	ano = ""
	curr_ann = 1;
	ann_ = ano*"$(curr_ann)"; if(isAnnotate == false) ann_ = false; annotate_center=false; end
	plot_and_annotate(tq1, q1, O = O, label = q1, textsize = textsize,num = num, ann = ann_,width = width);  
	for i in 1:(k)
		O = O + K1; 
		curr_ann+=1;
		ann_ = ano*"$(curr_ann)"; if(isAnnotate == false) ann_ = false; end
		plot_and_annotate(tq1, q1, O = O, label = q1, textsize = textsize,num = num, ann = ann_,width = width);  
	end
	O = O + K2; 
	ann_ = ano; if(isAnnotate == false) ann_ = false; end
	if(p1 == 0//1)  
		if(!isFlip) tp1 = tp1 .+ (-im); else tp1 = tp1 .+ (-im + 1); end
	end

	plot_and_annotate(tp1, p1, O = O, label = p1, textsize = textsize,num = num, ann = ann_,width = width);  

	O = O + K3; 

	for i in 1:num_iters;
		O = O + K4
		if(isDoubleDecomp)
			plot_and_annotate_alt_tile!(q1, annotate_center=isAnnotate,width = width, offSet = O+offSet, textsize = textsize, num = num+1, isAnnotate = isAnnotate)
		else
			plot_and_annotate_single_alt!(q1, num = num, offSet=O+offSet, width = width, textsize=textsize, isAnnotate=isAnnotate, ann = ano*"$(i)")
		end
	end

	return plot!(legend=false)
end


function plot_and_annotate_corrected_tile_even!(p::Rational{Int64}; isLong = false, isFlip = false, annotate_center=true,width = 1, offSet = 0, textsize = 10, num = 1, isAnnotate = true, isDoubleDecomp=false, isPlotParent=false)
	if(p == 0//1 || p == 1//1) return plot!(); end

	str = rational_approximation(p,9999)[1];

	W1 = count_full(str, '1')%2

	O1,O2,O3,O4 = getTileOffSets(p, str);
	p3,q3,p4,q4 = getFracs(str); vp1,vp2 = vf(p4); vq1,vq2 = vf(q4)

	##initial tile offset 
	if(isFlip)
		O = O4 + vp2; 
	else
		O = O3; 
	end

	k = count_end(str, '3')
	str = str[1:end-k]
	p0,q0,p1,q1 = getFracs(str)

	## use the first two
	tp1, tq1 = getTile(p1), getTile(q1)
	tp1,tq1 = tp1.+offSet, tq1 .+offSet
	vp1,vp2 = vf(p1); vq1,vq2 = vf(q1)



	##number times alternate repeats
	num_iters = 2*(k+1); 
	#if(length(str) > 0 && str[end] == '2' && !isLong) num_iters-=1; end
	if(!isLong) num_iters-=1; end

	K1,K2,K3,K4,A1,A2,A3,A4 = get_even_l_offsets(str, isFlip)


	ano = ""
	curr_ann = 1;
	ann_ = ano*"$(curr_ann)"; if(isAnnotate == false) ann_ = false; annotate_center=false; end

	plot_and_annotate(tp1, p1, O = O, label = p1, textsize = textsize,num = num, ann = ann_,width = width);  
	for i in 1:(k)
		O = O + K1; 
		curr_ann+=1;
		ann_ = ano*"$(curr_ann)"; if(isAnnotate == false) ann_ = false; end

		plot_and_annotate(tp1, p1, O = O, label = p1, textsize = textsize,num = num, ann = ann_,width = width);  
	end
	O = O + K2; 
	ann_ = ano*"$(curr_ann)"; if(isAnnotate == false) ann_ = false; end

	if(q1 == 1//1 && isFlip) tq1 = tq1 .+ (1+im); end
	plot_and_annotate(tq1, q1, O = O, label = q1, textsize = textsize,num = num, ann = ann_, width = width);  

	O = O + K3; 

	for i in 1:num_iters;
		O = O + K4
		if(isDoubleDecomp)
			plot_and_annotate_alt_tile!(p1, annotate_center=annotate_center,width = width, offSet = O+offSet, textsize = textsize, num = num, isAnnotate = isAnnotate)
		else
			plot_and_annotate_single_alt!(p1, num = num, offSet=O+offSet, width = width, textsize=textsize, isAnnotate=isAnnotate, ann = ano*"$(i)")
		end
	end

	return plot!(legend=false)
end



function plot_and_annotate_corrected_tile!(p::Rational{Int64};isLong = false, isFlip = false, annotate_center=true,width = 1, offSet = 0, textsize = 10, num = 1, isAnnotate = true, isDoubleDecomp=false, isPlotParent=false)
	if(isOdd(p))
		return plot_and_annotate_corrected_tile_odd!(p, isLong = isLong, isFlip = isFlip, annotate_center=annotate_center,width = width, offSet = offSet, textsize = textsize, num = num, isAnnotate = isAnnotate, isDoubleDecomp=isDoubleDecomp, isPlotParent=isPlotParent)
	else
		return plot_and_annotate_corrected_tile_even!(p, isLong = isLong, isFlip = isFlip, annotate_center=annotate_center,width = width, offSet = offSet, textsize = textsize, num = num, isAnnotate = isAnnotate, isDoubleDecomp=isDoubleDecomp, isPlotParent=isPlotParent)
	end
end

## AN L-correction overlaid on a standard tile
#str = "111112"
#p0,q0,p1,q1 = getFracs(str);
#p = p0; 
#isFlip = false;
#textsize=8
#plot_and_annotate_std_tile(p, width = 0.5, textsize = textsize,isAnnotate=false);
#plot_and_annotate_corrected_tile!(p, width = 1, isFlip = false, textsize = textsize)
#savefig(plot!(dpi=512), "odd_tile_$(str)_overlay")