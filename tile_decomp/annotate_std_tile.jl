using Plots, LaTeXStrings


function getOddEvenParents(p::Rational{Int64})
	p1,p2 = getP(p)
	if(isOdd(p1)) p1,q1 = p1,p2; else p1,q1 = p2,p1; end 
end

function plot_and_annotate(t::Array{Complex{Int64},1}, p::Rational{Int64}; O = 0, label = false, width = 1, ann = false, textsize = 10, isDouble = false, num = 1)
	if(isOdd(p)) color = "blue"; else color = "red"; end
	plot_t!(t.+O, p, label = label, color = color, width = width); 
	O = O + t[1]; ## shift to tile offset
	if(ann != false)
		v1,v2 = vf(p); if(isOdd(p)) v1 = 2*v1; end 
		if(isDouble == false)
			O = O + v1/2 + v2/2; 
		elseif(isDouble && isOdd(p))
			O = O + v1/2 + v2
		elseif(isDouble && isEven(p))
			O = O + v1 + v2/2
		end
		if(isDouble) T = "T^d"; else T = "T"; end
		if(num == false)
			if(isOdd(p)) p_s = "_p" else p_s ="_q"; end
		else
			if(isOdd(p)) p_s = "(p_{$(num)})" else p_s ="(q_{$(num)})"; end
		end
		p_s = p_s*"^{$(ann)}";
		#if(isOdd(p)) p_s = "(p_$(num))^{$(ann)}" else p_s = "(q_$(num))^{$(ann)}";  end

		ann = LaTeXString("\$"*T*p_s*"\$"); 
		annotate!(real(O),imag(O), text(ann, :black, :center, textsize))
	end
	return plot!()
end


#some of these arguments are not used
function plot_and_annotate_single_std!(p::Rational{Int64}; width = 1, offSet = 0, isDouble = false, textsize = 10, num = 1, isAnnotate = true, isRecurse = false)
	t0 = getTile(p); 
	plot_and_annotate(t0, p, O = offSet, label = p, textsize = textsize,num = num, width=width)  
end

function plot_and_annotate_single_std(p::Rational{Int64}; width = 1, offSet = 0, isDouble = false, textsize = 10, num = 1, isAnnotate = true, isRecurse = false)
	plot(aspect_ratio = true, axis = :off, ticks = :none)
	return plot_and_annotate_single_std!(p, width=width, offSet=offSet, isDouble=isDouble,  textsize=textsize, num=num, isAnnotate=isAnnotate, isRecurse=isRecurse)
end


function plot_and_annotate_std_tile!(p::Rational{Int64}; width = 1, offSet = 0, isDouble = false, textsize = 10, num = 1, isAnnotate = true, isRecurse = false)
	if(p == 0//1 || p == 1//1) return plot!(); end
	str = rational_approximation(p,9999)[1]; 
	t0 = getTile(p); 

	if(all2(str) || all3(str)) 
		return plot_and_annotate(t0, p, O = offSet, label = p, textsize = textsize,num = num, width=width)  
	end


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


	if(isEven(p)) #two cases are the same
		if(isLCorrect)
			#L corrections
			O_ = O3 + vp2+offSet
			plot_and_annotate_corrected_tile!(p1, offSet = O_, isFlip = true, annotate_center=annotate_center,width = width,  textsize = textsize, num = num+1, isAnnotate = isAnnotate, isDoubleDecomp=isDoubleDecomp, isPlotParent=isPlotParent)
			O_ = O4+offSet
			plot_and_annotate_corrected_tile!(p1, offSet = O_, isFlip = false,annotate_center=annotate_center,width = width,  textsize = textsize, num = num+1, isAnnotate = isAnnotate, isDoubleDecomp=isDoubleDecomp, isPlotParent=isPlotParent)

			#standard tiles
			plot_and_annotate(tp1, p1, O = O3, label = p1, textsize = textsize,num = num, ann = ano1, width=width)  
			plot_and_annotate(tp1, p1, O = O4+vp2, textsize = textsize,num = num, ann = anx2, width=width)  
		elseif(isDouble)
			plot_and_annotate(tp1d, p1, O = O3, label = p1, textsize = textsize,num = num, ann = ano, isDouble=true, width=width)    
			plot_and_annotate(tp1d, p1, O = O4, textsize = textsize,num = num, ann = anx, isDouble=true, width=width)  
		else
			plot_and_annotate(tp1, p1, O = O3, label = p1, textsize = textsize,num = num, ann = ano1, width=width)  
			plot_and_annotate(tp1, p1, O = O3+vp2, textsize = textsize,num = num, ann = ano2, width=width)  

			plot_and_annotate(tp1, p1, O = O4, textsize = textsize,num = num, ann = anx1, width=width)  
			plot_and_annotate(tp1, p1, O = O4+vp2, textsize = textsize,num = num, ann = anx2, width=width)  
		end
		plot_and_annotate(tq1, q1, O = O1, label = q1, textsize = textsize,num = num, ann = ano, width=width)    
		plot_and_annotate(tq1, q1, O = O2, label = q1, textsize = textsize,num = num, ann = anx, width=width)    

	else
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
				plot_and_annotate(tq1, q1, O = O1, label = q1, textsize = textsize,num = num, ann = ano, width=width)  
				plot_and_annotate(tq1, q1, O = O2+vq1, label = q1, textsize = textsize,num = num, ann = anx2, width=width)  

			else
				#L corrections
				O_ = O1 + offSet; 
				plot_and_annotate_corrected_tile!(q1, offSet = O_, isFlip = false,annotate_center=annotate_center,width = width,  textsize = textsize, num = num+1, isAnnotate = isAnnotate, isDoubleDecomp=isDoubleDecomp, isPlotParent=isPlotParent)
				O_ = O2 + vq1+offSet; 
				plot_and_annotate_corrected_tile!(q1, offSet = O_, isFlip = true,annotate_center=annotate_center,width = width,  textsize = textsize, num = num+1, isAnnotate = isAnnotate, isDoubleDecomp=isDoubleDecomp, isPlotParent=isPlotParent)

				#standard tiles
				plot_and_annotate(tq1, q1, O = O1+vq1, label = q1, textsize = textsize,num = num, ann = ano1, width=width)   
				plot_and_annotate(tq1, q1, O = O2, label = q1, textsize = textsize,num = num, ann = anx2, width=width)    

			end
		elseif(isDouble)
			plot_and_annotate(tq1d, q1, O = O1, label = q1, textsize = textsize,num = num, ann = ano, isDouble=true, width=width)  
			plot_and_annotate(tq1d, q1, O = O2, label = q1, textsize = textsize,num = num, ann = anx, isDouble = true, width=width)  
		else
			plot_and_annotate(tq1, q1, O = O1, label = q1, textsize = textsize,num = num, ann = ano1, width=width)  
			plot_and_annotate(tq1, q1, O = O1+vq1, textsize = textsize,num = num, ann = ano2, width=width)  

			plot_and_annotate(tq1, q1, O = O2, label = q1, textsize = textsize,num = num, ann = anx1, width=width)  
			plot_and_annotate(tq1, q1, O = O2+vq1, textsize = textsize,num = num, ann = anx2, width=width)  
		end

		plot_and_annotate(tp1, p1, O = O3, label = p1, textsize = textsize,num = num, ann = ano, width=width)  
		plot_and_annotate(tp1, p1, O = O4, textsize = textsize,num = num, ann = anx, width=width)  

	end


	return plot!()


end

function plot_and_annotate_std_tile(p::Rational{Int64}; width = 1, offSet = 0, isDouble = false, textsize = 10, num = 1, isAnnotate = true, isRecurse = false)
	plot(aspect_ratio = true, axis = :off, ticks = :none)
	return plot_and_annotate_std_tile!(p, width=width, offSet=offSet, isDouble=isDouble,  textsize=textsize, num=num, isAnnotate=isAnnotate, isRecurse=isRecurse)
end


# note this doesnt' work for when you have a correction 
## subscript num are the two parents of the even parent
## subscript num+1 are the two parents of the odd parent
function double_decomp_plot_and_annotate_std_tile!(p::Rational{Int64}; width = 1, offSet = 0, isDouble = false, textsize = 10,  num = 1, isRecurse = false, isAnnotate=false)
	if(p == 0//1 || p == 1//1) return plot!(); end
	str = rational_approximation(p,9999)[1]; # be careful
	if(length(str) ==1)
		return plot_and_annotate_std_tile!(p,  width = width, offSet = offSet, isDouble = isDouble, textsize = textsize, num = num+1)
	end
	p0,q0,p1,q1 = getFracs(str); 

	t0 = getTile(p)
	W1 = count_full(str, '1')%2

	vp1,vp2 = vf(p1)
	vq1,vq2 = vf(q1)


	O = offSet; 
	O1,O2,O3,O4 = getTileOffSets(p, str);


	if(isRecurse == true)
		f! = double_decomp_plot_and_annotate_std_tile!;
		textsize = textsize -1; 
		num = num + 1; 
		textsize = 1; 
	else
		f! = plot_and_annotate_std_tile!;
	end

	if(isEven(p))
		f!(q1,offSet = O + O1, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse, isDouble = isDouble, textsize = textsize, num = num)
		f!(q1,offSet = O + O2, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse, isDouble = isDouble, textsize = textsize, num = num)
		f!(p1,offSet = O + O3, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse, isDouble = isDouble, textsize = textsize, num = num+1)
		f!(p1,offSet = O + O3+vp2, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse, isDouble = isDouble, textsize = textsize, num = num+1)
		f!(p1,offSet = O + O4, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse, isDouble = isDouble, textsize = textsize, num = num+1)
		f!(p1,offSet = O + O4+vp2, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse, isDouble = isDouble, textsize = textsize, num = num+1)
	else
		f!(q1,offSet = O + O1, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse, isDouble = isDouble, textsize = textsize, num = num)
		f!(q1,offSet = O + O1+vq1, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse, isDouble = isDouble, textsize = textsize, num = num)
		f!(q1,offSet = O + O2, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse, isDouble = isDouble, textsize = textsize, num = num)
		f!(q1,offSet = O + O2+vq1, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse, isDouble = isDouble, textsize = textsize, num = num)
		f!(p1,offSet = O + O3, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse, isDouble = isDouble, textsize = textsize, num = num+1)
		f!(p1,offSet = O + O4, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse, isDouble = isDouble, textsize = textsize, num = num+1)
	end
	return plot!()


end

function double_decomp_plot_and_annotate_std_tile(p::Rational{Int64}; width = 1, offSet = 0, isDouble = false, textsize = 10,  num = 1, isRecurse = false, isAnnotate=false)
	plot(aspect_ratio = true, axis = :off, ticks = :none)
	return double_decomp_plot_and_annotate_std_tile!(p, width = width, offSet = offSet, isDouble = isDouble, textsize = textsize,  num = num, isRecurse = isRecurse, isAnnotate=isAnnotate)
end
