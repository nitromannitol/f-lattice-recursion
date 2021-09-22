## uses functions from annotate_std_tile.jl
## and get_alt_tile_offset.jl
## and tile_offsets



function plot_and_annotate_alt_tile!(p::Rational{Int64}; annotate_center=true,width = 1, offSet = 0, textsize = 10, num = 1, isAnnotate = true, isRecurse = false, isDoubleDecomp=false, isPlotParent=false)
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


	if(isEven(p))
		if(isLCorrect)

			O_ = O1+offSet;
			plot_and_annotate_corrected_tile!(q1, offSet = O_, isFlip = false, isLong = true, annotate_center=annotate_center,width = width,  textsize = textsize, num = num+1, isAnnotate = isAnnotate, isDoubleDecomp=isDoubleDecomp, isPlotParent=isPlotParent)

			O_ = O2+Od+offSet;
			plot_and_annotate_corrected_tile!(q1, offSet = O_, isFlip = true, isLong = true, annotate_center=annotate_center,width = width,  textsize = textsize, num = num+1, isAnnotate = isAnnotate, isDoubleDecomp=isDoubleDecomp, isPlotParent=isPlotParent)


			plot_and_annotate(tq1, q1, O = O1+(Od), label = q1, textsize = textsize,num = num, ann = ano2)  
			plot_and_annotate(tq1, q1, O = O2, label = q1, textsize = textsize,num = num, ann = anx1)  


		else
			plot_and_annotate(tq1, q1, O = O1, label = q1, textsize = textsize,num = num, ann = ano1)  
			plot_and_annotate(tq1, q1, O = O1+(Od), label = q1, textsize = textsize,num = num, ann = ano2)  
			plot_and_annotate(tq1, q1, O = O2, label = q1, textsize = textsize,num = num, ann = anx1)  
			plot_and_annotate(tq1, q1, O = O2+(Od), label = q1, textsize = textsize,num = num, ann = anx2)  
		end

		plot_and_annotate(tp1, p1, O = O3, label = p1, textsize = textsize,num = num, ann = ano)  
		plot_and_annotate(tp1, p1, O = O4, textsize = textsize,num = num, ann = anx)
	else
		plot_and_annotate(tq1, q1, O = O1, label = q1, textsize = textsize,num = num, ann = ano)  
		plot_and_annotate(tq1, q1, O = O2, label = q1, textsize = textsize,num = num, ann = anx)  

		if(isLCorrect)
			O_ = O3+Od+offSet;
			plot_and_annotate_corrected_tile!(p1, offSet = O_, isFlip = true, isLong = true, annotate_center=annotate_center,width = width,  textsize = textsize, num = num+1, isAnnotate = isAnnotate, isDoubleDecomp=isDoubleDecomp, isPlotParent=isPlotParent)

			O_ = O4+offSet;
			plot_and_annotate_corrected_tile!(p1, offSet = O_, isFlip = false, isLong = true, annotate_center=annotate_center,width = width,  textsize = textsize, num = num+1, isAnnotate = isAnnotate, isDoubleDecomp=isDoubleDecomp, isPlotParent=isPlotParent)


			plot_and_annotate(tp1, p1, O = O3, label = p1, textsize = textsize,num = num, ann = ano1)  
			plot_and_annotate(tp1, p1, O = O4+(Od), textsize = textsize,num = num, ann = anx2)  

		else

			plot_and_annotate(tp1, p1, O = O3, label = p1, textsize = textsize,num = num, ann = ano1)  
			plot_and_annotate(tp1, p1, O = O3+(Od), textsize = textsize,num = num, ann = ano2)  

			plot_and_annotate(tp1, p1, O = O4, textsize = textsize,num = num, ann = anx1)
			plot_and_annotate(tp1, p1, O = O4+(Od), textsize = textsize,num = num, ann = anx2)  


		end
	end


	#double decomposition
	if(isLCorrect) return plot!(); end #double overlap case

	##otherwise plot alterante
	if(isOdd(p))
		plot_and_annotate_single_alt!(q1, offSet=offSet+get_odd_offset(p1,q1,str), textsize=textsize, width=width,isAnnotate=annotate_center, num = num)
	else
		plot_and_annotate_single_alt!(p1, offSet=offSet+get_even_offset(p1,q1,str), textsize=textsize, width=width,isAnnotate=annotate_center, num = num)
	end

	return plot!()


end



function plot_and_annotate_alt_tile(p::Rational{Int64}; annotate_center=true,width = 1, offSet = 0, textsize = 10, num = 1, isAnnotate = true, isRecurse = false)
	plot(aspect_ratio = true, axis = :off, ticks = :none)
	return plot_and_annotate_alt_tile!(p,annotate_center=annotate_center,width = width, offSet = offSet, textsize = textsize, num = num, isAnnotate = isAnnotate, isRecurse = isRecurse)
end




## subscript num are the two parents of the even parent
## subscript num+1 are the two parents of the odd parent
function double_decomp_plot_and_annotate_alt_tile!(p::Rational{Int64}; isDoubleDecomp = false, width = 1, offSet = 0, textsize = 10,  num = 1, isRecurse = false, isAnnotate=false)
	if(p == 0//1 || p == 1//1) return plot!(); end
	str = rational_approximation(p,9999)[1]; # be careful
	if(length(str) ==1)
		return plot_and_annotate_std_tile!(p,  width = width, offSet = offSet,  textsize = textsize, num = num+1)
	end
	p0,q0,p1,q1 = getFracs(str); 



	if(isDoubleDecomp)


		O1,O2,O3,O4,Od = getAltTileOffSets(p, str);
		if(isRecurse == true)
			f! = double_decomp_plot_and_annotate_std_tile!;
			textsize = textsize -1; 
			num = num + 1; 
			textsize = 1; 
		else
			f! = plot_and_annotate_std_tile!;
		end
		O = offSet;

		if(isEven(p))
			f!(q1,offSet = O + O1, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse,  textsize = textsize, num = num)
			f!(q1,offSet = O + O1+Od, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse,  textsize = textsize, num = num)
			f!(q1,offSet = O + O2, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse,  textsize = textsize, num = num)
			f!(q1,offSet = O + O2+Od, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse,  textsize = textsize, num = num)
			f!(p1,offSet = O + O3, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse,  textsize = textsize, num = num+1)
			f!(p1,offSet = O + O4, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse,  textsize = textsize, num = num+1)

		else
			f!(q1,offSet = O + O1, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse,  textsize = textsize, num = num)
			f!(q1,offSet = O + O2, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse,  textsize = textsize, num = num)
			f!(p1,offSet = O + O3, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse,  textsize = textsize, num = num+1)
			f!(p1,offSet = O + O3+Od, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse,  textsize = textsize, num = num+1)
			f!(p1,offSet = O + O4, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse,  textsize = textsize, num = num+1)
			f!(p1,offSet = O + O4+Od, width = width, isAnnotate = isAnnotate, isRecurse = isRecurse,  textsize = textsize, num = num+1)

		end
	else
		plot_and_annotate_alt_tile!(p, annotate_center=false,offSet =offSet, width = width,textsize = textsize, num = num, isAnnotate = isAnnotate, isRecurse = isRecurse)
	end

	num+=1; #fix this later to have more descriptive labels


	#avoid this case for now
	if(str[end] == '2' && isOdd(p) || str[end] == '3' && isEven(p)) return plot!(); end #double overlap case


	p_alt = p1; if(isOdd(p)) p_alt = q1; end
	plot_and_annotate_alt_tile!(p_alt, offSet =offSet + get_alt_offset(p, str), width = width,textsize = textsize, num = num+1, isAnnotate = isAnnotate, isRecurse = isRecurse)





	return plot!()


end

function double_decomp_plot_and_annotate_alt_tile(p::Rational{Int64}; isDoubleDecomp = false, width = 1, offSet = 0, isDouble = false, textsize = 10,  num = 1, isRecurse = false, isAnnotate=false)
	plot(aspect_ratio = true, axis = :off, ticks = :none)
	return double_decomp_plot_and_annotate_alt_tile!(p, isDoubleDecomp = isDoubleDecomp, width = width, offSet = offSet, textsize = textsize,  num = num, isRecurse = isRecurse, isAnnotate=isAnnotate)
end
