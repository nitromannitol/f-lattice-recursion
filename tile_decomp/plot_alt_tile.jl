## requires alt_tile word code


################################
#### comptue the center of an alternate tile
################################


function compute_corner_offset_odd(p1::Rational{Int64},q1::Rational{Int64}, str::String)
	vp1,vp2 = vf(p1)
	vq1,vq2 = vf(q1)
	W1 = count_full(str, '1')%2;
	if(W1 == 0)
		return 0
	else
		return -vq1;
	end
end

function compute_corner_offset_even(p1::Rational{Int64},q1::Rational{Int64}, str::String)
	vp1,vp2 = vf(p1)
	vq1,vq2 = vf(q1)
	W1 = count_full(str, '1')%2;
	if(W1 == 0)
		return -vq1
	else
		return 0

	end
end

function compute_alttile_center(p::Rational{Int64})
	if(p == 0//1 || p == 1//1) error("alt tiles don't exist for 0/1 and 1/1"); end
	str = rational_approximation(p,9999)[1]
	p0,q0,p1,q1 = getFracs(str); 


	vp1,vp2 = vf(p1)
	vq1,vq2 = vf(q1)

	tile_x = 2vp1 + 2vq1; 
	tile_y = vq2 + 2vp2; 

	W1 = count_full(str, '1')%2;

	if( (p.num+p.den +W1)%2 == 1)
		corner_offset = 0
	else
		corner_offset = -vq1;
	end
	return corner_offset + (tile_x)/2 + (tile_y)/2; 

end


function plot_and_annotate_single_alt!(p::Rational{Int64}; num = 1, offSet=0, width = 1, textsize=5, isAnnotate=true, ann = 1)
	if(p == 0//1 || p == 1//1) error("alt tiles don't exist for 0/1 and 1/1"); end

	str = rational_approximation(p,9999)[1]
	t0 = get_alt_tile(str,p);
	t0 = t0 .+offSet;
	plot_t!(t0, label = "$(p) alt", color = "gray", width = width)
	
	
	if(isAnnotate) ##annotate
		O = t0[1] + compute_alttile_center(p)
		if(isOdd(p)) 
			ann = LaTeXString("\$\\hat{T}(p_$(num))\$"); 
		else
			ann = LaTeXString("\$\\hat{T}(q_$(num))\$"); 
		end
		annotate!(real(O),imag(O), text(ann, :black, :center, textsize))
	end
	return plot!()

end