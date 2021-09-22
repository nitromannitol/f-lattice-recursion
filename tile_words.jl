plot_t!(t0,p=0; color = "red", label = false, width = 2, linealpha = 1)  = plot!(real(t0),imag(t0),label=label,color = color, legend = :bottomleft, width = width, linealpha = linealpha)


##builds the bottom edge of a tile
function get_bottom_edge(p)  
	if(p <= 0 || p >= 1) return []; end 
	str,_=rational_approximation(p,9999)
	p1,q1,p0,q0 = getFracs(str)
	k = length(str)
	#only 2
	if(all2(str))
		if(isOdd(p))
			w1 = [1]
			for i in 1:(2*k+1)
				w1 = [w1; 1; im]
			end
			w1 = [w1; 1; 1]
		else
			w1 = [1]
			for i in 1:2*k
				w1 = [w1; 1; im]
			end
			w1 = [w1; 1; 1]
		end


	#only 3 
	elseif(all3(str))
		if(isOdd(p))
			w1 = Array{Complex{Int64},1}()
			for i in 1:k
				w1 = [w1; 1; 1]
			end
			w1 = [w1; 1; im]
			for i in 1:(k+1)
				w1 = [w1; 1; 1]
			end
			w1 = [w1; 1]

		else
			w1 = [1]
			for i in 1:(k+1)
				w1 = [w1; 1; 1]
			end
		end
	else
		w1,v1 = get_bottom_edge(p0), get_bottom_edge(q0); 

		if(isOdd(p))
			W1 = [w1; im; v1; im; v1]
		else
			W1 = [w1; im; v1]
		end

		if(count_full(str, '1')%2 == 0)
			if(isOdd(p))
				W1 = [v1; im; v1; im; w1]
			else
				W1 = [v1; im; w1]
			end

		end
		w1 = W1;
	end

	return w1;

end

function getWords(p)
	w1 = get_bottom_edge(p)
	w2 = im .* get_bottom_edge(getFlippedRep(p));
	return w1,w2;
end

function getTile(p)
	if(p == 0//1)
		t0 = cumsum([0; 1; im; im; -1; -im; -im])
	elseif(p == 1//1)
		t0 = cumsum([0; 1; im; -1; -im])
	end
	if(p == 0//1 || p == 1//1) return t0; end

	w1,w2 = getWords(p)
	t0 = [w1; w2; involution(reverse(w1)); involution(reverse(w2))];
	t0 = cumsum([0;t0])
	return t0; 

end






### ALT TILES


function get_alt_tile(str, p) 
	kk = count_full(str, '1')%2; 

	if(p == 0//1 || p ==1//1)
		return getTile(p)
	end

	if(all2(str) || all3(str))
		return get_base_alt_tile(str,p)
	end

	p0,q0,p1,q1 = getFracs(str);
	wp1, wp2 = getWords(p1);
	wq1, wq2 = getWords(q1); 
	w1,w2,w3, w4 = map(z->Array{Complex{Int64},1}(),[1,2,3,4]);

	if(!isOdd(p))
		if(kk == 1)
			w1 = [wp1; im; wq1];
			w2 = [wp2; -1];
			w3 = [im; wq1]
			w4 = [wq2; -1; wp2]
		else
	
			w1 = [wq1; im; wp1]
			w2 = [wp2; -1; wq2];
			w3 = [involution(reverse(wq1)); -im];
			w4 = [-1; wp2]

		end
	else
		if(kk==1)
			w1 = [wp1; im; wq1;]
			w2 = [wq2; -1; wp2];
			w3 = [involution(reverse(wq1)); -im]
			w4 = [-1; wp2]
		else
			w1 = [wq1; im; wp1];
			w2 = [wp2; -1]
			w3 = [im; wq1]
			w4 = [wp2; -1; wq2]
		end
	end


	t0 = [w1; w2; w3; w4; involution(reverse(w1)); involution(reverse(w2)); involution(reverse(w3)); involution(reverse(w4))]
	t0 = cumsum([0; t0]);
	return t0;

end


function get_base_alt_tile(str, p)
	if( p != 1//3 && all2(str)) # so that 1/3 goes to the correct one

		if(isOdd(p))


			q = getFlippedRep(p); 
			d = q.den; 
			w1 = [1; 1] 
			for i in 1:(d-3)
				w1 = [w1; im; 1]
			end
			w1 = [w1; 1; 1]
			w2 = [im for i in 1:d-1]; 
			w2 = [w2; -1]; 
			for i in 1:d; w2 = [w2; im]; end

			t0 = [w1; w2; involution(reverse(w1)); involution(reverse(w2))];
			t0 = cumsum([0;t0])


		else
			p = getFlippedRep(p); 
			d = p.den; 

			w1 = [1; 1]
			for i in 1:(d-2) 
				w1 = [w1; im; 1]
			end
			w1 = [w1; 1]
			w2 = [im for i in 1:d]
			w2 = [w2; -1; -1]
			for i in 1:(d-1)
				w2 = [w2; im]
			end

			t0 = [w1; w2; involution(reverse(w1)); involution(reverse(w2))]
			t0 = cumsum([0;t0])


	
		end
		return t0; 
	elseif(p!= 1//2 && all3(str))
		if(isOdd(p))

		d = p.den

		w1 = [1 for i in 1:d]
		w1 = [w1; im; im;]
		for i in 1:(d-1)
			w1 = [w1; 1]
		end

		w2 =[im; im]
		for i in 1:(d-2)
			w2 = [w2; -1; im]
		end
		w2 = [w2; im]

		t0 = [w1; w2; involution(reverse(w1)); involution(reverse(w2))]
		t0 = cumsum([0;t0])
		

		else
			 d = p.den; 
			 #oriented at correct corner
			w1 = [1 for i in 1:d]
			w2 = [im; im] 
			for i in 1:(d-3)
				w2 = [w2; -1; im]
			end
			w2 = [w2; im; im]
			w3 = [-1 for i in 1:d-1];
			w3 = [w3; -im]


			 t0 = [w1; w2; w3; involution(reverse(w1)); involution(reverse(w2)); involution(reverse(w3))];
			 t0 = cumsum([0;t0])


			#oriented at lower left corner
			# w1 = [1 for i in 1:d-1]
			# w1 = [w1; im]
			# for i in 1:d  w1 = [w1; 1]; end

			# w2 = [im; im] 
			# for i in 1:(d-3)
			# 	w2 = [w2; -1; im]
			# end
			# w2 = [w2; im; im]
			# t0 = [w1; w2; involution(reverse(w1)); involution(reverse(w2))];
			# t0 = cumsum([0;t0])
	

		end
		return t0; 
	end
end



##plots a tile and draws the pattern
function full_plot_std(A,p; isAnnotate=false, width = 1, isTile = false, isDecomp = false)
	A = unpad(pad(A)[1])[1]
	if(isTile)
		C0,B = get_tiling(A,p)
		C0 = unpad(C0)[1]
		x,y = get_lower_left_vertex(C0);
		drawPattern(C0)
		O = x + y*im; 
		v1,v2 = vf(p); if(isOdd(p)) v1 = v1*2; end
		for i in 0:2, j in 0:2
			if(isDecomp)
				plot_and_annotate_std_tile!(p, offSet = O+v1*i+v2*j, isAnnotate=isAnnotate, width = width); 
			else
				plot_and_annotate_single_std!(p, offSet = O+v1*i+v2*j, isAnnotate=isAnnotate, width = width); 

			end
		end
	else
		drawPattern(A,p=p, isTile = isTile)
		x,y = get_lower_left_vertex(A); O = x+y*im;
		if(isDecomp)
			plot_and_annotate_std_tile!(p, offSet = O, isAnnotate=isAnnotate, width = width); 
		else
			plot_and_annotate_single_std!(p, offSet = O, isAnnotate=isAnnotate, width = width); 

		end
	end
	return plot!(legend=false)

end

function full_plot_alt(B,q; isAnnotate=false, width = 1, isTile = false, isDecomp = true)
	B = unpad(pad(B)[1])[1]

	if(isTile)
		C0,B = get_tiling(B,q)
		C0 = unpad(C0)[1]
		x,y = get_lower_left_vertex(C0);
 		O = x+y*im; O = O + translate_alts(q);
 		drawPattern(C0); 
		v1,v2 = vf(q); if(isOdd(q)) v1 = v1*2; end
		for i in 0:2, j in 0:2
			if(isDecomp)
 				plot_and_annotate_alt_tile!(q, offSet = O+v1*i+v2*j, isAnnotate=isAnnotate, width = width);
 			else
 				plot_and_annotate_single_alt!(q, offSet = O+v1*i+v2*j, isAnnotate=isAnnotate, width = width);
 			end
		end

	else
		drawPattern(B,p=q)
		x,y = get_lower_left_vertex(B); O = x+y*im; 
		O = O + translate_alts(q);
		if(isDecomp)
 			plot_and_annotate_alt_tile!(q, offSet = O, isAnnotate=isAnnotate, width = width);
 		else
 			plot_and_annotate_single_alt!(q, offSet = O, isAnnotate=isAnnotate, width = width);
 		end
 	end
	return plot!(legend=false)

end
