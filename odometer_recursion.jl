#s is the type of child
function iterate_mat(mats, alt_mats, str, s)
	Gp1,Gq1, Gp0,Gq0, Gp1_,Gq1_,Gp0_,Gq0_  = center_mats(mats...,alt_mats...)
	p1,q1,p0,q0 = getFracs(str); 

	str = str*s

	## pick the two parents 
	if(str[end]=='1')
		Gp1,Gq1,Gp1_,Gq1_ = Gp1,Gq1,Gp1_,Gq1_
		p1,q1 = p1,q1
	elseif(str[end] == '2')
		Gp1,Gq1,Gp1_,Gq1_ = Gp1,Gq0,Gp1_,Gq0_
		p1,q1 = p1,q0
	elseif(str[end] == '3')
		Gp1,Gq1,Gp1_,Gq1_ = Gp0,Gq1,Gp0_,Gq1_
		p1,q1 = p0,q1
	end


	Gp1,Gq1,Gp1_,Gq1_ = pad(Gp1,Gq1,Gp1_,Gq1_)
	vp1,vp2 = vf(p1); vq1,vq2 = vf(q1); vp1 = 2vp1;
	ap1,ap2 = af(p1); aq1,aq2 = af(q1); ap1 = 2ap1;

	mats = Gp1,Gq1, Gp1_, Gq1_
	vs = vp1,vp2, vq1, vq2; 
	as = ap1,ap2, aq1, aq2;


	return mats, str

end

#get the affine factors of the parents
function getAffine(str)
	p0,q0,p1,q1 = getFracs(str); 
	vp1,vp2 = vf(p1); vq1,vq2 = vf(q1); vp1 = 2vp1;
	ap1,ap2 = af(p1); aq1,aq2 = af(q1); ap1 = 2ap1;
	return vp1,vp2, vq1, vq2, ap1, ap2, aq1, aq2;
end





## produces standard children when no correction is needed
function get_std_child_unc(mats, str; getOdd = true)
	W1 = (sum([s == '1' for s in str]))%2
	Gp1,Gq1, Gp1_, Gq1_ = mats

	if(getOdd)
		vs,as = getOddTileFull(str)
		curr_mats = Gq1,Gq1,Gq1,Gq1,Gp1,Gp1
	else
		vs,as = getEvenTileFull(str)
		curr_mats = Gq1,Gq1,Gp1,Gp1,Gp1,Gp1
	end

	C0,C1,C2,C3,C4,C5 = map(translate_va, curr_mats, vs, as)

	## set it up so things overlap
	if(getOdd && W1 == 1) 
		#C2,C3 = C3,C2; 
		C2,C3,C4,C5 = C4,C5,C2,C3
	elseif(getOdd && W1 == 0 )
		C2,C3,C4,C5 = C4,C5,C2,C3
	end 
	if(!getOdd) C1,C2,C3 = C2,C3,C1; end

	for C1 in [C1,C2,C3,C4,C5]
		C0 = min_nan.(C0, C1 .- minimum(nan_ign, (C1 - C0))); 
	end
	C0 = unpad(C0)[1]
	return C0
end


function translate_alts(mats, str)
	## check if a translation is needed
	Gp1,Gq1, Gp1_, Gq1_ = mats

	#check Gq1_
	q = getFracs(str)[4]; 
	if(q != 1//1)
		str_q, fracs_q = rational_approximation(q, 999);
		p2,q2,p3,q3 = fracs_q

		if(q2 != q) error("prob_q"); end
		W_q = (sum([s == '1' for s in str_q]))%2
		if(W_q == 1) # no correction needed

		else
			v1,v2 = vf(p3); v1 = v1*2; 
			v3,v4 = vf(q3); 
			O = v3-v2; 
			if(imag(O) < 0) #no correction needed

			else
				if(all3(str_q)) O = O + im; end
				Gq1_ = translate_va(Gq1_,-O,0)
			end

		end
	end


	p = getFracs(str)[3]; 
	if(p != 0//1)
		str_p, fracs_p = rational_approximation(p, 999);
		p2,q2,p3,q3 = fracs_p
		if(p2 != p) error("prob_p"); end
		W_p = (sum([s == '1' for s in str_p]))%2
		if(W_p == 0) #no correction needed

		else
			v1,v2 = vf(p3); v1 = v1*2; 
			v3,v4 = vf(q3); 
			O = v3-v2
			if(imag(O) < 0) 
			else
				Gp1_ = translate_va(Gp1_, -O, 0)
			end

		end
	end

	return Gp1_, Gq1_;



end



## produces alternate children when no correction is needed
function get_alt_child_unc(mats, str; getOdd=true)
	W1 = (sum([s == '1' for s in str]))%2
	Gp1,Gq1, Gp1_, Gq1_ = mats
	Gp1_,Gq1_ = translate_alts(mats, str);

	if(getOdd)
		vs,as = getAltOddTileFull(str)
		curr_mats = Gq1,Gq1,Gp1,Gp1,Gp1,Gp1
	else
		vs,as = getAltEvenTileFull(str)
		curr_mats = Gq1,Gq1,Gq1,Gq1,Gp1,Gp1
	end

	C0,C1,C2,C3,C4,C5 = map(translate_va, curr_mats, vs, as)

	if(getOdd)
		G_ = Gq1_
		v_, a_ = get_alt_full_offset(getFracs(str)[1],str)
		C1,C2,C3 = C2,C3,C1
	else
		G_ = Gp1_
		v_, a_ = get_alt_full_offset(getFracs(str)[2],str)
		C2,C3,C4 = C4,C2,C3
	end


	C_ = translate_va(G_, v_,a_)
	for C1 in [C_, C1,C2,C3,C4,C5]
		C0 = min_nan.(C0, C1 .- minimum(nan_ign, (C1 - C0))); 
	end

	C0 = unpad(C0)[1]
	return C0

end


## plot a quadruple
function plot_mats(mats, str; fullPlot = false)
	fracs= getFracs(str); fracs = getFracs(str)
	if(fullPlot == false)
		plt1 = drawPattern(mats[1])
		plt2 = drawPattern(mats[2])
		plt3 = drawPattern(mats[3])
		plt4 = drawPattern(mats[4])
		plts = plt1,plt2,plt3,plt4
	else
		f(A,p) = drawPattern(A,p=p);
		plts = map(f, mats, fracs);
	end
	return plot(plts...)
end




## get corrected decompositions

function get_std_child_c(mats, l_mats, str; getOdd = false)
	W1 = (sum([s == '1' for s in str]))%2
	Gp1,Gq1, Gp1_, Gq1_ = mats


	Gl,Gr = l_mats
	#translate l corrections
	if(!getOdd)
		x2,y2 = get_upper_right_vertex(Gp1);
	else
		x2,y2 = get_upper_right_vertex(Gq1);
	end
	if(str[end] != '3')
		if(W1 == 1)
			x1,y1 = get_upper_right_vertex(Gr)
			O = (x2-x1) + (y2-y1)*im;
			Gr = translate_va(Gr, O, 0)
		else
			x1,y1 = get_upper_right_vertex(Gl)
			O = (x2-x1) + (y2-y1)*im;
			Gl = translate_va(Gl, O, 0)
		end
	else
		x1,y1 = get_upper_right_vertex(Gl)
		O = (x2-x1) + (y2-y1)*im;
		Gl = translate_va(Gl, O, 0)
	end


	if(getOdd)
		vs,as = getOddTileFull(str)
		if(W1 == 0) curr_mats = Gq1,Gr,Gl,Gq1,Gp1,Gp1; end
		if(W1 == 1) curr_mats = Gr,Gq1,Gq1,Gl,Gp1,Gp1; end #overlap
	else
		vs,as = getEvenTileFull(str)
		curr_mats = Gq1,Gq1,Gp1,Gl,Gr,Gp1
	end

	C0,C1,C2,C3,C4,C5 = map(translate_va, curr_mats, vs, as)

	## set it up so things overlap
	if(getOdd)
		C2,C3,C4,C5 = C4,C5,C2,C3
	end 
	if(!getOdd) C1,C2,C3 = C2,C3,C1; end

	for C1 in [C1,C2,C3,C4,C5]
	
		C0 = min_nan.(C0, C1 .- minimum(nan_ign, (C1 - C0))); 
	end
	C0 = unpad(C0)[1];

	return C0


end


function get_alt_child_c(mats, l_mats, str; getOdd=true)
	W1 = (sum([s == '1' for s in str]))%2
	Gp1,Gq1, Gp1_, Gq1_ = mats
	#Gl, Gr = correct_l_offset(mats, l_mats, str)
	Gl,Gr = l_mats


	#translate l corrections
	if(getOdd)
		x2,y2 = get_upper_right_vertex(Gp1);

	else
		x2,y2 = get_upper_right_vertex(Gq1);
	end

	if(str[end] != '3')
		if(W1 == 1)
			x1,y1 = get_upper_right_vertex(Gr)
			O = (x2-x1) + (y2-y1)*im;
			Gr = translate_va(Gr, O, 0)
		else
			x1,y1 = get_upper_right_vertex(Gl)
			O = (x2-x1) + (y2-y1)*im;
			Gl = translate_va(Gl, O, 0)
		end
	else
		x1,y1 = get_upper_right_vertex(Gl)
		O = (x2-x1) + (y2-y1)*im;
		Gl = translate_va(Gl, O, 0)
	end



	#pick subtiles
	if(getOdd)
		vs,as = getAltOddTileFull(str)
		curr_mats = Gq1,Gq1,Gp1,Gl,Gr,Gp1
	else
		vs,as = getAltEvenTileFull(str)
		curr_mats = Gr,Gq1,Gq1,Gl,Gp1,Gp1;
	end
	#translate
	C0,C1,C2,C3,C4,C5 = map(translate_va, curr_mats, vs, as)


	#change order so things overlap
	C1,C2,C3,C4,C5 = C2,C3,C4,C5, C1;


	for C1 in [C1,C2,C3,C4,C5]
		C0 = min_nan.(C0, C1 .- minimum(nan_ign, (C1 - C0))); 
	end

	C0 = unpad(C0)[1]

	return C0

end


## l corrections

## this produce sthe l corrected version of a standard tile

#we use the parent Farey pair of str
#so guranteed that length(str) > 0
function get_l_correct(mats, alt_mats, ancestor_str, curr_str; getOdd=true, isFlip = false, isLong = false)
	k = length(curr_str)-length(ancestor_str)

	Gp1,Gq1, Gp2,Gq2, Gp1_,Gq1_,Gp2_,Gq2_  = center_mats(mats...,alt_mats...)
	

	p1,q1,p0,q0 = getFracs(ancestor_str); 
	vp1,vp2, vq1, vq2, ap1, ap2, aq1, aq2 = getAffine(ancestor_str); 
	num_iters = 2*(k+1); 

	Gp1,Gq1, Gp2,Gq2, Gp1_,Gq1_,Gp2_,Gq2_ = pad(Gp1,Gq1, Gp2,Gq2, Gp1_,Gq1_,Gp2_,Gq2_);

	curr_mats = Gp2,Gq2, Gp2_, Gq2_;
	Gp2_,Gq2_ = translate_alts(curr_mats, curr_str);


	if(!getOdd)
		K1,K2,K3,K4,A1,A2,A3,A4 = get_even_l_offsets(ancestor_str, isFlip)
		G1,G2,G3 = pad(Gp2, Gq2, Gp2_); 
		#if(ancestor_str[end] == '2' && !isLong) num_iters-=1; end
		if(!isLong) num_iters-=1; end
	else
		K1,K2,K3,K4,A1,A2,A3,A4 = get_odd_l_offsets(ancestor_str, isFlip)
		G1,G2,G3 = pad(Gq2, Gp2, Gq2_); 
		#if(ancestor_str[end] == '3' && !isLong) num_iters-=1; end
		if(!isLong) num_iters-=1; end #testing 
	end


	Ov = 0; Oa = 0;
	## initial offset
	if(getOdd)
		(v1,v2,v3,v4,v5,v6),(a1,a2,a3,a4,a5,a6) = getOddTileFull(curr_str) 
	else
		(v1,v2,v3,v4,v5,v6),(a1,a2,a3,a4,a5,a6) = getEvenTileFull(curr_str) 
	end
	if(isFlip && getOdd)
		Ov = Ov + v4
		Oa = Oa + a4
	elseif(isFlip && !getOdd)
		Ov = Ov + v6;
		Oa = Oa + a6; 
	end
	C0 = translate_va(G1,Ov, Oa); 
	for i in 1:k
		Ov = Ov + K1; Oa = Oa + A1; 
		C1 = translate_va(G1, Ov, Oa);
		C0 = min_nan.(C0, C1 .- minimum(nan_ign, (C1 - C0))); 
	end
	Ov = Ov + K2; Oa = Oa + A2; 
	C1 = translate_va(G2, Ov, Oa);

	## shifting the base case tiles
	if(all3(ancestor_str) && getOdd && !isFlip) C1 = translate_va(C1, -im,0); end;
	if(all3(ancestor_str) && getOdd && isFlip) C1 = translate_va(C1, -im+1,0) end; 
	if(all2(ancestor_str) && !getOdd && isFlip) C1 = translate_va(C1, 1+im,0); end;
	C0 = min_nan.(C0, C1 .- minimum(nan_ign, (C1 - C0))); 

	Ov = Ov + K3; Oa = Oa + A3; 
	for i in 1:num_iters
		Ov = Ov + K4; Oa = Oa + A4; 
		C1 = translate_va(G3, Ov, Oa);
		C0 = min_nan.(C0, C1 .- minimum(nan_ign, (C1 - C0))); 
	end
	C0 = unpad(C0)[1];
	C0 = C0 .+ minimum(nan_ign, -C0);
	return C0; 

end





function get_corrections(curr_str, anc_data2,anc_data3, isLong = false)
	Gp1l = get_l_correct(anc_data2..., curr_str, getOdd = true, isFlip = true, isLong = isLong); 
	Gp1r = get_l_correct(anc_data2..., curr_str, getOdd = true, isFlip = false, isLong = isLong); 

	Gq1l = get_l_correct(anc_data3..., curr_str, getOdd = false, isFlip = true, isLong = isLong); 
	Gq1r = get_l_correct(anc_data3..., curr_str, getOdd = false, isFlip = false, isLong = isLong); 


	l_mats = Gp1l, Gp1r, Gq1l,Gq1r; 
	return l_mats; 
end



#rec_data = mats,alt_mats, l_mats, l_mats_l, anc_data2, anc_data3;
function get_next_mat(rec_data, curr_str, s)
	mats,alt_mats, l_mats, l_mats_l, anc_data2, anc_data3 = rec_data

	Gp1,Gq1, Gp0,Gq0, Gp1_,Gq1_,Gp0_,Gq0_, Gp1l, Gp1r, Gq1l,Gq1r, Gp1l_, Gp1r_, Gq1l_,Gq1r_  = center_mats(mats...,alt_mats..., l_mats..., l_mats_l...)
	p1,q1,p0,q0 = getFracs(curr_str); 


	## pick the two parents 
	if(s=='1')
		Gp1,Gq1,Gp1_,Gq1_ = Gp1,Gq1,Gp1_,Gq1_
		Gl,Gr,Gl_,Gr_ = Gp1,Gq1,Gp1_,Gq1_ 
		p1,q1 = p1,q1
	elseif(s == '2')
		Gp1,Gq1,Gp1_,Gq1_ = Gp1,Gq0,Gp1_,Gq0_
		Gl,Gr = Gp1l, Gp1r; 
		Gl_,Gr_ = Gp1l_,Gp1r_
		p1,q1 = p1,q0
	elseif(s == '3')
		Gp1,Gq1,Gp1_,Gq1_ = Gp0,Gq1,Gp0_,Gq1_
		Gl,Gr = Gq1l, Gq1r; 
		Gl_,Gr_ = Gq1l_,Gq1r_
		p1,q1 = p0,q1
	end
	Gp1,Gq1,Gp1_,Gq1_, Gl, Gr, Gl_,Gr_ = pad(Gp1,Gq1,Gp1_,Gq1_, Gl, Gr, Gl_,Gr_)
	mats = Gp1,Gq1, Gp1_, Gq1_;
	if(s != '1') l_mats = Gl,Gr; l_mats_l = Gl_,Gr_; end; #otherwise no correction needed


	# create the alternate and standard child
	if(s == '1')
		Gp2 = get_std_child_unc(mats, curr_str, getOdd = true);
		Gq2 = get_std_child_unc(mats, curr_str, getOdd = false);
		Gp2_ = get_alt_child_unc(mats, curr_str, getOdd = true);
		Gq2_ = get_alt_child_unc(mats, curr_str, getOdd = false);
	elseif(s == '2')
		Gp2 = get_std_child_unc(mats, curr_str, getOdd = true);
		Gq2 = get_std_child_c(mats, l_mats, curr_str, getOdd = false);
		Gp2_ = get_alt_child_c(mats, l_mats_l, curr_str, getOdd = true);
		Gq2_ = get_alt_child_unc(mats,curr_str, getOdd = false);
	elseif(s == '3')
		Gp2 = get_std_child_c(mats, l_mats, curr_str, getOdd = true);
		Gq2 = get_std_child_unc(mats, curr_str, getOdd = false);
		Gp2_ = get_alt_child_unc(mats, curr_str, getOdd = true);
		Gq2_ = get_alt_child_c(mats,l_mats_l, curr_str, getOdd = false);
	end


	#update the quadruple
	mats = Gp2,Gq2, Gp1, Gq1
	alt_mats = Gp2_, Gq2_, Gp1_, Gq1_



	####
	#update ancestor data for corrections
	####
	if(s == '1')
		anc_data2 = mats,alt_mats,curr_str
		anc_data3 = mats,alt_mats,curr_str
	elseif(s == '2')
		anc_data3 = mats,alt_mats,curr_str
	elseif(s == '3')
		anc_data2 = mats,alt_mats,curr_str
	end

	######
	#build_l corrections
	#####
	l_mats = get_corrections(curr_str, anc_data2,anc_data3, false)
	l_mats_l = get_corrections(curr_str, anc_data2, anc_data3, true)


	return mats,alt_mats, l_mats, l_mats_l, anc_data2, anc_data3

end





## given the helper codes in odometer_recursion 
## produce odometer for given string
function get_odometers(str)
	curr_str = get_first_ancestor(str)
	mats, alt_mats = get_base_quadruple(curr_str)
	curr_ind = length(curr_str)+1;

	anc_data2 = mats,alt_mats,curr_str;  #only update if there is a 3 or 1
	anc_data3 = mats,alt_mats,curr_str;  #only update if there is a 2 or 1


	## lcorrections for the child pair
	#one of these is not used in each step (and is wrong in the first step lol)
	if(length(curr_str) > 0)
		l_mats = get_corrections(curr_str, anc_data2, anc_data3, false);
		l_mats_l = get_corrections(curr_str, anc_data2, anc_data3, true);

	else
		l_mats = mats; 
		l_mats_l = mats;
	end
	rec_data = mats,alt_mats, l_mats, l_mats_l, anc_data2, anc_data3;

	i = curr_ind
	for i in curr_ind:length(str)
		curr_str = str[1:i];
		s = curr_str[end];
		rec_data = get_next_mat(rec_data, curr_str, s)
	end
	mats,alt_mats, l_mats, l_mats_l, anc_data2, anc_data3 = rec_data

	return mats, alt_mats; 

end