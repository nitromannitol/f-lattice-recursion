## helper functions to compute standard tile offsets 


################################################
################################################
### 			STANDARDS  			########
################################################
################################################


## in the order T(q_1)^+ T(q_1)^- T^d(p_1)^+ T^d(p_1)^-
function getEvenTileOffSets(str::String)
	p0,q0,p1,q1 = getFracs(str); vp1,vp2 = vf(p1); vq1,vq2 = vf(q1);
	W1 = count_full(str, '1')%2;
	if(W1 == 0)
		return 0, 2vp1 +2vp2, vq1, vq2; 

	else
		return 2vp1, 2vp2, 0, vq2+vq1;
	end
end

function getEvenTileAffine(str::String)
	p0,q0,p1,q1 = getFracs(str); ap1,ap2 = af(p1); aq1,aq2 = af(q1);
	W1 = count_full(str, '1')%2;
	if(W1 == 0)
		return 0, ap1 +2ap2, aq1, aq2; 

	else
		return 2ap1, 2ap2, 0, aq2+aq1;
	end
end

function getEvenTileFull(str::String)
	p0,q0,p1,q1 = getFracs(str); vp1,vp2 = vf(p1); vq1,vq2 = vf(q1); vp1 = 2vp1; 
	ap1,ap2 = af(p1); aq1,aq2 = af(q1); ap1 = 2ap1; 
	
	v1,v2,v3,v4 = getEvenTileOffSets(str)
	a1,a2,a3,a4 = getEvenTileAffine(str)

	v1,v2,v3,v4,v5,v6 = v1, v2, v3, v3+vp2, v4, v4+vp2
	a1,a2,a3,a4,a5,a6 = a1, a2, a3, a3+ap2, a4, a4+ap2

	return (v1,v2,v3,v4,v5,v6),(a1,a2,a3,a4,a5,a6)
end


## in the order T^d(q_1)^+ T^d(q_1)^- T(p_1)^+ T(p_1)^-
function getOddTileOffSets(str::String)
	p0,q0,p1,q1 = getFracs(str); vp1,vp2 = vf(p1); vq1,vq2 = vf(q1);
	W1 = count_full(str, '1')%2;
	if(W1 == 0)
		return 0, 2vp1+vp2, 2vq1, vq2;
	else
		return 2vp1, vp2, 0, 2vq1 + vq2;
	end
end

function getOddTileAffine(str::String)
	p0,q0,p1,q1 = getFracs(str); ap1,ap2 = af(p1); aq1,aq2 = af(q1);
	W1 = count_full(str, '1')%2;
	if(W1 == 0)
		return 0, 2ap1+ap2, 2aq1, aq2;
	else
		return 2ap1, ap2, 0, 2aq1 + aq2;
	end
end

function getOddTileFull(str::String)
	p0,q0,p1,q1 = getFracs(str); vp1,vp2 = vf(p1); vq1,vq2 = vf(q1); vp1 = 2vp1; 
	ap1,ap2 = af(p1); aq1,aq2 = af(q1); ap1 = 2ap1; 
	
	v1,v2,v3,v4 = getOddTileOffSets(str)
	a1,a2,a3,a4 = getOddTileAffine(str)

	v1,v2,v3,v4,v5,v6 = v1, v1+vq1, v2,v2+vq1, v3,v4
	a1,a2,a3,a4,a5,a6 = a1, a1+aq1, a2,a2+aq1, a3,a4

	return (v1,v2,v3,v4,v5,v6),(a1,a2,a3,a4,a5,a6)
end


function getTileOffSets(p::Rational{Int64}, str::String)
	if(isEven(p))
		return getEvenTileOffSets(str);
	else
		return getOddTileOffSets(str);
	end
end



################################################
################################################
### 			ALTERNATES  			########
################################################
################################################

## in the order T(q_1)^+,T(q_1)^-, T(p_1)^{+,1} T(p_1)^{-,1}
## last term offset between T(p_1)^{,1} and T(p_1)^{,2}
function getAltOddTileOffSets(str::String)
	p0,q0,p1,q1 = getFracs(str); vp1,vp2 = vf(p1); vq1,vq2 = vf(q1);
	W1 = count_full(str, '1')%2;
	if(W1 == 0)
		return  0,2vp1+2vp2+vq1,vq1,vq2,(vq1+vp2)
	else
		return 2vp1,2vp2-vq1,0,vq1+vq2, (-vq1+vp2);
	end
end


function getAltOddTileAffine(str::String)
	p0,q0,p1,q1 = getFracs(str); ap1,ap2 = af(p1); aq1,aq2 = af(q1);
	W1 = count_full(str, '1')%2;
	if(W1 == 0)
		return  0,2ap1+2ap2+aq1,aq1,aq2,(aq1+ap2)
	else
		return 2ap1,2ap2-aq1,0,aq1+aq2, (-aq1+ap2);
	end
end

function getAltOddTileFull(str::String)
	v1,v2,v3,v4,vs = getAltOddTileOffSets(str)
	a1,a2,a3,a4,as = getAltOddTileAffine(str)

	v1,v2,v3,v4,v5,v6 = v1,v2,v3,v3+vs,v4,v4+vs;
	a1,a2,a3,a4,a5,a6 = a1,a2,a3,a3+as,a4,a4+as

	return (v1,v2,v3,v4,v5,v6),(a1,a2,a3,a4,a5,a6)

end



## in the order T(q_1)^+,1  T(q_1)^-,1 T(p_1)^+ T(p_1)^- 
## with final term being the offset between T(q_1)^{,1} and T(q_1)^{,2}
function getAltEvenTileOffSets(str::String)
	p0,q0,p1,q1 = getFracs(str); vp1,vp2 = vf(p1); vq1,vq2 = vf(q1);
	W1 = count_full(str, '1')%2;
	if(W1 == 0)
		return 0, 2vp1+vp2, vq1,vp2-vq1+vq2, (-vq1+vp2);
	else
		return  2vp1, vp2, 0,vp2+vq2+2vq1, (vq1+vp2);

	end
end

function getAltEvenTileAffine(str::String)
	p0,q0,p1,q1 = getFracs(str); ap1,ap2 = af(p1); aq1,aq2 = af(q1);
	W1 = count_full(str, '1')%2;
	if(W1 == 0)
		return 0, 2ap1+ap2, aq1,ap2-aq1+aq2, (-aq1+ap2);
	else
		return  2ap1, ap2, 0,ap2+aq2+2aq1, (aq1+ap2);

	end
end


function getAltEvenTileFull(str::String)
	v1,v2,v3,v4,vs = getAltEvenTileOffSets(str)
	a1,a2,a3,a4,as = getAltEvenTileAffine(str)

	v1,v2,v3,v4,v5,v6 = v1,v1+vs,v2,v2+vs,v3,v4
	a1,a2,a3,a4,a5,a6 = a1,a1+as,a2,a2+as,a3,a4

	return (v1,v2,v3,v4,v5,v6),(a1,a2,a3,a4,a5,a6)

end




function getAltTileOffSets(p::Rational{Int64}, str::String)
	if(isEven(p))
		return getAltEvenTileOffSets(str);
	else
		return getAltOddTileOffSets(str);
	end
end



################################
#### get offSets to the alternate tile
################################

# maybe not correct y offset for the double overlaps

## get the offset to the corner of the alt tile in the even case
function get_even_offset(p1::Rational{Int64},q1::Rational{Int64}, str::String)
	vp1,vp2 = vf(p1)
	vq1,vq2 = vf(q1)
	ss = str[end];W1 = count_full(str, '1')%2;

	#compute the offset to the center
	if(W1 == 0)
		if(ss == '1')
			return  vp2 + (2vp1-vq1)
		elseif(ss =='2')
			return vq2
		elseif(ss == '3')
			return (vp2 + vq1-2*(vq1-2vp1)); #double overlap 
		end
	else
		if(ss == '1')
			return vp2 + vq1

		elseif(ss =='2')
			return 2vq1+vq2

		elseif(ss == '3')
			return vq1+vp2; #double overlap
		end
	end


end



## get the offset to the corner of the alt tile in the even case
function get_even_affine(p1::Rational{Int64},q1::Rational{Int64}, str::String)
	ap1,ap2 = af(p1)
	aq1,aq2 = af(q1)
	ss = str[end];W1 = count_full(str, '1')%2;

	#compute the offset to the center
	if(W1 == 0)
		if(ss == '1')
			return  ap2 + (2ap1-aq1)
		elseif(ss =='2')
			return aq2
		elseif(ss == '3')
			return (ap2 + aq1-2*(aq1-2ap1)); #double overlap 
		end
	else
		if(ss == '1')
			return ap2 + aq1

		elseif(ss =='2')
			return 2aq1+aq2

		elseif(ss == '3')
			return aq1+ap2; #double overlap
		end
	end


end



## get the offset to the corner of the alt tile in the odd case
function get_odd_offset(p1::Rational{Int64},q1::Rational{Int64}, str::String)
	vp1,vp2 = vf(p1)
	vq1,vq2 = vf(q1)
	ss = str[end];W1 = count_full(str, '1')%2;

	#compute the offset to the center
	if(W1 == 0)
		if(ss == '1')
			return  vq1+vp2
		elseif(ss =='2')
			return (vp2+vq1); #double overlap
		elseif(ss == '3')
			return (vq1+vp2); 
		end
	else
		if(ss == '1')
			return vp2 + 2vp1-vq1

		elseif(ss =='2')
			return vp2; #double overlap

		elseif(ss == '3')
			return vp2+2vp1-vq1; 
		end
	end


end



## get the offset to the corner of the alt tile in the odd case
function get_odd_affine(p1::Rational{Int64},q1::Rational{Int64}, str::String)
	ap1,ap2 = af(p1)
	aq1,aq2 = af(q1)
	ss = str[end];W1 = count_full(str, '1')%2;

	#compute the offset to the center
	if(W1 == 0)
		if(ss == '1')
			return  aq1+ap2
		elseif(ss =='2')
			return (ap2+aq1); #double overlap
		elseif(ss == '3')
			return (aq1+ap2); 
		end
	else
		if(ss == '1')
			return ap2 + 2ap1-aq1

		elseif(ss =='2')
			return ap2; #double overlap

		elseif(ss == '3')
			return ap2+2ap1-aq1; 
		end
	end


end


#requires str to be nonempty
function get_alt_offset(p::Rational{Int64}, str::String)
	p0,q0,p1,q1 = getFracs(str)
	if(isOdd(p))
		return get_odd_offset(p1,q1,str);
	else
		return get_even_offset(p1,q1,str);
	end
end


function get_alt_affine(p::Rational{Int64}, str::String)
	p0,q0,p1,q1 = getFracs(str)
	if(isOdd(p))
		return get_odd_affine(p1,q1,str);
	else
		return get_even_affine(p1,q1,str);
	end
end


function get_alt_full_offset(p::Rational{Int64}, str::String)
	return get_alt_offset(p,str), get_alt_affine(p,str)

end




## goes from lower left corner 
## of an alt tile to the correct corner
function translate_alts(p)
	alt_offset = 0; 
	if(isOdd(p))
		if(p != 0//1)
			str_p, fracs_p = rational_approximation(p, 999);
			p2,q2,p3,q3 = fracs_p
			if(p2 != p) error("prob_p"); end
			W_p = (sum([s == '1' for s in str_p]))%2
			if(W_p == 1) 
				v1,v2 = vf(p3); v1 = v1*2; 
				v3,v4 = vf(q3); 
				O = v3-v2
				if(imag(O) >= 0) 
					alt_offset = O; 
				end
			end
		end
	else
		q = p;
		if(q != 1//1)
			str_q, fracs_q = rational_approximation(q, 999);
			p2,q2,p3,q3 = fracs_q
			if(q2 != q) error("prob_q"); end
			W_q = (sum([s == '1' for s in str_q]))%2
			if(W_q == 0) 
				v1,v2 = vf(p3); v1 = v1*2; 
				v3,v4 = vf(q3); 
				O = v3-v2; 
				if(imag(O) >= 0) 
					if(all3(str_q)) O = O + im; end
					alt_offset = O; 
				end

			end
		end
	end
	return alt_offset; 
end





#######################
## L CORRECTION OFFSETS
#######################

function get_even_l_offsets(str, isFlip)
	p1,q1,p2,q2 = getFracs(str)
	vp1,vp2 = vf(p2); vq1,vq2 = vf(q2)
	ap1,ap2 = af(p2); aq1,aq2 = af(q2)
	W1 = (sum([s == '1' for s in str]))%2
	if(W1 == 1 && isFlip == false)
		K1 = 2vp1;  A1 = 2ap1;
		K2 = 2vp1; A2 = 2ap1
		if(str[end] == '2')
			K3 = (vq1-2vp1) + vq2-vp2+vq1
			A3 = (aq1-2ap1) + aq2-ap2+aq1
		else
			K3 = (vq1-2vp1); A3 = (aq1-2ap1);
		end
		K4 = vp2; A4 = (ap2); 
	elseif(W1 == 0 && isFlip == false)
		K1 = -2vp1; A1 = (-2ap1);
		K2 = -vq1; A2 = (-aq1);
		if(str[end] == '2')
			K3 =  vq2-vp2; A3 = (aq2-ap2); 
		else
			K3 = 2vp1-vq1; A3 = (2ap1-aq1);
		end
		K4 = vp2;  A4 = ap2; 
	elseif(W1 == 1 && isFlip == true)
		K1 = -2vp1; A1 = (-2ap1); 
		K2 = -vq1+vp2-vq2; #change orientation
		A2 = -aq1+ap2-aq2;
		if(str[end] == '2')
			K3 = vq2-vp2+vq1
			A3 = aq2-ap2+aq1
		else
			K3 = 0
			A3 = 0
		end
		K4 = -vp2; A4 = -ap2; 
	elseif(W1 == 0 && isFlip == true)
		K1 = 2vp1;  A1 = 2ap1;
		K2 = 2vp1+vp2-vq2; #change orientation
		A2 = 2ap1+ap2-aq2;
		if(str[end] == '2')
			K3 = vq2-vp2-2vp1+vq1
			A3 = aq2-ap2-2ap1+aq1
		else
			K3 = 0
			A3 = 0
		end
		K4 = -vp2; A4 = -ap2 
	end
	return K1,K2,K3,K4,A1,A2,A3,A4
end


function get_odd_l_offsets(str, isFlip)
	p1,q1,p2,q2 = getFracs(str)
	vp1,vp2 = vf(p2); vq1,vq2 = vf(q2)
	ap1,ap2 = af(p2); aq1,aq2 = af(q2)
	W1 = (sum([s == '1' for s in str]))%2

	if(W1 == 1 && isFlip == false)
		K1 = vq2; A1 = aq2;
		K2 = vq2+(vq1-2vp1); A2 = aq2 + (aq1-2ap1);
		K3 = (-vq2+vp2)+(2vp1-vq1); A3 = (-aq2+ap2)+(2ap1-aq1);
		K4 = -vq1; A4 = -aq1; 
	elseif(W1 == 0 && isFlip == false)
		K1 = vq2; A1 = aq2;
		K2 = vq2; A2 = aq2;
		K3 = (-vq2+vp2); A3 = (-aq2+ap2);
		K4 = vq1; A4 = aq1;
	elseif(W1 == 1 && isFlip == true)
		K1 = -vq2; A1 = -aq2;
		K2 = -vp2; A2 = -ap2; 
		K3 = 2vp1-vq1; A3 = 2ap1-aq1;
		K4 = vq1; A4 = aq1;
	elseif(W1 == 0 && isFlip == true)
		K1 = -vq2; A1 = -aq2;
		K2 = -vp2-(2vp1-vq1); A2 = -ap2-(2ap1-aq1);
		K3 = 0; A3 = 0;
		K4 = -vq1; A4 = -aq1;
	end
	return K1,K2,K3,K4,A1,A2,A3,A4

end

