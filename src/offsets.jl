# Offsets used to place the parent tiles when building a child tile.
_offset_zero() = Complex{BigInt}(0, 0)
_exact_fraction(p::Rational) = BigInt(numerator(p)) // BigInt(denominator(p))
_exact_fraction(p::Integer) = BigInt(p) // BigInt(1)
count_full(s::AbstractString, c::Char) = count(==(c), s)
all2(s::AbstractString) = all(==('2'), s)
all3(s::AbstractString) = all(==('3'), s)
isOdd(p::Rational) = isodd(BigInt(numerator(p)) + BigInt(denominator(p)))
isEven(p::Rational) = !isOdd(p)

function vf(n::Integer, d::Integer)
    N,D = BigInt(n),BigInt(d)
    return Complex{BigInt}(D,N),Complex{BigInt}(N-D,N+D)
end
vf(p::Rational) = vf(numerator(p),denominator(p))
# m(n,d)*v1 = 0 and m(n,d)*v2 = (n,-d)
function af(n::Integer,d::Integer)
    return _offset_zero(),Complex{BigInt}(BigInt(n),-BigInt(d))
end
af(p::Rational) = af(numerator(p),denominator(p))
function getChildren(o::Rational,e::Rational)
    O,E = _exact_fraction(o),_exact_fraction(e)
    return (numerator(E)+numerator(O))//(denominator(E)+denominator(O)),
           (numerator(E)+2*numerator(O))//(denominator(E)+2*denominator(O))
end
function getNextFrac(fracs, c::Char)
    p,q,a,b = map(_exact_fraction,fracs)
    if c=='1'; a,b=p,q
    elseif c=='2'; a,b=p,b
    elseif c=='3'; a,b=a,q
    else; throw(ArgumentError("Farey word must contain only 1, 2, 3"))
    end
    pp,qq=getChildren(a,b)
    return pp,qq,a,b
end
function getFracs(s::AbstractString)
    f=(big(1)//big(2),big(1)//big(3),big(0)//big(1),big(1)//big(1))
    for c in s; f=getNextFrac(f,c); end
    return f
end
_between(x,a,b) = min(a,b)<=x<=max(a,b)
"""Return the word at which x first appears, and the fractions for that word.
x must be strictly between 0 and 1."""
function rational_approximation(x::Rational)
    X=_exact_fraction(x)
    0<X<1 || throw(DomainError(X,"x must be strictly between 0 and 1"))
    f=getFracs(""); word=IOBuffer()
    while X!=f[1] && X!=f[2]
        p,q,a,b=f
        c = _between(X,p,q) ? '1' : _between(X,p,b) ? '2' : '3'
        c=='3' && !_between(X,a,q) && error("invalid Farey interval")
        write(word,c); f=getNextFrac(f,c)
    end
    return String(take!(word)),f
end
function getModFareyParents(p::Rational)
    _,f=rational_approximation(p)
    return f[3],f[4]
end
function getParents(o::Rational,e::Rational)
    O,E=_exact_fraction(o),_exact_fraction(e)
    O==0 && throw(DomainError(O,"0 has no parents"))
    return (numerator(E)-numerator(O))//(denominator(E)-denominator(O)),
           (2*numerator(O)-numerator(E))//(2*denominator(O)-denominator(E))
end
# Rotations by 90 degrees clockwise and counterclockwise.
Rot(z::Complex{<:Integer}) = Complex{BigInt}(BigInt(imag(z)),-BigInt(real(z)))
rotccw(z::Complex{<:Integer}) = Complex{BigInt}(-BigInt(imag(z)),BigInt(real(z)))
arithmetic_dual(p::Rational) = (big(1)-_exact_fraction(p))/(big(1)+_exact_fraction(p))
function getAffine(s::AbstractString)
    _,_,p,q=getFracs(s)
    p1,p2=vf(p);q1,q2=vf(q);a1,a2=af(p);b1,b2=af(q)
    return 2*p1,p2,q1,q2,2*a1,a2,b1,b2
end

## helper functions to compute standard tile offsets 


################################################
################################################
### 			STANDARDS  			########
################################################
################################################


## in the order T(q_1)^+ T(q_1)^- T^d(p_1)^+ T^d(p_1)^-
function getEvenTileOffSets(str::AbstractString)
	p0,q0,p1,q1 = getFracs(str); vp1,vp2 = vf(p1); vq1,vq2 = vf(q1);
	W1 = count_full(str, '1')%2;
	if(W1 == 0)
		return _offset_zero(), 2vp1 +2vp2, vq1, vq2; 

	else
		return 2vp1, 2vp2, _offset_zero(), vq2+vq1;
	end
end

function getEvenTileAffine(str::AbstractString)
	p0,q0,p1,q1 = getFracs(str); ap1,ap2 = af(p1); aq1,aq2 = af(q1);
	W1 = count_full(str, '1')%2;
	if(W1 == 0)
		return _offset_zero(), ap1 +2ap2, aq1, aq2; 

	else
		return 2ap1, 2ap2, _offset_zero(), aq2+aq1;
	end
end

function getEvenTileFull(str::AbstractString)
	p0,q0,p1,q1 = getFracs(str); vp1,vp2 = vf(p1); vq1,vq2 = vf(q1); vp1 = 2vp1; 
	ap1,ap2 = af(p1); aq1,aq2 = af(q1); ap1 = 2ap1; 
	
	v1,v2,v3,v4 = getEvenTileOffSets(str)
	a1,a2,a3,a4 = getEvenTileAffine(str)

	v1,v2,v3,v4,v5,v6 = v1, v2, v3, v3+vp2, v4, v4+vp2
	a1,a2,a3,a4,a5,a6 = a1, a2, a3, a3+ap2, a4, a4+ap2

	return (v1,v2,v3,v4,v5,v6),(a1,a2,a3,a4,a5,a6)
end


## in the order T^d(q_1)^+ T^d(q_1)^- T(p_1)^+ T(p_1)^-
function getOddTileOffSets(str::AbstractString)
	p0,q0,p1,q1 = getFracs(str); vp1,vp2 = vf(p1); vq1,vq2 = vf(q1);
	W1 = count_full(str, '1')%2;
	if(W1 == 0)
		return _offset_zero(), 2vp1+vp2, 2vq1, vq2;
	else
		return 2vp1, vp2, _offset_zero(), 2vq1 + vq2;
	end
end

function getOddTileAffine(str::AbstractString)
	p0,q0,p1,q1 = getFracs(str); ap1,ap2 = af(p1); aq1,aq2 = af(q1);
	W1 = count_full(str, '1')%2;
	if(W1 == 0)
		return _offset_zero(), 2ap1+ap2, 2aq1, aq2;
	else
		return 2ap1, ap2, _offset_zero(), 2aq1 + aq2;
	end
end

function getOddTileFull(str::AbstractString)
	p0,q0,p1,q1 = getFracs(str); vp1,vp2 = vf(p1); vq1,vq2 = vf(q1); vp1 = 2vp1; 
	ap1,ap2 = af(p1); aq1,aq2 = af(q1); ap1 = 2ap1; 
	
	v1,v2,v3,v4 = getOddTileOffSets(str)
	a1,a2,a3,a4 = getOddTileAffine(str)

	v1,v2,v3,v4,v5,v6 = v1, v1+vq1, v2,v2+vq1, v3,v4
	a1,a2,a3,a4,a5,a6 = a1, a1+aq1, a2,a2+aq1, a3,a4

	return (v1,v2,v3,v4,v5,v6),(a1,a2,a3,a4,a5,a6)
end


function getTileOffSets(p::Rational{BigInt}, str::AbstractString)
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
function getAltOddTileOffSets(str::AbstractString)
	p0,q0,p1,q1 = getFracs(str); vp1,vp2 = vf(p1); vq1,vq2 = vf(q1);
	W1 = count_full(str, '1')%2;
	if(W1 == 0)
		return _offset_zero(),2vp1+2vp2+vq1,vq1,vq2,(vq1+vp2)
	else
		return 2vp1,2vp2-vq1,_offset_zero(),vq1+vq2, (-vq1+vp2);
	end
end


function getAltOddTileAffine(str::AbstractString)
	p0,q0,p1,q1 = getFracs(str); ap1,ap2 = af(p1); aq1,aq2 = af(q1);
	W1 = count_full(str, '1')%2;
	if(W1 == 0)
		return _offset_zero(),2ap1+2ap2+aq1,aq1,aq2,(aq1+ap2)
	else
		return 2ap1,2ap2-aq1,_offset_zero(),aq1+aq2, (-aq1+ap2);
	end
end

function getAltOddTileFull(str::AbstractString)
	v1,v2,v3,v4,vs = getAltOddTileOffSets(str)
	a1,a2,a3,a4,as = getAltOddTileAffine(str)

	v1,v2,v3,v4,v5,v6 = v1,v2,v3,v3+vs,v4,v4+vs;
	a1,a2,a3,a4,a5,a6 = a1,a2,a3,a3+as,a4,a4+as

	return (v1,v2,v3,v4,v5,v6),(a1,a2,a3,a4,a5,a6)

end



## in the order T(q_1)^+,1  T(q_1)^-,1 T(p_1)^+ T(p_1)^- 
## with final term being the offset between T(q_1)^{,1} and T(q_1)^{,2}
function getAltEvenTileOffSets(str::AbstractString)
	p0,q0,p1,q1 = getFracs(str); vp1,vp2 = vf(p1); vq1,vq2 = vf(q1);
	W1 = count_full(str, '1')%2;
	if(W1 == 0)
		return _offset_zero(), 2vp1+vp2, vq1,vp2-vq1+vq2, (-vq1+vp2);
	else
		return  2vp1, vp2, _offset_zero(),vp2+vq2+2vq1, (vq1+vp2);

	end
end

function getAltEvenTileAffine(str::AbstractString)
	p0,q0,p1,q1 = getFracs(str); ap1,ap2 = af(p1); aq1,aq2 = af(q1);
	W1 = count_full(str, '1')%2;
	if(W1 == 0)
		return _offset_zero(), 2ap1+ap2, aq1,ap2-aq1+aq2, (-aq1+ap2);
	else
		return  2ap1, ap2, _offset_zero(),ap2+aq2+2aq1, (aq1+ap2);

	end
end


function getAltEvenTileFull(str::AbstractString)
	v1,v2,v3,v4,vs = getAltEvenTileOffSets(str)
	a1,a2,a3,a4,as = getAltEvenTileAffine(str)

	v1,v2,v3,v4,v5,v6 = v1,v1+vs,v2,v2+vs,v3,v4
	a1,a2,a3,a4,a5,a6 = a1,a1+as,a2,a2+as,a3,a4

	return (v1,v2,v3,v4,v5,v6),(a1,a2,a3,a4,a5,a6)

end




function getAltTileOffSets(p::Rational{BigInt}, str::AbstractString)
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
function get_even_offset(p1::Rational{BigInt},q1::Rational{BigInt}, str::AbstractString)
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
function get_even_affine(p1::Rational{BigInt},q1::Rational{BigInt}, str::AbstractString)
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
function get_odd_offset(p1::Rational{BigInt},q1::Rational{BigInt}, str::AbstractString)
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
function get_odd_affine(p1::Rational{BigInt},q1::Rational{BigInt}, str::AbstractString)
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
function get_alt_offset(p::Rational{BigInt}, str::AbstractString)
	p0,q0,p1,q1 = getFracs(str)
	if(isOdd(p))
		return get_odd_offset(p1,q1,str);
	else
		return get_even_offset(p1,q1,str);
	end
end


function get_alt_affine(p::Rational{BigInt}, str::AbstractString)
	p0,q0,p1,q1 = getFracs(str)
	if(isOdd(p))
		return get_odd_affine(p1,q1,str);
	else
		return get_even_affine(p1,q1,str);
	end
end


function get_alt_full_offset(p::Rational{BigInt}, str::AbstractString)
	return get_alt_offset(p,str), get_alt_affine(p,str)

end




## goes from lower left corner 
## of an alt tile to the correct corner
function translate_alts(p)
	alt_offset = _offset_zero(); 
	if(isOdd(p))
		if(p != 0//1)
			str_p, fracs_p = rational_approximation(p);
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
			str_q, fracs_q = rational_approximation(q);
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
			K3 = _offset_zero()
			A3 = _offset_zero()
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
			K3 = _offset_zero()
			A3 = _offset_zero()
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
		K3 = _offset_zero(); A3 = _offset_zero();
		K4 = -vq1; A4 = -aq1;
	end
	return K1,K2,K3,K4,A1,A2,A3,A4

end

