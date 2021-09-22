using Plots 



#############################################
#ignores NaN
#############################################
function nan_ign(x)
	if( isnan(x))
		return Inf
	else 
		return x
	end
end
#############################################
## minimum of two numbers sending NaN to Inf
#############################################
function min_nan(x,y)
	if(isnan(x))
		return y
	elseif(isnan(y))
		return x
	else
		return min(x,y)
	end
end


#############################################
## compute the discrete Laplacian with 0 boundary conditions
#############################################
function laplacian(G, offSet = 0)
	A = fill(NaN, size(G))
	for x in 2:(size(G,1)-1), y in 2:(size(G,2)-1)
		if((x+y)%2 == offSet)
			A[x,y] = -2*G[x,y] + G[x+1,y] + G[x-1,y];
		else
			A[x,y] = -2*G[x,y] + G[x,y+1] + G[x,y-1];
		end
	end
	return A
end


#############################################
## quintuple the size of a partial odometer by 
## padding it with NaNs
#############################################

function pad(A...)
	N1,N2 = 0,0;
	for A1 in A
		N1, N2 = max(N1, size(A1,1)), max(N2, size(A1,2));
	end
	N1,N2 = ceil(Int, 4*N1), ceil(Int, 4*N2);
	M1 = ceil(Int, N1*1/4); M2 = ceil(Int, N2*1/4);

	G_ = map(A1->fill(NaN,N1,N2),A)
	for i in 1:length(G_)
		G1 = A[i]; 
		G_[i][M1:M1+size(G1,1)-1, M2:M2+size(G1,2)-1] .= G1;
	end
	return G_
end

#############################################
#remove nanpadding from G
#############################################
function unpad(A...)
	G = fill(0, size(A[1]));
	for A1 in A
		G = G + isfinite.(A1);
	end
	xs = findall(vec(sum(G,dims=2)).>0)
	x1 = xs[1]-1; x2 = xs[end]+1;
	ys = findall(vec(sum(G,dims=1)).>0)
	y1 = ys[1]-1; y2 = ys[end]+1;
	

	x1,y1 = max(x1, 1), max(y1,1);
	x2,y2 = min(x2, size(G,1)), min(y2, size(G,2));

	return map(A1->A1[x1:x2,y1:y2],A)
end

###########################
### centers two matrices so that
### the lower left vertices match
###########################
function center_mats(Bp,Bq)
	Bp,Bq = pad(Bp,Bq);
	x1,y1 = get_lower_left_vertex(Bp)
	x2,y2 = get_lower_left_vertex(Bq)
	O = (x2-x1) + im*(y2-y1);
	Bp = translate_va(Bp, O, 0)
	return unpad(Bp,Bq)
end

###########################
### centers arbitrary number of matrices
### so that the lower left vertices match
###########################
function center_mats(A...)
	A = pad(A...)
	x1,y1 = get_lower_left_vertex(A[1])
	O = fill(0im, length(A))
	for i in 1:length(A)
		x2,y2 = get_lower_left_vertex(A[i])
		O[i] = (x1-x2) + im*(y1-y2);
	end
	f(A,O) = translate_va(A,O,0);
	A = map(f, A,O);
	return unpad(A...)
end




#############################################
## translates a partial odometer by vector v
#############################################

function translate_v(G,v)
	## translate by v 
	G_ = fill(NaN, size(G))
	for x in 1:size(G,1), y in 1:size(G,2)
		xx = x + real(v); yy = y + imag(v);
		if(1 <= xx <= size(G,1) && 1 <= yy <= size(G,2))
			G_[xx,yy] = G[x,y]
		end
	end	
	return G_
end

#############################################
#Translate a partial odometer by vector v and slope a
#############################################

function translate_va(G,v,a)
  G = translate_v(G, v)
  Y1 = [x for x in 0:size(G,1)-1, y in 1:size(G,2)]
  Y2 = [y for x in 1:size(G,1), y in 0:size(G,2)-1]

  Z = Y1 + im .* Y2
  G = G + real(a'*Z)
end




#############################################
## AFFINE FUNCTIONS
#############################################

t(n,d) = (d^2 + 2*d*n - n^2)
m(n,d) = 1//(t(n,d))*[-n^2 d*n; d*n -d^2]
ve(n,d) = [d; n], [n-d; n+d]
vf(n,d) = (d +  n*im), ((n-d) + (n+d)*im)
vf(p) = vf(p.num, p.den);
af(p) = af(p.num, p.den)
function af(n,d)
	v1,v2 = ve(n,d);
	mat = m(n,d)
	a1,a2 = mat*v1, mat*v2;
	return a1[1] + im*a1[2], a2[1] + im*a2[2];
end


function getChildren(o, e)
	c1 = (e.num + o.num)//(e.den+o.den)
	c2 = (e.num+2*o.num)//(e.den+2*o.den)
	return c1,c2
end

function getParents(o,e)
	if(o == 0//1)
		error("too_small")
	end
	e2 = (2*o.num-e.num)//(2*o.den-e.den)
	o2 = (e.num-o.num)//(e.den-o.den)
	return o2,e2
end

function switch(p,q)
	pp = p; qq = q; 
	if((pp.num+pp.den)%2 == 0)
		p = qq;
		q = pp;
	end
	return p,q
end


#is p an ancestor of q
function isAncestor(p,q)
	if(q == 0//1 || q == 1//1) return false; end
	if(q == p) return true; end

	q1,q2 = getParents(q)
	bool1 = isAncestor(p,q1)
	bool2 = isAncestor(p, q2)
	return (bool1 || bool2)
end

####################################
### get the tiling of a tile odometer
####################################
function get_tiling(A, p)
		v1,v2 = vf(p); a1,a2 = af(p)
		if(isOdd(p)) v1 = 2v1; end
		C0, B = pad(A,A);
		if(size(A,1) < 50) C0,B = pad(C0,B); end 
		for (i,j) in [(0,0), (-1,0), (1,0), (0,1), (0,-1), (1,1), (1,-1), (-1,1), (-1,-1) ]
			v = v1*i + v2*j; a = a1*i + a2*j; 
			C1 = translate_va(B, v,a);
			z =  minimum(nan_ign, (C1 - C0));
			C1 = C1 .-z
			C0 = min_nan.(C0, C1); 
		end
		return C0,B;
end


#############################################
## draws the pattern associated with the odometer
## by drawing the tiling and then truncating
#############################################
function drawPattern(G; II= 0, p = false, title = "", isTile=false)
	if(p != false)
		C0,B = get_tiling(G,p)
		


		#center at the lower left vertex
		x,y = get_lower_left_vertex(C0);
		II = (II+x+y)%2; 

		S = laplacian(C0,II); if(length(unique(S)) > 3) II=(II+1)%2; end


		#two way of doing it
		S = laplacian(C0,II); if(!isTile) S[isnan.(B)].=NaN; end; S = unpad(S)[1];

	else

		x,y = get_lower_left_vertex(G);
		II = (II+x+y)%2; 
		
		S = laplacian(G,II); if(length(unique(S)) > 3) II=(II+1)%2; end; 
		S = laplacian(G,II);
		S[isnan.(S).*isfinite.(G)] .= 0; 
	end
	if(title == "" && p != false) title = p; end
 	return heatmap(S', title = title, legend = false, aspect_ratio = true, axis = false, ticks = false, c = :grays)
end



####################################
### FUNCTIONS FOR FAREY CONSTRUCTION
####################################

function getFarey(n1,d1,n2,d2)
	p = (n1+n2)//(d1+d2);
	n = p.num;
	d = p.den 
	return n,d
end

function findFareyParents(p)
	nU,dU,nL,dL = 1,1,0,1
	if(p <= 0 || p >= 1) return -1,-1,-1,-1; end;
	while(true)
		println( (nL, dL), " ", (nU, dU))
		n,d = getFarey(nU,dU,nL,dL)
		if(p < n//d)
			nU = n; dU = d
		elseif(p > n//d)
			nL = n; dL = d; 
		elseif(p == n//d)
			println(p)
			println("Done!")
			break;
		end
	end
end

function getFareyParents(n,d)
	p = n//d
	if(p < 0 || p > 1) return -1,-1,-1,-1; end;
	nU,dU,nL,dL = 1,1,0,1
	while(true)
		n,d = getFarey(nU,dU,nL,dL)
		if(p < n//d)
			nU = n; dU = d
		elseif(p > n//d)
			nL = n; dL = d; 
		elseif(p == n//d)
			return nL,dL, nU,dU
		end
	end
end



function getFlippedRep(p)
	n = p.num; d = p.den; 
	return (d-n)//(n+d)
end

function printFlippedRep(p)
	n = p.num; d = p.den; 

	println((n//d), " ", (d-n)//(n+d))

end

## relic of old code 
function init_mats()
	mats = G_1_2,G_1_3, G_0_1,G_1_1;
	alt_mats = G_1_2_alt, G_1_3_alt, G_0_1_alt, G_1_1_alt;
	fracs = 1//2,1//3,0//1,1//1;
	rot_mats = G_1_2_r, G_1_3_r, G_0_1_r, G_1_1_r


	## save offsets for each pair 
	p1_q1 = 0,0,0,2;
	p1_q2 = -2,-2,2,-2;
	p2_q1 = -im,-im,im,2-im; 
	offSets = p1_q1, p1_q2, p2_q1
	reps = 1,0,0;
	tsum = 1,0,0,1; 
	return mats, alt_mats, rot_mats, fracs, offSets,reps,tsum
end


function isOdd(p)
	if((p.num+p.den)%2 == 1) return true; end
	return false; 
end
function isEven(p) return !isOdd(p); end





####################################
### FUNCTIONS FOR WORDS
####################################

### word helper functions
function count_full(str, c)
	return sum([s==c for s in str])
end

function get_first_ancestor(str)
	z = findfirst('1',str)
	if(z == nothing)

	else
		str = str[1:z-1]
	end

	z2 = findfirst('2',str)
	z3 = findfirst('3',str)
	if(z2 == nothing || z3 == nothing)
		return str
	else
		z = max(z2,z3)
		return str[1:z-1]
	end
end


function involution(w)
	return .- reverse(w)
end

function count_end(str,ss)
	str = reverse(str)
	k = 0
	for i in 1:length(str)
		if(str[i] == ss)
			k+=1
		else
			break;
		end
	end
	return k 
end



################################################
## assuming a background filled with NaN
## and a finite figure, gets a vertex
################################################


function get_lower_right_vertex(D)
	for x in reverse(1:size(D,1))
		for y in 1:size(D,2)
			if(isfinite(D[x,y]))
				return x,y
			end
		end
	end
end

function get_upper_left_vertex(D)
	for x in 1:size(D,1)
		for y in reverse(1:size(D,2))
			if(isfinite(D[x,y]))
				return x,y
			end
		end
	end
end

function get_lower_left_vertex(D)
	for y in 1:size(D,2)
		for x in 1:size(D,1)
			if(isfinite(D[x,y]))
				return x,y
			end
		end
	end
end


function get_upper_right_vertex(D)
	for y in reverse(1:size(D,2))
		for x in reverse(1:size(D,1))
			if(isfinite(D[x,y]))
				return x,y
			end
		end
	end
end


