## contains explicit formulas to construct the base cases

function all3(str) return length(str) == 0 || (!occursin('1', str) && !occursin('2', str)); end
function all2(str) return length(str) == 0 || (!occursin('1', str) && !occursin('3', str)); end


### ##################
### STANDARD STAIRCASES 
### ##################

## 3^k 
function get_staircase_even(str, isFull=true)
	fracs = getFracs(str);
	q = fracs[2].den;
	d = q
	A = fill(NaN, 2*d-1, d+2) 

	for i in 1:size(A,1), j in 1:size(A,2)
		x = i-(d-1); y = j-1; ## lower left vertex is 0,0 
		z1 = -1//2*(y)*(y+1);
		z2 = y+min(x-1,0)
		if(2-d <= x <= d && 0 <= y <= d+1 && 1 <= (x+y) <=  d+2)
			A[i,j] = z1+z2
		end
		if(isFull)
			if((x,y) == (0,0) || (x,y) == (2, d+1) )
				A[i,j] = z1+z2+1;
			end
		end
	end
	return A; 

end


## 2^k
function get_staircase_odd(str, isFull = true)
	fracs = getFracs(str);
	p = fracs[1];
	d = getFlippedRep(p).den;

	A = fill(NaN, d+2,2*d-1) 

	for i in 1:size(A,1), j in 1:size(A,2)
		x = i-1; y = j; #lower left vertex is 0,1
		z1 = -floor(Int, (x-y)^2//4)
		z2 = min(d-y,0)

		if(0 <= x <= d+1 && 1 <= y <= 2*d-1 && -1 <= (y-x) <=  d)
			A[i,j] =z1+z2
		end

		if(isFull)
			#fill the corners ADD 1 
			if( (x,y) == (0,d+1) || (x,y) == (d+1,d-1)) 
				A[i,j] = z1+z2+1
			end
		end
	end


	return A

end



### ##################
### ALTERNATE STAIRCASES 
### ##################

## 3^k 
function get_staircase_even_alt(str, isFull = true)
	fracs = getFracs(str);
	q = fracs[2].den;
	d = q
	A = fill(NaN, 3*d-3, d+3) 
	C = fill((NaN,NaN), size(A));

	for i in 1:size(A,1), j in 1:size(A,2)
		x = (i-1)-(2*d-3-1); #middle corner is (0,0)
		y = (j-1)-1


		z1 = -1//2*(y)*(y+1);
		z2 = min(0, 2-x)
		z3 = min(0, d-2+x)
		z4 = max((x+y-2),0);

		if(-d+2 <= (x+y) <= d+2 && 0 <= y <= d && -(2*d-4) <= x<=d)
			A[i,j]=z1+z2+z3+z4;
			C[i,j] = (x,y)
		end

		if(isFull)

			#fill bottom 
			if( (y == -1 && -(d-1) <= x <= 0) || (y == d+1 &&  3-(d-2)-1 <= x <= 3)) 
				A[i,j]=z1+z2+z3+z4;
				C[i,j] = (x,y)

			end

			#correct diagonal 
			if( (x+y)==2 && 0 < y < d)
				A[i,j]=z1+z2+z3+z4+1;
			end

			#fill corners
			if((x,y) == (-(d-1),0) || (x,y) == (3,d))
				A[i,j]=z1+z2+z3+z4+1
				C[i,j] = (x,y)

			end
			#fill extra two corner
			if((x,y) == (-(d-1),-1) || (x,y) == (3,d+1))
				A[i,j]=z1+z2+z3+z4+1
				C[i,j] = (x,y)

			end

		end
	end
	return A;
end


## 2^k
function get_staircase_odd_alt(str, isFull = true)
	fracs = getFracs(str);
	p = fracs[1];
	d = getFlippedRep(p).den;
	A = fill(NaN, d+3, 3*d-3);  #right tile
	C = fill((NaN,NaN), size(A));
	D = fill(NaN, size(A)); 
	for i in 1:size(A,1), j in 1:size(A,2)
		x = i-1-1; y = j;  #lower left corner is 0,1

		z1 = -floor(Int, (x-y)^2//4)
		z2 = min(d-y-1,0)
		z3 = min((2*d-1)-y,0)
		z4 = max((y-x)-d+1,0)

		if(0 <= x <= d && 1 <= y <= 3*d-3 && -1 <= (y-x) <= 2*d-1)
			A[i,j] = z1+z2+z3+z4
			C[i,j] = (x,y)

		end

		if(y-x == d -1 && 0 < x < d)
			A[i,j] = z1+z2+z3+z4+1; #correct diagonal
		end


		if(isFull)


			#fill bottom 
			if( (x == -1 && d+1 <= y <= 2*d) || (x == d+1 &&  d-2 <= y <= 2*d-3))
				A[i,j]=z1+z2+z3+z4;
				C[i,j] = (x,y)

			end

			#fill corners
			if( (x,y) == (0,2*d) || (x,y) == (d,d-2))
				A[i,j] =  z1+z2+z3+z4+1
				C[i,j] = (x,y)

			end


			#fill extra two corners
			if((x,y) == (-1,2*d) || (x,y) == (d+1,d-2))
				A[i,j]=z1+z2+z3+z4+1
				C[i,j] = (x,y)

			end


		
		end


		D[i,j] = z4

	end
	return A;
end



### ##################
### STANDARD DOUBLED STAIRCASES 
### ##################



function get_doubled_ladder_odd(str)
	p = getFracs(str)[1]
	d = p.den
	A = fill(NaN, 3*d-1, d+3)
	C = [(NaN, NaN) for x in 1:size(A,1), y in 1:size(A,2)]
	D = fill(NaN, size(A));
	for i in 1:size(A,1), j in 1:size(A,2)
		x = i-d+1# lower left corner is (0,0)
		y = j -1; 
		D[i,j] = x+y

		z1 = -y*(y+1)//2
		z2 = y+min(0, x-1)
		z3 = max((y+x)-d-2, 0)
		z4 = min(d+1-x,0)
		if(1 <= (x+y) <= 2*d+3 && -d+2 <= x <= 2*d && 0 <= y <= d+2)
			if( (y == 0 && x >= d) || (y == d+2 && x <= 2)) #remove the strips

			else
				A[i,j] = z1+z2+z3+z4
				C[i,j] = (x,y)
			end
		end

		## add the corners
		if( (x,y) == (0,0) || (x,y) == (d+2,d+2))
			A[i,j] = z1+z2+z3+z4+1;
			C[i,j] = (x,y)
		end

		#increment diagonal 
		if( (x+y)==d+2 && 1 < y < d+1)
			A[i,j] = z1+z2+z3+z4+1
			C[i,j] = (x,y)
		end


	end
	return A
end




## 2^k
function get_doubled_ladder_even(str)
	p = getFlippedRep(getFracs(str)[2])
	d = p.den
	A = fill(NaN, d+3, 3*d-1)
	C = [(NaN, NaN) for x in 1:size(A,1), y in 1:size(A,2)]
	D = fill(NaN, size(A));

	for i in 1:size(A,1), j in 1:size(A,2)
		x = i-2 #lower left corner is (0,1)
		y = j; 
		
		z1 = -floor(Int, (x-y)^2//4)
		z2 = min(0, d-y)
		z3 = max(0, y-x -d); 
		z4 = min(0, 2*d-y)


		if(-1 <= x <= d+1 && 1 <= y <= 3*d-1 && 2*d+1 >= y-x >= -1)
			if( (x == -1 && y <= d+1)  || (x == (d+1) && y >= 2*d-1))
				#remove strips
			else
				A[i,j] = z1+z2+z3+z4
			end
		end  

		## add corners 
		if( (x,y) == (-1, 2*d+1) || (x,y) == (d+1, d-1))
			A[i,j]=z1+z2+z3+z4+1;
		end

		#fix diagonal 
		if( (y-x) == d && 0 < x < d)
			A[i,j]=z1+z2+z3+z4+1;
		end


		D[i,j] = x-y
		C[i,j] = (x,y)

	end
	return A;
 
end

### ##################
### ALTERNATE DOUBLED STAIRCASES 
### ##################


## 3^k 
function get_doubled_ladder_odd_alt(str)
	p = getFracs(str)[1]
	d = p.den
	B = fill(NaN, d*3-2, d+4);
	C = fill((NaN,NaN), size(B)); 
	D = fill(NaN, size(B))
	for i in 1:size(B,1), j in 1:size(B,2)
		x = i - d; y = j-2; #bottom left corner is -1,-1
		C[i,j] = (x,y)

		z1 = -y*(y+1)//2;
		z2 = min(x,0)
		z3 = max((x+y-d), 0)
		z4 = min(d-x-1,0)


		if(-1 <= (x+y) <= 2*d+1 &&  -d+1 <= x <= 2*d-2 && -1 <= y <= d+2)

			if( (-1 <= y <= 0 && x >= d) || d+1 <= y <= d+2 && x <= -1)
				#avodi 

			else
				B[i,j] = z1+z2+z3+z4
				D[i,j]=x+y

				if((x+y)==d && 0 < x < d-1)
					B[i,j]=z1+z2+z3+1
				end
			end
		end

		#fill corners 
		if( (x,y) == (-1,-1) || (x,y) == (d,d+2))
			B[i,j] = z1 + z2 + z3 + z4 + 1;
		end

	end
	return B; 

end

## 2^k
function get_doubled_ladder_even_alt(str, isOffset=false)
	fracs = getFracs(str)

	p = getFlippedRep(fracs[2]); 
	d = p.den; 


	B = fill(NaN, d+4, d*3-2);
	C = fill((NaN,NaN), size(B)); 
	D = fill(NaN, size(B))
	for i in 1:size(B,1), j in 1:size(B,2)
		x = i-3; y = j; #bottom left corner is 0,1
		C[i,j] = (x,y)

		z1 = -floor(Int, (x-y)^2//4)
		z2 = min(0, d-y)
		z4 = min(2d-1-y,0)
		z3 = max(y-x-d,0);

		if( 2d+1 >= (y-x)  >= -1 && -2 <= x <= d+1 && 1 <= y <= 3d-2)
			if(  (x <= -1 && y <=d-1 ) || ( x >= d && y >= 2d ) )

				#ignore strips
			else
				B[i,j] = z1+z2+z3+z4

				if((y-x)==d && 0 < x < d-1)
					B[i,j]=z1+z2+z3+z4+1
				end
			end
		end

		# add corners
		if((x,y) == (-2, 2d) || (x,y) ==  (d+1, d-1))
			B[i,j]=z1+z2+z3+z4+1
		end

		D[i,j]=y

	end
	return B;
end


#### #### #### #### #### #### #### #### #### #### #### 
#### GET THE BASE QUADRUPLE
#### BASE QUADRUPLE CODE
#### #### #### #### #### #### #### #### #### #### #### 


#given a nan matrix, creates the 0/1 tile
function get_0_1(O)
	G_0_1 = fill(NaN, size(O))
	x,y = get_lower_left_vertex(O)
	G_0_1[x:x+1,y:y+2] .= 
	[
	   0.0    0.0    -1
	   0.0    0.0    -1
	]
	return G_0_1;
end
#given nan, creates the 1/1 tile
function get_1_1(O)
	G_1_1 = fill(NaN, size(O))
	x,y = get_lower_left_vertex(O)
	G_1_1[x:x+1,y:y+1] .= 
	[
	   0.0    -1    
	   0.0    0.0   
	]
	return G_1_1;
end

## note that this gives the smaller 1/2 and 1/3
## but that's OK
function get_base_quadruple(str)
	if(length(str) == 0)
		#standards
		Op1 = get_staircase_odd(str)
		Oq1 = get_staircase_even(str)
		Op0 = get_0_1(Op1)
		Oq0 = get_1_1(Oq1)

		#alternates 
		Op1_ = get_staircase_odd_alt(str)
		Oq1_ = get_staircase_even_alt(str)
		Op0_ = get_0_1(Op1)
		Oq0_ = get_1_1(Oq1)


	elseif(all2(str)) #all 2
		#standards
		Op1 = get_staircase_odd(str)
		Oq1 = get_doubled_ladder_even(str)
		Op0 = get_staircase_odd(str[1:end-1])
		Oq0 = get_1_1(Oq1)


		#alternates
		Op1_ = get_staircase_odd_alt(str)
		Oq1_ = get_doubled_ladder_even_alt(str)
		Op0_ = get_staircase_odd_alt(str[1:end-1])
		Oq0_ = get_1_1(Oq1)


	elseif(all3(str)) #all 3
		#standards
		Op1 = get_doubled_ladder_odd(str)
		Oq1 = get_staircase_even(str)
		Op0 = get_0_1(Op1)
		Oq0 = get_staircase_even(str[1:end-1])

		#alternates
		Op1_ = get_doubled_ladder_odd_alt(str)
		Oq1_ = get_staircase_even_alt(str)
		Op0_ = get_0_1(Op1)
		Oq0_ = get_staircase_even_alt(str[1:end-1])


	
	else
		error("not a base case")
	end
	mats = Op1,Oq1,Op0,Oq0
	alt_mats = Op1_,Oq1_,Op0_,Oq0_

	return mats, alt_mats
end
