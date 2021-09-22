################################################
## functions which compute the rational recursion 
################################################

## get rational approximation of x 
## using the modified Farey recursion
function rational_approximation(x, depth)
	str = ""
	p1,q1,p2,q2 = 1//2, 1//3, 0//1, 1//1;
	for i in 1:depth
		if(isapprox(x,q1)|| isapprox(x,p1)) break; end
		if(p1 <= x <= q1 || p1 >= x >= q1)
			z = 1; 
			p2,q2 = p1,q1;
		elseif(p1 <= x <= q2 || p1 >= x >= q2)
			z = 2; 
			p2,q2 = p1,q2;
		elseif(p2 <= x <= q1 || p2 >= x >= q1)
			z = 3;
			p2,q2 = p2,q1;
		end
		p1,q1 = getChildren(p2,q2)
		str = "$(str)$(z)"
	end
	fracs = p1,q1,p2,q2
	return str, fracs
end


function getFracs(str)
	p1,q1,p2,q2 = 1//2, 1//3, 0//1, 1//1;
	for i in 1:length(str)
		ss = str[i]
		p1,q1,p2,q2 = getNextFrac((p1,q1,p2,q2),ss)
	end
	return p1,q1,p2,q2
end


function getNextFrac(fracs,ss)
	p1,q1,p2,q2 = fracs
	if(ss == '1')
		p2,q2 = p1,q1; 
	elseif(ss == '2')
		p2,q2 = p1,q2;
	elseif(ss == '3')
		p2,q2 = p2,q1;
	else
		error("something went wrong")
	end
	return getChildren(p2,q2)..., p2,q2
end 




function getModFareyParents(p::Rational{Int64})
	str, fracs = rational_approximation(p, 999)
	p1,q1,p2,q2 = fracs

	if(isOdd(p)) 
		if(p1 != p) error("not found"); end 
	else
		if(q1 != p) error("not found"); end
	end

	return p2,q2
end

