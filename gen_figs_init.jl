## explore the tree up to a given depth recursively 
## and produce pics

global dir_str = "test_figs"



function save_curr_fig(mats, alt_mats, str::String)
	Ap, Aq = mats[1],mats[2]
	Bp,Bq = alt_mats[1], alt_mats[2]
	p,q,p1,q1 = getFracs(str);


	#comment in if you want to save the full plots
	isTile = false; isDecomp = true;
	plt1 = full_plot_std(Ap, p, isTile = isTile, isDecomp = isDecomp); plot!(title = "std -$(p1)"); 
	plt2 = full_plot_std(Aq, q, isTile = isTile, isDecomp = isDecomp); plot!(title = "std -$(q1)");
	plt3 = full_plot_alt(Bp, p, isTile = isTile, isDecomp = isDecomp); plot!(title = "alt -$(p1)");
	plt4 = full_plot_alt(Bq, q, isTile = isTile, isDecomp = isDecomp); plot!(title = "alt -$(q1)");
	#make sure the full directory exists
	savefig(plot(dpi=1024, plt1,plt2,plt3,plt4), "$(dir_str)/full_plots/_$(str)") 
	return; 

	A_c,A_cf = get_fig_tiling(Ap, p)
	B_c,B_cf = get_fig_tiling(Bp, p)
	save("$(dir_str)/odd/$(str)_.png", rotl90(A_c,3))
	save("$(dir_str)/tiling/odd/$(str)_.png", rotl90(A_cf,3))


	save("$(dir_str)/odd/$(str)_alt.png", rotl90(B_c,3))
	save("$(dir_str)/tiling/odd/$(str)_alt.png", rotl90(B_cf,3))


	A_c,A_cf = get_fig_tiling(Aq, q)
	B_c,B_cf = get_fig_tiling(Bq, q)
	save("$(dir_str)/even/$(str)_.png", rotl90(A_c,3))
	save("$(dir_str)/tiling/even/$(str)_.png", rotl90(A_cf,3))

	save("$(dir_str)/even/$(str)_alt.png", rotl90(B_c,3))
	save("$(dir_str)/tiling/even/$(str)_alt.png", rotl90(B_cf,3))


end

function get_fig_tiling(A, p)
	v1,v2 = vf(p); a1,a2 = af(p)
	if(isOdd(p)) v1 = 2v1; end
	C0, B = pad(A,A);
	if(size(A,1) < 50) C0,B = pad(C0,B); end
	for (i,j) in [(0,0), (-1,0), (1,0), (0,1), (0,-1), (1,1), (1,-1), (-1,1), (-1,-1) ]
		v = v1*i + v2*j; a = a1*i + a2*j; 
		C1 = translate_va(B, v,a);
		z =  minimum(nan_ign, (C1 - C0));
		if(!isfinite(z)) z = 0; end
		C1 = C1 .-z
		C0 = min_nan.(C0, C1); 
	end
	#center at the lower left vertex
	x,y = get_lower_left_vertex(B);
	II = (x+y)%2; 

	S = laplacian(C0,II); if(length(unique(S)) > 3) S = laplacian(C0,(II+1)%2); end
	S_f = copy(S);

	S[isnan.(B)].=NaN;


	S = unpad(S)[1];
	S_f = unpad(S_f)[1];
	S = S[2:end-1,2:end-1]; #remove the border of nans 
	S_f = S_f[2:end-1,2:end-1];

	S[isnan.(S)].=1;
	S_f[isnan.(S_f)].=1;

	S = S .+2; 
	S_f = S_f .+2; 

	c0 = colorant"black"
	c1 = colorant"lightgray"
	c2 = colorant"white"
	c_p = [c0, c1, c2]
	

	S = Int.(S); 	
	S_f = Int.(S_f); 

	S_f[S_f.>=3].=3; S_f[S_f.<=0].=3; 
	S[S.>=3].=3; S[S.<=0].=3; 
	
	A_c = IndirectArray(S, c_p)
	A_cf = IndirectArray(S_f, c_p); 
	return A_c, A_cf; 

end

function gen_figs(curr_depth::Int,max_lattice_size::Int, max_depth::Int, curr_str::String, rec_data)
	mats,alt_mats, l_mats, l_mats_l, anc_data2, anc_data3 = rec_data;
	save_curr_fig(mats, alt_mats, curr_str)
	println("saved: _$(curr_str)")
	if(curr_depth >= max_depth) return; end
	lattice_size = max(compute_odd_lattice_size(curr_str), compute_even_lattice_size(curr_str));
	if(lattice_size >= max_lattice_size) return; end

	for s in ['1', '2', '3']
		next_str = curr_str*s; 
		if(all3(next_str) || all2(next_str)) continue; end
		gen_figs(curr_depth+1, max_lattice_size, max_depth, next_str, get_next_mat(rec_data, next_str, s))
	end
end



function initialize_data(curr_str)
	if(!all2(curr_str) && !all3(curr_str))
		error("not a base case")
	end
	mats, alt_mats = get_base_quadruple(curr_str)
	anc_data2 = mats,alt_mats,curr_str;  #only update if there is a 3 or 1
	anc_data3 = mats,alt_mats,curr_str;  #only update if there is a 2 or 1
	if(length(curr_str) > 0)
		l_mats = get_corrections(curr_str, anc_data2, anc_data3, false);
		l_mats_l = get_corrections(curr_str, anc_data2, anc_data3, true);

	else
		l_mats = mats; 
		l_mats_l = mats;
	end
	rec_data = mats,alt_mats, l_mats, l_mats_l, anc_data2, anc_data3;
	return rec_data;
end

using LinearAlgebra
function compute_odd_lattice_size(str)
	p0,q0,p1,q1 = getFracs(str); 
	v1, v2 = vf(p0); v1 = v1*2; 
	M = [real(v1) real(v2); imag(v1) imag(v2)]
	return det(M)
end
function compute_even_lattice_size(str)
	p0,q0,p1,q1 = getFracs(str); 
	v1, v2 = vf(q0);
	M = [real(v1) real(v2); imag(v1) imag(v2)]
	return det(M)
end
function gen_all_figs(max_k, max_depth, max_lattice_size = 999) 
	#iterate over all the base cases
	for i in 0:max_k
		for s in ['2', '3']
			if(s == '3' && i == 0) continue; end #don't double dip
			curr_str = ""; for j in 1:i; curr_str = curr_str*s; end

			lattice_size = max(compute_odd_lattice_size(curr_str), compute_even_lattice_size(curr_str));
			if(lattice_size > max_lattice_size) continue; end

			rec_data = initialize_data(curr_str)
			println("the ancestor is: _$(curr_str)")
			curr_depth = 0;
			gen_figs(curr_depth, max_lattice_size, max_depth, curr_str, rec_data)
		end
	end
end

## code to recursively generate the appendix
#max_k = 8; max_depth = 4;
#gen_all_figs(2, 2)


