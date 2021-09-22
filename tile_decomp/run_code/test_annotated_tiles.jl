## code for looking at pictures of alt tiles

str = "11113"
str = "11111111"

p0,q0,p1,q1 = getFracs(str); 


## EVEN
p = q0; 

v1,v2 = vf(p); if(isOdd(p)) v1 = v1*2; end 
textsize = 7; width = 1


p, q = getParents(p); if(isOdd(q)) p,q = q,p; end 
vp1,vp2 = vf(p); vp1=2vp1; vq1,vq2 = vf(q); 



plot_and_annotate_alt_tile(q0, textsize=textsize, width = width);
for v in [v1+vq1, -v2, -v2+v1+vq1]
#for v in [v1,-v2, -v2+v1]
	plot_and_annotate_alt_tile!(q0, offSet = v, textsize = textsize, width = width);
end
plt0 = plot!(legend=false);



plot_and_annotate_alt_tile(q0, textsize=textsize, width = width);
for v in [v1, -v2, -v2+v1]
	plot_and_annotate_alt_tile!(q0, offSet = v, textsize = textsize, width = width);
end
plt1 = plot!(legend=false);


plot(plt0,plt1)


#p = p0; 
#double_decomp_plot_and_annotate_alt_tile(p, isDoubleDecomp=true);
#plt1 = plot!(legend=false);
#lot_and_annotate_alt_tile(p);
#plt2 = plot!(legend=false);
#double_decomp_plot_and_annotate_std_tile(q1); plt3 = plot!(legend=false);
#plot(plt1,plt2,plt3)



## compare against pattern
#mats, alt_mats, rot_mats, fracs, offSets,reps,tsum = iterate_str(str, init_mats()...)
#Bp1,Bq1,p1,q1 = get_alts(mats, alt_mats, rot_mats, fracs, offSets,reps,tsum);
#Ap1,Aq1 = mats[1],mats[2];

#plt1 = full_plot_alt(Bp1,p1)
#plt2 = full_plot(mats[4],fracs[4])
#double_decomp_plot_and_annotate_alt_tile(p, isDoubleDecomp=true);plt3 = plot!(legend=false);

savefig(plot(plt0, dpi = 512), "alt_tile_tiling_1.png")
savefig(plot(plt1, dpi = 512), "alt_tile_tiling_2.png")
