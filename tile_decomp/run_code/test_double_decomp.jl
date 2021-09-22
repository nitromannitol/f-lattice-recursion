## code for generating double decomposition 

### S2 CASE



str = "1111112"

p0,q0,p1,q1 = getFracs(str); 
isFlip = false; 
isAnnotate=false
textsize = 10; 
p = q0; 

plot(aspect_ratio = true, axis = :off, ticks = :none, legend=false)
plt1 = plot_and_annotate_corrected_tile!(p, isFlip = isFlip, textsize = textsize);
#plot(plt1,legend)


textsize=6
plot(aspect_ratio = true, axis = :off, ticks = :none, legend=false)

offSet=0
k = count_end(str, '3')
str = str[1:end-k]
p0,q0,p1,q1 = getFracs(str)

## use the first two
tp1, tq1 = getTile(p1), getTile(q1)
tp1,tq1 = tp1.+offSet, tq1 .+offSet
vp1,vp2 = vf(p1); vq1,vq2 = vf(q1)
O=0;
num=2
#plot_and_annotate(tp1, p1, O = O, label = p1, textsize = textsize,num = num, ann = "");
plot_and_annotate_std_tile!(p1,width = width, offSet = offSet, isDouble = false, textsize = textsize, num = num+1, isAnnotate = true, isRecurse = false)

O = O -vq1;
plot_and_annotate(tq1, q1, O = O, label = q1, textsize = textsize,num = num, ann = "");
O = O + vq2
plot_and_annotate_alt_tile!(p1, annotate_center=true,width = width, offSet = O+offSet, textsize = textsize, num = num+1, isAnnotate = true)
#else
#plot_and_annotate_single_alt!(p1, num = num, offSet=O+offSet, width = 1, textsize=textsize, isAnnotate=isAnnotate, ann = "")


plt2 = plot!(legend=false)

plot(plt1, plt2)

#savefig(plot(plt1,dpi=512), "l_correct_s2")
#savefig(plot(plt2,dpi=512),"l_correct_s2_double_decomp")



### S1 CASE

str = "11111111"

p0,q0,p1,q1 = getFracs(str); 
isFlip = false; 
isAnnotate=false
textsize = 10; 
p = q0; 

plot(aspect_ratio = true, axis = :off, ticks = :none, legend=false)
plt1 = plot_and_annotate_corrected_tile!(p, isFlip = isFlip, textsize = textsize);
#plot(plt1,legend)


textsize=6
plot(aspect_ratio = true, axis = :off, ticks = :none, legend=false)

offSet=0
k = count_end(str, '3')
str = str[1:end-k]
p0,q0,p1,q1 = getFracs(str)

## use the first two
tp1, tq1 = getTile(p1), getTile(q1)
tp1,tq1 = tp1.+offSet, tq1 .+offSet
vp1,vp2 = vf(p1); vq1,vq2 = vf(q1)
O=0;
num=2
plot_and_annotate_std_tile!(p1,width = width, offSet = O, isDouble = false, textsize = textsize, num = num+1, isAnnotate = true, isRecurse = false)

O = O -vq1;
plot_and_annotate_std_tile!(q1,width = width, offSet = O, isDouble = false, textsize = textsize, num = num+1, isAnnotate = true, isRecurse = false)
O = O + vp2+2vp1-vq1
plot_and_annotate_alt_tile!(p1, annotate_center=true,width = width, offSet = O+offSet, textsize = textsize, num = num+1, isAnnotate = true)
#else
#plot_and_annotate_single_alt!(p1, num = num, offSet=O+offSet, width = 1, textsize=textsize, isAnnotate=isAnnotate, ann = "")


plt2 = plot!(legend=false)

plot(plt1, plt2)


savefig(plot(plt1,dpi=512), "Figures/l_correct/s1")
savefig(plot(plt2,dpi=512),"Figures/l_correct/s1_double")
