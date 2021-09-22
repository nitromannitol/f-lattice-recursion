## here we check that l corrections overlay properly on the standard tile



## 4 cases
str = "11111111"; 

tt = 5

global bool1 = true ;
#toggle to get correct double decomp for plt1 -true or plt2 -false must be true for plt3,plt4



p1,q1,p0,q0 = getFracs(str); 
plt1 = plot_and_annotate_std_tile(q1, isAnnotate=true,textsize = tt); plot!(legend=false);
#plt1 = plot_and_annotate_single_std(q1, isAnnotate=true,textsize = tt); plot!(legend=false);

plot_and_annotate_corrected_tile!(q1, width = 1, isFlip = false, textsize = tt)



p1,q1,p0,q0 = getFracs(str*"2"); 
#plt2 = plot_and_annotate_std_tile(q1, isAnnotate=true,textsize = tt); plot!(legend=false);
plt2 = plot_and_annotate_single_std(q1, isAnnotate=true,textsize = tt); plot!(legend=false);

plot_and_annotate_corrected_tile!(q1, width = 1, isFlip = false, textsize = tt)

p1,q1,p0,q0 = getFracs(str*"3"); 
#plt3 = plot_and_annotate_std_tile(q1, isAnnotate=true,textsize = tt); plot!(legend=false);
plt3 = plot_and_annotate_single_std(q1, isAnnotate=true,textsize = tt); plot!(legend=false);

plot_and_annotate_corrected_tile!(q1, width = 1, isFlip = false, textsize = tt)

p1,q1,p0,q0 = getFracs(str*"23"); 
#plt4 = plot_and_annotate_std_tile(q1, isAnnotate=true,textsize = tt); plot!(legend=false);
plt4 = plot_and_annotate_single_std(q1, isAnnotate=true,textsize = tt); plot!(legend=false);

plot_and_annotate_corrected_tile!(q1, width = 1, isFlip = false, textsize = tt)



k=3;ll = true
p1,q1,p0,q0 = getFracs(str*"3"^k); 
plt5 = plot_and_annotate_std_tile(q1, isAnnotate=true,textsize = tt); plot!(legend=false);
#plt5 = plot_and_annotate_single_std(q1, isAnnotate=true,textsize = tt); plot!(legend=false);

plot_and_annotate_corrected_tile!(q1, width = 1, isFlip = false, textsize = tt, isLong = ll)

p1,q1,p0,q0 = getFracs(str*"23"*"3"^(k-1)); 
plt6 = plot_and_annotate_std_tile(q1, isAnnotate=true,textsize = tt); plot!(legend=false);
#plt6 = plot_and_annotate_single_std(q1, isAnnotate=true,textsize = tt); plot!(legend=false);

plot_and_annotate_corrected_tile!(q1, width = 1, isFlip = false, textsize = tt, isLong = ll)


