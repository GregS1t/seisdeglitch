% Can run only after fit_glitch.m exec
% to generate clusters_of_glitches variable

h =[];
f = 1;
for i=2:length(clusters_of_glitches)
        h(f) = length(clusters_of_glitches{i});
        f=f+1;
end

figure (1000)
histogram(h)
xlabel("Size of the clusters of glitches")
ylabel("Number of clusters")
title("Histogram of the size of the cluster of glitches")
grid on


