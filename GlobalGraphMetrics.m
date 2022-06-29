==================================================
%% Compute the global graph metrics, clustering coefficient and shortest path length 
%% across all subjects for their weighted/binarized network, using different 
%% proportional thresholds. Generate 50 randomized networks with 10 iterations across 
%% all subjects. Compute normalized graph metrics across all subjects.  
==================================================

% Weighted networks 
% Let's compute L while setting negative weights to zero and inverting the weights
Lw = zeros(size(newconn,3),1);

for i = 1:size(newconn,3)
W = newconn(:,:,i);
W(W<0) = 0;
avgw (i) = mean(nonzeros(W));
W(W>0) = 1./W(W>0);
Dw = distance_wei(W);
Lw(i) = mean(squareform(Dw));
end

% Let's compute C while setting negative weights to zero and scaling the connection weights
% by the maximum occurring weight before computing the C.
Cw = zeros(size(newconn,3),1);

for i = 1:size(newconn,3)
    W = newconn(:,:,i);
    W(W<0) = 0;
    Wscaled = W/max(W(:));
    Cw(i) = mean(clustering_coef_wu(Wscaled));
end

%% Randomized networks
for i = 1:size(newconn,3)
    for j = 1:50
        W = newconn(:,:,i);
        W(W<0) = 0;       
        Wscaled = W/max(W(:));
        B = randmio_und(Wscaled,10);
        Crandomw(i,j) = mean(clustering_coef_wu(B));
        
        W(W>0) = 1./W(W>0);
        Wgroup = randmio_und(W,10);
        D = distance_wei(Wgroup);
        Lrandomw(i,j) = mean(squareform(D));
    end
end

% Let's compute one random network per subject
for i = 1:size(newconn,3)
    Lrandomwmain(i) = mean(Lrandomw(i),2);
    Crandomwmain(i) = mean(Crandomw(i),2);
end

%% Normalized metrics
     Lnormalizedw = Lw./mean(Lrandomwmain);
     Cnormalizedw = Cw./mean(Crandomwmain);
     
==================================================
% Binarized networks (here we bring examples with PT% density of 0.4)
% Let's create matrices with density of 0.4 
m = size(newconn,3);
d = zeros(m,1);
L = zeros(m,1);
C = zeros(m,1);

for i = 1:m
    thr_m = threshold_proportional(newconn(:,:,i),0.4);
    d(i) = density_und(thr_m);
    D = distance_bin(thr_m);
    L(i) = mean(squareform(D));
    C(i) = mean(clustering_coef_bu(thr_m));
end

%% Let's compute randomized networks for all subjects with density 0.4
Lrandom = zeros(m,50);
Crandom = zeros(m,50);

for i = 1:m
    thr_m = threshold_proportional(newconn(:,:,i),0.4);
    for j = 1:50
        B = randmio_und(thr_m,10);
        Drandom = distance_bin(B);
        Lrandom(i,j) = mean(squareform(Drandom));
        Crandom(i,j) = mean(clustering_coef_bu(B));
    end
end

% Let's compute one random network per subject
for i = 1:m
    Lrandommain(i) = mean(Lrandom(i),2);
    Crandommain(i) = mean(Crandom(i),2);
end

%% Normalized metrics
% Let's create normalizations of L and C when density is 0.4
     Lnormalized = L./mean(Lrandommain);
     Cnormalized = C./mean(Crandommain);
   
% if needed set any NaN/Inf values to 0
   
