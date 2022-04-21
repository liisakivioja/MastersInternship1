==================================================
%% Compute the global graph metrics clustering coefficient and shortest 
%% path length across all subjects using different thresholds. Generate 
%% 100 randomized networks across all subjects. Compute normalized graph
%% metrics across all subjects. Finish by computing the small-world index
%% across all subjects.
==================================================

%% Let's compute graph metrics across all subjects with only positive weights
m = size(newconn,3);
d = zeros(m,1);
L = zeros(m,1);
C = zeros(m,1);

for i = 1:m
    A = double(newconn(:,:,i)>0);
    A(A < 0) = 0;
    d(i) = density_und(A);
    D = distance_bin(A);
    L(i) = mean(squareform(D));
    C(i) = mean(clustering_coef_bu(A));
end

L(isnan(L)) = 0;

%% Let's compute graph metrics across all subjects with weights exceeding threshold 0.15
d2 = zeros(m,1);
L2 = zeros(m,1);
C2 = zeros(m,1);

for i = 1:m
    A2 = double(newconn(:,:,i)>0.15);
    A2(A2 < 0) = 0;
    d2(i) = density_und(A2);
    D2 = distance_bin(A2);
    L2(i) = mean(squareform(D2));
    C2(i) = mean(clustering_coef_bu(A2));
end

L2(isnan(L2)) = 0;

%% Let's create matrices with density of 0.3
d3 = zeros(m,1);
L3 = zeros(m,1);
C3 = zeros(m,1);

for i = 1:m
    thr_m = threshold_proportional(newconn(:,:,i),0.3);
    d3(i) = density_und(thr_m);
    D3 = distance_bin(thr_m);
    L3(i) = mean(squareform(D3));
    C3(i) = mean(clustering_coef_bu(thr_m));
end

%% Let's compute randomized networks for all subjects with postive weights only
Lrandom = zeros(m,100);
Crandom = zeros(m,100);

for i = 1:m
    A = double(newconn(:,:,i)>0);
    for j = 1:100
        B = randmio_und(A,10);
        Drandom = distance_bin(B);
        Lrandom(i,j) = mean(squareform(Drandom));
        Crandom(i,j) = mean(clustering_coef_bu(B));
    end
end

%% Let's compute randomized networks for all subjects with weights exceeding 0.15
Lrandom2 = zeros(m,100);
Crandom2 = zeros(m,100);

for i = 1:m
    A2 = double(newconn(:,:,i)>0.15);
    for j = 1:100
        B2 = randmio_und(A2,10);
        Drandom2 = distance_bin(B2);
        Lrandom2(i,j) = mean(squareform(Drandom2));
        Crandom2(i,j) = mean(clustering_coef_bu(B2));
    end
end

%% Let's compute randomized networks for all subjects with density 0.3
Lrandom3 = zeros(m,100);
Crandom3 = zeros(m,100);

for i = 1:m
    thr_m2 = threshold_proportional(newconn(:,:,i),0.3);
    for j = 1:100
        B3 = randmio_und(thr_m2,10);
        Drandom3 = distance_bin(B3);
        Lrandom3(i,j) = mean(squareform(Drandom3));
        Crandom3(i,j) = mean(clustering_coef_bu(B3));
    end
end

%% Let's compute L- and Crandom across all subjects
%When only positive weights are present
for i = 1:m
    Lrandommain(i) = mean(Lrandom(i),2);
    Crandommain(i) = mean(Crandom(i),2);
end

%When weights exceeding 0.15 are present
for i = 1:m
    Lrandommain2(i) = mean(Lrandom2(i),2);
    Crandommain2(i) = mean(Crandom2(i),2);
end

%When density is 0.3
for i = 1:m
    Lrandommain3(i) = mean(Lrandom3(i),2);
    Crandommain3(i) = mean(Crandom3(i),2);
end

%% Let's create normalizations of L and C 
%When only positive weights are present
 for i = 1:m
     Lnormalized(i) = L(i)/mean(Lrandommain(i));
     Cnormalized(i) = C(i)/mean(Crandommain(i));
 end 
 
%When weights exceeding 0.15 are present
  for i = 1:m
     Lnormalized2(i) = L2(i)/mean(Lrandommain2(i));
     Cnormalized2(i) = C2(i)/mean(Crandommain2(i));
  end 

%When density is 0.3
 for i = 1:m
     Lnormalized3(i) = L3(i)/mean(Lrandommain3(i));
     Cnormalized3(i) = C3(i)/mean(Crandommain3(i));
 end 
 
%% Let's compute SW index across subjects
%When only positive weights are present
for i = 1:m
    SW(i) = Cnormalized(i)/Lnormalized(i);
end

%When weights exceeding 0.15 are present
for i = 1:m
    SW2(i) = Cnormalized2(i)/Lnormalized2(i);
end

%When density is 0.3
for i = 1:m
    SW3(i) = Cnormalized3(i)/Lnormalized3(i);
end
