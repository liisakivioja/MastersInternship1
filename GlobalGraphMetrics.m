==================================================
%% Compute the global graph metrics clustering coefficient and shortest 
%% path length across all subjects for their binarized network, using 
%% different proportional thresholds. Generate 50 randomized networks 
%% across all subjects. Compute normalized graph metrics across all 
%% subjects. Finish by computing the small-world index across all subjects.
==================================================

%% Let's create matrices with density of 0.45
m = size(newconn,3);
d = zeros(m,1);
L = zeros(m,1);
C = zeros(m,1);

for i = 1:m
    thr_m = threshold_proportional(newconn(:,:,i),0.45);
    d(i) = density_und(thr_m);
    D = distance_bin(thr_m);
    L(i) = mean(squareform(D));
    C(i) = mean(clustering_coef_bu(thr_m));
end

%% Let's create matrices with density of 0.4
d2 = zeros(m,1);
L2 = zeros(m,1);
C2 = zeros(m,1);

for i = 1:m
    thr_m = threshold_proportional(newconn(:,:,i),0.4);
    d2(i) = density_und(thr_m);
    D2 = distance_bin(thr_m);
    L2(i) = mean(squareform(D2));
    C2(i) = mean(clustering_coef_bu(thr_m));
end

%% Let's create matrices with density of 0.35
d3 = zeros(m,1);
L3 = zeros(m,1);
C3 = zeros(m,1);

for i = 1:m
    thr_m = threshold_proportional(newconn(:,:,i),0.35);
    d3(i) = density_und(thr_m);
    D3 = distance_bin(thr_m);
    L3(i) = mean(squareform(D3));
    C3(i) = mean(clustering_coef_bu(thr_m));
end

%% Let's create matrices with density of 0.3
d4 = zeros(m,1);
L4 = zeros(m,1);
C4 = zeros(m,1);

for i = 1:m
    thr_m = threshold_proportional(newconn(:,:,i),0.3);
    d4(i) = density_und(thr_m);
    D4 = distance_bin(thr_m);
    L4(i) = mean(squareform(D4));
    C4(i) = mean(clustering_coef_bu(thr_m));
end

%% Let's create matrices with density of 0.25
d5 = zeros(m,1);
L5 = zeros(m,1);
C5 = zeros(m,1);

for i = 1:m
    thr_m = threshold_proportional(newconn(:,:,i),0.25);
    d5(i) = density_und(thr_m);
    D5 = distance_bin(thr_m);
    L5(i) = mean(squareform(D5));
    C5(i) = mean(clustering_coef_bu(thr_m));
end

%% Let's create matrices with density of 0.2
d6 = zeros(m,1);
L6 = zeros(m,1);
C6 = zeros(m,1);

for i = 1:m
    thr_m = threshold_proportional(newconn(:,:,i),0.2);
    d6(i) = density_und(thr_m);
    D6 = distance_bin(thr_m);
    L6(i) = mean(squareform(D6));
    C6(i) = mean(clustering_coef_bu(thr_m));
end

%% Let's create matrices with density of 0.15
d7 = zeros(m,1);
L7 = zeros(m,1);
C7 = zeros(m,1);

for i = 1:m
    thr_m = threshold_proportional(newconn(:,:,i),0.15);
    d7(i) = density_und(thr_m);
    D7 = distance_bin(thr_m);
    L7(i) = mean(squareform(D7));
    C7(i) = mean(clustering_coef_bu(thr_m));
end

=================================================

%% Let's compute randomized networks for all subjects with density 0.45
Lrandom = zeros(m,50);
Crandom = zeros(m,50);

for i = 1:m
    thr_m = threshold_proportional(newconn(:,:,i),0.45);
    for j = 1:50
        B = randmio_und(thr_m,10);
        Drandom = distance_bin(B);
        Lrandom(i,j) = mean(squareform(Drandom));
        Crandom(i,j) = mean(clustering_coef_bu(B));
    end
end

%% Let's compute randomized networks for all subjects with density 0.40
Lrandom2 = zeros(m,50);
Crandom2 = zeros(m,50);

for i = 1:m
    thr_m = threshold_proportional(newconn(:,:,i),0.4);
    for j = 1:50
        B2 = randmio_und(thr_m,10);
        Drandom2 = distance_bin(B2);
        Lrandom2(i,j) = mean(squareform(Drandom2));
        Crandom2(i,j) = mean(clustering_coef_bu(B2));
    end
end

%% Let's compute randomized networks for all subjects with density 0.35
Lrandom3 = zeros(m,50);
Crandom3 = zeros(m,50);

for i = 1:m
    thr_m = threshold_proportional(newconn(:,:,i),0.35);
    for j = 1:50
        B3 = randmio_und(thr_m,10);
        Drandom3 = distance_bin(B3);
        Lrandom3(i,j) = mean(squareform(Drandom3));
        Crandom3(i,j) = mean(clustering_coef_bu(B3));
    end
end

%% Let's compute randomized networks for all subjects with density 0.3
Lrandom4 = zeros(m,50);
Crandom4 = zeros(m,50);

for i = 1:m
    thr_m = threshold_proportional(newconn(:,:,i),0.3);
    for j = 1:50
        B4 = randmio_und(thr_m,10);
        Drandom4 = distance_bin(B4);
        Lrandom4(i,j) = mean(squareform(Drandom4));
        Crandom4(i,j) = mean(clustering_coef_bu(B4));
    end
end

%% Let's compute randomized networks for all subjects with density 0.25
Lrandom5 = zeros(m,50);
Crandom5 = zeros(m,50);

for i = 1:m
    thr_m = threshold_proportional(newconn(:,:,i),0.25);
    for j = 1:50
        B5 = randmio_und(thr_m,10);
        Drandom5 = distance_bin(B5);
        Lrandom5(i,j) = mean(squareform(Drandom5));
        Crandom5(i,j) = mean(clustering_coef_bu(B5));
    end
end

%% Let's compute randomized networks for all subjects with density 0.2
Lrandom6 = zeros(m,50);
Crandom6 = zeros(m,50);

for i = 1:m
    thr_m = threshold_proportional(newconn(:,:,i),0.2);
    for j = 1:50
        B6 = randmio_und(thr_m,10);
        Drandom6 = distance_bin(B6);
        Lrandom6(i,j) = mean(squareform(Drandom6));
        Crandom6(i,j) = mean(clustering_coef_bu(B6));
    end
end

%% Let's compute randomized networks for all subjects with density 0.15
Lrandom7 = zeros(m,50);
Crandom7 = zeros(m,50);

for i = 1:m
    thr_m = threshold_proportional(newconn(:,:,i),0.15);
    for j = 1:50
        B7 = randmio_und(thr_m,10);
        Drandom7 = distance_bin(B7);
        Lrandom7(i,j) = mean(squareform(Drandom7));
        Crandom7(i,j) = mean(clustering_coef_bu(B7));
    end
end

=======================================================

%% Let's compute L- and Crandom across all subjects
%When density is 0.45
for i = 1:m
    Lrandommain(i) = mean(Lrandom(i),2);
    Crandommain(i) = mean(Crandom(i),2);
end

%When density is 0.4
for i = 1:m
    Lrandommain2(i) = mean(Lrandom2(i),2);
    Crandommain2(i) = mean(Crandom2(i),2);
end

%When density is 0.35
for i = 1:m
    Lrandommain3(i) = mean(Lrandom3(i),2);
    Crandommain3(i) = mean(Crandom3(i),2);
end

%When density is 0.3
for i = 1:m
    Lrandommain4(i) = mean(Lrandom4(i),2);
    Crandommain4(i) = mean(Crandom4(i),2);
end

%When density is 0.25
for i = 1:m
    Lrandommain5(i) = mean(Lrandom5(i),2);
    Crandommain5(i) = mean(Crandom5(i),2);
end

%When density is 0.2
for i = 1:m
    Lrandommain6(i) = mean(Lrandom6(i),2);
    Crandommain6(i) = mean(Crandom6(i),2);
end

%When density is 0.15
for i = 1:m
    Lrandommain7(i) = mean(Lrandom7(i),2);
    Crandommain7(i) = mean(Crandom7(i),2);
end

================================================
%% Let's create normalizations of L and C 
%When density is 0.45
 for i = 1:m
     Lnormalized(i) = L(i)/mean(Lrandommain(i));
     Cnormalized(i) = C(i)/mean(Crandommain(i));
 end 
 
%When density is 0.4
  for i = 1:m
     Lnormalized2(i) = L2(i)/mean(Lrandommain2(i));
     Cnormalized2(i) = C2(i)/mean(Crandommain2(i));
  end 
  
%When density is 0.35
  for i = 1:m
     Lnormalized3(i) = L3(i)/mean(Lrandommain3(i));
     Cnormalized3(i) = C3(i)/mean(Crandommain3(i));
  end 

%When density is 0.3
 for i = 1:m
     Lnormalized4(i) = L4(i)/mean(Lrandommain4(i));
     Cnormalized4(i) = C4(i)/mean(Crandommain4(i));
 end 
 
%When density is 0.25
 for i = 1:m
     Lnormalized5(i) = L5(i)/mean(Lrandommain5(i));
     Cnormalized5(i) = C5(i)/mean(Crandommain5(i));
 end 
 
 %When density is 0.2
 for i = 1:m
     Lnormalized6(i) = L6(i)/mean(Lrandommain6(i));
     Cnormalized6(i) = C6(i)/mean(Crandommain6(i));
 end 

  %When density is 0.15
 for i = 1:m
     Lnormalized7(i) = L7(i)/mean(Lrandommain7(i));
     Cnormalized7(i) = C7(i)/mean(Crandommain7(i));
 end
 
 % if needed set any NaN/Inf values to 0
   
===============================================
%% Let's compute SW index across subjects
  %When density is 0.45
for i = 1:m
    SW(i) = Cnormalized(i)/Lnormalized(i);
end


  %When density is 0.4
for i = 1:m
    SW2(i) = Cnormalized2(i)/Lnormalized2(i);
end

  %When density is 0.35
for i = 1:m
    SW3(i) = Cnormalized3(i)/Lnormalized3(i);
end

%When density is 0.3
for i = 1:m
    SW4(i) = Cnormalized4(i)/Lnormalized4(i);
end

%When density is 0.25
for i = 1:m
    SW5(i) = Cnormalized5(i)/Lnormalized5(i);
end

%When density is 0.2
for i = 1:m
    SW6(i) = Cnormalized6(i)/Lnormalized6(i);
end

%When density is 0.15
for i = 1:m
    SW7(i) = Cnormalized7(i)/Lnormalized7(i);
end

 % if needed set any NaN/Inf values to 0
