addpath('C:\Users\Liisa\Documents\MATLAB_masters1\hcp_fc_liisa/hcp');
addpath('C:\Users\Liisa\Documents\MATLAB_masters1\2019_03_03_BCT');
load('connectivity_fmri_gmean_0.01-0.1_lausanne120.mat');
load('region_properties_lausanne120_8seeds.mat');

%% Filter out unsuitable matrices in subjects
cmatrix = squeeze(connectivity(:,:,1,:));
subjects = any(cmatrix,[1 2]);
SUBJECTS = cmatrix(:,:,subjects==1); %disregard connectivity matrices
VOLS = regionProperties(:,4,subjects==1); %disregard volumes

%% Let's compute the total volume for each subject
m = size(SUBJECTS, 3);
vol = zeros(m,1);

for i = 1:m
    temp = sum(VOLS(:,:,i));
    vol(i,:) = temp;
end

%% Let's visualize the data
figure; imagesc(connectivity(:,:,1,1)); colorbar 
hold on
xlabel('Nodes'); ylabel('Nodes'); %title('FC matrix');

%Let's now create a weighted group-averaged FC matrix
Wfc = mean(SUBJECTS(:,:,:),3);
figure; imagesc(Wfc); colorbar
hold on
xlabel('Nodes'); ylabel('Nodes'); %title('Weighted group averaged FC matrix');

%% Clean the dataset from outliers
% Average functional connectivity per subject
for i = 1:size(SUBJECTS,3)
    sl = (squareform(SUBJECTS(:,:,i)));
    avg = mean(sl,'all'); 
    fconn(i,:) = avg;
end

s = [1:size(SUBJECTS,3)].'; %create subjects list
fconn = [fconn s];

Q1 = quantile(fconn,[0.25]); %0.0101
Q2 = quantile(fconn,[0.5]); %0.0155
Q3 = quantile(fconn,[0.75]); %0.0222
IQR = iqr(fconn); %0.0122
afc = fconn(fconn(:,1)<(Q1-1.5*IQR)|fconn(:,1)>(Q3+1.5*IQR),:);
% In total 21 subs. Visually, subs 91,118,363 and 434 stand out.

%% Also, let's find subjects with strongly deviating connectivity patterns
T = zeros(m);
J = SUBJECTS(:,:,:);
J = permute(reshape(J, 114, 114, 470), [3 1 2]);
T = corr(J(:,:)'); 

%Let's visualize T
figure; imagesc(T); colorbar

%% Get a single score per subject
score = mean(T,2);
subs = [1:size(score)].';
scatter(subs,score);

subscore = []; 
subscore = [subs score];

Q1s = quantile(score,[0.25]); %0.0.5435
Q2s = quantile(score,[0.5]); %0.5673
Q3s = quantile(score,[0.75]); %0.5896
IQR2 = iqr(score); %0.0461
scorepersub = subscore(subscore(:,2)<(Q1s-1.5*IQR2)|subscore(:,2)>(Q3s+1.5*IQR2),:);
% In total 10 subs. Visually, subs 1,253,411 and 434 stand out.
% Subs 118 and 434 are identified outliers in both techniques.

%% Let's disregard these subjects
SUBJECTS(:,:,118) = [];
SUBJECTS(:,:,434) = [];
vol(118) = [];
vol(434) = [];
subs = [1:size(SUBJECTS,3)].';

%% Let's compute graph metrics across all subjects with only positive weights
m = size(SUBJECTS,3);
d = zeros(m,1);
L = zeros(m,1);
C = zeros(m,1);

for i = 1:m
    A = double(SUBJECTS(:,:,i)>0);
    A(A < 0) = 0;
    d(i) = density_und(A);
    D = distance_bin(A);
    L(i) = mean(squareform(D));
    C(i) = mean(clustering_coef_bu(A));
end

L(isnan(L)) = 0;

%% Let's compute graph metrics across all subjects with weights exceeding threshold 0.15
m = size(SUBJECTS,3);
d2 = zeros(m,1);
L2 = zeros(m,1);
C2 = zeros(m,1);

for i = 1:m
    A2 = double(SUBJECTS(:,:,i)>0.15);
    A2(A2 < 0) = 0;
    d2(i) = density_und(A2);
    D2 = distance_bin(A2);
    L2(i) = mean(squareform(D2));
    C2(i) = mean(clustering_coef_bu(A2));
end

L2(isnan(L2)) = 0;

%% Let's also create matrices with density of 0.3
for i = 1:m
    thr_m = threshold_proportional(SUBJECTS(:,:,i),0.3);
    d3(i) = density_und(thr_m);
    D3 = distance_bin(thr_m);
    L3(i) = mean(squareform(D3));
    C3(i) = mean(clustering_coef_bu(thr_m));
end

%% Let's compute randomized networks for all subjects with postive weights only
Lrandom = zeros(m,50);
Crandom = zeros(m,50);

for i = 1:m
    A = SUBJECTS(:,:,i)>0;
    for j = 1:10
        B = randmio_und(A,10);
        Drandom = distance_bin(B);
        Lrandom(i,j) = mean(squareform(Drandom));
        Crandom(i,j) = mean(clustering_coef_bu(B));
    end
end

%% Let's compute randomized networks for all subjects with weights exceeding 0.15
Lrandom2 = zeros(m,50);
Crandom2 = zeros(m,50);

for i = 1:m
    A2 = SUBJECTS(:,:,i)>0.15;
    for j = 1:10
        B2 = randmio_und(A2,10);
        Drandom2 = distance_bin(B2);
        Lrandom2(i,j) = mean(squareform(Drandom2));
        Crandom2(i,j) = mean(clustering_coef_bu(B2));
    end
end

%% Let's compute randomized networks for all subjects with density 0.3
Lrandom3 = zeros(m,50);
Crandom3 = zeros(m,50);

for i = 1:m
    thr_m2 = threshold_proportional(SUBJECTS(:,:,i),0.3);
    for j = 1:10
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
 
