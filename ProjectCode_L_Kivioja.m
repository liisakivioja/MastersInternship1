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

%% Create arrays with subject nr, L, C and their brain vol
info = [];
info = [subs L C vol];
info2 = [];
infor2 = [subs L2 C2 vol];

