=============================================================================
%% Clean the functional connectivity dataset by first disregarding unsuitable
%% subjects and then applying the outlier criterion on subject's total FC and
%% on intersubject correlation coefficients between functional connectivity
%% matrices. 
%%
%% Start by loading the functional connectivity matrix and the volume matrix.
=============================================================================
  
%% Filter out unsuitable matrices in subjects
cmatrix = squeeze(connectivity(:,:,1,:));
subjects = any(cmatrix,[1 2]);

newconn = cmatrix(:,:,subjects==1); %disregard their connectivity matrices
VOLS = regionProperties(:,4,subjects==1); %disregard their volumes

%% Let's compute the total volume for each subject
m = size(newconn, 3);
vol = zeros(m,1);
for i = 1:m
    temp = sum(VOLS(:,:,i));
    vol(i,:) = temp;
end

%% Clean the dataset from outliers
%Average functional connectivity per subject
FC = zeros(m,1);
for i = 1:size(newconn,3)
    temp = newconn(:,:,i);
    index = find(temp>0);
    FC(i,:) = mean(temp(index));
end

Q1 = quantile(FC,[0.25]); 
Q2 = quantile(FC,[0.5]); 
Q3 = quantile(FC,[0.75]);
IQR = iqr(FC); 
outl1 = fconn(FC(:,1)<(Q1-1.5*IQR)|FC(:,1)>(Q3+1.5*IQR),:);
clear Q1 Q2 Q3 IQR
  
%% Find subjects with strongly deviating connectivity patterns
T = zeros(m);
J = newconn(:,:,:);
J = permute(reshape(J, 114, 114, 470), [3 1 2]); %114 is the amount of nodes, 470 is the amount of subjects
T = corr(J(:,:)'); 

%% Get a single score per subject
score = mean(T,2);

Q1 = quantile(score,[0.25]); 
Q2 = quantile(score,[0.5]); 
Q3 = quantile(score,[0.75]); 
IQR = iqr(score);
outl2 = score(score<(Q1-1.5*IQR)|score>(Q3+1.5*IQR),:); 

%Let's disregard the subjects identified as potential outliers in both methods.
%For example subject 118 was identified as an outlier
newconn(:,:,1) = [];
vol(1) = [];
VOLS(:,:,1) = [];
FC(1) = [];

subs = [1:size(newconn,3)].'; %update the amount of included subjects        
