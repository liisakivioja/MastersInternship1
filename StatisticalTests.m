=============================================================================
%% Lets evaluate our results statistically. First, let's find the correlation between 
%% cortical volume and overall functional connectivity, and cortical volume and graph metrics 
%% across the whole HCP dataset using Pearson correlation. Then, let's assess the statistical 
%% differences in chosen metrics using student's t-tests. This is done by grouping students 
%% into top n = 200, top n = 100, top n = 50 subjects with higher and lower cortical volumes, 
%% and then performing the t-test between the chosen groups.
=============================================================================

%% Let's test the correlations across the whole HCP dataset
% First, the correlation between cortical volume and functional connectivity
[R,P] = corrcoef(vol,FC);

%Then, the correlation between cortical volume and shortest path length
[R,P] = corrcoef(vol,Lw);

%Lastly, the correlation between cortical volume and clustering coefficient
[R,P] = corrcoef(vol,Cw);

=============================================================================

%% Let's perform the t-test on the top n = 50, n = 100 and n = 200 subjects
topw = [Lw Cw vol];
topw = sortrows(topw,3);

% top 50
top50l = topw(1:50,:);
top50h = topw(406:end,:);
[h,p,ci,stats] = ttest2(top50h(:,1),top50l(:,1)); 
[h,p,ci,stats] = ttest2(top50h(:,2),top50l(:,2));

% top 100
top100l = topw(1:100,:);
top100h = topw(356:end,:);
[h,p,ci,stats] = ttest2(top100h(:,1),top100l(:,1));
[h,p,ci,stats] = ttest2(top100h(:,2),top100l(:,2));

% top 200
top200l = topw(1:200,:);
top200h = topw(256:end,:);
[h,p,ci,stats] = ttest2(top200h(:,1),top200l(:,1));
[h,p,ci,stats] = ttest2(top200h(:,2),top200l(:,2));
