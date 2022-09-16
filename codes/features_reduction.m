clear;clc;%%%%%%%%selected delta-FC features based on univariate pearson correlation

load('E:\IFC\upload_github\RT.mat')%%load observed RT
load('E:\IFC\upload_github\pcor_diff.mat')%%load delta-FC across training
load('E:\IFC\upload_github\headmotion.mat')%%load mean frame-wise displacement parameters

for subject = 1:length(RT)
    RT(subject,:) = [];
    pcor_cell(subject,:) = [];
    for i = 1:size(pcor_cell,2)
        [r,p] = corr(RT,pcor_cell(:,i));
        p_cell(:,i) = p;
    end
    ccc = find(p_cell < 0.01);
    RT = []; pcor_cell = [];
    load('E:\IFC\upload_github\RT.mat')%%load observed RT
    load('E:\IFC\upload_github\pcor_diff.mat')%%load delta-FC across training
    selected_TFC = pcor_cell(:,ccc);
    for j = 1:length(selected_TFC)
        [r,p] = corr(headmotion,selected_TFC(:,j),'type','Spearman');
        p_headcell(:,j) = p;
        index = find(p_headcell < 0.05);
    end
    selected_TFC(:,index) = []; ccc(:,index) = [];%%%%%%%%exclude features related to headmotion
    clear p_headcell index
    cell_selected_TFC{subject} = selected_TFC; features_index{subject} = ccc;
    clear p_cell ccc selected_TFC  
end
%%%
%%
clear;clc; % find common features index across each LOOCV
load('E:\IFC\upload_github\feature_index.mat')

for subjects  = 1:length(features_index)
    FC_number = 1:60031;%%%%%% total link number
    bbbb      = features_index{subjects};
        for j = 1:length(bbbb)
            FC_number(FC_number==bbbb(:,j))=NaN;
        end
    features_index_subjects(subjects,:) = isnan(FC_number);
    clear FC_number bbbb 
end

features_index_subjects = double(features_index_subjects);
mean_features           = mean(features_index_subjects); 
common_index            = find(mean_features==1); 
common_index            = common_index';
