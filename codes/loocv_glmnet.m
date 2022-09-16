clear;clc; %%% nested LOOCV prediction
options       = glmnetSet;
all_alpha     = [];
aaa           = [0:0.01:1];%%% grid of alpha value

load('E:\IFC\upload_github\common_index.mat')%%load common features'index
load('E:\IFC\upload_github\feature_index.mat')%%load features' index across each LOOCV
load('E:\IFC\upload_github\RT.mat')%%load observed RT
load('E:\IFC\upload_github\pcor001.mat')%%load features across each LOOCV

for alphaa        = aaa
    options.alpha = alphaa;
    for subject = 1:length(RT)
        select_TFC_lesion = cell_selected_TFC{subject};
        newx = select_TFC_lesion(subject,:);
        RT(subject,:) = []; select_TFC_lesion(subject,:) = []; result(100,1) = 0;
        for iterate = 1:100  
            CVerr = cvglmnet(select_TFC_lesion,RT,'gaussian',options,'mae',3);
            ncoef = cvglmnetCoef(CVerr,'lambda_min'); ncoef(1,:) = []; ncoef = ncoef';
            nomrlized_ncoef = ncoef.*(std(select_TFC_lesion)./std(RT));
            [C,ia,ib] = intersect(features_index{subject}',common_index,'sorted');
            selected_index(iterate,:) = nomrlized_ncoef(:,ia);
            result(iterate,:) = glmnetPredict(CVerr.glmnet_fit,newx,CVerr.lambda_min,'link');
            CVerr = []; ncoef = []; nomrlized_ncoef = []; C = []; ia = []; ib = []; 
        end
        index = find(aaa==alphaa);
        predictive_result{index}(subject,:) = mean(result);
        beta_result{index}(subject,:) = mean(selected_index); RT = []; select_TFC_lesion = []; 
        load('E:\IFC\upload_github\RT.mat')%%load observed RT
        load('E:\IFC\upload_github\pcor001.mat')%%load features across each LOOCV
    end
    disp(['alphaa1--' num2str(index), ' is done'])
end

cd('E:\IFC\upload_github')
save(['PCOR_005_RESULTS.mat'],...
    'predictive_result');
save(['PCOR_005_beta.mat'],...
    'beta_result');




%%%%%%%%%%%%%%%%%% permutation test
%%
clear;clc;
options                   = glmnetSet;
options.alpha             = 0.99; %%% optimal alpha value
permutation_results(48,1) = 0;

load('E:\IFC\upload_github\RT.mat')%%load observed RT
load('E:\IFC\upload_github\pcor001.mat')%%load features across each LOOCV

for subject = 1:length(RT)
    select_TFC_lesion = cell_selected_TFC{subject};
    newx = select_TFC_lesion(subject,:);
    RT(subject,:) = [];  select_TFC_lesion(subject,:) = []; 
        for interted = 1:10000
            sort_RT = randperm(size(RT,1));
            random_RT = RT(sort_RT,:);
            CVerr = cvglmnet(select_TFC_lesion,random_RT,'gaussian',options,'mae',3);
            result = glmnetPredict(CVerr.glmnet_fit,newx,CVerr.lambda_min,'link');  
            CVerr = [];
            permutation_results(subject,interted) = result;
            clear sort_RT random_RT
        end
    disp(['subject' num2str(subject), ' is done']); RT = []; select_TFC_lesion = []; 
    load('E:\IFC\upload_github\RT.mat')%%load observed RT
    load('E:\IFC\upload_github\pcor001.mat')%%load features across each LOOCV
end

cd('E:\IFC\upload_github')
save(['permutation_results.mat'],...
    'permutation_results');
