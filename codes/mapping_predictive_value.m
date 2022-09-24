clear;clc;
load('F:\HBM_paper\real_leave_one_subjects_out\james_shine\beta\PCOR_001_beta.mat')%%%%%% load selected beta weights
results = abs(final_index{100})';%%%%%% optimal alpha value

aaa   =  1:60031;%%%%%% total link number
m     =  length(aaa);
n     =  (1+sqrt(1+8*m))/2 ; 
AA    =  zeros(n,n);
index = 1;
for i = 1:n
    for j = i+1:n
        AA(i,j) = aaa(index);  
        index = index+1;
    end
end
AA                   = triu(AA);
AA                   = AA+AA';
AA(eye(size(AA))~=0) = 0;

load('E:\IFC\upload_github\common_index.mat')%%load common features'index
for i = 1:length(common_index)
    [x,y] = find(common_index(i)==AA);
    x_cell(i) = x(1); 
    y_cell(i) = y(1);
end
combine = [x_cell',y_cell'];

for j = 1:48 %%subject index
    for i = 1:347 %%brain index
        [x,y] = find(combine==i);
        regional_subject(i,j) = sum(results(x,j));
        clear x y 
    end
end
    
%%%%%%%
%% mapp value to brain %add rp_writefile and rp_readfile functions to matlab path

[A,B,C] = rp_readfile('E:\HBM_paper_revision\Parcels-19cwpgu\Parcels\Parcels_MNI_333.nii');%%%%%% load cortical mask
for j = 1:48
    selected_subject = regional_subject(:,j);
    for i = 1:length(selected_subject)
        A(A==i) = selected_subject(i,:);
    end
    filedir = 'F:\HBM_paper\real_leave_one_subjects_out\james_shine\beta\pcor005\control_headmotion\correct_error';%%%%%% output dir
    filename = [filedir filesep,'subject',num2str(j),'.nii'];
    rp_writefile(A,filename,size(A),[3,3,3],C,'single');
    clear A B C selected_subject
    [A,B,C] = rp_readfile('E:\HBM_paper_revision\Parcels-19cwpgu\Parcels\Parcels_MNI_333.nii');%%%%%% load cortical mask
end

[A,B,C] = rp_readfile('F:\HBM_paper\real_leave_one_subjects_out\james_shine\beta\pcor001\control_headmotion\correct_error\subcortical.nii');%%%%%% load subcortical mask
index   = [336 337 338 340]; %%%%%% nonzero weights in subcortical
for j = 1:48
    selected_subject = regional_subject(:,j);
    for i = index
        A(A==i)=selected_subject(i,:);
    end
    filedir = 'F:\HBM_paper\real_leave_one_subjects_out\james_shine\beta\pcor001\control_headmotion\correct_error\subcortical';%%%%%% output dir
    filename = [filedir filesep,'subject',num2str(j),'.nii'];
    rp_writefile(A,filename,size(A),[3,3,3],C,'single');
    clear A B C selected_subject
    [A,B,C] = rp_readfile('F:\HBM_paper\real_leave_one_subjects_out\james_shine\beta\pcor001\control_headmotion\correct_error\subcortical.nii');%%%%%% load subcortical mask
end




