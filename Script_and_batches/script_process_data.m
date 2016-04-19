%% SCRIPT
%
% Script written to process data from Sophie Gillain.
% Data are organised into subdirectories named IDxxx (xxx = subjct's ID)
% and its subdirectories for FLAIR image, mask and MPM's.
% 1/ Reoganize the data into some sudirectories
% 2/ create the mask images with the LST tool for all the subjects.
%__________________________________________________________________________
% Copyright (C) 2016 Cyclotron Research Centre

% Written by Christophe Phillips
% University of Liège, Belgium

% Find list of subjects' directorie
DATA_rootpth = 'C:\IRM\MRI\MRI_Gabi\MRI_newprotocole_subjectmaps_preprocessingKB';
[files,dirs] = spm_select('List',DATA_rootpth,'^ID.*');
Ndirs = size(dirs,1)

% Rearrange T1w files
for ii=1:Ndirs
    newDir = fullfile(DATA_rootpth,dirs(ii,:),'classified_data','T1_LST');
    mkdir(newDir)
    pth_T1w = fullfile(DATA_rootpth,dirs(ii,:),'classified_data','MP','MT');
    fn_T1w = spm_select('FPList',pth_T1w,'^s.*_T1w\.nii$');
    copyfile(fn_T1w,newDir)
end

% Build T1 & FLAIR image lists
fn_T1 = cell(Ndirs,1);
fn_FLAIR = cell(Ndirs,1);
for ii=1:Ndirs
    dr_T1 = fullfile(DATA_rootpth,dirs(ii,:),'classified_data','T1_LST');
    fn_T1{ii} = spm_select('FPList',dr_T1,'^s.*_T1w\.nii$');
    dr_FLAIR = fullfile(DATA_rootpth,dirs(ii,:),'classified_data','Flair');
    fn_FLAIR{ii} = spm_select('FPList',dr_FLAIR,'^s.*\.nii$');
end

% Build LST batch
fn_batch = fullfile('C:\IRM\MRI\MRI_Gabi\batch_cp','batch_LSTproces_empty.mat');
load(fn_batch)
matlabbatch{1}.spm.tools.LST.lesiongrow.data_T1 = fn_T1;
matlabbatch{1}.spm.tools.LST.lesiongrow.data_FLAIR = fn_FLAIR;
fn_batch_filled = fullfile('C:\IRM\MRI\MRI_Gabi\batch_cp','batch_proces_filled.mat');
save(fn_batch_filled,'matlabbatch')


