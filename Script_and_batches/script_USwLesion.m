%% SCRIPT
%
% Script written for the processing of Sophie Gillain's data.
% Data are organised into subdirectories named IDxxx (xxx = subjct's ID)
% and its subdirectories for FLAIR image, mask and MPM's.
% It applies the whole USwLesion processing: segmentation, BD smoothing and
% parameter extraction
%__________________________________________________________________________
% Copyright (C) 2016 Cyclotron Research Centre

% Written by Christophe Phillips
% University of Liège, Belgium


% Find list of subjects' directorie
DATA_rootpth = 'C:\IRM\MRI\MRI_Gabi\MRI_newprotocole_subjectmaps_preprocessingKB';
[files,dirs] = spm_select('List',DATA_rootpth,'^ID.*');
Ndirs = size(dirs,1)

% fn_batch_empty = fullfile(spm('dir'),'toolbox','USwLesion','batch_MSprocess_empty');
fn_batch_empty = 'batch_MSprocess_empty';

%% Loop over all subjects
% -> select data
% -> fill empty batch
% -> process subject
for ii=1:Ndirs
    % Select data
    dr_FLAIR = fullfile(DATA_rootpth,dirs(ii,:),'classified_data','Flair');
    fn_FLAIR = spm_select('FPList',dr_FLAIR,'^s.*\.nii$');
    dr_MPM = fullfile(DATA_rootpth,dirs(ii,:),'classified_data','MP','MT');
    fn_MT  = spm_select('FPList',dr_MPM,'^s.*_MT\.nii$');
    fn_A   = spm_select('FPList',dr_MPM,'^s.*_A\.nii$');
    fn_R1  = spm_select('FPList',dr_MPM,'^s.*_R1\.nii$');
    fn_R2s = spm_select('FPList',dr_MPM,'^s.*_R2s\.nii$');
    dr_MSK = fullfile(DATA_rootpth,dirs(ii,:),'classified_data','T1_LST');
    fn_msk = spm_select('FPList',dr_MSK,'^lesion_.*\.nii$');
    % Fill & save batch
    eval(fn_batch_empty)
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {{fn_msk}};
    matlabbatch{2}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = { ...
        {fn_MT}, {fn_A}, {fn_R1}, {fn_R2s} };
    matlabbatch{3}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {{fn_FLAIR}};
    save(['batch_USwL_',dirs(ii,:)],'matlabbatch')
%     % Execute batch
%      spm_jobman('run',matlabbatch);
end

%% Next step:
