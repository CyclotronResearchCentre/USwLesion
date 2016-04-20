%% SCRIPT: MS data
%
% Script to rearrange, i.e. copy, the subset of patient data for the whole
% processing (segmentation, smoothing, parameter extraction).
% 
% Steps:
% - Take the list from pm_data
% - create subdirectory for each patient
% - copy on the the images needed, i.e. 4 MPMs, FLAIR and lesionmask
% 
% In the mean time check that the right number of files are copied from one
% place to the other...
%__________________________________________________________________________
% Copyright (C) 2016 Cyclotron Research Centre

% Written by Christophe Phillips
% University of Liège, Belgium

% Predefined directories
rootdir = '/Users/emilielommers/Desktop/IRM_cyclo';
targetdir = 'SAFEfromCP';
maskdir = fullfile(rootdir,'LESION_MASK');

% Get list of subjects
opt = struct('group',1,'study','transverse');
data = pm_data(opt);
Nd = numel(data);
% Nd = 5

% Expected number/types of files
tFiles = {'FLAIR','lesionmask','MPMs'};
nFiles = [2 2 4]; % .img/hdr for FLAIR & lesionmask, .nii for MPMs
check_copy = zeros(Nd,3);
MPMs_suffix = {'A','MT','R1','R2s'};

% Loop through the patients
for ii=1:Nd
    targetdir_ii = fullfile(rootdir,targetdir,data(ii).id);
    
    % Create target directory
    mkdir(targetdir_ii)
    
    % Find FLAIR and copy
    fnFLAIR = spm_select('FPList',data(ii).FLAIR,'^s.*\.[ihn][mdi][gri]$'); % Select .img/.hdr/.nii
    for jj=1:size(fnFLAIR,1)
        status_copy = copyfile(deblank(fnFLAIR(jj,:)),targetdir_ii);
        if status_copy
            check_copy(ii,1) = check_copy(ii,1)+1;
        end
    end
    
    % Find lesionmask and copy
    mask_filt = sprintf('^%s_lesionmask%s$',data(ii).id2,'\.[ihn][mdi][gri]');
    fnMASK = spm_select('FPList',maskdir,mask_filt); % Select .img/.hdr/.nii
    for jj=1:size(fnMASK,1)
        status_copy = copyfile(deblank(fnMASK(jj,:)),targetdir_ii);
        if status_copy
            check_copy(ii,2) = check_copy(ii,2)+1;
        end
    end    
    
    % Find MPMs and copy
    for jj=1:numel(MPMs_suffix)
        MPMjj_filt = sprintf('^s[0123456789].*_%s%s',MPMs_suffix{jj},'\.nii$');
        % Should start with 's' followed by a number, then the map type and
        % finish with .nii
        fnMPMjj = spm_select('FPList',data(ii).MP,MPMjj_filt);
        status_copy = copyfile(fnMPMjj,targetdir_ii);
        if status_copy
            check_copy(ii,3) = check_copy(ii,3)+1;
        end
    end
end

% Check things went all right
l_prob = cell(1,3);
for ii=1:3 % checking 3 types of files
    l_prob{ii} = find(check_copy(:,ii)~=nFiles(ii));
    if ~isempty(l_prob{ii})
        fprintf('\n Problem copying %s image for subject(s):\n',tFiles{ii});
        for jj = l_prob{ii}'
            fprintf('\t %s\n',data(jj).id);
        end
    end
end

