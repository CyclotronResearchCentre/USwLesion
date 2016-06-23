%% Batch to process all the MSchallenge data.
%
% This works subject by subject and performs the followign operation to
% each of them:
% 1. find the filenames of all useful images: T1, T2, FLAIR, lesion_mask
%    (reference & estimated, i.e. lesMsk_ref & lesMsk_est)
% 2. launch USwithLesion with T1, T2, FLAIR, lesMsk_est
% 3. binarize (and cleanup) the estimated lesion posterior probability
%    map including some ICV masking, i.e. c3 -> bc3 and cbc3
% 4. with lesMsk_ref as reference for each of them, assess the quality of
%    the estimated lesion volumes of lesMsk_est, bc3 and cbc3.
%
% NOTE
% - the "estimated lesion masks" , as calculated by the LST package [1],
%   have already been produced from the T1 and FLAIR images. This is the
%   1st % approximation for the segmentation, thus is also incuded in the
%   comparison. The goal is to show that USwL improves the solution
%   compared % to this 1st estimate. :-)
% - since we used the LST package, the FLAIR image has been bias-corrected.
%   This one is thus used instead of the original one (rm*_FLAIR.nii)
%
% REF:
% [1] Schmidt et al., Neuroimage, 2012, "An automated tool for detection of
%     FLAIR-hyperintense white-matter lesions in Multiple Sclerosis"
%     http://dx.doi.org/10.1016/j.neuroimage.2011.11.032
%_______________________________________________________________________
% Copyright (C) 2016 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

%% Define & intialize stuff
% rootDir = 'D:\ccc_DATA\MSchallenge';
rootDir = 'D:\ccc_DATA\MSchallenge\TESTING\';
dn = {'CHB','UNC'};
fn_results = fullfile(rootDir,'MSchal_USwL_results');
res = struct('mJ',[],'mHd',[],'overlap',[]);
fn_batch = 'batch_segm_MSchallenge_empty.mat';

% Set option for binarizing
opt_bin = struct(... % DBinarization option
    'bin_classic', true, ...
    'bin_ML', true, ...
    'bin_mltr', false, ...
    'singleImg', false, ...
    'ListOut', 3, ...
    'thr', .2);
% Set option for cleaning
clean_k = 64; % must be at leat 8mm^3, i.e. 4x4x4 voxels
clean_prefix = 'c';
% Set option for image overlap
opt_ImgOv.mask = ''; % specified on the fly

%% build subject directories list, based on dn
dname = char(0);
for ii=1:numel(dn)
    [dirs] = spm_select('FPList',rootDir,'dir',['^',dn{ii},'.*']);
    if ~isempty(dirs)
        dname = char(dname,dirs);
    end
end
dname(1,:) = [];
nDir = size(dname,1);
res(nDir,3) = res;

%%
MBempty = load(fullfile(rootDir,fn_batch));
prblm = false(nDir,1);
%% Loop over the subjects
for idir = 1:nDir
    try
        subjDir = dname(idir,:);
        % 1. Collect filenames
        fn_T1 = spm_select('FPList',subjDir,'^.*_T1.nii$');
        fn_T2 = spm_select('FPList',subjDir,'^.*_T2.nii$');
        fn_FLAIR = spm_select('FPList',subjDir,'^rm.*_FLAIR.nii$');
        fn_lesmsk_ref = spm_select('FPList',subjDir,'^.*_lesion.nii$');
        fn_lesmsk_est = spm_select('FPList',subjDir,'^ples.*_FLAIR.nii$');
        % 2. Launch USwithLesion
        matlabbatch = MBempty.matlabbatch;
        matlabbatch{1}.spm.tools.USwLtools.uswl.imgMsk = {fn_lesmsk_est};
        matlabbatch{1}.spm.tools.USwLtools.uswl.imgRef = {fn_T1};
        matlabbatch{1}.spm.tools.USwLtools.uswl.imgMPM = cellstr(char(fn_T1, fn_T2));
        matlabbatch{1}.spm.tools.USwLtools.uswl.imgOth = {fn_FLAIR};
        spm_jobman('run', matlabbatch);
        % get segmented images and ICVmsk
        fn_c1234 = spm_select('FPList',subjDir,'^c[1234]k.*_T1.nii$');
        fn_ICVmsk = spm_select('FPList',subjDir,'^icv_k.*_T1.nii$');
        % 3. Binarize and cleanup posterior lesion map
        [fn_les_bin,fn_nc] = crc_binarize_segm(fn_c1234,fn_ICVmsk,opt_bin);
        fn_les_clbin = crc_lesion_cleanup(fn_les_bin,clean_k,clean_prefix);
        % 4. Apply comparison function & collect results
        fn_ref = fn_lesmsk_ref;
        fn_test = cellstr(char(fn_lesmsk_est,fn_les_bin,fn_les_clbin));
        opt_ImgOv.mask = fn_ICVmsk;
        for jj=1:numel(fn_test)
            [mJ,mHd,overlap] = image_overlap(fn_ref,fn_test{jj},opt_ImgOv);
            res(idir,jj).mJ = mJ;
            res(idir,jj).mHd = mHd;
            res(idir,jj).overlap = overlap;
        end
    catch
        prblm(idir) = true;
    end
end
% save things!
save(fn_results,'res')

%% Collect all results in arrays for further analysis
% methods names
methods_label = {'LST','bin_USwL','cbin_USwL'};
stat_labels = {'modified Jaccard','mean Hausdorff dist','cluster','voxel'};
% results
modJaccard = zeros(nDir,3);
mHausdDist = zeros(nDir,3);
cl_stat = zeros(nDir,3,2);
cl_label = {'tp','fp'};
vx_stat = zeros(nDir,3,6);
vx_label = {'tp','fp','tn','fn','mcc','CK'};
confMat = zeros(nDir,3,4);
for idir = 1:nDir
    for jj=1:3
        modJaccard(idir,jj) = res(idir,jj).mJ;
        mHausdDist(idir,jj) = res(idir,jj).mHd;
        cl_stat(idir,jj,1) = res(idir,jj).overlap.cluster.tp;
        cl_stat(idir,jj,2) = res(idir,jj).overlap.cluster.fp;
        vx_stat(idir,jj,1) = res(idir,jj).overlap.voxel.tp;
        vx_stat(idir,jj,2) = res(idir,jj).overlap.voxel.fp;
        vx_stat(idir,jj,3) = res(idir,jj).overlap.voxel.tn;
        vx_stat(idir,jj,4) = res(idir,jj).overlap.voxel.fn;
        vx_stat(idir,jj,5) = res(idir,jj).overlap.voxel.mcc;
        vx_stat(idir,jj,6) = res(idir,jj).overlap.voxel.CK;
        confMat(idir,jj,:) = res(idir,jj).overlap.voxel.cm(:);
    end
end

%% Display these results
% Jaccard
ms_modJaccard = [mean(modJaccard) ; std(modJaccard)];
figure, hold on
bar(1:3,ms_modJaccard(1,:))
errorbar(1:3,ms_modJaccard(1,:),ms_modJaccard(2,:),'.')
title(stat_labels{1})
% Hausdorff
ms_mHausdDist = [mean(mHausdDist) ; std(mHausdDist)];
figure, hold on
bar(1:3,ms_mHausdDist(1,:))
errorbar(1:3,ms_mHausdDist(1,:),ms_mHausdDist(2,:),'.')
title(stat_labels{2})
% Cluster
ms_cluster = [mean(cl_stat) ; std(cl_stat)];
figure,
subplot(2,1,1), hold on
bar(1:3,ms_cluster(:,1))
errorbar(1:3,ms_mHausdDist(:,1,1),ms_mHausdDist(:,2,1),'.')
subplot(2,1,2), hold on
bar(1:3,ms_cluster(:,2))
errorbar(1:3,ms_mHausdDist(:,2,1),ms_mHausdDist(:,2,2),'.')


% end
