% script_process_MSchallengeData
%
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
rootDir = 'D:\ccc_DATA\MSchallenge';
% rootDir = 'D:\ccc_DATA\MSchallenge\TESTING\';
dn = {'CHB','UNC'};
fn_results = fullfile(rootDir,'MSchal_USwL_results');
res = struct('mJ',[],'mHd',[],'overlap',[]);
fn_batch = 'batch_segm_MSchallenge_empty.mat';

% Choose what to do: 'Segment','BinAndClean' and/or 'Compare'
% operation = {'Segment','BinAndClean','Compare'};
operation = {'Compare'};

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

%% Getting started
MBempty = load(fullfile(rootDir,fn_batch));
prblm = false(nDir,1);

%% Loop over the subjects
fprintf('\n\n')
for idir = 1:nDir
    try
        subjDir = dname(idir,:);
        fprintf('Processing subject %d out of %d in : %s\n',idir,nDir,subjDir)
        % 1. Collect filenames
        fn_T1 = spm_select('FPList',subjDir,'^.*_T1.nii$');
        fn_T2 = spm_select('FPList',subjDir,'^.*_T2.nii$');
        fn_FLAIR = spm_select('FPList',subjDir,'^rm.*_FLAIR.nii$');
        fn_lesmsk_ref = spm_select('FPList',subjDir,'^.*_lesion.nii$');
        fn_lesmsk_est = spm_select('FPList',subjDir,'^ples.*_FLAIR.nii$');
        % 2. Launch USwithLesion
        if any(strcmp('Segment',operation))
            matlabbatch = MBempty.matlabbatch;
            matlabbatch{1}.spm.tools.USwLtools.uswl.imgMsk = {fn_lesmsk_est};
            matlabbatch{1}.spm.tools.USwLtools.uswl.imgRef = {fn_T1};
            matlabbatch{1}.spm.tools.USwLtools.uswl.imgMPM = cellstr(char(fn_T1, fn_T2));
            matlabbatch{1}.spm.tools.USwLtools.uswl.imgOth = {fn_FLAIR};
            spm_jobman('run', matlabbatch);
        end
        % get segmented images and ICVmsk
        fn_c1234 = spm_select('FPList',subjDir,'^c[1234]k.*_T1.nii$');
        fn_ICVmsk = spm_select('FPList',subjDir,'^icv_k.*_T1.nii$');
        % 3. Binarize and cleanup posterior lesion map
        if any(strcmp('BinAndClean',operation))
            [fn_les_bin,fn_nc] = crc_binarize_segm(fn_c1234,fn_ICVmsk,opt_bin);
            fn_les_clbin = crc_lesion_cleanup(fn_les_bin,clean_k,clean_prefix);
        else
            fn_les_bin = spm_select('FPList',subjDir,'^bc3k.*_T1.nii$');
            fn_nc = spm_select('FPList',subjDir,'^ncc1k.*_T1.nii$');
            fn_les_clbin = spm_select('FPList',subjDir,'^cbc3k.*_T1.nii$');
        end
        % 4. Apply comparison function & collect results
        if any(strcmp('Compare',operation))
            fn_ref = fn_lesmsk_ref;
            fn_test = cellstr(char(fn_lesmsk_est,fn_les_bin,fn_les_clbin));
            opt_ImgOv.mask = fn_ICVmsk;
            for jj=1:numel(fn_test)
                [mJ,mHd,overlap] = image_overlap(fn_test{jj},fn_ref,opt_ImgOv);
                res(idir,jj).mJ = mJ;
                res(idir,jj).mHd = mHd;
                res(idir,jj).overlap = overlap;
            end
        end
    catch
        prblm(idir) = true;
    end
end
% save things!
save(fn_results,'res')

% List of problems
if any(prblm)
    fprintf('Problems in :\n')
    for idir = find(prblm')
        fprintf('\t %s\n',dname(idir,:))
    end
end

%% Collect all results in arrays for further analysis
% methods names
methods_label = {'LST','bin_USwL','cbin_USwL'};
stat_labels = {'modified Jaccard','mean Hausdorff dist','cluster','voxel'};

% list of subjects to consider!
% A few subjects have very small lesions and LST does NOT find them -> poor
% starting estimate, -> USwL cannot help
l_exclude = [2 3 10 11 15 18 19];
l_keep = 1:nDir; l_keep(l_exclude) = [];
nRes = numel(l_keep);

% results
modJaccard = zeros(nRes,3);
mHausdDist = zeros(nRes,3);
cl_stat_tp = zeros(nRes,3);
cl_stat_fp = zeros(nRes,3);
% cl_stat = zeros(nRes,3,2);
% cl_label = {'tp','fp'};
% vx_stat = zeros(nRes,3,6);
% vx_label = {'tp','fp','tn','fn','mcc','CK'};
vx_stat_tp = zeros(nRes,3);
vx_stat_fp = zeros(nRes,3);
vx_stat_tn = zeros(nRes,3);
vx_stat_fn = zeros(nRes,3);
vx_stat_mcc = zeros(nRes,3);
vx_stat_CK = zeros(nRes,3);
confMat = zeros(nRes,3,4);
for ii = 1:nRes
    idir = l_keep(ii);
    for jj=1:3
        modJaccard(ii,jj) = res(idir,jj).mJ;
        mHausdDist(ii,jj) = res(idir,jj).mHd;
%         cl_stat(ii,jj,1) = res(idir,jj).overlap.cluster.tp;
%         cl_stat(ii,jj,2) = res(idir,jj).overlap.cluster.fp;
%         vx_stat(ii,jj,1) = res(idir,jj).overlap.voxel.tp;
%         vx_stat(ii,jj,2) = res(idir,jj).overlap.voxel.fp;
%         vx_stat(ii,jj,3) = res(idir,jj).overlap.voxel.tn;
%         vx_stat(ii,jj,4) = res(idir,jj).overlap.voxel.fn;
%         vx_stat(ii,jj,5) = res(idir,jj).overlap.voxel.mcc;
%         vx_stat(ii,jj,6) = res(idir,jj).overlap.voxel.CK;
        cl_stat_tp(ii,jj) = res(idir,jj).overlap.cluster.tp;
        cl_stat_fp(ii,jj) = res(idir,jj).overlap.cluster.fp;
        vx_stat_tp(ii,jj) = res(idir,jj).overlap.voxel.tp;
        vx_stat_fp(ii,jj) = res(idir,jj).overlap.voxel.fp;
        vx_stat_tn(ii,jj) = res(idir,jj).overlap.voxel.tn;
        vx_stat_fn(ii,jj) = res(idir,jj).overlap.voxel.fn;
        vx_stat_mcc(ii,jj) = res(idir,jj).overlap.voxel.mcc;
        vx_stat_CK(ii,jj) = res(idir,jj).overlap.voxel.CK;
        confMat(ii,jj,:) = res(idir,jj).overlap.voxel.cm(:);
    end
end

%% Display these results
% Jaccard & Hausdorff distance
crc_MSchalResults(cat(3,modJaccard,mHausdDist),char('modified Jaccard','mean Hausdorff dist'));

% Cluster
crc_MSchalResults(cat(3,cl_stat_tp,cl_stat_fp),char('clust_TP','clust_FP'));

% Voxels tp/fp/fn/tn
crc_MSchalResults( ...
    cat(3,vx_stat_tp,vx_stat_fp,vx_stat_fn,vx_stat_tn) , ...
    char('vox.TP','vox.FP','vox.FN','vox.TN'));

% Voxels mcc/CK
crc_MSchalResults(cat(3,vx_stat_mcc,vx_stat_CK),char('vox.Mcc','vox.CK'));

return
%%
idir = 15;
subjDir = dname(idir,:)
fn_images = char(...
    spm_select('FPList',subjDir,'^k.*_T1.nii$'), ...
    spm_select('FPList',subjDir,'^.*_lesion.nii$'), ...
    spm_select('FPList',subjDir,'^ples.*_FLAIR.nii$'), ...
    spm_select('FPList',subjDir,'^c3.*_T1.nii$'), ...
    spm_select('FPList',subjDir,'^bc3.*_T1.nii$'), ...
    spm_select('FPList',subjDir,'^cbc3.*_T1.nii$') ...
    );
spm_check_registration(fn_images)

% end
