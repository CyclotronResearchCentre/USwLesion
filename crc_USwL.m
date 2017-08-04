function fn_out = crc_USwL(fn_in,options)
% 
% Function doing all the work of "Unified segmentation with lesion".
% 
% INPUT
% - fn_in   : cell array (1x4) of input filenames
%       {1} : lesion mask image
%       {2} : reference structural image, used for the 1st segmentation
%       {3} : multi-parametric maps 
%       {4} : other structural images
% - options :
%       imgTpm  : tissue probability maps (just one 4D file)
%       img4US  : flag -> images to use for the 2nd segmentation 
%                   (0, ref only; 1, MPMs only; 2, MPMs+others)
%                   MPMs
%       biasreg : bias regularisation value
%       biasfwhm: bias FWHM value (Inf = no bias correction) 
%       biaswr  : flag -> save bias corrected and/or bias field ([0/1 0/1])
%       NbGaussian : number of Gaussians per tissue class, incl. lesion.
%       tpm4lesion : flag -> TPM(s) affected by the lesion
%                      (0, GM; 1, WM; 2, GM+WM; 3, GM+WM+CSF) 
%       thrMPM  : preliminary thresholding the MPMs, 0/1
%       ICVmsk  : mask the MPMs & Other images by created ICV-mask, 0/1
%       mrf     : MRF parameter
%       cleanup : Standard post-US clean-up, should most likely not be 
%                 applied (because there is an extra lesion-class) 
%                 -> better set it to 0
%       thrLesion: lesion mask extent thresholding 
%                   (0, no thresh; k, keep size>k; Inf, keep only largest)
% 
% OUTPUT
% - fn_out :
%       wstruct     : warped structural reference
%       fxMPM_i    : 'fixed' MPMs, i = 1, 2,... #MPM
%       fxMPMmsk_i : map of fixed voxels for i^th MPM
%       ICVmsk      : intra-cranial volume mask, as generated from str-ref
%       kMPM_i      : masked i^th MPM
%       kOth_i      : masked i^th Other image
%       wMPM_i      : warped (masked) i^th MPM
%       wOth_i      : warped (masked) i^th other image
%       TPMl        : subject specific TPM with lesion
%       segmImg     : structure with posterior tissue probabilities
%           c(i)    : class #i in subject space
%           wc(i)   : class #i in MNI space
%           mwc(i)  : modulated class #i in MNI space
% 
% OPERATIONS
% Here are the main steps:
%   0. Clean up of the MPM images, based on each maps value range
%   1. "Trim 'n grow" the mask image : -> t_Msk / dt_Msk
%       - remov the "small" MS patches using a simple criteria: number of
%         voxels in patch must be > minNR    -> t_Msk
%       - then grow volume by 1 voxel -> dt_Msk
%   2. Apply the mask on the reference structural images -> k_sRef
%   3. Segment the masked structural (k_sRef), normalize the cleaned up
%      mask (t_Msk) and smooth it -> new TPM for the lesion.
%   4. Update the TPMs to include a 7th tissue class -> TPMms
%   Note that the lesion is inserted in *3rd position*, between WM and CSF!
%   5. Do the segmentation with the new TPM_ms
%       img4US = 0 -> Structural reference only
%              = 1 -> all MPMs
%              = 2 -> all MPMs + others
%   6. Apply the deformation onto the MPMs -> warped MPMs
%   7. Collect all the image filenames created
%
% Check the Readme file for further processing details.
%_______________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium
% Cyril Pernet updated few bits to work with no MPM images + added structural
% normalization and N Gaussians - Edinburgh Imaging, The University of Edinburgh

%% Input data
% fn_in{1} = Mask image
% fn_in{2} = structural reference
% fn_in{3} = All MPM's
% fn_in{4} = Other images
if isempty(fn_in{3})
    nMPM = 0;
else
    nMPM = size(fn_in{3},1);
end
if isempty(fn_in{4})
    nOth = 0;
else
    nOth = size(fn_in{4},1);
end
pth = spm_file(fn_in{2},'path');

% By default no cleanup after US-with-lesion segmentation!
if ~isfield(options,'cleanup') || ...
        isempty(options.cleanup) || strcmp(options.cleanup,'<UNDEFINED>')
    options.cleanup = 0;
end

% Check #Gaussians and #TPMs for the USwLesion segmentation.
NbGaussian = options.NbGaussian;
fn_tpm_USwL = options.imgTpm{1};
fn_tpm_USwL = spm_file(fn_tpm_USwL,'number','');
Vtpm_USwL = spm_vol(fn_tpm_USwL);
if numel(NbGaussian)~=(numel(Vtpm_USwL)+1)
    error('There are %d tpm (incl. lesion) but only %d #Gaussians provided', ...
        numel(Vtpm_USwL)+1,numel(NbGaussian));
end

%% Define processing parameters for the creation of the updated TPM
opt = crc_USwL_get_defaults('uTPM');

%% 0. Clean up of the MPM images!
% Need to know the order of the images, ideally MT, A, R1, R2s and should
% check with their filename? based on '_MT', '_A', '_R1', '_R2s'?
fn_in_3_orig = fn_in{3};
if options.thrMPM && nMPM ~= 0
    % Fix MPMs
    fn_in{3} = crc_fix_MPMintens(fn_in{3});
    
    % Check if reference structural is one of the MPMs 
    % -> if match, use the fixed version!
    is_StrucRef_MPM = strcmp(fn_in{2}, cellstr(fn_in_3_orig));
    
    % Update the filenames with fixed MPMs
    if any(is_StrucRef_MPM)
        fn_in{2} = fn_in{3}(find(is_StrucRef_MPM),:);
    end
end

%% 1. "Trim 'n grow" the mask image : -> t_Msk / dt_Msk
% - remov the "small" lesion patches using a simple criteria: volume of
%   lesion patch must be > minVol -> creates on the drive t_Msk
% - then grow volume by nDilate voxel(s) -> creates on the drive dt_Msk
[fn_tMsk,fn_dtMsk] = mask_trimNgrow(fn_in{1},opt.minVol,opt.nDilate);

%% 2. Apply the lesion mask on the reference structural images -> k_sRef
fn_kMTw = spm_file(fn_in{2},'prefix','k');
Vi(1) = spm_vol(fn_in{2});
Vi(2) = spm_vol(fn_dtMsk);
Vo = Vi(1);
Vo.fname = fn_kMTw;
Vo = spm_imcalc(Vi,Vo,'i1.*(((i2>.5)-1)./((i2>.5)-1))');

%% 3. Segment the masked structural (k_sRef), normalize the cleaned up mask
% (t_Msk) and smooth it -> new TPM for the lesion.
% Then create an ICV mask (for MPM's ICV masking)
% Here use the standard TPMs, as this is a reference structural image
% masked for the lesion, i.e. "cost function masking" approach.
clear matlabbatch
fn_tpm_msksegm = crc_USwL_get_defaults('msksegm.imgTpm');
matlabbatch = batch_normalize_smooth( ...
    fn_kMTw, ... % masked structural image used for the warping estimation
    fn_tMsk, ... % cleaned up lesion mask to be warped into MNI
    fn_tpm_msksegm{1}, ... % filename of tissue probability map
    opt.smoKern); % smoothing applied on the normalized lesion mask -> new prior
spm_jobman('run', matlabbatch);
fn_swtMsk = spm_file(fn_tMsk,'prefix','sw'); % smooth normalized lesion mask
fn_wtMsk = spm_file(fn_tMsk,'prefix','w'); %#ok<*NASGU> % normalized lesion mask

% Build ICV, with some cleaning up
fn_TCin = char( ...
    spm_file(fn_kMTw,'prefix','c1'),... % GM
    spm_file(fn_kMTw,'prefix','c2'),... % WM
    spm_file(fn_kMTw,'prefix','c3'),... % CSF
    fn_dtMsk); % dilated lesion mask, as masked struct-ref
opt_ICV = struct( ...
    'fn_ref', fn_kMTw, ...
    'fn_warp', spm_file(fn_kMTw,'prefix','y_'), ...
    'fn_iwarp', spm_file(fn_kMTw,'prefix','iy_'), ...
    'smoK', 4);
fn_icv_out = crc_build_ICVmsk(fn_TCin,opt_ICV);
fn_ICV = deblank(fn_icv_out(1,:));
fn_swICV = deblank(fn_icv_out(3,:));

% Delete unnecesary files
to_delete = char( fn_TCin(1:3,:), opt_ICV.fn_warp, opt_ICV.fn_iwarp, ...
    spm_file(fn_kMTw,'suffix','_seg8','ext','mat') );
for ii=1:size(to_delete,1), delete(deblank(to_delete(ii,:))); end

% Apply the mask -> this possibly overwrites the masked struct reference
if options.ICVmsk && nMPM ~= 0 % ICV-mask the MPMs & others
    fn_tmp = [];
    for ii=1:nMPM
        fn_MPM_ii = deblank(fn_in{3}(ii,:));
        fn_tmp = char(fn_tmp,mask_img(fn_MPM_ii,fn_ICV,'k'));
    end
    fn_in{3} = fn_tmp(2:end,:);
    % Mask other images too!
    fn_tmp = [];
    for ii=1:nOth
        fn_Oth_ii = deblank(fn_in{4}(ii,:));
        fn_tmp = char(fn_tmp,mask_img(fn_Oth_ii,fn_ICV,'k'));
    end
    fn_in{4} = fn_tmp(2:end,:);
end

%% 4. Update the TPMs to include an extra tissue class -> TPMms
% Note that the lesion is inserted in *3rd position*, between WM and CSF!
opt_tpm = struct(...
    'tpm4lesion', options.tpm4lesion, ... % tissues to be modified for lesion (0/1/2/3) for GM / WM / GM+WM / GM+WM+CSF
    'fn_tpm', fn_tpm_USwL, ... % tpm file name
    'tpm_ratio', opt.tpm_ratio, ... % ratio between healthy and lesion tissue
    'min_tpm_icv', opt.min_tpm_icv, ... % minimum value in intracranial volume
    'min_tpm', opt.min_tpm); % minum value overall
if options.ICVmsk && nMPM ~= 0 % ICV-mask the TPMs
    opt_tpm.fn_swICV = fn_swICV; % smoothed-warped ICV mask to apply on TPMs
end
fn_TPMl = update_TPM_with_lesion(opt_tpm, fn_swtMsk); % that creates the new tissue class tmp

%% 5. Do the segmentation with the new TPM_ms
% img4US = 0 -> Structural reference only
%        = 1 -> all MPMs
%        = 2 -> all MPMs + others

% the segmentation options are changed for the nb of gaussians + because we
% have now lesions we increase the clean up (Markov = 2) but set cleanup to
% none

switch options.img4US
    case 0
        fn_Img2segm = fn_in{2}; %#ok<*CCAT1>
    case 1
        %         if isempty(fn_in{3}) % if no MPM
        %             if ~nOth % and no others
        %                 fn_Img2segm = char(fn_in{2}); % use struct
        %             else
        %                 fn_Img2segm = char(fn_in{2}, fn_in{4}); % otherwise use struct and others
        %             end
        %         else
        fn_Img2segm = fn_in{3}; % else as requested use MPM only
        %         end
    case 2
        fn_Img2segm = char(fn_in{3} , fn_in{4});
        %         if isempty(fn_in{3}) % if no MPM
        %             if ~nOth % and no others
        %                 fn_Img2segm = char(fn_in{2}); % use struct
        %             else
        %                 fn_Img2segm = char(fn_in{2}, fn_in{4}); % otherwise use struct and others
        %             end
        %         else
        %             fn_Img2segm = char(fn_in{3}, fn_in{4}); % else as requested use MPM, and all others
        %         end
end

opt_segm = struct( ...
    'b_param', [options.biasreg options.biasfwhm], ...
    'b_write', options.biaswr, ...
    'nGauss', NbGaussian, ...
    'mrf', options.mrf, ...
    'cleanup', options.cleanup);

clear matlabbatch
[matlabbatch] = batch_segment_l(fn_Img2segm, fn_TPMl, opt_segm);
spm_jobman('run', matlabbatch);

fn_Cimg   = spm_select('FPList',pth, ...  % native space
    '^c[0-9].*',spm_file(fn_Img2segm(1,:),'basename'),'\.nii$'); 
fn_rCimg  = spm_select('FPList',pth, ...  % native dartel imported
    '^rc[0-9].*',spm_file(fn_Img2segm(1,:),'basename'),'\.nii$');
fn_wCimg  = spm_select('FPList',pth, ...  % warped
    '^wc[0-9].*',spm_file(fn_Img2segm(1,:),'basename'),'\.nii$');
fn_mwCimg = spm_select('FPList',pth, ...  % modulated warped
    '^mwc[0-9].*',spm_file(fn_Img2segm(1,:),'basename'),'\.nii$');
% When using BG, i.e. extended TPM, then recombine GM with BG
%  -> add c8 onto c1 -> only 1 image (c1) with GM + c8 with BG.
if numel(NbGaussian)==8
    add_2_images(fn_Cimg([1 end],:),  fn_Cimg(1,:), 2); % uint8
    add_2_images(fn_rCimg([1 end],:), fn_rCimg(1,:), 16); % float32
    add_2_images(fn_wCimg([1 end],:), fn_wCimg(1,:), 2); % uint8
    add_2_images(fn_mwCimg([1 end],:), fn_mwCimg(1,:), 16); % float32
end

% Rebuild ICV mask from latest segmentation, with some cleaning up
fn_TCin = fn_Cimg(1:4,:);  % GM, WM, Lesion, CSF
opt_ICV = struct( ...
    'fn_ref', spm_file(fn_in{3}(1,:)), ...
    'fn_warp', spm_file(fn_in{3}(1,:),'prefix','y_'), ...
    'smoK', 0);
fn_icv_out = crc_build_ICVmsk(fn_TCin,opt_ICV);
fn_ICV = deblank(fn_icv_out(1,:));
fn_wICV = deblank(fn_icv_out(2,:));

% Apply mask on GM, WM, Lesion and CSF (+BG if there)
for ii=1:4
    mask_img(fn_Cimg(ii,:),fn_ICV,'');
    mask_img(fn_rCimg(ii,:),fn_ICV,'');
    mask_img(fn_wCimg(ii,:),fn_wICV,'');
    mask_img(fn_mwCimg(ii,:),fn_wICV,'');
end

if numel(NbGaussian)==8
    mask_img(fn_Cimg(end,:),fn_ICV,'');
    mask_img(fn_rCimg(end,:),fn_ICV,'');
    mask_img(fn_wCimg(end,:),fn_wICV,'');
    mask_img(fn_mwCimg(end,:),fn_wICV,'');
end

% Apply the final ICV mask on MPM & Others
if options.ICVmsk && nMPM ~= 0 % ICV-mask the MPMs & others
    fn_tmp = [];
    for ii=1:nMPM
        fn_MPM_ii = deblank(fn_in{3}(ii,:));
        fn_tmp = char(fn_tmp,mask_img(fn_MPM_ii,fn_ICV,'k'));
    end
    fn_in{3} = fn_tmp(2:end,:);
    % Mask other images too!
    fn_tmp = [];
    for ii=1:nOth
        fn_Oth_ii = deblank(fn_in{4}(ii,:));
        fn_tmp = char(fn_tmp,mask_img(fn_Oth_ii,fn_ICV,'k'));
    end
    fn_in{4} = fn_tmp(2:end,:);
end

%% 6. Apply the deformation onto the MPMs/Other -> warped MPMs/Other
% Re-create an ICV mask

fn_warp = spm_file(fn_Img2segm(1,:),'prefix','y_');
% Apply on all images: struct ref + MPM + others, if available

fn_img2warp = {fn_in{2}}; % Always include struct-ref to be warped
% Note that the warps are applied on the possible thresholded-masked
% MPMs/Other, depending on options chosen.
if ~isempty(fn_in{3})
    fn_img2warp = {char(fn_img2warp{1}, fn_in{3})};
end
if ~isempty(fn_in{4})
    fn_img2warp = {char(fn_img2warp{1} , fn_in{4})};
end
clear matlabbatch

if options.biaswr(2) % bias image corrected image screated
    fn_img2warp = {spm_file(fn_img2warp{1},'prefix','m')};
end

[matlabbatch] = batch_normalize_MPM(fn_img2warp,fn_warp);
spm_jobman('run', matlabbatch);

fn_warped_struct = spm_file(fn_in{2},'prefix','w');
if ~isempty(fn_in{3})
    fn_warped_MPM = spm_file(fn_in{3},'prefix','w');
else
    fn_warped_MPM = '';
end
if ~isempty(fn_in{4})
    fn_warped_Oth = spm_file(fn_in{4},'prefix','w');
else
    fn_warped_Oth = '';
end
% fn_mwTC = char( ...
%     spm_file(fn_in{3}(1,:),'prefix','smwc1'), ...
%     spm_file(fn_in{3}(1,:),'prefix','smwc2'), ...
%     spm_file(fn_in{3}(1,:),'prefix','smwc3') ); %#ok<*NASGU>

%% 7. Collect all the image filenames created

% There must always be a struct-ref -> warped one
fn_out.wstruct = {deblank(fn_warped_struct)};

if options.thrMPM && nMPM ~= 0
    for ii=1:nMPM
        fn_out.(sprintf('fxMPM_%d',ii)) = ...
            {spm_file(deblank(fn_in_3_orig(ii,:)),'prefix','t')};
        fn_out.(sprintf('fxMPMmsk_%d',ii)) = ...
            {spm_file(fn_in_3_orig(ii,:),'prefix','fx_')};
    end
end

fn_out.ICVmsk = {fn_ICV};

if options.ICVmsk && nMPM ~= 0;
    for ii=1:nMPM
        fn_out.(sprintf('kMPM_%d',ii)) = {deblank(fn_in{3}(ii,:))};
    end
    for ii=1:nOth
        fn_out.(sprintf('kOth_%d',ii)) = {deblank(fn_in{4}(ii,:))};
    end
end

if ~isempty(fn_warped_MPM) % warped MPMs
    for ii=1:size(fn_warped_MPM,1)
        fn_out.(sprintf('wMPM_%d',ii)) = {deblank(fn_warped_MPM(ii,:))};
    end
end

if ~isempty(fn_warped_Oth)
    for ii=1:size(fn_warped_Oth,1) % warped Others
        fn_out.(sprintf('wOth_%d',ii)) = {deblank(fn_warped_Oth(ii,:))};
    end
end

% subject specific TPM with lesion
fn_out.TPMl = {fn_TPMl};

% segmented tissues
fn_out.segmImg.c1 = {deblank(fn_Cimg(1,:))}; % GM
fn_out.segmImg.c2 = {deblank(fn_Cimg(2,:))}; % WM
fn_out.segmImg.c3 = {deblank(fn_Cimg(3,:))}; % Lesion
fn_out.segmImg.c4 = {deblank(fn_Cimg(4,:))}; % CSF
% warped segmented tissues
fn_out.segmImg.wc1 = {deblank(fn_wCimg(1,:))}; % warped GM
fn_out.segmImg.wc2 = {deblank(fn_wCimg(2,:))}; % warped WM
fn_out.segmImg.wc3 = {deblank(fn_wCimg(3,:))}; % warped Lesion
fn_out.segmImg.wc4 = {deblank(fn_wCimg(4,:))}; % warped CSF
% modulated warped segmented tissues
fn_out.segmImg.mwc1 = {deblank(fn_mwCimg(1,:))}; % modulated warped GM
fn_out.segmImg.mwc2 = {deblank(fn_mwCimg(2,:))}; % modulated warped WM
fn_out.segmImg.mwc3 = {deblank(fn_mwCimg(3,:))}; % modulated warped Lesion
fn_out.segmImg.mwc4 = {deblank(fn_mwCimg(4,:))}; % modulated warped CSF

if options.thrLesion ~= 0
    crc_lesion_cleanup(fn_out.segmImg.c3,options.thrLesion);
end

end

% =======================================================================
%% SUBFUNCTIONS
% =======================================================================

%% STEP 1: 
% Removing small lesion patches from mask
function [fn_tMsk,fn_dtMsk] = mask_trimNgrow(P_in,minVol,nDilate)
% 1) Trim a mask image by removing bits that would be too small to really
%    matter according to medical criteria
%   For example cf. E. Lommers and MS patients:
%    "Lesions will ordinarily be larger than 3 mm in cross section"
%    A cube of 2x2x2 mm^3 has a diagonal of sqrt(12)~3.4mm and
%    and a volume of 8 mm^3 -> minVol = 8 [DEF]
%   -> fn_tMsk used for the new TPM_ms
% 2) Then grow the volume by nDilate voxels [2, DEF]
%   -> fn_dtMsk used for the masking for the 1st warping

if nargin<3
    nDilate = 2;
end
if nargin<2
    minNr = 8;
end

% 1) Load things
V = spm_vol(P_in);
[Msk,XYZ] = spm_read_vols(V);

% Ensures values are [0 1], in case scaling was wrong, e.g. [0 255], or
% there are some tiny negative values, e.g. if mask was resampled
if max(Msk(:))>1 || min(Msk(:))<0
    fprintf('WARNING: some bad values in the lesion mask!\n')
    fprintf('\tValue range : [%1.4f %1.4f] -> Setting it to [0 1]\n', ...
        min(Msk(:)), max(Msk(:)))
    Msk(Msk>1e-6) = 1; % non-zero values, as in >1e-6, set to 1
    Msk(Msk<0) = 0; % anything below zero set to zero.
end

XYZvx = V.mat\[XYZ ; ones(1,size(XYZ,2))];
vx_vol = abs(det(V.mat(1:3,1:3)));
minNr = minVol/vx_vol;

% 2) Clean up
lMsk  = find(Msk(:)>0);
lXYZvx = XYZvx(1:3,(lMsk));
vMsk = Msk(lMsk);
[Ncl,Zcl,Mcl,Acl,XYZcl] = spm_max(vMsk,lXYZvx); %#ok<*NASGU,*ASGLU>

nrA = max(Acl);
l_all = ones(length(lMsk),1);
n_rem = 0;
for ii=1:nrA
    % deal with regions, one by one
    l_ii = find(Acl == ii);
    n_ii = length(l_ii);
    if n_ii<minNr
        l_all(l_ii) = 0;
        n_rem = n_rem+1;
    end
end
Msk_nM = Msk;
l_rem = find(l_all==0);
Msk_nM(lMsk(l_rem)) = 0; %#ok<*FNDSB>

% 4) Save 1st image fn_tMsk
V_nM = V;
V_nM.fname = spm_file(V.fname,'prefix','t');
V_nM = spm_create_vol(V_nM);
V_nM = spm_write_vol(V_nM,Msk_nM);
fn_tMsk = V_nM.fname;

% 5) dilate mask
if nDilate
    dMsk_nM = imdilate(~~Msk_nM,ones(3,3,3));
end
if nDilate>1
    for ii=1:nDilate-1
        dMsk_nM = imdilate(dMsk_nM,ones(3,3,3));
    end
end

% 6) Save 2nd image fn_dtMsk
V_nM = V;
V_nM.fname = spm_file(V.fname,'prefix','dt');
V_nM = spm_create_vol(V_nM);
V_nM = spm_write_vol(V_nM,dMsk_nM);

fn_dtMsk = V_nM.fname;

end

%% STEP 3: 
% Creating the normalization batch to 
%   + normalize the masked structural image
%   + generate an ICV mask

function matlabbatch = batch_normalize_smooth(fn_kRef,fn_tMsk,fn_TPM,smoKern)
% matlabbatch = batch_normalize_smooth(fn_kRef,fn_tMsk,fn_TPM,smoKern)
% This batch includes:
% [1,2] defining inputs
% [3] segmentation of the masked structural
% [4,5] writing out + smoothing the normalized lesion mask
%
% INPUT:
% - fn_kRef : masked structural image used for the warping estimation
% - fn_tMsk : cleaned up lesion mask to be warped into MNI
% - fn_TPM  : filename of tissue probability map
% - smoKern : smoothing applied on the normalized lesion mask -> new prior
%
% OUTPUT:
% - matlabbatch : operation batch
% - fn_ICV : file name to ICV mask created

pth_img = spm_file(fn_tMsk,'path');
% fn_ICV = spm_file(spm_file(fn_kRef,'prefix','icv_'),'filename');

% get the defaults
segm_def = crc_USwL_get_defaults('msksegm');

% Build the batch structure
matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'LesionMask';
matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {{fn_tMsk}};
matlabbatch{2}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'MaskedRefStruct';
matlabbatch{2}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {{fn_kRef}};
matlabbatch{3}.spm.spatial.preproc.channel.vols(1) = cfg_dep('Named File Selector: MaskedRefStruct(1) - Files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
matlabbatch{3}.spm.spatial.preproc.channel.biasreg = segm_def.biasreg;
matlabbatch{3}.spm.spatial.preproc.channel.biasfwhm = segm_def.biasfwhm;
matlabbatch{3}.spm.spatial.preproc.channel.write = segm_def.biaswr;
for ii=1:6
    matlabbatch{3}.spm.spatial.preproc.tissue(ii).tpm = {spm_file(fn_TPM,'number',ii)};
    matlabbatch{3}.spm.spatial.preproc.tissue(ii).ngaus = segm_def.NbGaussian(ii);
    matlabbatch{3}.spm.spatial.preproc.tissue(ii).native = segm_def.native(ii,:);
    matlabbatch{3}.spm.spatial.preproc.tissue(ii).warped = [0 0];
end
matlabbatch{3}.spm.spatial.preproc.warp.mrf = segm_def.mrf;
matlabbatch{3}.spm.spatial.preproc.warp.cleanup = segm_def.cleanup; %% the cleanup is ad-hoc by default leave 1
matlabbatch{3}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{3}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{3}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{3}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{3}.spm.spatial.preproc.warp.write = [1 1];
matlabbatch{4}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Named File Selector: LesionMask(1) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
matlabbatch{4}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72 ; 90 90 108];
matlabbatch{4}.spm.spatial.normalise.write.woptions.vox = [1.5 1.5 1.5];
matlabbatch{4}.spm.spatial.normalise.write.woptions.interp = 1; % Only trilinear interpolation
matlabbatch{5}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{5}.spm.spatial.smooth.fwhm = smoKern*[1 1 1];
matlabbatch{5}.spm.spatial.smooth.dtype = 16;
matlabbatch{5}.spm.spatial.smooth.im = 0;
matlabbatch{5}.spm.spatial.smooth.prefix = 's';

end

%% STEP 4: 
% Updating the TPM with an extra class, the lesion
% Note that in the resulting TPM the lesion is inserted in *3rd position*,
% i.e. between WM and CSF!
function fn_TPMl = update_TPM_with_lesion(opt, fn_swtMsk)
% fn_TPMl = update_TPM_with_lesion(opt, fn_swtMsk)
%
% INPUT
% - opt : structure with a few parameters
%     .tpm4lesion : tissues to be modified for lesion (0/1/2/3) for
%                       GM / WM / GM+WM / GM+WM+CSF
%     .fn_tpm : tpm file name
%     .tpm_ratio : ratio between healthy and lesion tissue
%     .min_tpm_icv : minimum value in intracranial volume
%     .min_tpm : minum value overall
%     .fn_swICV : [optional] smoothed-warped ICV mask to apply on TPMs
% - fn_swtMsk : filename of smoothed normalized cleaned lesion mask, to be
%               used to create the lesion tissue class

% 0) select TPM and load
[pth,fnam,ext,num] = spm_fileparts(opt.fn_tpm);
fn_TPM = fullfile(pth,[fnam,ext]); % ensuring I load all TPMs together.
Vtpm     = spm_vol(fn_TPM);
tpm_orig = spm_read_vols(Vtpm);
Ntpm_o = numel(Vtpm);
if Ntpm_o==6
    % Standard TPM
    tpm_std = true;
elseif Ntpm_o==7
    % Special TPM for MPM (with BG)
    tpm_std = false;
else
    error('Wrong number of TPM''s.');
end
if tpm_std
    tpm_GM   = squeeze(tpm_orig(:,:,:,1));
    tpm_WM   = squeeze(tpm_orig(:,:,:,2)); % used later on to define ICV
    tpm_CSF  = squeeze(tpm_orig(:,:,:,3));
else
    tpm_GM   = squeeze(sum(tpm_orig(:,:,:,[1 7]),4));
    alpha    = squeeze(tpm_orig(:,:,:,1)./tpm_orig(:,:,:,7)); % ratio GM/BG
    tpm_WM   = squeeze(tpm_orig(:,:,:,2)); % used later on to define ICV
    tpm_CSF  = squeeze(tpm_orig(:,:,:,3));
end

% Read in the healthy tissue prob map
% + define some masks for where lesion could also be and ICV area
msk_ICV = tpm_WM>=opt.min_tpm_icv;
switch opt.tpm4lesion
    case 0 % GM only
        tpm_healthy = tpm_GM;
        msk_les_possible = (tpm_GM>=tpm_WM) & (tpm_GM>=tpm_CSF);
        
    case 1 % WM only
        tpm_healthy = tpm_WM;
        msk_les_possible = (tpm_WM>=tpm_GM) & (tpm_WM>=tpm_CSF);
    case 2 % WM+GM
        tpm_healthy = tpm_GM+tpm_WM;
        msk_les_possible = tpm_WM+tpm_GM>=tpm_CSF;
    case 3 % WM+GM+CSF
        tpm_healthy = tpm_GM+tpm_WM+tpm_CSF;
        msk_les_possible = tpm_healthy>=opt.min_tpm_icv;
    otherwise
        error('Wrong tissue flag');
end
% load smooth lesion mask = tentative lesion tpm
Vl = spm_vol(fn_swtMsk);
tpm_l = spm_read_vols(Vl);

% Ensures values are [0 1], in case scaling was wrong, e.g. [0 255], or
% there are some tiny negative values, e.g. if mask was resampled
if max(tpm_l(:))>1 || min(tpm_l(:))<0
    fprintf('WARNING: some bad values in the lesion prior map!\n')
    fprintf('\tValue range : [%1.4f %1.4f] -> Setting it to [0 1]\n', ...
        min(tpm_l(:)), max(tpm_l(:)))
    tpm_l = tpm_l/max(tpm_l(:)); % rescale to 1
    tpm_l(tpm_l<0) = 0; % anything below zero set to zero.
end

% define where lesion could also be as smoothed version of msk_les_possible
prob_l_possible = uint8(zeros(size(msk_les_possible)));
fwhm = 4./sqrt(sum(Vtpm(1).mat(1:3,1:3).^2)); % 4mm expressed in voxels
spm_smooth(uint8(msk_les_possible*255),prob_l_possible,fwhm);
prob_l_possible = double(prob_l_possible)/255;

% 1) scale lesion tpm and adjust with healthy tissue prob map in ICV
% 2) ensure minium value all over
% 3) concatenate by setting lesion at last position & adjust 'other' class
tpm_Lu = (1-1/opt.tpm_ratio)*tpm_l.*tpm_healthy; % update lesion tpm
l_les_possible = find( (prob_l_possible(:)*opt.min_tpm_icv>tpm_Lu(:)) & (prob_l_possible(:)>0) );
tpm_Lu(l_les_possible) = prob_l_possible(l_les_possible)*opt.min_tpm_icv; % min_tpm_icv x prob of possible lesion
tpm_Lu(tpm_Lu<opt.min_tpm) = opt.min_tpm; % at least min_tpm everywhere
tpm_ext = cat(4,tpm_orig,tpm_Lu); % put lesion at the end

switch opt.tpm4lesion % update healthy tissue classes
    case 0 % GM only
        tpm_GMu = tpm_healthy - tpm_Lu;
        % equiv. to tpm_GMu = tpm_GM .* (1 - (1-1/opt.tpm_ratio) * tpm_l);
        if tpm_std
            tpm_GMu(tpm_GMu<opt.min_tpm) = opt.min_tpm;
            tpm_GMu(msk_ICV & tpm_GMu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
            tpm_ext(:,:,:,1) = tpm_GMu; % update GM
        else
            tpm_Gu = tpm_GMu./(1+alpha);
            tpm_Pu = tpm_GMu - tpm_Gu;
            tpm_Gu(tpm_Gu<opt.min_tpm) = opt.min_tpm;
            tpm_Gu(msk_ICV & tpm_Gu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
            tpm_Pu(tpm_Gu<opt.min_tpm) = opt.min_tpm;
            tpm_Pu(msk_ICV & tpm_Pu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
            tpm_ext(:,:,:,1) = tpm_Gu; % update GM
            tpm_ext(:,:,:,7) = tpm_Pu; % update BG
        end
    case 1 % WM only
        tpm_WMu = tpm_healthy - tpm_Lu;
        % equiv. to tpm_WMu = tpm_WM .* (1 - (1-1/opt.tpm_ratio) * tpm_l);
        tpm_WMu(tpm_WMu<opt.min_tpm) = opt.min_tpm;
        tpm_WMu(msk_ICV & tpm_WMu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
        tpm_ext(:,:,:,2) = tpm_WMu; % update WM
    case 2 % WM+GM
        tpm_WMu = tpm_WM .* (1 - (1-1/opt.tpm_ratio) * tpm_l);
        tpm_WMu(tpm_WMu<opt.min_tpm) = opt.min_tpm;
        tpm_WMu(msk_ICV & tpm_WMu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
        tpm_GMu = tpm_GM .* (1 - (1-1/opt.tpm_ratio) * tpm_l);
        if tpm_std
            tpm_GMu(tpm_GMu<opt.min_tpm) = opt.min_tpm;
            tpm_GMu(msk_ICV & tpm_GMu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
            tpm_ext(:,:,:,1) = tpm_GMu; % update GM
            tpm_ext(:,:,:,2) = tpm_WMu; % update WM
        else
            tpm_Gu = tpm_GMu./(1+alpha);
            tpm_Pu = tpm_GMu - tpm_Gu;
            tpm_Gu(tpm_Gu<opt.min_tpm) = opt.min_tpm;
            tpm_Gu(msk_ICV & tpm_Gu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
            tpm_Pu(tpm_Gu<opt.min_tpm) = opt.min_tpm;
            tpm_Pu(msk_ICV & tpm_Pu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
            tpm_ext(:,:,:,1) = tpm_Gu; % update GM
            tpm_ext(:,:,:,7) = tpm_Pu; % update BG
            tpm_ext(:,:,:,2) = tpm_WMu; % update WM
        end
    case 3 % WM+GM+CSF
        tpm_WMu = tpm_WM .* (1 - (1-1/opt.tpm_ratio) * tpm_l);
        tpm_WMu(tpm_WMu<opt.min_tpm) = opt.min_tpm;
        tpm_WMu(msk_ICV & tpm_WMu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
        tpm_CSFu = tpm_CSF .* (1 - (1-1/opt.tpm_ratio) * tpm_l);
        tpm_CSFu(tpm_CSFu<opt.min_tpm) = opt.min_tpm;
        tpm_CSFu(msk_ICV & tpm_CSFu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
        tpm_GMu = tpm_GM .* (1 - (1-1/opt.tpm_ratio) * tpm_l);
        if tpm_std
            tpm_GMu(tpm_GMu<opt.min_tpm) = opt.min_tpm;
            tpm_GMu(msk_ICV & tpm_GMu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
            tpm_ext(:,:,:,1) = tpm_GMu; % update GM
            tpm_ext(:,:,:,2) = tpm_WMu; % update WM
            tpm_ext(:,:,:,3) = tpm_CSFu; % update CSF
        else
            tpm_Gu = tpm_GMu./(1+alpha);
            tpm_Pu = tpm_GMu - tpm_Gu;
            tpm_Gu(tpm_Gu<opt.min_tpm) = opt.min_tpm;
            tpm_Gu(msk_ICV & tpm_Gu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
            tpm_Pu(tpm_Gu<opt.min_tpm) = opt.min_tpm;
            tpm_Pu(msk_ICV & tpm_Pu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
            tpm_ext(:,:,:,1) = tpm_Gu; % update GM
            tpm_ext(:,:,:,7) = tpm_Pu; % update BG
            tpm_ext(:,:,:,2) = tpm_WMu; % update WM
            tpm_ext(:,:,:,3) = tpm_CSFu; % update CSF
        end
    otherwise
        error('Wrong tissue flag');
end

% list of tissue, without air
ltpm = [1:5 7:Ntpm_o+1];

% Mask out with swICV, if provided
if isfield(opt,'fn_swICV');
    V_swICV = spm_vol(opt.fn_swICV);
    swICV = spm_read_vols(V_swICV);
    for ii = ltpm
        tmp = tpm_ext(:,:,:,ii).*swICV;
        tmp(tmp<opt.min_tpm) = opt.min_tpm;
        tpm_ext(:,:,:,ii) = tmp;
    end
end

% Update 'other', which is in 6th position of original tpm
tpm_ext(:,:,:,6) = 1 - sum(tpm_ext(:,:,:,ltpm),4);

% 4) save the TPMl, with lesion between WM & CSF, in subject's data dir.
fn_TPMl = fullfile(spm_file(fn_swtMsk,'path'), ...
    spm_file(spm_file(fn_TPM,'filename'),'suffix','_les'));
Vtpm_l = Vtpm;
% adding 1 tpm
Vtpm_l(Ntpm_o+1) = Vtpm(Ntpm_o);
mem_sz = Vtpm(2).pinfo(3)-Vtpm(1).pinfo(3);
Vtpm_l(Ntpm_o+1).pinfo(3) = Vtpm_l(Ntpm_o+1).pinfo(3) + mem_sz;
Vtpm_l(Ntpm_o+1).n(1) = Ntpm_o+1;
if tpm_std
    tc_order = [1 2 7 3 4 5 6]; % the lesion class is inserted in 3rd position!
else
    tc_order = [1 2 8 3 4 5 6 7]; % the lesion class is inserted in 3rd position!
end
for ii=1:Ntpm_o+1
    Vtpm_l(ii).fname = fn_TPMl;
    Vtpm_l(ii) = spm_create_vol(Vtpm_l(ii));
    Vtpm_l(ii) = spm_write_vol(Vtpm_l(ii),tpm_ext(:,:,:,tc_order(ii)));
end

end

%% STEP 5: 
% Creating the segmentatin batch with 7 (or 8) tissue clasess
% + smoothing of modulated warped tissue classes
function [matlabbatch] = batch_segment_l(P,Ptpm_l,opt)
% [matlabbatch] = batch_segment_l(P,Ptpm,param)
%
% INPUT:
% - P     : cell array of structural image filenames, e.g. the MT image.
%           If multiple images are passed, then they enter as different
%           channels
% - Ptpm  : tissue probability maps, inlcuding the MS lesion
% - opt   : structure with some options
%   . b_param : bias correction parameters [regularisation fwhm]
%   . b_write : write out bias corrected
%   . nGauss  : Number of Gaussians per tissue class
%   . mrf     : mrf parameter
%   . cleanup : cleanup parameter

% Multiple channels & number of tissue classes
nP = size(P,1);
nG = numel(opt.nGauss);

b_param = opt.b_param;
b_write = opt.b_write;

if size(b_param,1)<nP
    b_param = ones(nP,1)*b_param(1,:) ;
    % ensure 1 set of bias correction param per channel
end
if size(b_write,1)<nP
    b_write = ones(nP,1)*b_write(1,:) ;
    % ensure 1 set of bias correction param per channel
end

matlabbatch{1}.spm.spatial.preproc.channel.vols = {deblank(P(1,:))};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = b_param(1,1);
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = b_param(1,2);
matlabbatch{1}.spm.spatial.preproc.channel.write = b_write(1,:);
% Adding other channels, if provided
if nP>1
    for ii=2:nP
        matlabbatch{1}.spm.spatial.preproc.channel(ii) = ...
            matlabbatch{1}.spm.spatial.preproc.channel(ii-1);
        matlabbatch{1}.spm.spatial.preproc.channel(ii).vols = {deblank(P(ii,:))};
        matlabbatch{1}.spm.spatial.preproc.channel(ii).biasreg = b_param(ii,1);
        matlabbatch{1}.spm.spatial.preproc.channel(ii).biasfwhm = b_param(ii,2);
        matlabbatch{1}.spm.spatial.preproc.channel(ii).write = b_write(ii,:);
    end
end
% Define TPM's
if nG==7
    cr_native = [1 1 ; 1 1 ; 1 1 ; 1 0 ; 1 0 ; 1 0 ; 0 0 ];
    cr_warped = [1 1 ; 1 1 ; 1 1 ; 1 1 ; 0 0 ; 0 0 ; 0 0 ];
else % Using TPM with BG
    cr_native = [1 1 ; 1 1 ; 1 1 ; 1 0 ; 1 0 ; 1 0 ; 0 0 ; 1 1 ];
    cr_warped = [1 1 ; 1 1 ; 1 1 ; 1 1 ; 0 0 ; 0 0 ; 0 0 ; 1 1 ];
end
for ii = 1:nG
    matlabbatch{1}.spm.spatial.preproc.tissue(ii).tpm = {[Ptpm_l,',',num2str(ii)]};
    matlabbatch{1}.spm.spatial.preproc.tissue(ii).ngaus = opt.nGauss(ii);
    matlabbatch{1}.spm.spatial.preproc.tissue(ii).native = cr_native(ii,:);
    matlabbatch{1}.spm.spatial.preproc.tissue(ii).warped = cr_warped(ii,:);
end
% Define other parameters
matlabbatch{1}.spm.spatial.preproc.warp.mrf = opt.mrf;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = opt.cleanup;
% matlabbatch{1}.spm.spatial.preproc.warp.mrf = 0;
% matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 0;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
% % Smoothing a bit
% matlabbatch{2}.spm.spatial.smooth.data(1) = ...
%     cfg_dep('Segment: mwc1 Images', substruct('.','val', '{}',{1}, ...
%     '.','val', '{}',{1}, '.','val', '{}',{1}), ...
%     substruct('.','tiss', '()',{1}, '.','mwc', '()',{':'}));
% matlabbatch{2}.spm.spatial.smooth.data(2) = ...
%     cfg_dep('Segment: mwc2 Images', substruct('.','val', '{}',{1}, ...
%     '.','val', '{}',{1}, '.','val', '{}',{1}), ...
%     substruct('.','tiss', '()',{2}, '.','mwc', '()',{':'}));
% matlabbatch{2}.spm.spatial.smooth.data(3) = ...
%     cfg_dep('Segment: mwc3 Images', substruct('.','val', '{}',{1}, ...
%     '.','val', '{}',{1}, '.','val', '{}',{1}), ...
%     substruct('.','tiss', '()',{3}, '.','mwc', '()',{':'}));
% matlabbatch{2}.spm.spatial.smooth.data(4) = ...
%     cfg_dep('Segment: mwc8 Images', substruct('.','val', '{}',{1}, ...
%     '.','val', '{}',{1}, '.','val', '{}',{1}), ...
%     substruct('.','tiss', '()',{8}, '.','mwc', '()',{':'}));
% matlabbatch{2}.spm.spatial.smooth.fwhm = [2 2 2];
% matlabbatch{2}.spm.spatial.smooth.dtype = 0;
% matlabbatch{2}.spm.spatial.smooth.im = 0;
% matlabbatch{2}.spm.spatial.smooth.prefix = 's';

end

%% STEP 6: 
% Creating the normalization batch for the MPM
function [matlabbatch] = batch_normalize_MPM(fn_img2warp,fn_warp)
% [malabbatch] = batch_normalize_MPM(fn_img2warp,fn_warp)
%
% INPUT
% - fn_2warp : cell array of filenames of images to warp
% - fn_wapr  : file name of warping image

matlabbatch{1}.spm.spatial.normalise.write.subj.def = {fn_warp};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(char(fn_img2warp));
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
    78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;

end

%% ADDING 2 IMAGES TOGETHER
function add_2_images(fn_in,fn_out,dtype)
% Adding 2 images together
% No checking whatsoever!
Vi = spm_vol(fn_in);
Vo = spm_vol(fn_out);
fl.dtype = dtype;
Vo = spm_imcalc(Vi, Vo, 'i1+i2' ,fl);

end

%% MASKING AN IMAGE
function fn_kimg = mask_img(fn_img,fn_msk,prefix)
% Masking an image by a mask, adding a prefix to the masked image.
if nargin<2, prefix = 'k'; end
fn_kimg = spm_file(fn_img,'prefix',prefix);
Vi = spm_vol(char(fn_img,fn_msk));
Vo = Vi(1);
Vo.fname = fn_kimg;
Vo = spm_imcalc(Vi,Vo,'i1.*(i2>0)');

end
