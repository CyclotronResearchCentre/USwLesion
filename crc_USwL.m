function fn_out = crc_USwL(fn_in,options)
% 
% Function doing all the work of "Unified segmentation with lesion".
% 
% INPUT
% - fn_in   : cell array (1x4) of input filenames
%       {1} : lesion mask image, where 1 is lesion, 0 is healthy. Must be
%             provided.
%       {2} : reference structural image, used for the 1st US with a "cost
%             function masking" (CFM) approach. Must be provided.
%       {3} : structural image(s), used for the multi-channel "US with 
%             Lesion". If none, then the ref-struct is used.
%       {4} : other structural images to be warped along
% - options :
%       imgTpm  : tissue probability maps (just one 4D file)
%       biasreg : bias regularisation value
%       biasfwhm: bias FWHM value (Inf = no bias correction) 
%       biaswr  : flag -> save bias corrected and/or bias field ([0/1 0/1])
%       NbGaussian : number of Gaussians per tissue class, incl. lesion.
%       tpm4lesion : flag -> TPM(s) affected by the lesion
%                      (0, GM; 1, WM; 2, GM+WM; 3, GM+WM+CSF) 
%       ICVmsk  : mask the images images by created ICV-mask, [0/1], in
%                 native and warped space
%       mrf     : MRF parameter
% 
% OUTPUT
% - fn_out :
%       wstruct     : warped structural reference, 1st US with CFM
%       ICVmsk      : intra-cranial volume mask, as generated from str-ref
%       kStruc_i    : masked i^th structural images
%       kOth_i      : masked i^th other image
%       wStruc_i    : warped (masked) i^th structural image
%       wOth_i      : warped (masked) i^th other image
%       TPMl        : subject specific TPM with lesion
%       segmImg     : structure with posterior tissue probabilities
%           c(i)    : class #i in subject space
%           wc(i)   : class #i in MNI space
%           mwc(i)  : modulated class #i in MNI space
% 
% OPERATIONS
% Here are the main steps:
%   1. "Trim 'n grow" the mask image {1} : -> t_Msk / dt_Msk
%       - remov the "small" MS patches using a simple criteria: volume of 
%         patch must be >= minVol (mm3) -> t_Msk image
%       - then grow volume by nDilate voxel(s) -> dt_Msk image
%      minVol & nDilate are set in the crc_USwL_defaults file!
%   2. Apply the mask on the reference structural images {2} -> k_sRef
%   3. Segment the masked structural (k_sRef), normalize the cleaned up
%      mask (t_Msk) and smooth it 
%      -> new tissue probability map for the lesion
%   4. Update the TPMs to include a 7th tissue class -> TPMms
%   Note that the lesion is inserted in *3rd position*, between WM and CSF!
%   5. Do the segmentation with the new TPM_ms on the structural images {3}
%   6. Apply the deformation onto the structural & other images {3,4} 
%       -> create warped images
%   7. Collect all the image filenames created
% 
% NOTES
% - Check the Readme file for further processing details.
% - some parameters/constants are set in the crc_USwL_defaults file and NOT
%   presented in the batch GUI:
%   * minVol    -> minimum volume (mm3) of lesion patch in mask
%   * nDilate   -> #steps for the dilation process of the lesion mask
% - the lesion probability maps is inserted in *3rd position*, i.e. between
%   the WM and CSF tissue probability maps!
%_______________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium
% Cyril Pernet updated few bits to work with no MPM images + added structural
% normalization and N Gaussians - Edinburgh Imaging, The University of Edinburgh

% 
% % TESTING parfor
% nG = options.NbGaussian;
% opt = crc_USwL_get_defaults('uTPM');
% scDefReg = crc_USwL_get_defaults('msksegm.scDefReg');
% 
% fprintf('\nSubject %s: nDil = %d, nG_c3 = %d, ScalReg = %1.1g\n', ...
%     spm_file(fn_in{1},'basename'),opt.nDilate, nG(3), scDefReg);
% fn_out = {};
% 
% return
% 

%% Input data
% fn_in{1} = Mask image
% fn_in{2} = structural reference, for 1st CFM segmentation
% fn_in{3} = All structurals for USwL
% fn_in{4} = Other images (warped along)

% Deal with Struct for USwL
if isempty(fn_in{3})
    fn_in{3} = fn_in{2};
end
nStruc = size(fn_in{3},1);

% Deal with others
if isempty(fn_in{4})
    nOth = 0;
else
    nOth = size(fn_in{4},1);
end
pth = spm_file(fn_in{2},'path');

% Deal with options
if nargin<2, options = struct; end
options_def = struct( ...
	'imgTpm', {crc_USwL_get_defaults('segment.imgTpm')}, ...
    'biasreg', crc_USwL_get_defaults('segment.biasreg'), ...
    'biasfwhm', crc_USwL_get_defaults('segment.biasfwhm'), ...
    'biaswr',crc_USwL_get_defaults('segment.biaswr'), ...
    'NbGaussian', crc_USwL_get_defaults('segment.NbGaussian'), ... 
    'tpm4lesion', crc_USwL_get_defaults('segment.tpm4lesion'), ... 
    'ICVmsk', crc_USwL_get_defaults('segment.ICVmsk'), ...
    'mrf', crc_USwL_get_defaults('segment.mrf') );
options = crc_check_flag(options_def,options); % Check and pad filter structure


% By default no cleanup after US-with-lesion segmentation!
if ~isfield(options,'cleanup') || ...
        isempty(options.cleanup) || strcmp(options.cleanup,'<UNDEFINED>')
    options.cleanup = 0;
end

% Check #Gaussians and #TPMs for the USwLesion segmentation.
NbGaussian = options.NbGaussian;
fn_tpm_USwL = options.imgTpm{1};
fn_tpm_USwL = spm_file(fn_tpm_USwL,'number',''); % remove any number from SPM select
Vtpm_USwL = spm_vol(fn_tpm_USwL);
if numel(NbGaussian)~=(numel(Vtpm_USwL)+1)
    error('There are %d tpm''s (incl. lesion) but only %d #Gaussians provided', ...
        numel(Vtpm_USwL)+1,numel(NbGaussian));
end

%% Define processing parameters for the creation of the updated TPM
opt = crc_USwL_get_defaults('uTPM');

%% 1. "Trim 'n grow" the mask image : -> t_Msk / dt_Msk
% - remov the "small" lesion patches using a simple criteria: volume of
%   lesion patch must be > minVol -> creates on the drive t_Msk
% - then grow volume by nDilate voxel(s) -> creates on the drive dt_Msk
[fn_tMsk,fn_dtMsk] = mask_trimNgrow(fn_in{1},opt.minVol,opt.nDilate);

%% 2. Apply the lesion mask on the reference structural images -> k_sRef
fn_kRef = spm_file(fn_in{2},'prefix','k');
Vi(1) = spm_vol(fn_in{2});
Vi(2) = spm_vol(fn_dtMsk);
Vo = Vi(1);
Vo.fname = fn_kRef;
Vo = spm_imcalc(Vi,Vo,'i1.*(((i2>.5)-1)./((i2>.5)-1))');

%% 3. Segment the masked reference (k_sRef), normalize the cleaned up mask
% (t_Msk) and smooth it -> new TPM for the lesion.
% Then create an ICV mask (for Struct's ICV masking)
% Here use the standard TPMs, as this is a reference structural image
% masked for the lesion, i.e. use the "cost function masking" approach.
clear matlabbatch
fn_tpm_msksegm = crc_USwL_get_defaults('msksegm.imgTpm');
scDefReg = crc_USwL_get_defaults('msksegm.scDefReg');

% Create a matlabbbatch with normalize and smooth operations, pass some key
% input then rely on default values in crc_USwL_defaults.m ('msksegm' bit) 
matlabbatch = crt_batch_normalize_smooth( ...
    fn_kRef, ... % masked structural image used for the warping estimation
    fn_tMsk, ... % cleaned up lesion mask to be warped into MNI
    fn_tpm_msksegm{1}, ... % filename of tissue probability map
    opt.smoKern , ... % smoothing applied on the normalized lesion mask -> new prior
    scDefReg); % scaling of warp regularisation
spm_jobman('run', matlabbatch);
fn_swtMsk = spm_file(fn_tMsk,'prefix','sw'); % smooth normalized lesion mask
fn_wtMsk = spm_file(fn_tMsk,'prefix','w'); %#ok<*NASGU> % normalized lesion mask

% Build ICV, with some cleaning up
fn_TCin = char( ...
    spm_file(fn_kRef,'prefix','c1'),... % GM
    spm_file(fn_kRef,'prefix','c2'),... % WM
    spm_file(fn_kRef,'prefix','c3'),... % CSF
    fn_dtMsk); % dilated lesion mask, as masked struct-ref
opt_ICV = struct( ...
    'fn_ref', fn_kRef, ...
    'fn_warp', spm_file(fn_kRef,'prefix','y_'), ...
    'fn_iwarp', spm_file(fn_kRef,'prefix','iy_'), ...
    'smoK', 4);
fn_icv_out = crc_build_ICVmsk(fn_TCin,opt_ICV);
fn_ICV = deblank(fn_icv_out(1,:)); % -> in subject space
fn_swICV = deblank(fn_icv_out(3,:)); % -> in MNI space, smoothed

% Delete files from this CFM-US step that are not useful anymore:
% c1/2/3 images + (inverse) warps and *_seg8.mat files
to_delete = char( fn_TCin(1:3,:), opt_ICV.fn_warp, opt_ICV.fn_iwarp, ...
    spm_file(fn_kRef,'suffix','_seg8','ext','mat') );
 % for ii=1:size(to_delete,1), delete(deblank(to_delete(ii,:))); end

% Apply the mask (using a sub-function)
% -> this possibly overwrites the masked struct reference
% -> replace names in the fn_in celle array
if options.ICVmsk % ICV-mask the MPMs & others
    fn_tmp = [];
    for ii=1:nStruc
        fn_MPM_ii = deblank(fn_in{3}(ii,:));
        fn_tmp = char(fn_tmp, ...
            mask_img(fn_MPM_ii,fn_ICV,'k'));
    end
    fn_in{3} = fn_tmp(2:end,:);
    % Mask other images too!
    fn_tmp = [];
    for ii=1:nOth
        fn_Oth_ii = deblank(fn_in{4}(ii,:));
        fn_tmp = char(fn_tmp, ...
            mask_img(fn_Oth_ii,fn_ICV,'k'));
    end
    if nOth
        fn_in{4} = fn_tmp(2:end,:);
    end
end

%% 4. Update the TPMs to include an extra tissue class -> TPMms
% Note that the lesion is inserted in *3rd position*, between WM and CSF!
opt_tpm = struct(...
    'tpm4lesion', options.tpm4lesion, ... % tissues to be modified for lesion (0/1/2/3) for GM/WM/GM+WM/ GM+WM+CSF
    'fn_tpm', fn_tpm_USwL, ... % tpm file name
    'tpm_ratio', opt.tpm_ratio, ... % ratio between healthy and lesion tissue
    'min_tpm_icv', opt.min_tpm_icv, ... % minimum value in intracranial volume
    'min_tpm', opt.min_tpm); % minum value overall
if options.ICVmsk && nStruc ~= 0 % ICV-mask the TPMs
    opt_tpm.fn_swICV = fn_swICV; % smoothed-warped ICV mask to apply on TPMs
end
fn_TPMl = update_TPM_with_lesion(opt_tpm, fn_swtMsk); % that creates the new tissue class tmp

%% 5. Do the segmentation with the new TPM_ms & create final ICV mask
fn_Img2segm = fn_in{3};

scDefReg = crc_USwL_get_defaults('segment.scDefReg');
% scaling of warping regularisation
opt_segm = struct( ...
    'b_param', [options.biasreg options.biasfwhm], ...
    'b_write', options.biaswr, ...
    'nGauss', NbGaussian, ...
    'mrf', options.mrf, ...
    'cleanup', options.cleanup, ...
    'scDefReg', scDefReg);

% Create the matlabbatch & run it
clear matlabbatch
[matlabbatch] = crt_batch_USwL(fn_Img2segm, fn_TPMl, opt_segm);
spm_jobman('run', matlabbatch);

fn_Cimg   = spm_select('FPList',pth, ...  % native space
    ['^c[0-9].*',spm_file(fn_Img2segm(1,:),'basename'),'\.nii$']); 
fn_rCimg  = spm_select('FPList',pth, ...  % native dartel imported
    ['^rc[0-9].*',spm_file(fn_Img2segm(1,:),'basename'),'\.nii$']);
fn_wCimg  = spm_select('FPList',pth, ...  % warped
    ['^wc[0-9].*',spm_file(fn_Img2segm(1,:),'basename'),'\.nii$']);
fn_mwCimg = spm_select('FPList',pth, ...  % modulated warped
    ['^mwc[0-9].*',spm_file(fn_Img2segm(1,:),'basename'),'\.nii$']);

% When using seperate GM for BG, i.e. extended TPM, then it is usefull to 
% recombine GM-w/o-BG with GM-BG
%  -> add c8 onto c1 -> only 1 image (c1) with GM + c8 with GM-BG.
% Use a specific sub-function to preserve informations!
if numel(NbGaussian)==8
    add_2_images(fn_Cimg([1 end],:));
    add_2_images(fn_rCimg([1 end],:));
    add_2_images(fn_wCimg([1 end],:));
    add_2_images(fn_mwCimg([1 end],:));
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

if options.ICVmsk % ICV-mask GM, WM, Lesion and CSF (+BG if there)
    for ii=1:4
        mask_img(fn_Cimg(ii,:),fn_ICV,'');
        mask_img(fn_wCimg(ii,:),fn_wICV,'');
        mask_img(fn_mwCimg(ii,:),fn_wICV,'');
    end
    % Case of GM-BG specific tissue class
    if numel(NbGaussian)==8
        mask_img(fn_Cimg(end,:),fn_ICV,'');
        mask_img(fn_wCimg(end,:),fn_wICV,'');
        mask_img(fn_mwCimg(end,:),fn_wICV,'');
    end
end

% Pick up bias corrected images, if they were produced.
fn_mStruc = spm_file(fn_in{3},'prefix','m');
exist_mStruc = false(nStruc,1);
for ii=1:nStruc
    if exist(fn_mStruc(ii,:),'file')
        exist_mStruc(ii) = true;
    end
end
if all(exist_mStruc)
    fn_in{3} = fn_mStruc;
end

% Apply the final ICV mask on struct & Others
if options.ICVmsk % ICV-mask the Struct & others
    fn_tmp = [];
    for ii=1:nStruc
        fn_Struc_ii = deblank(fn_in{3}(ii,:));
        fn_tmp = char(fn_tmp, ...
            mask_img(fn_Struc_ii,fn_ICV,'k'));
    end
    fn_ICV_in3 = fn_tmp(2:end,:);
    % Mask other images too!
    fn_tmp = [];
    for ii=1:nOth
        fn_Oth_ii = deblank(fn_in{4}(ii,:));
        fn_tmp = char(fn_tmp, ...
            mask_img(fn_Oth_ii,fn_ICV,'k'));
    end
    if nOth
        fn_ICV_in4 = fn_tmp(2:end,:);
    end
end

%% 6. Apply the deformation onto the Stru/Other -> warped Struc/Other
% Re-create an ICV mask

fn_warp = spm_file(fn_Img2segm(1,:),'prefix','y_');
% Apply on all images: structs + others, if available

fn_img2warp = {fn_in{3}}; % Use Struc images, Ref image not necessary 
if ~isempty(fn_in{4})
    fn_img2warp = {char(fn_img2warp{1} , fn_in{4})};
end

clear matlabbatch
[matlabbatch] = crt_batch_normalize_StrucOth(fn_img2warp,fn_warp);
spm_jobman('run', matlabbatch);

fn_warped_struct = spm_file(fn_in{2},'prefix','w');
if ~isempty(fn_in{3})
    fn_warped_Struc = spm_file(fn_in{3},'prefix','w');
else
    fn_warped_Struc = '';
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

if options.thrMPM && nStruc ~= 0
    for ii=1:nStruc
        fn_out.(sprintf('fxMPM_%d',ii)) = ...
            {spm_file(deblank(fn_in_3_orig(ii,:)),'prefix','t')};
        fn_out.(sprintf('fxMPMmsk_%d',ii)) = ...
            {spm_file(fn_in_3_orig(ii,:),'prefix','fx_')};
    end
end

fn_out.ICVmsk = {fn_ICV};

if options.ICVmsk && nStruc ~= 0;
    for ii=1:nStruc
        fn_out.(sprintf('kMPM_%d',ii)) = {deblank(fn_in{3}(ii,:))};
    end
    for ii=1:nOth
        fn_out.(sprintf('kOth_%d',ii)) = {deblank(fn_in{4}(ii,:))};
    end
end

if ~isempty(fn_warped_Struc) % warped MPMs
    for ii=1:size(fn_warped_Struc,1)
        fn_out.(sprintf('wStruc_%d',ii)) = {deblank(fn_warped_Struc(ii,:))};
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

% if options.thrLesion ~= 0
%     crc_lesion_cleanup(fn_out.segmImg.c3,options.thrLesion);
% end
% fn_out = {};
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
    Msk(Msk>.5) = 1; % non-zero values, as in >1e-6, set to 1
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
else
    dMsk_nM = ~~Msk_nM;
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

function matlabbatch = crt_batch_normalize_smooth(fn_kRef,fn_tMsk,fn_TPM,smoKern,scDefReg)
% matlabbatch = crt_batch_normalize_smooth(fn_kRef,fn_tMsk,fn_TPM,smoKern,scDefReg)
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
% - scDefReg: scaling of default warping regularisation
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
matlabbatch{3}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2]*scDefReg; 
matlabbatch{3}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{3}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{3}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{3}.spm.spatial.preproc.warp.write = [1 1];
matlabbatch{4}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Named File Selector: LesionMask(1) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
matlabbatch{4}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72 ; 90 90 108];
matlabbatch{4}.spm.spatial.normalise.write.woptions.vox = [1.5 1.5 1.5];
matlabbatch{4}.spm.spatial.normalise.write.woptions.interp = 1; % 1 -> trilinear, 0 -> NN
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
% Implicit smoothing is 4mm!!!
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
    % Special TPM for MPM (with separate GM class for subcortical areas)
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
l_les_possible = find( (prob_l_possible(:)*opt.min_tpm_icv>tpm_Lu(:)) ...
                     & (prob_l_possible(:)>0) );
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
function [matlabbatch] = crt_batch_USwL(P,Ptpm_l,opt)
% [matlabbatch] = crt_batch_USwL(P,Ptpm,param)
%
% INPUT:
% - P     : cell array of structural image filenames, e.g. the MT image.
%           If multiple images are passed, then they enter as different
%           channels
% - Ptpm  : tissue probability maps, including one for the lesion (3rd pos)
% - opt   : structure with some options
%   . b_param : bias correction parameters [regularisation fwhm]
%   . b_write : write out bias corrected
%   . nGauss  : Number of Gaussians per tissue class
%   . mrf     : mrf parameter
%   . cleanup : cleanup parameter
%   . scDefReg: scaling of default warping regularisation

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
elseif nG==8 % Using TPM with specific GM of basal Ganglia
    cr_native = [1 1 ; 1 1 ; 1 1 ; 1 0 ; 1 0 ; 1 0 ; 0 0 ; 1 1 ];
    cr_warped = [1 1 ; 1 1 ; 1 1 ; 1 1 ; 0 0 ; 0 0 ; 0 0 ; 1 1 ];
else
    error('USwL:multichanSegm','Wrong number of Gaussians/tissue classes.');
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
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2]*opt.scDefReg;
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
function [matlabbatch] = crt_batch_normalize_StrucOth(fn_img2warp,fn_warp)
% [malabbatch] = crt_batch_normalize_StrucOth(fn_img2warp,fn_warp)
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
function add_2_images(fn_in)
% Adding 2 images together
% No checking whatsoever!
% Deal with sum through the nifti object -> save in 1st one
% Why ?
% -> preserves crucial orientation info in rc* files for Dartel
% -> much faster than ImCalc
V1 = spm_vol(fn_in(1,:));
V2 = spm_vol(fn_in(2,:));
V1.private.dat(:,:,:) = V1.private.dat(:,:,:) + V2.private.dat(:,:,:);

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
