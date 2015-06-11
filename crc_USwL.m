function fn_out = crc_USwL(job)

testing = true;

if ~testing
    %% Define defaults processing parameters
    opt = struct( ...
        'fn_tpm', 'nwTPM_sl2.nii', ... % filename of tpm to use (could also be 'TPM.nii')
        'minNr', 8, ...    % #voxels in lesion patch must be > minNr
        'nDilate', 2, ...  % # of dilation step
        'smoKern', 2, ... % smoothing (in mm) of the warped lesion mask
        'tpm_ratio', 100, ... % ratio of lesion/tpm
        'min_tpm', 1e-6, ... % minimum value of tpm
        'min_tpm_icv', 1e-3, ... % minimum value of tpm in intracranial volume
        'b_param', [.00001 Inf], ... % no bias correction needed
        'b_write', [0 0] ... % not writing bias corrected images
        );
    %     'img2segm', 'MPM+FLAIR', ... % use all images
    %     'fwhm', 8 ... % smoothing for the preprocessed maps
    
    %% Collect input -> to fit into previously written code. :-)
    fn_in{1} = spm_file(job.imgMsk{1},'number','');
    fn_in{2} = spm_file(job.imgRef{1},'number','');
    fn_in{3} = char(spm_file(job.imgMPM,'number',''));
    fn_in{4} = char(spm_file(job.imgOth,'number',''));
    
    %% 1. "Trim 'n grow" the mask image : -> t_Msk / dt_Msk
    % - remov the "small" MS patches using a simple criteria: number of voxels
    %   in patch must be > minNR    -> t_Msk
    % - then grow volume by 1 voxel -> dt_Msk
    [fn_tMsk,fn_dtMsk] = mask_trimNgrow(fn_in{1},opt.minNr,opt.nDilate);
    
    %% 2. Apply the mask on the reference structural images -> k_sRef
    fn_kMTw = spm_file(fn_in{2},'prefix','k');
    Vi(1) = spm_vol(fn_in{2});
    Vi(2) = spm_vol(fn_dtMsk);
    Vo = Vi(1);
    Vo.fname = fn_kMTw;
    Vo = spm_imcalc(Vi,Vo,'i1.*(((i2>.5)-1)./((i2>.5)-1))');
    pth = spm_file(fn_in{1},'path');
    
    %% 3. Segment the masked structural (k_sRef), normalize the cleaned up mask
    % (t_Msk) and smooth it -> new TPM for the lesion.
    clear matlabbatch
    [matlabbatch] = batch_normalize_smooth(fn_kMTw,fn_tMsk,opt.smoKern);
    spm_jobman('run', matlabbatch);
    fn_swtMsk = spm_file(fn_tMsk,'prefix','sw');
    fn_wtMsk = spm_file(fn_tMsk,'prefix','w');
    
    %% 4. Update the TPMs to include a 7th tissue class -> TPMms
    % Note that the lesion is inserted in *3rd position*, between WM and CSF!
    opt_tpm = struct(...
        'tpm4lesion', job.options.tpm4lesion, ...
        'fn_tpm', opt.fn_tpm, ...
        'tpm_ratio', opt.tpm_ratio, ...
        'min_tpm_icv', opt.min_tpm_icv, ...
        'min_tpm', opt.min_tpm);
    
    fn_TPMl = update_TPM_with_lesion(opt_tpm, fn_swtMsk);
    % fn_TPMl = fullfile(spm_file(fn_swtMsk,'path'), ...
    %               spm_file(opt.fn_tpm,'suffix','_l'));
    
    %% 5. Do the segmentation with the new TPM_ms
    % img4US = 0 -> Structural reference only
    %        = 1 -> all MPMs
    %        = 2 -> all MPMs + others
    
    switch job.options.img4US
        case 0
            fn_Img2segm = fn_in{2}; %#ok<*CCAT1>
        case 1
            fn_Img2segm = fn_in{3};
        case 2
            fn_Img2segm = char(fn_in{3} , fn_in{4});
    end
    opt_segm  =struct( ...
        'b_param', opt.b_param, ...
        'b_write', opt.b_write);
    clear matlabbatch
    [matlabbatch] = batch_segment_l(fn_Img2segm, fn_TPMl, opt_segm);
    spm_jobman('run', matlabbatch);
    
    %% 6. Apply the deformation onto the MPMs -> warped MPMs
    
    fn_warp = spm_file(fn_Img2segm(1,:),'prefix','y_');
    % Apply on all images: MPM + others
    fn_img2warp = {char(fn_in{3} , fn_in{4})};
    clear matlabbatch
    [matlabbatch] = batch_normalize_MPM(fn_img2warp,fn_warp);
    spm_jobman('run', matlabbatch);
    fn_warped_MPM = spm_file(fn_in{3},'prefix','w');
    fn_warped_Oth = spm_file(fn_in{4},'prefix','w');
    fn_mwTC = char( ...
        spm_file(fn_in{3}(1,:),'prefix','smwc1'), ...
        spm_file(fn_in{3}(1,:),'prefix','smwc2'), ...
        spm_file(fn_in{3}(1,:),'prefix','smwc3') );
    
    %% 8. Collect all the image filenames created    
    if ~isempty(fn_warped_MPM) % warped MPMs
        for ii=1:size(fn_warped_MPM,1)
          fn_out.(['wMPM',num2str(ii)]) = {deblank(fn_warped_MPM(ii,:))};
        end
    end
    if ~isempty(fn_warped_Oth)
        for ii=1:size(fn_warped_MPM,1) % warped Others
          fn_out.(['wOth',num2str(ii)]) = {deblank(fn_warped_Oth(ii,:))};
        end
    end
    tmp = spm_select('FPList',pth,'^c[0-9].*\.nii$'); % segmented tissues
    fn_out.segmImg.c1 = {deblank(tmp(1,:))}; % GM
    fn_out.segmImg.c2 = {deblank(tmp(2,:))}; % WM
    fn_out.segmImg.c3 = {deblank(tmp(3,:))}; % Lesion
    fn_out.segmImg.c4 = {deblank(tmp(4,:))}; % CSF
    tmp = spm_select('FPList',pth,'^wc[0-9].*\.nii$'); % segmented tissues
    fn_out.segmImg.wc1 = {deblank(tmp(1,:))}; % warped GM
    fn_out.segmImg.wc2 = {deblank(tmp(2,:))}; % warped WM
    fn_out.segmImg.wc3 = {deblank(tmp(3,:))}; % warped Lesion
    fn_out.segmImg.wc4 = {deblank(tmp(4,:))}; % warped CSF
    tmp = spm_select('FPList',pth,'^mwc[0-9].*\.nii$'); % segmented tissues
    fn_out.segmImg.mwc1 = {deblank(tmp(1,:))}; % modulated warped GM
    fn_out.segmImg.mwc2 = {deblank(tmp(2,:))}; % modulated warped WM
    fn_out.segmImg.mwc3 = {deblank(tmp(3,:))}; % modulated warped Lesion
    fn_out.segmImg.mwc4 = {deblank(tmp(4,:))}; % modulated warped CSF
    fn_out.TPMl = fn_TPMl;
else
    fn_in{1} = spm_file(job.imgMsk{1},'number','');
    fn_in{2} = spm_file(job.imgRef{1},'number','');
    fn_in{3} = char(spm_file(job.imgMPM,'number',''));
    fn_in{4} = char(spm_file(job.imgOth,'number',''));
    pth = spm_file(fn_in{1},'path');
    if ~isempty(fn_in{3})
        fn_warped_MPM = spm_file(fn_in{3},'prefix','w');
        for ii=1:size(fn_warped_MPM,1)
          fn_out.(['wMPM',num2str(ii)]) = {deblank(fn_warped_MPM(ii,:))};
        end
    end
    if ~isempty(fn_in{4})
        fn_warped_Oth = spm_file(fn_in{4},'prefix','w');
        for ii=1:size(fn_warped_MPM,1)
          fn_out.(['wOth',num2str(ii)]) = {deblank(fn_warped_Oth(ii,:))};
        end
    end
    tmp = spm_select('FPList',pth,'^c[0-9].*\.nii$'); % segmented tissues
    fn_out.segmImg.c1 = {deblank(tmp(1,:))}; % GM
    fn_out.segmImg.c2 = {deblank(tmp(2,:))}; % WM
    fn_out.segmImg.c3 = {deblank(tmp(3,:))}; % Lesion
    fn_out.segmImg.c4 = {deblank(tmp(4,:))}; % CSF
    tmp = spm_select('FPList',pth,'^wc[0-9].*\.nii$'); % segmented tissues
    fn_out.segmImg.wc1 = {deblank(tmp(1,:))}; % warped GM
    fn_out.segmImg.wc2 = {deblank(tmp(2,:))}; % warped WM
    fn_out.segmImg.wc3 = {deblank(tmp(3,:))}; % warped Lesion
    fn_out.segmImg.wc4 = {deblank(tmp(4,:))}; % warped CSF
    tmp = spm_select('FPList',pth,'^mwc[0-9].*\.nii$'); % segmented tissues
    fn_out.segmImg.mwc1 = {deblank(tmp(1,:))}; % modulated warped GM
    fn_out.segmImg.mwc2 = {deblank(tmp(2,:))}; % modulated warped WM
    fn_out.segmImg.mwc3 = {deblank(tmp(3,:))}; % modulated warped Lesion
    fn_out.segmImg.mwc4 = {deblank(tmp(4,:))}; % modulated warped CSF
    fn_TPMl = spm_select('FPList',pth,'TPM.*\.nii$');
    fn_out.TPMl = {fn_TPMl};
end
end

%% =======================================================================
%% SUBFUNCTIONS
%% =======================================================================

%% STEP 1: Removing small lesion patches from mask
function [fn_tMsk,fn_dtMsk] = mask_trimNgrow(P_in,minNr,nDilate)
% 1) Trim a mask image by removing bits that would be to small to really
%    matter according to medical criteria (cf. E. Lommers):
%    "Lesions will ordinarily be larger than 3 mm in cross section"
%    With 1x1x1mm^3 voxels, a cube of 2x2x2 voxels has a diagonal of
%    sqrt(12)~3.4mm and counts 8 voxels -> minNr = 8
% -> fn_tMsk used for the new TPM_ms
% 2) Then grow the volume by 2 voxels
% -> fn_dtMsk used for the masking for the 1st warping

if nargin<3
    nDilate = 2;
end
if nargin<2
    minNr = 8;
end

% 1) Load things
V = spm_vol(P_in);
[Msk,XYZ] = spm_read_vols(V);
XYZvx = V.mat\[XYZ ; ones(1,size(XYZ,2))];

% 2) Clean up
lMsk  = find(Msk(:)>0);
lXYZvx = XYZvx(1:3,(lMsk));
vMsk = Msk(lMsk);

[Ncl,Zcl,Mcl,Acl,XYZcl] = spm_max(vMsk,lXYZvx); %#ok<*NASGU,*ASGLU>
% figure, hist(unique(Ncl))
% figure, hist(unique(Ncl(Ncl<100)))

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

%% STEP 3: Creating the normalization batch for the masked structural image
function [matlabbatch] = batch_normalize_smooth(fn_kRef,fn_tMsk,smoKern)
% [matlabbatch] = batch_normalize_smooth(fn_kRef,fn_tMsk,smoKern)
%
% INPUT:
% - fn_kRef : masked structural image used for the warping estimation
% - fn_tMsk : cleaned up lesion mask to be warped into MNI
% - smoKern : smoothing applied on the normalized lesion mask -> new prior
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {fn_kRef};
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = {fn_tMsk};
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-90 -126 -72
    90 90 108];
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [1.5 1.5 1.5];
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 1; % trilinear
matlabbatch{2}.spm.spatial.smooth.data(1) = ...
    cfg_dep('Normalise: Estimate & Write: Normalised Images (Subj 1)', ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, ...
    '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{2}.spm.spatial.smooth.fwhm = smoKern*[1 1 1];
matlabbatch{2}.spm.spatial.smooth.dtype = 16;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's';

end

%% STEP 4: Updating the TPM with a 7th class, the lesion
% Note that the lesion is inserted in *3rd position*, between WM and CSF!
function fn_TPMl = update_TPM_with_lesion(opt, fn_swtMsk)
% fn_TPMl = update_TPM_with_lesion(opt, fn_swtMsk)
%
% INPUT
% - opt : structure with a few parameters
%     .tpm4lesion : tissues to be modified for lesion
%     .fn_tpm : tpm file name
%     .tpm_ratio : ration between WM and lesion
%     .min_tpm_icv : minimum value in intracranial volume
%     .min_tpm : minum value overall
% - fn_swtMsk : filename of smoothed normalized cleaned lesion mask, to be
%               used as the 7th tissue class

% NOTE:
% Only supporting WM lesion at the moment!!!
if opt.tpm4lesion~=0
    error('LESIONS IN GM NOT SUPPORTED YET!!!');
end
% 0) select TPM and load
fn_TPM = fullfile(spm('dir'),'tpm',opt.fn_tpm);
Vtpm = spm_vol(fn_TPM);
tpm_orig = spm_read_vols(Vtpm);
tpm_WM = squeeze(tpm_orig(:,:,:,2));
Vl = spm_vol(fn_swtMsk);
tpm_l = spm_read_vols(Vl);

% 1) scale MS and adjust WM in intracranial volume
tpm_l = (1-1/opt.tpm_ratio)*tpm_l.*tpm_WM;
ll = find(tpm_WM>=1e-3 & tpm_l<opt.min_tpm_icv); %#ok<*BDSCI> % find intracranial volume from WM tpm
tpm_WM = tpm_WM - tpm_l;
tpm_l(ll) = opt.min_tpm_icv;

% 2) ensure minium value all over
tpm_WM(tpm_WM<opt.min_tpm) = opt.min_tpm;
tpm_l(tpm_l<opt.min_tpm) = opt.min_tpm;

% 3) concatenate by setting lesion at #7 & adjust 'other' class
tpm_ext = cat(4,tpm_orig,tpm_l);
tpm_ext(:,:,:,2) = tpm_WM; % update WM
tpm_ext(:,:,:,6) = 1-sum(tpm_ext(:,:,:,[1:5 7]),4); % update 'other'

% 4) save the TPMl, with lesion in #3
fn_TPMl = fullfile(spm_file(fn_swtMsk,'path'),spm_file(opt.fn_tpm,'suffix','_l'));
Vtpm_l = Vtpm;
Vtpm_l(7) = Vtpm(6);
mem_sz = Vtpm(2).pinfo(3)-Vtpm(1).pinfo(3);
Vtpm_l(7).pinfo(3) = Vtpm_l(7).pinfo(3) + mem_sz;
Vtpm_l(7).n(1) = 7;
tc_order = [1 2 7 3 4 5 6]; % the lesion class is inserted in 3rd position!
for ii=1:7
    Vtpm_l(ii).fname = fn_TPMl;
    Vtpm_l(ii) = spm_create_vol(Vtpm_l(ii));
    Vtpm_l(ii) = spm_write_vol(Vtpm_l(ii),tpm_ext(:,:,:,tc_order(ii)));
end

end

%% STEP 5: Creating the segmentatin batch with 7 tissue clasess
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

% Multiple channels?
nP = size(P,1);

if nargin<3
    opt = struct(...
        'b_param',ones(nP,1)*[.00001 Inf],...
        'b_write',zeros(nP,2));
    % By default no bias
end

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
% Define 7 TPM's
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[Ptpm_l,',1']}; % GM
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 2; % increased from default value 1, could be Inf for nonparametric approach
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[Ptpm_l,',2']}; % WM
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 2; % increased from default value 1
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[Ptpm_l,',3']}; % Lesion TPM!
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[Ptpm_l,',4']}; % CSF
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[Ptpm_l,',5']}; % bones
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[Ptpm_l,',6']}; % soft tissues
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(7).tpm = {[Ptpm_l,',7']}; % others
matlabbatch{1}.spm.spatial.preproc.tissue(7).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(7).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(7).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
matlabbatch{2}.spm.spatial.smooth.data(1) = ...
    cfg_dep('Segment: mwc1 Images', substruct('.','val', '{}',{1}, ...
    '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','tiss', '()',{1}, '.','mwc', '()',{':'}));
matlabbatch{2}.spm.spatial.smooth.data(2) = ...
    cfg_dep('Segment: mwc2 Images', substruct('.','val', '{}',{1}, ...
    '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','tiss', '()',{2}, '.','mwc', '()',{':'}));
matlabbatch{2}.spm.spatial.smooth.data(3) = ...
    cfg_dep('Segment: mwc3 Images', substruct('.','val', '{}',{1}, ...
    '.','val', '{}',{1}, '.','val', '{}',{1}), ...
    substruct('.','tiss', '()',{3}, '.','mwc', '()',{':'}));
matlabbatch{2}.spm.spatial.smooth.fwhm = [2 2 2];
matlabbatch{2}.spm.spatial.smooth.dtype = 0;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's';

end

%% STEP 6: Creating the normalization batch for the MPM
function [matlabbatch] = batch_normalize_MPM(fn_img2warp,fn_warp)
% [malabbatch] = batch_normalize_MPM(fn_img2warp,fn_warp)
%
% INPUT
% - fn_2warp : cell array of filenames of images to warp
% - fn_wapr  : file name of warping image

%-----------------------------------------------------------------------
% Job saved on 10-Mar-2015 18:50:37 by cfg_util (rev $Rev: 6134 $)
% spm SPM - SPM12 (12.0)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {fn_warp};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(char(fn_img2warp));
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
    78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;

end



