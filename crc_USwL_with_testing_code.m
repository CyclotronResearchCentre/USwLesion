function fn_out = crc_USwL(job)
% Doing all the work of "Unified segmentation with lesion".
% Here are the main steps:
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


testing = true;
% testing = false;

%% Collect input -> to fit into previously written code. :-)
fn_in{1} = spm_file(job.imgMsk{1},'number',''); % Mask image
fn_in{2} = spm_file(job.imgRef{1},'number',''); % structural reference
fn_in{3} = char(spm_file(job.imgMPM,'number','')); % All MPM's
fn_in{4} = char(spm_file(job.imgOth,'number','')); % Other images
nMPM = size(fn_in{3},1);

if ~testing
    %% Define defaults processing parameters
    opt = struct( ...
        'minNr', 8, ...    % #voxels in lesion patch must be > minNr
        'nDilate', 2, ...  % # of dilation step
        'smoKern', 2, ... % smoothing (in mm) of the warped lesion mask
        'tpm_ratio', 100, ... % ratio of lesion/tpm
        'min_tpm', 1e-6, ... % minimum value of tpm overall
        'min_tpm_icv', 1e-3, ... % minimum value of tpm in intracranial volume
        'b_param', [.00001 Inf], ... % no bias correction needed
        'b_write', [0 0] ... % not writing bias corrected images
        );
    %     'thrMPM', true, ...    % threshold MPM images to avoid unsually large/negative values
    %     'ICVmskMPM', true, ... % mask the MPMs to keep the ICV = skull strip
    
    
    %% 0. Clean up of the MPM images!
    % Need to know the order of the images, ideally MT, A, R1, R2 and should
    % check with their filename? based on '_MT', '_A', '_R1', '_R2'?
    if job.options.thrMPM
        strMPM = {'_A', '_MT', '_R1', '_R2'}; nSt = numel(strMPM);
        thrMPM = [200 5 5 100]; % Thresholds for A, MT, R1 & R2.
        fn_tmp = [];
        for ii=1:nMPM % Loop over MPM files
            mtch = zeros(nSt,1);
            for jj=1:nSt
                tmp = strfind(spm_file(fn_in{3}(ii,:),'filename'),strMPM{jj});
                if ~isempty(tmp), mtch(jj) = tmp(end); end % pick last index if many
            end
            [~,p_mtch] = max(mtch);
            if p_mtch
                fn_tmp = char( fn_tmp , ...
                    fix_MPMintens(deblank(fn_in{3}(ii,:)),thrMPM(p_mtch)));
            else
                fprintf('\nCould not fix file : %s',fn_in{3}(ii,:))
                fn_tmp = char( fn_tmp , deblank(fn_in{3}(ii,:)));
            end
        end
        fn_in{3} = fn_tmp(2:end,:);
    end
    
    %% 1. "Trim 'n grow" the mask image : -> t_Msk / dt_Msk
    % - remov the "small" lesion patches using a simple criteria: number of
    %   voxels in patch must be > minNR    -> t_Msk
    % - then grow volume by 1 voxel     -> dt_Msk
    [fn_tMsk,fn_dtMsk] = mask_trimNgrow(fn_in{1},opt.minNr,opt.nDilate);
    
    %% 2. Apply the mask on the reference structural images -> k_sRef
    fn_kMTw = spm_file(fn_in{2},'prefix','k');
    Vi(1) = spm_vol(fn_in{2});
    Vi(2) = spm_vol(fn_dtMsk);
    Vo = Vi(1);
    Vo.fname = fn_kMTw;
    Vo = spm_imcalc(Vi,Vo,'i1.*(((i2>.5)-1)./((i2>.5)-1))');
    pth = spm_file(fn_in{2},'path');
    
    %% 3. Segment the masked structural (k_sRef), normalize the cleaned up mask
    % (t_Msk) and smooth it -> new TPM for the lesion.
    % Then create an ICV mask for MPM's ICV masking
    clear matlabbatch
    [matlabbatch,fn_ICV] = batch_normalize_smooth(fn_kMTw,fn_tMsk,job.options.imgTpm{1},opt.smoKern);
    spm_jobman('run', matlabbatch);
    fn_swtMsk = spm_file(fn_tMsk,'prefix','sw'); % smooth normalized lesion mask
    fn_wtMsk = spm_file(fn_tMsk,'prefix','w'); %#ok<*NASGU> % normalized lesion mask
    
    if job.options.ICVmsk % ICV-mask the MPMs
        fn_in_3_orig = fn_in{3};
        fn_tmp = [];
        for ii=1:nMPM
            fn_MPM_ii = deblank(fn_in{3}(ii,:));
            Vi(1) = spm_vol(fn_MPM_ii);
            Vi(2) = spm_vol(fn_ICV);
            Vo = Vi(1);
            Vo.fname = spm_file(fn_MPM_ii,'prefix','k');
            Vo = spm_imcalc(Vi,Vo,'i1.*i2');
            fn_tmp = char(fn_tmp,Vo.fname);
        end
        fn_in{3} = fn_tmp(2:end,:);
    end
    
    %% 4. Update the TPMs to include a 7th tissue class -> TPMms
    % Note that the lesion is inserted in *3rd position*, between WM and CSF!
    opt_tpm = struct(...
        'tpm4lesion', job.options.tpm4lesion, ...
        'fn_tpm', job.options.imgTpm, ...
        'tpm_ratio', opt.tpm_ratio, ...
        'min_tpm_icv', opt.min_tpm_icv, ...
        'min_tpm', opt.min_tpm);
    
    fn_TPMl = update_TPM_with_lesion(opt_tpm, fn_swtMsk);
    
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
        spm_file(fn_in{3}(1,:),'prefix','smwc3') ); %#ok<*NASGU>
    
    
    %% 7. Collect all the image filenames created
    if job.options.thrMPM
        for ii=1:nMPM
            fn_out.(sprintf('thrMPM%d',ii)) = ...
                {spm_file(deblank(fn_in_3_orig(ii,:)),'prefix','t')};
            fn_out.(sprintf('thrMPMmsk%d',ii)) = ...
                {spm_file(fn_in_3_orig(ii,:),'prefix','msk_')};
        end
    end
    fn_out.ICVmsk = {fn_ICV};
    if job.options.ICVmsk
        for ii=1:nMPM
            fn_out.(sprintf('kMPM%d',ii)) = {deblank(fn_in{3}(ii,:))};
        end
    end
    if ~isempty(fn_warped_MPM) % warped MPMs
        for ii=1:size(fn_warped_MPM,1)
            fn_out.(sprintf('wMPM%d',ii)) = {deblank(fn_warped_MPM(ii,:))};
        end
    end
    if ~isempty(fn_warped_Oth)
        for ii=1:size(fn_warped_Oth,1) % warped Others
            fn_out.(sprintf('wOth%d',ii)) = {deblank(fn_warped_Oth(ii,:))};
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
    fn_out.TPMl = {fn_TPMl};
else
    pth = spm_file(fn_in{1},'path');
    fn_ICV = spm_file(fn_in{2},'prefix','icv_k');

    if ~isempty(fn_in{3})
        fn_warped_MPM = spm_file(fn_in{3},'prefix','wkt');
        for ii=1:size(fn_warped_MPM,1)
            fn_out.(['wMPM',num2str(ii)]) = {deblank(fn_warped_MPM(ii,:))};
        end
    end
    if ~isempty(fn_in{4})
        fn_warped_Oth = spm_file(fn_in{4},'prefix','w');
        for ii=1:size(fn_warped_Oth,1)
            fn_out.(['wOth',num2str(ii)]) = {deblank(fn_warped_Oth(ii,:))};
        end
    end
    if job.options.thrMPM
        for ii=1:nMPM
            fn_out.(sprintf('thrMPM%d',ii)) = ...
                {spm_file(deblank(fn_in{3}(ii,:)),'prefix','t')};
            fn_out.(sprintf('thrMPMmsk%d',ii)) = ...
                {spm_file(fn_in{3}(ii,:),'prefix','msk_')};
        end
    end
    fn_out.options.ICVmsk = {fn_ICV};
    if job.options.ICVmsk
        for ii=1:nMPM
            fn_out.(sprintf('kMPM%d',ii)) = ...
                {spm_file(deblank(fn_in{3}(ii,:)),'prefix','kt')};
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

%% STEP 0: Fixing intensities of MPM images
function fn_out = fix_MPMintens(fn_in,thrMPM)
% Make sure that MPM intensities are within [0 thrMPM] by capping the
% values. The resulting image is written out with the prefix 't'.
% On top, create a binary mask of voxels that were "fixed", with a value
% of 1 if the voxel value was <0, or 2 if >thrMPM.

crt_mask = true;

V = spm_vol(fn_in);
dd = spm_read_vols(V);
sz_dd = size(dd); dd = dd(:);

if crt_mask
    ll_fix = (dd<0) + (dd>thrMPM)*2;
    Vf = V;
    Vf.dt(1) = 2; % uint8
    Vf.fname = spm_file(V.fname,'prefix','msk_');
    Vf.descrip = 'fixed voxels, 1 if <0 and 2 if >thrMPM';
    Vf = spm_create_vol(Vf);
    Vf = spm_write_vol(Vf,reshape(ll_fix,sz_dd));
end

dd = abs(dd);
NaboveThr = sum(dd>thrMPM);
dd(dd>thrMPM) = thrMPM * (1 + randn(NaboveThr,1)*1e-3);

dd = reshape(dd,sz_dd);
Vc = V;
Vc.fname = spm_file(V.fname,'prefix','t');
Vc = spm_create_vol(Vc);
Vc = spm_write_vol(Vc,dd);
fn_out = Vc.fname;


end

%% STEP 1: Removing small lesion patches from mask
function [fn_tMsk,fn_dtMsk] = mask_trimNgrow(P_in,minNr,nDilate)
% 1) Trim a mask image by removing bits that would be to small to really
%    matter according to medical criteria (cf. E. Lommers):
%    "Lesions will ordinarily be larger than 3 mm in cross section"
%    With 1x1x1mm^3 voxels, a cube of 2x2x2 voxels has a diagonal of
%    sqrt(12)~3.4mm and counts 8 voxels -> minNr = 8 [DEF]
%   -> fn_tMsk used for the new TPM_ms
% 2) Then grow the volume by 2 voxels [DEF]
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
XYZvx = V.mat\[XYZ ; ones(1,size(XYZ,2))];

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

%% STEP 3: Creating the normalization batch for the masked structural image
function [matlabbatch,fn_ICV] = batch_normalize_smooth(fn_kRef,fn_tMsk,fn_TPM,smoKern)
% [matlabbatch,fn_ICV] = batch_normalize_smooth(fn_kRef,fn_tMsk,fn_TPM,smoKern)
% This includes:
% - segmentation of the masked structural
% - writing out + smoothing the normalized lesion mask
% - creating the ICV-mask from c1/c2/c3/lesion-mask
% - deleting temporary files
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
fn_ICV = spm_file(fn_kRef,'prefix','icv_');

matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'LesionMask';
matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {{fn_tMsk}};
matlabbatch{2}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'MaskedRefStruct';
matlabbatch{2}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {{fn_kRef}};
matlabbatch{3}.spm.spatial.preproc.channel.vols(1) = cfg_dep('Named File Selector: MaskedRefStruct(1) - Files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
matlabbatch{3}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{3}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{3}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(1).tpm = {spm_file(fn_TPM,'number',1)};
matlabbatch{3}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{3}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(2).tpm = {spm_file(fn_TPM,'number',2)};
matlabbatch{3}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{3}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(3).tpm = {spm_file(fn_TPM,'number',3)};
matlabbatch{3}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{3}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(4).tpm = {spm_file(fn_TPM,'number',4)};
matlabbatch{3}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{3}.spm.spatial.preproc.tissue(4).native = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(5).tpm = {spm_file(fn_TPM,'number',5)};
matlabbatch{3}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{3}.spm.spatial.preproc.tissue(5).native = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(6).tpm = {spm_file(fn_TPM,'number',6)};
matlabbatch{3}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{3}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{3}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{3}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{3}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{3}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{3}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{3}.spm.spatial.preproc.warp.write = [0 1];
matlabbatch{4}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Named File Selector: LesionMask(1) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
matlabbatch{4}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72 ; 90 90 108];
matlabbatch{4}.spm.spatial.normalise.write.woptions.vox = [1.5 1.5 1.5];
matlabbatch{4}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{5}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{5}.spm.spatial.smooth.fwhm = smoKern*[1 1 1];
matlabbatch{5}.spm.spatial.smooth.dtype = 16;
matlabbatch{5}.spm.spatial.smooth.im = 0;
matlabbatch{5}.spm.spatial.smooth.prefix = 's';
matlabbatch{6}.spm.util.imcalc.input(1) = cfg_dep('Named File Selector: LesionMask(1) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
matlabbatch{6}.spm.util.imcalc.input(2) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
matlabbatch{6}.spm.util.imcalc.input(3) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
matlabbatch{6}.spm.util.imcalc.input(4) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
matlabbatch{6}.spm.util.imcalc.output = 'tmp.nii';
matlabbatch{6}.spm.util.imcalc.outdir = {pth_img};
matlabbatch{6}.spm.util.imcalc.expression = 'sum(X)';
matlabbatch{6}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{6}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{6}.spm.util.imcalc.options.mask = 0;
matlabbatch{6}.spm.util.imcalc.options.interp = 1;
matlabbatch{6}.spm.util.imcalc.options.dtype = 2;
matlabbatch{7}.spm.spatial.smooth.data(1) = cfg_dep('Image Calculator: ImCalc Computed Image: tmp.nii', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{7}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{7}.spm.spatial.smooth.dtype = 0;
matlabbatch{7}.spm.spatial.smooth.im = 0;
matlabbatch{7}.spm.spatial.smooth.prefix = 's';
matlabbatch{8}.spm.util.imcalc.input(1) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{8}.spm.util.imcalc.output = fn_ICV;
matlabbatch{8}.spm.util.imcalc.outdir = {pth_img};
matlabbatch{8}.spm.util.imcalc.expression = 'i1>.3';
matlabbatch{8}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{8}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{8}.spm.util.imcalc.options.mask = 0;
matlabbatch{8}.spm.util.imcalc.options.interp = 1;
matlabbatch{8}.spm.util.imcalc.options.dtype = 2;
matlabbatch{9}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('Segment: Seg Params', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','param', '()',{':'}));
matlabbatch{9}.cfg_basicio.file_dir.file_ops.file_move.files(2) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
matlabbatch{9}.cfg_basicio.file_dir.file_ops.file_move.files(3) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
matlabbatch{9}.cfg_basicio.file_dir.file_ops.file_move.files(4) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
matlabbatch{9}.cfg_basicio.file_dir.file_ops.file_move.files(5) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{9}.cfg_basicio.file_dir.file_ops.file_move.files(6) = cfg_dep('Image Calculator: ImCalc Computed Image: tmp.nii', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{9}.cfg_basicio.file_dir.file_ops.file_move.files(7) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{9}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;

end

%% STEP 4: Updating the TPM with a 7th class, the lesion
% Note that the lesion is inserted in *3rd position*, between WM and CSF!
function fn_TPMl = update_TPM_with_lesion(opt, fn_swtMsk)
% fn_TPMl = update_TPM_with_lesion(opt, fn_swtMsk)
%
% INPUT
% - opt : structure with a few parameters
%     .tpm4lesion : tissues to be modified for lesion (0/1/2/3) for
%                       GM / WM / GM+WM / GM+WM+CSF
%     .fn_tpm : tpm file name
%     .tpm_ratio : ration between WM and lesion
%     .min_tpm_icv : minimum value in intracranial volume
%     .min_tpm : minum value overall
% - fn_swtMsk : filename of smoothed normalized cleaned lesion mask, to be
%               used to create the lesion tissue class

% 0) select TPM and load
[pth,fnam,ext,num] = spm_fileparts(opt.fn_tpm);
fn_TPM = fullfile(pth,[fnam,ext]); % ensuring I load all 6 TPMs together.
Vtpm     = spm_vol(fn_TPM);
tpm_orig = spm_read_vols(Vtpm);
tpm_GM   = squeeze(tpm_orig(:,:,:,1));
tpm_WM   = squeeze(tpm_orig(:,:,:,2)); % used later on to define ICV
tpm_CSF  = squeeze(tpm_orig(:,:,:,3));
switch opt.tpm4lesion % Read in the healthy tissue prob map.
    case 0 % GM only
        tpm_healthy = tpm_GM;
    case 1 % WM only
        tpm_healthy = tpm_WM;
    case 2 % WM+GM
        tpm_healthy = tpm_GM+tpm_WM;
    case 3 % WM+GM+CSF
        tpm_healthy = tpm_GM+tpm_WM+tpm_CSF;
    otherwise
        error('Wrong tissue flag');
end
Vl = spm_vol(fn_swtMsk);
tpm_l = spm_read_vols(Vl);

% 1) scale lesion tpm and adjust healthy tissue prob map in ICV
% 2) ensure minium value all over
% 3) concatenate by setting lesion at #7 & adjust 'other' class
tpm_Lu = (1-1/opt.tpm_ratio)*tpm_l.*tpm_healthy; % update lesion tpm
tpm_Lu(tpm_WM>=opt.min_tpm_icv & tpm_Lu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
tpm_ext = cat(4,tpm_orig,tpm_Lu);
switch opt.tpm4lesion % update healthy tissues
    case 0 % GM only
        tpm_GMu = tpm_healthy - tpm_Lu;
        % equiv. to tpm_GMu = tpm_GM .* (1 - (1-1/opt.tpm_ratio) * tpm_l);
        tpm_GMu(tpm_GMu<opt.min_tpm) = opt.min_tpm;
        tpm_GMu(tpm_WM>=opt.min_tpm_icv & tpm_GMu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
        tpm_ext(:,:,:,1) = tpm_GMu; % update GM
    case 1 % WM only
        tpm_WMu = tpm_healthy - tpm_Lu;
        % equiv. to tpm_WMu = tpm_WM .* (1 - (1-1/opt.tpm_ratio) * tpm_l);
        tpm_WMu(tpm_WMu<opt.min_tpm) = opt.min_tpm;
        tpm_WMu(tpm_WM>=opt.min_tpm_icv & tpm_WMu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
        tpm_ext(:,:,:,2) = tpm_WMu; % update WM
    case 2 % WM+GM
        tpm_WMu = tpm_WM .* (1 - (1-1/opt.tpm_ratio) * tpm_l);
        tpm_WMu(tpm_WMu<opt.min_tpm) = opt.min_tpm;
        tpm_WMu(tpm_WM>=opt.min_tpm_icv & tpm_WMu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
        tpm_GMu = tpm_GM .* (1 - (1-1/opt.tpm_ratio) * tpm_l);
        tpm_GMu(tpm_GMu<opt.min_tpm) = opt.min_tpm;
        tpm_GMu(tpm_WM>=opt.min_tpm_icv & tpm_GMu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
        tpm_ext(:,:,:,1) = tpm_GMu; % update GM
        tpm_ext(:,:,:,2) = tpm_WMu; % update WM
    case 3 % WM+GM+CSF
        tpm_WMu = tpm_WM .* (1 - (1-1/opt.tpm_ratio) * tpm_l);
        tpm_WMu(tpm_WMu<opt.min_tpm) = opt.min_tpm;
        tpm_WMu(tpm_WM>=opt.min_tpm_icv & tpm_WMu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
        tpm_GMu = tpm_GM .* (1 - (1-1/opt.tpm_ratio) * tpm_l);
        tpm_GMu(tpm_GMu<opt.min_tpm) = opt.min_tpm;
        tpm_GMu(tpm_WM>=opt.min_tpm_icv & tpm_GMu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
        tpm_CSFu = tpm_CSF .* (1 - (1-1/opt.tpm_ratio) * tpm_l);
        tpm_CSFu(tpm_CSFu<opt.min_tpm) = opt.min_tpm;
        tpm_CSFu(tpm_WM>=opt.min_tpm_icv & tpm_CSFu<opt.min_tpm_icv) = opt.min_tpm_icv; % at least min_tpm_icv in ICV
        tpm_ext(:,:,:,1) = tpm_GMu; % update GM
        tpm_ext(:,:,:,2) = tpm_WMu; % update WM
        tpm_ext(:,:,:,3) = tpm_CSFu; % update CSF
    otherwise
        error('Wrong tissue flag');
end
tpm_ext(:,:,:,6) = 1 - sum(tpm_ext(:,:,:,[1:5 7]),4); % update 'other'

% 4) save the TPMl, with lesion in #3, in subject's data directory.
fn_TPMl = fullfile(spm_file(fn_swtMsk,'path'), ...
    spm_file(spm_file(fn_TPM,'filename'),'suffix','_les'));
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
nGauss = [2 2 2 2 3 4 4]; % Note: GM & WM 1->2,  Lesion -> 2, other 2->4.
cr_native = [1 1 ; 1 1 ; 1 1 ; 1 0 ; 1 0 ; 1 0 ; 0 0 ];
cr_warped = [1 1 ; 1 1 ; 1 1 ; 1 1 ; 0 0 ; 0 0 ; 0 0 ];
for ii = 1:7
    matlabbatch{1}.spm.spatial.preproc.tissue(ii).tpm = {[Ptpm_l,',',num2str(ii)]};
    matlabbatch{1}.spm.spatial.preproc.tissue(ii).ngaus = nGauss(ii);
    matlabbatch{1}.spm.spatial.preproc.tissue(ii).native = cr_native(ii,:);
    matlabbatch{1}.spm.spatial.preproc.tissue(ii).warped = cr_warped(ii,:);
end
% Define other parameters
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
% matlabbatch{1}.spm.spatial.preproc.warp.mrf = 0;
% matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 0;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
% Smoothing a bit
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

matlabbatch{1}.spm.spatial.normalise.write.subj.def = {fn_warp};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(char(fn_img2warp));
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
    78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;

end
