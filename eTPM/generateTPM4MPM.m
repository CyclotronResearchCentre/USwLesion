%% Script to generate some MPM specific TPM
%
% WHY:
% From some data check with several (elderly) subjects, the basal ganglia 
% (BG) surrounding shows very different voxel intensities from those of the
% rest of the GM (due to % iron deposit with age). This is particularly 
% visible in the A (much darker voxels) and R2s (much brighter voxels) 
% images. There is no obvious difference in the R1 neither MT imges.
% This change in intensities is most probably linked to normal ageing as
% observed in the Callaghan et al., 2014, paper. This seems to be the main
% reason why segmentation is usually performed on a MT image as the
% intensities do not seem to be (too affected) by ageing...
% In practice though, when attempting a multi-channel segmentation (with 
% lesion or not), i.e. using the 4 images (R1/MT/A/R2s) for segmentation,
% then segmentation will *fail* for the "pallidum area". Those voxel end up
% being classified as skin or anything but GM.
% In practice Callaghan et al. report the following areas to be affected by
% ageing in term of R2* values: pallidum, caudate, putamen, and thalami. 
% This approximately covers the basal ganglia (BG) area.
% 
% HOW:
% One solution is to improve the unified-segmentation (US) model. In this 
% case, one could include a specific tissue class for the BG in order to 
% catch its very specific image intensities. Given the US flexibility at
% handling several channels and tissue classes, that should work smoothly.
% The only requirement is to update the TPM with a new prior tissue class 
% (for the BG area ) and to update the other priors, moslty the GM one.
% There is one atlas distributed with SPM12, let's use it!
% 
% In practice this messes up a bit the general flow of images as the GM is
% now "split" into 2 different images. Moreover some part of code are
% hard-coded with the first 3 segmented images to be GM-WM-CSF. One way out
% of this is to include the BG class at the end (7th position), then after 
% the segmentation c1 (GM minus BG) and c7 (BG only) images can
% be added together to produce a "new" c1 image with the full GM posterior 
% probability map.
% 
% NOTE:
% The same kind of problem occurs when a lesion class is introduced between
% the WM and CSF, so it is important to check any SPM code dealing with
% segmented tissue classes (e.g. cleaning!)
% 
% REFERENCE:
% Callaghan et al, 2014, Widespread age-related differences in the human 
% brain microstructure revealed by quantitative magnetic resonance imaging.
% Neurobiology of Aging, 35:1862-1872
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

%% DEFINE a few filenames & parameters
fn_TPM = 'eTPM.nii';
fn_atlas = 'labels_Neuromorphometrics.nii';
% fn_bg = 'msk_basalganglia.nii';
fn_brainp = 'msk_BrainParts.nii';
% where #2 = cortical GM+WM, #3 = subcortex, #4 = cerebellum, #5 =
% subcortex without BG, #6 = Basal Ganglia only.
dr_TPM = fullfile(spm('dir'),'tpm');
dr_TPMuswl = fullfile(spm_file(which('tbx_cfg_USwLesion.m'),'path'),'eTPM');
mask_smoothing = [1 1 1]*4; % smoothing (in mm of FWHM) for the BG mask

%% GET TPMs
fn_TPMs = spm_select('ExtFPList',dr_TPM,fn_TPM,1:10);
nTPMs = size(fn_TPMs,1);
Vtpm = spm_vol(fn_TPMs);

val_tpm = spm_read_vols(Vtpm);
SZ = size(val_tpm);
vx_sz = sqrt(sum(Vtpm(1).mat(1:3,1:3).^2));

%% Some checks & scaling
% vv = reshape(val_tpm,[prod(SZ(1:3)) SZ(4)])';
% sc_tpm = sum(vv);
% vv = vv./(ones(SZ(4),1)*sc_tpm);
% val_tpm = reshape(vv',SZ);
% 
% fplot(sum(vv)-1)
% sum(sum(vv)>1)

%% GET Basal Ganglia part + a little bit of smoothing
fn_brp = fullfile(dr_TPMuswl,fn_brainp);
Vbgm = spm_vol(spm_file(fn_brp,'number',6)); % basal ganglia mask
val_bgm = spm_read_vols(Vbgm);
sval_bgm = zeros(SZ(1:3));
spm_smooth(val_bgm,sval_bgm,mask_smoothing./vx_sz); % extend a bit by smoothing

%% Update TPM by introducing a new tpm between GM and WM
min_tpm = crc_USwL_get_defaults('uTPM.min_tpm');
min_tpm_icv = crc_USwL_get_defaults('uTPM.min_tpm_icv');

vval_tpm = reshape(val_tpm , [prod(SZ(1:3)) SZ(4)]) - min_tpm;
vval_utpm = zeros(prod(SZ(1:3)),SZ(4)+1);

vval_utpm(:,7) = (vval_tpm(:,1)- min_tpm_icv).*sval_bgm(:) ; % BG mask
vval_utpm(:,1) = vval_tpm(:,1) - vval_utpm(:,7); % GM minus BG
vval_utpm(:,2:SZ(4)) = vval_tpm(:,2:SZ(4)); % rest

vval_utpm = vval_utpm*(1+SZ(4)*min_tpm)/(1+(SZ(4)+1)*min_tpm)+ min_tpm ; % non-zero everywhere
val_utpm = reshape(vval_utpm,[SZ(1:3) SZ(4)+1]);

%% Save into file
fn_TPMu = fullfile(dr_TPMuswl, spm_file(fn_TPM,'suffix','_wBG'));
Vtpm_u = Vtpm;
Vtpm_u(7) = Vtpm(6);
mem_sz = Vtpm(2).pinfo(3)-Vtpm(1).pinfo(3);
Vtpm_u(7).pinfo(3) = Vtpm_u(7).pinfo(3) + mem_sz;
Vtpm_u(7).n(1) = 7;
for ii=1:7
    Vtpm_u(ii).fname = fn_TPMu;
    Vtpm_u(ii) = spm_create_vol(Vtpm_u(ii));
    Vtpm_u(ii) = spm_write_vol(Vtpm_u(ii),val_utpm(:,:,:,ii));
end

% % Some checks
% zz = sum(vval_utpm,2)';
% fplot(zz-1)
% zz = sum(vval_tpm,2)';
% fplot(zz)
% 
% zz = min(vval_utpm,[],2)';
% fplot(zz)
% any(zz<0)