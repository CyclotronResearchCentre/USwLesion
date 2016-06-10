%% Script to generate some MPM specific TPM
%
% WHY:
% From some data check with several (elderly) subjects, the pallidum shows 
% very different voxel intensities from those of the rest of the GM.
% This is particularly visible in the A (much darker voxels) and R2s (much
% brighter voxels) images. There is no obvious difference in the R1 neither
% MT imges.
% This change in intensities is most probably linked to normal ageing as
% observed in the Calaghan et al., 2014, paper. This seems to be the main
% reason why segmentation is usually performed on a MT image as the
% intensities do not seem to be (too affected) by ageing...
% In practice though, when attempting a multi-channel segmentation (with 
% lesion or not), i.e. using the 4 images (R1/MT/A/R2s) for segmentation,
% then segmentation will *fail* for the pallidum area. Those voxel end up
% being classified as skin or anything but GM.
% 
% HOW:
% One solution is to improve the unified-segmentation (US) model. In this 
% case, one could include a specific tissue class for the pallidum in order
% to catch its very specific image intensities. Given the US flexibility at
% handling several channels and tissue classes, that should work smoothly.
% The only requirement is to update the TPM with a new prior tissue class 
% (for the pallidum) and to update the other priors, moslty the GM one.
% There is one atlas distributed with SPM12, let's use it!

%% DEFINE a few filenames
fn_TPM = 'unwTPM_sl2.nii';
fn_atlas = 'labels_Neuromorphometrics.nii';
fn_pallidum = 'b_pallidum.nii';
dr_TPM = fullfile(spm('dir'),'tpm');
dr_SB = fullfile(spm('dir'),'toolbox','USwLesion','Script_and_batches');

%% GET TPMs
fn_TPMs = spm_select('ExtFPList',dr_TPM,fn_TPM,1:10);
nTPMs = size(fn_TPMs,1);
Vtpm = spm_vol(fn_TPMs);

val_tpm = spm_read_vols(Vtpm);
SZ = size(val_tpm);

%% Some checks
vv = reshape(val_tpm,[prod(SZ(1:3)) SZ(4)])';
sc_tpm = sum(vv);
vv = vv./(ones(SZ(4),1)*sc_tpm);
val_tpm = reshape(vv',SZ);

% fplot(sum(vv))
% sum(sum(vv)>1)

%% EXTRACT Pallidum from atlas, binarized and smoothed (4mm FWHM)

% The L/R pallidum are coded as #55 and #56 in the atlas
clear matlabbatch
matlabbatch{1}.spm.util.imcalc.input = {fullfile(dr_TPM,fn_atlas)};
matlabbatch{1}.spm.util.imcalc.output = fn_pallidum;
matlabbatch{1}.spm.util.imcalc.outdir = {dr_SB};
matlabbatch{1}.spm.util.imcalc.expression = '(i1==55) + (i1==56)';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 2;
matlabbatch{2}.spm.spatial.smooth.data(1) = ...
    cfg_dep('Image Calculator: ImCalc Computed Image: b_pallidum.nii', ...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, ...
                      '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{2}.spm.spatial.smooth.fwhm = [4 4 4];
matlabbatch{2}.spm.spatial.smooth.dtype = 16;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's';
spm_jobman('run', matlabbatch);

% load pallidum tpm
fn_pallidum = ['s',fn_pallidum];
Vpal = spm_vol(fullfile(dr_SB,fn_pallidum));
val_pal = spm_read_vols(Vpal);

%% Update TPM by introducing a new tpm between GM and WM
min_tpm = crc_USwL_get_defaults('uTPM.min_tpm');
min_tpm_icv = crc_USwL_get_defaults('uTPM.min_tpm_icv');
vval_tpm = reshape(val_tpm , [prod(SZ(1:3)) SZ(4)]) - min_tpm;
vval_utpm = zeros(prod(SZ(1:3)),SZ(4)+1);

vval_utpm(:,2) = vval_tpm(:,1).*val_pal(:) - min_tpm_icv; % Pallidum
vval_utpm(:,1) = vval_tpm(:,1) - vval_utpm(:,2); % GM
vval_utpm(:,3:SZ(4)+1) = vval_tpm(:,2:SZ(4)); % rest
vval_utpm = vval_utpm*(1+SZ(4)*min_tpm)/(1+(SZ(4)+1)*min_tpm)+ min_tpm ; % non-zero everywhere
val_utpm = reshape(vval_utpm,[SZ(1:3) SZ(4)+1]);

% Save into file
fn_TPMu = fullfile(dr_SB, spm_file(fn_TPM,'suffix','_uMPM'));
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


% zz = sum(vval_utpm,2)';
% fplot(zz)
% zz = sum(vval_tpm,2)';
% fplot(zz)

