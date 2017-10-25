%% Script to generate some masks for brain split
% 
% WHY:
% For the MPMs analysis of MS patients, it's not unusual to look at
% mean/std values over parts of the GM , e.g. the cortex.
% Here we just split the brain GM into 3 main parts: cortex, cerrebellum
% and subcortical.
%
% HOW
% We simply use the neuromorphometric atals provided in SPM12 and combine
% the various 'patches' together into those 3 different masks.
% 
% OUTPUT
% A single file is generated, msk_BrainParts.nii, with 4 different volumes 
% covering the brain: 
% - cortex without WM, 
% - cortex with WM, 
% - subcortex (pons + basal ganglia), and
% - cerebellum.
% 
% NOTE
% The mask generated include *both* GM ande WM! It is up to the user to
% also provide a segmented GM or WM map, e.g. c1*.nii and c2*.nii, in order
% to decide which tissue type to work with.
% 
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

%% DEFINE a few filenames
fn_atlas = 'labels_Neuromorphometrics.nii';
fn_parts = 'msk_BrainParts.nii';
dr_TPM = fullfile(spm('dir'),'tpm'); % SPM's tpm folder
dr_TPMuswl = fullfile(spm_file(which('tbx_cfg_USwLesion.m'),'path'), ...
    'eTPM'); % USwL's eTPM folder
mask_smoothing = 4*ones(1,3); % in mm
sval_thr = .1; % threshold for smoothed volume

%% LOAD the atlas
Vatlas = spm_vol(fullfile(dr_TPM,fn_atlas));
vx_sz = sqrt(sum(Vatlas.mat(1:3,1:3).^2));

val_atl = spm_read_vols(Vatlas);
SZ = size(val_atl);

l_Rinds =  unique(val_atl(:));

%% CORTICAL part without WM
% all regions with ind>99 are cortical bits
val_cort = double(val_atl>99);
sval_cort = zeros(SZ);
spm_smooth(val_cort,sval_cort,mask_smoothing./vx_sz); % extend a bit by smoothing

%% CORTICAL part with WM
% all regions with ind>99 are cortical bits, #44&45 and the cortical WM
val_cort_wWM = double(val_atl>99) +  double(val_atl == 44)+  double(val_atl == 45) ;
sval_cort_wWM = zeros(SZ);
spm_smooth(val_cort_wWM,sval_cort_wWM,mask_smoothing./vx_sz); % extend a bit by smoothing

%% SUBCORTICAL part
l_rois = [ 23 30 31 32 35 36 37 47 48 55 56 57 58 59 60 61 62 75 76 ];
val_subc = zeros(SZ);
sval_subc = zeros(SZ);
for ii=l_rois
    val_subc = val_subc +  double(val_atl == ii);
end
spm_smooth(val_subc,sval_subc,mask_smoothing./vx_sz); % extend a bit by smoothing

%% CEREBELLUM part
% List of L/R regions and coresponding indexes in the atlas from SPM12
% - cerebellum exterior, #38 and #39
% - cerebellar vermal lobules, #71 72 73
% - cerebellar WM, #40 41
l_rois = [38 39 71 72 73 40 41];

val_cereb = zeros(SZ);
sval_cereb = zeros(SZ);
for ii=l_rois
    val_cereb = val_cereb +  double(val_atl==ii);
end
spm_smooth(val_cereb,sval_cereb,mask_smoothing./vx_sz); % extend a bit by smoothing

%% Build the 4D volume
% Keep a voxel in a mask if sval>.1 and larger than the other 2
val_parts  = zeros([SZ 3]);
val_parts(:,:,:,1) = sval_cort>sval_thr & ...
    sval_cort>sval_cereb & sval_cort>sval_subc; % -> cortical part
val_parts(:,:,:,3) = sval_subc>sval_thr & ...
    sval_subc>sval_cereb & sval_subc>sval_cort; % -> sub-cortical part
val_parts(:,:,:,4) = sval_cereb>sval_thr & ...
    sval_cereb>sval_cort & sval_cereb>sval_subc; % -> cerebellum part

val_parts(:,:,:,2) = sval_cort_wWM>sval_thr & ...
    ~val_parts(:,:,:,3) & ~val_parts(:,:,:,3); % -> cortical part with WM

%% save the 4D volume
fn_mask = fullfile(dr_TPMuswl,fn_parts);
% assuming that Vatlas is a uint8 file -> 1 voxel = 1 byte
descrip = {'Cotical GM mask', 'Cotical GM+WM mask', ...
            'Subcortical GM mask','Cerebellar GM mask'};
for ii=1:4
    V_mask(ii) = Vatlas; %#ok<*SAGROW>
    V_mask(ii).pinfo(1) = 1/255;
    V_mask(ii).pinfo(3) = V_mask(ii).pinfo(3) + (ii-1)*prod(SZ);
    V_mask(ii).fname = fn_mask;
    V_mask(ii).n(1) = ii;
    V_mask(ii).descrip = descrip{ii};
    V_mask(ii) = spm_create_vol(V_mask(ii));
    V_mask(ii) = spm_write_vol(V_mask(ii),val_parts(:,:,:,ii));
end


% %% Bits to create images -> easier to figure out what's in/out
% % Build function for imcalc function.
% func_imcalc = sprintf('(i1==%d)',l_rois(1));
% for ii=2:numel(l_rois)
%     func_imcalc = [func_imcalc, sprintf(' + (i1==%d)',l_rois(ii))]; %#ok<*AGROW>
% end
% % func_imcalc = 'i1>99';
% 
% clear matlabbatch
% matlabbatch{1}.spm.util.imcalc.input = {fullfile(dr_TPM,fn_atlas)};
% matlabbatch{1}.spm.util.imcalc.output = 'tmp.nii';
% matlabbatch{1}.spm.util.imcalc.outdir = {dr_TPMuswl};
% matlabbatch{1}.spm.util.imcalc.expression = func_imcalc;
% matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
% matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
% matlabbatch{1}.spm.util.imcalc.options.mask = 0;
% matlabbatch{1}.spm.util.imcalc.options.interp = 1;
% matlabbatch{1}.spm.util.imcalc.options.dtype = 2;
% % spm_jobman('run', matlabbatch);
% matlabbatch{2}.spm.spatial.smooth.data(1) = ...
%     cfg_dep(['Image Calculator: ImCalc Computed Image: ', 'tmp.nii'], ...
%             substruct('.','val', '{}',{1}, '.','val', '{}',{1}, ...
%                       '.','val', '{}',{1}), substruct('.','files'));
% matlabbatch{2}.spm.spatial.smooth.fwhm = [1 1 1]*mask_smoothing;
% matlabbatch{2}.spm.spatial.smooth.dtype = 16;
% matlabbatch{2}.spm.spatial.smooth.im = 0;
% matlabbatch{2}.spm.spatial.smooth.prefix = 's';
% spm_jobman('run', matlabbatch);
