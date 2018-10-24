%% Script to generate some masks for brain splitting
% 
% WHY:
% For the MPMs analysis of MS patients, it's not unusual to look at
% mean/std values over sub-parts of the GM or WM, e.g. the cortex.
% Here we just split the brain GM into a few parts, see list here under.
%
% HOW
% We simply use the neuromorphometric atlas provided in SPM12 and combine
% the various 'patches' together into those 6 different masks.
% 
% OUTPUT
% A single file is generated, msk_BrainParts.nii, with 6 different volumes 
% covering the brain: 
% 1. cortex without WM, 
% 2. cortex with WM, 
% 3. subcortex (Basal Ganglia, Amygdala, brainstem),
% 4. cerebellum.
% 5. subcortex without the Basal Ganglia
% 6. Basal Ganglia only (pallidum, caudate, putamen, nucleaus accumbens,
%   thalamus proper)
% 
% NOTE
% For an individual subject, the masks generated likely include *both* GM
% and WM voxels!
% It is up to the user to also provide a segmented GM or WM map, e.g. 
% wc1*.nii and wc2*.nii images, in order to specific voxel values in the
% brain part(s) of interest.
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

%% DEFINE a few filenames & parameters
fn_atlas = 'labels_Neuromorphometrics.nii';
fn_parts = 'msk_BrainParts.nii';
dr_TPM = fullfile(spm('dir'),'tpm'); % SPM's tpm folder
dr_TPMuswl = fullfile(spm_file(which('tbx_cfg_USwLesion.m'),'path'), ...
    'eTPM'); % USwL's eTPM folder
mask_smoothing = ones(1,3)*4; % kernel size (in mm) for cortex & cerebellum
mask_smoothing2 = mask_smoothing*2; % twice for subcortical
sval_thr = .1; % threshold for smoothed volume
nVols = 6;

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

%% SUBCORTICAL part: Basal Ganglial, Amygdala, brainstem,...
l_rois = [ 23 30 31 32 35 36 37 47 48 55 56 57 58 59 60 61 62 75 76 ];
val_subc = zeros(SZ);
sval_subc = zeros(SZ);
for ii=l_rois
    val_subc = val_subc +  double(val_atl == ii);
end
spm_smooth(val_subc,sval_subc,mask_smoothing2./vx_sz); % extend a bit by smoothing

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

%% SUBCORTICAL part, without Basal Ganglia
l_rois = [ 31 32 35 47 48 61 62 75 76 ];
val_subc_woBG = zeros(SZ);
sval_subc_woBG = zeros(SZ);
for ii=l_rois
    val_subc_woBG = val_subc_woBG +  double(val_atl == ii);
end
spm_smooth(val_subc_woBG,sval_subc_woBG,mask_smoothing2./vx_sz); % extend a bit by smoothing

%% SUBCORTICAL part, Basal Ganglia only
% BG area includes these bits
% - pallidum, #55 and #56
% - caudate, #36 and #37
% - putamen, #57 and #58
% - nucleaus accumbens, #23 and #30
% - thamus proper, #59 and #60
l_rois = [55 56 36 37 57 58 23 30 59 60];
val_subc_BG = zeros(SZ);
sval_subc_BG = zeros(SZ);
for ii=l_rois
    val_subc_BG = val_subc_BG +  double(val_atl == ii);
end
spm_smooth(val_subc_BG,sval_subc_BG,mask_smoothing2./vx_sz); % extend a bit by smoothing

%% Build the 4D volume
% Keep a voxel in a mask if sval>.1 and larger than the other 2
val_parts  = zeros([SZ nVols]);
val_parts(:,:,:,1) = sval_cort>sval_thr & ...
    sval_cort>sval_cereb & sval_cort>sval_subc; % -> cortical part
val_parts(:,:,:,3) = sval_subc>sval_thr & ...
    sval_subc>=sval_cereb & sval_subc>=sval_cort; % -> sub-cortical part
val_parts(:,:,:,4) = sval_cereb>sval_thr & ...
    sval_cereb>=sval_cort & sval_cereb>=sval_subc; % -> cerebellum part

% Keep a voxel if sval>.1 and not part of cerebelum or sub-cortex
val_parts(:,:,:,2) = sval_cort_wWM>sval_thr & ...
    ~val_parts(:,:,:,3) & ~val_parts(:,:,:,4); % -> cortical part with WM

% Keep a voxel if sval>.1 and in sub-cortex but not in the other one
val_parts(:,:,:,5) = sval_subc_woBG>sval_thr & ...
    sval_subc_woBG>=sval_subc_BG & val_parts(:,:,:,3); % -> sub-cortical wo BG
val_parts(:,:,:,6) = sval_subc_BG>sval_thr & ...
    sval_subc_BG>=sval_subc_woBG & val_parts(:,:,:,3); % -> BG only

%% save the 4D volume
fn_mask = fullfile(dr_TPMuswl,fn_parts);
% assuming that Vatlas is a uint8 file -> 1 voxel = 1 byte
descrip = {'Cotical GM mask', 'Cotical GM+WM mask', ...
            'Subcortical GM mask','Cerebellar GM mask', ...
            'Subcortical wo BG GM mask', 'Basal Ganglia GM mask'};
for ii=1:nVols
    V_mask(ii) = Vatlas; %#ok<*SAGROW>
    V_mask(ii).pinfo(1) = 1/255;
    V_mask(ii).pinfo(3) = V_mask(ii).pinfo(3) + (ii-1)*prod(SZ);
    V_mask(ii).fname = fn_mask;
    V_mask(ii).n(1) = ii;
    V_mask(ii).descrip = descrip{ii};
    V_mask(ii) = spm_create_vol(V_mask(ii));
    V_mask(ii) = spm_write_vol(V_mask(ii),val_parts(:,:,:,ii));
end
