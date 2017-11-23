function fn_out = crc_ExtractParam_MPMs(fn_img,opt)
% Extracting parameters from the segmented images.
%
% Which parameters should extracted?
% Here are a few things done so far:
%   1. tICV -> reference volume
%   2. match between mask and segmented lesion volume
%   3. See list here under
%   4. MPMvalues in the follwoing tissue types: cortical GM, deep GM (no
%      cerebellum nor brainstem), cortical NA, lesion
%      Note: cortical = brain - (central, cerebellum, brainstem), as
%      defined from the brainparts mask image.
%   5. some summary statistics: min, max , mean, median, std, skewness,
%      kurtosis, p10, & p90 of the values extracted at 4., for each MPM and
%      tissue types.
%
% INPUT
% - fn_img : structure with necessary filenames
%       .fn_cImg   : tissue classes (GM, WM, Les, CSF) or (GM, WM, CSF)
%       .fn_MPM    : quantitative MRIs (typically A, MT, R1, R2s)
%       .fn_seg8   : warping & segmentation parameters 
%       .fn_BrPart : brain parts mask for value estraction
%       .fn_MskLes : original lesion mask [optional]
% - opt    : structure with some options
%       .outdir    : output directory
%       .thrICV    : threshold for the creation of the tICV mask [def .5]
%       .thrLesion : lesion map threshold, for match with lesion mask
%       .thrTC     :
%       .lBrParts  :
%
% OUTPUT
% - fn_out  : filename of output matlab file with structure containing the
%             following fields
%
%
% REFS:
% tICV : http://www.sciencedirect.com/science/article/pii/S1053811914007769
%__________________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by C. Phillips, 2015.
% Cyclotron Research Centre, University of Liege, Belgium

fn_cImg = fn_img.fn_cImg;
fn_MPM  = fn_img.fn_MPM;
nMPM = size(fn_MPM,1); % 4 -> [GM, WM, Lesion, CSF] or 3 -> [GM, WM, CSF]
fn_seg8 = fn_img.fn_seg8;
fn_iwarp = fn_img.fn_iwarp;

if isfield(fn_img,'fn_MskLes')
    fn_MskLes  = fn_img.fn_MskLes;
    MskLes_match = true;
else
    fn_MskLes  = '';
    MskLes_match = false;
end

[pth,fn] = spm_fileparts(fn_MPM(1,:));
if isempty(opt.outdir)
    pth_out = pth;
else
    pth_out = opt.outdir;
end

%% 1. build ICV from the sum of GM, WM, CSF and lesion
fn_icv = spm_file(spm_file(fn_cImg(1,:),'filename'),'prefix','tICV_');
matlabbatch = create_mask(fn_cImg(1:end,:),fn_icv,pth,opt.thrICV);
spm_jobman('run', matlabbatch);
res.Picv = fn_icv;

%% 2. match between lesion mask and segmented tissue
if MskLes_match
    % Dice coefficient: http://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient
    % Jaccard index: http://en.wikipedia.org/wiki/Jaccard_index
    
    % In native space
    Vtmsk = spm_vol(fn_MskLes);
    Vc3 = spm_vol(fn_cImg(3,:));
    v_tmsk = spm_read_vols(Vtmsk);
    v_tmsk = v_tmsk(:) >.5; % thresholding at .5, as it should be a binary img
    v_c3 = spm_read_vols(Vc3);
    v_c3 = v_c3(:) > opt.thrLesion;
    
    N_tmsk = sum(v_tmsk);
    N_c3 = sum(v_c3);
    N_intersec = sum(v_tmsk.*v_c3);
    N_union = sum( (v_tmsk+v_c3)>0 );
    
    % Dice coefficient
    DC = 2*N_intersec / (N_tmsk + N_c3) ;
    % Jaccard index
    JI = N_intersec/N_union ;
    
    res.match = struct( ...
        'DiceC', DC, ...
        'JaccardI', JI, ...
        'N_tmsk', N_tmsk, ...
        'N_c3', N_c3, ...
        'N_intersec', N_intersec, ...
        'N_union', N_union);
end

%% 3. Volumes
% Volumes or fractions to extract:
% - TIV
% - GMV = Cortical and Deep Grey matter
% - WMV = NAWM + lesions
% - Brain Parenchymal Fraction= (GMV+WMV)/TIV
% - Grey Matter Fraction = GMV/TIV
% - White Matter Fraction =  WMV/TIV
% - Lesion Fraction = Lesion volume/TIV (more standard than relative than
%   fraction of WMV)
% - Lesion Fraction WM = Lesion volume/WMV

clear matlabbatch
% fn_seg8 = spm_select('FPList',pth,'^kt.*_seg8\.mat$');
matlabbatch{1}.spm.util.tvol.matfiles = {fn_seg8};
matlabbatch{1}.spm.util.tvol.tmax = nMPM;
matlabbatch{1}.spm.util.tvol.mask = {fullfile(spm('dir'),'tpm','mask_ICV.nii,1')};
matlabbatch{1}.spm.util.tvol.outf = fullfile(pth_out,['TV_',fn]);
spm_jobman('run', matlabbatch);
tmp = load(fn_seg8);
volumes = tmp.volumes;
tiv = sum(volumes.litres);
gmv = volumes.litres(1);
if nMPM ==4
    wmv = sum(volumes.litres([2 3])); % WM + lesion volume!
    lesv = volumes.litres(3);
else
    wmv = volumes.litres(2);
end

BrParenchFrac = (gmv+wmv)/tiv*100; % expressed in percentage
GMFrac = gmv/tiv*100; % expressed in percentage
WMFrac = wmv/tiv*100; % expressed in percentage
if nMPM ==4
    LesFrac = volumes.litres(3)/tiv*100; % expressed in percentage
    LesFracWM = volumes.litres(3)/wmv*100; % expressed in percentage
end
      
lesionVol = struct(...
    'volumes', volumes.litres, ...
    'tiv', tiv, ...
    'gmv', gmv, ...
    'wmv', wmv, ...
    'BrParenchFrac', BrParenchFrac, ...
    'GMFrac', GMFrac, ...
    'WMFrac', WMFrac);
if nMPM ==4
    lesionVol.lesv = lesv;
    lesionVol.LesFrac = LesFrac;
    lesionVol.LesFracWM = LesFracWM;
end

% Store in main results structure
res.lesionVol = lesionVol;

%% 4. extraction of MPM values for the GM/WM(/lesion)
% Accounting for the brain parts selected.

% First copy Brain Parts mask into subject temporary subdir
fn_brP = fullfile(spm_file(mfilename('fullpath'),'path'), ...
    'eTPM','msk_BrainParts.nii');
dr_tmp = fullfile(pth,'tmp');
if ~exist(dr_tmp,'dir'), mkdir(dr_tmp); end
copyfile(fn_brP,dr_tmp);
nBrP = 6; % Use all six images of 4D volume
fn_brP_loc = '';
for ii=1:nBrP
    fn_brP_loc = char(fn_brP_loc, ...
        spm_file(fn_brP,'path',dr_tmp,'number',ii));

end
fn_brP_loc(1,:) = [];

% Create MatlabBatch
clear matlabbatch
matlabbatch = create_MB(fn_iwarp,fn_cImg(1,:),fn_brP_loc);
spm_jobman('run', matlabbatch);








nMPM = size(fn_MPM,1);
Vmpm = spm_vol(fn_MPM);
if ~isempty(fn_MPMmsk)
    nMsk = size(fn_MPMmsk,1);
    if nMsk ~= nMPM
        error('USwL:ExParam','Num MPM images (%d) ~= Num Msk images (%d)', nMPM,nMsk);
    end
    Vmsk = spm_vol(fn_MPMmsk);
else
    nMsk = 0;
end

Vtc = spm_vol(fn_cImg);
v_tc123 = spm_read_vols(Vtc(1:3));
vt_tc123 = v_tc123 > opt.thrTC ;

vMPM = cell(3,1); % #tissue classes x #MPM
if nMsk
    vMPMmsk = cell(3,1);
    v_msk = spm_read_vols(Vmsk);
    v_msk = ~sum(v_msk,4);
end
for ii=1:nMPM
    v_mpm = spm_read_vols(Vmpm(ii));
    v_mpm = v_mpm(:);
    for jj=1:3 % tc
        tmp = vt_tc123(:,:,:,jj);
        vMPM{jj}(:,ii) = v_mpm(tmp(:));
        if nMsk
            tmp = tmp & v_msk;
            vMPMmsk{jj}(:,ii) = v_mpm(tmp(:));
        end
    end
end
res.vMPM = vMPM; %#ok<*STRNU>
if nMsk
    res.vMPMmsk = vMPMmsk; %#ok<*STRNU>
end

%% 5. some stats from the MPM values
% Find min-max, mean, median, std, skewness ,kurtosis for each MPM
mM = zeros(nMPM,2); mM(:,1) = Inf; mM(:,2) = -Inf;
meanVal = zeros(3,nMPM); medVal = zeros(3,nMPM); stdVal = zeros(3,nMPM);
skewVal = zeros(3,nMPM); kurtVal = zeros(3,nMPM);
if nMsk
    v_use = res.vMPMmsk;
else
    v_use = res.vMPM;
end
for ii=1:3 % tc
    for jj=1:nMPM % mpm
        % min/max total
        tmp_m = min(v_use{ii}(:,jj));
        if tmp_m<mM(jj,1), mM(jj,1) = tmp_m; end
        tmp_M = max(v_use{ii}(:,jj));
        if tmp_M>mM(jj,2), mM(jj,2) = tmp_M; end
        % mean/meadian/std/skewness/kurtosis
        meanVal(ii,jj) = mean(v_use{ii}(:,jj));
        medVal(ii,jj) = median(v_use{ii}(:,jj));
        stdVal(ii,jj) = std(v_use{ii}(:,jj));
        skewVal(ii,jj) = skewness(v_use{ii}(:,jj));
        kurtVal(ii,jj) = kurtosis(v_use{ii}(:,jj))-3;
    end
end
res.vMPMstats = struct('mM',mM,'meanVal',meanVal,'medVal',medVal,...
    'stdVal',stdVal,'skewVal',skewVal,'kurtVal',kurtVal);

%% 6. save things and pass out fn_out
fn_out.fn_ExParam = fullfile(pth_out,['ExP_',fn]);
save(fn_out.fn_ExParam,'res');


end

% =======================================================================
%% SUBFUNCTION
% =======================================================================
function matlabbatch = create_mask(fn_in,fn_out,pth,thr,smK)
% Batch to create a mask from a few segmented tissue classes.
% It works in 3 steps:
% 1. sum up the input images
% 2. smooth the sum
% 3. threshold the result

if nargin<5, smK = 8; end
if nargin<4, thr = .5; end

matlabbatch{1}.spm.util.imcalc.input = cellstr(fn_in);
matlabbatch{1}.spm.util.imcalc.output = 'tmp.nii';
matlabbatch{1}.spm.util.imcalc.outdir = {pth};
matlabbatch{1}.spm.util.imcalc.expression = 'sum(X)';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 2;
matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('Image Calculator: ImCalc Computed Image: output', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{2}.spm.spatial.smooth.fwhm = [1 1 1]*smK;
matlabbatch{2}.spm.spatial.smooth.dtype = 0;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's';
matlabbatch{3}.spm.util.imcalc.input(1) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{3}.spm.util.imcalc.output = fn_out;
matlabbatch{3}.spm.util.imcalc.outdir = {pth};
matlabbatch{3}.spm.util.imcalc.expression = sprintf('i1>%f',thr);
matlabbatch{3}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{3}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{3}.spm.util.imcalc.options.mask = 0;
matlabbatch{3}.spm.util.imcalc.options.interp = 1;
matlabbatch{3}.spm.util.imcalc.options.dtype = 2;
matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('Image Calculator: ImCalc Computed Image: output', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_move.files(2) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;

end

% =======================================================================
function matlabbatch = create_MB(fn_iwarp,fn_subj,fn_brP)
% Building matlabbatch for the normalize-write operation, muche easier than
% building the deformation and applying it manually.
% Bring SPM-ICV into subject space -> need to properly define the latter!

prec_round = 1e6;

V_icv = spm_vol(fn_subj);
% voxel size rounded to some precision
vx_sz = round(sqrt(sum(V_icv.mat(1:3,1:3).^2))*prec_round)/prec_round;
% defining BB
p_min = -V_icv.mat\[0 0 0 1]' ;
p_max = (V_icv.dim.*vx_sz)' + p_min(1:3) ;
img_bb = [-abs(p_min(1:3)') ; abs(p_max(1:3)')];

% % Bring in SPM-ICV into subject space
% fn_icvSPM = fullfile(spm('dir'),'tpm','mask_ICV.nii');
% fn_icvSPM_loc = fullfile(spm_file(fn_iwarp,'path'),'icv_SPM.nii');
% copyfile(fn_icvSPM,fn_icvSPM_loc)

matlabbatch{1}.spm.spatial.normalise.write.subj.def(1) = {spm_file(fn_iwarp,'number','')};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {fn_brP};
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = img_bb;
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = vx_sz;
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;

end
