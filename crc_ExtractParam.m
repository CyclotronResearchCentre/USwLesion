function fn_out = crc_ExtractParam(job)
% Extracting parameters from the segmented images.
%
% Which parameters should extracted?
% Here are a few things done so far:
%   1. tICV -> reference volume
%   2. match between mask and segmented lesion volume
%   3. %lesion in tICV, %lesion in WM volume
%   4. MPMvalues in GM, WM and lesion
%        -> check their histogram
% REFS:
% http://www.sciencedirect.com/science/article/pii/S1053811914007769
% Dice coefficient: http://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient
% Jaccard index: http://en.wikipedia.org/wiki/Jaccard_index
%__________________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by C. Phillips, 2015.
% Cyclotron Research Centre, University of Liege, Belgium

fn_cImg = char(job.cImg);
fn_MskLes  = char(job.imgMsk);
fn_MPM  = char(job.imgMPM);
fn_MPMmsk = char(job.imgMPMmsk);

[pth,fn] = spm_fileparts(fn_MPM(1,:));
if isempty(job.outdir{1})
    pth_out = pth;
else
    pth_out = job.outdir{1};
end
opt = job.opt;

%% 1. build ICV from the sum of GM, WM, CSF and lesion
fn_icv = spm_file(spm_file(fn_cImg(1,:),'filename'),'prefix','icv_');
matlabbatch = create_mask(fn_cImg(1:end,:),fn_icv,pth,opt.thrICV);
spm_jobman('run', matlabbatch);
res.Picv = fn_icv;

%% 2. match between lesion mask and segmented tissue
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

%% 3. Volumes
clear matlabbatch
fn_seg8 = spm_select('FPList',pth,'^.*_seg8\.mat$');
matlabbatch{1}.spm.util.tvol.matfiles = {fn_seg8};
matlabbatch{1}.spm.util.tvol.tmax = 4;
matlabbatch{1}.spm.util.tvol.mask = {fullfile(spm('dir'),'tpm','mask_ICV.nii,1')};
% matlabbatch{1}.spm.util.tvol.outf = 'test_volume';
matlabbatch{1}.spm.util.tvol.outf = fullfile(pth_out,['TV_',fn]);
spm_jobman('run', matlabbatch);
tmp = load(fn_seg8);
volumes = tmp.volumes;
tiv = sum(volumes.litres);
wmv = sum(volumes.litres([2 3])); % WM + lesion volume!
lesionV_per_tiv = volumes.litres(3)/tiv*100 ; % expresse in percentage
lesionV_per_wmv = volumes.litres(3)/wmv*100 ; % expresse in percentage
GM_per_tiv = volumes.litres(1)/tiv*100 ; %
WM_per_wmv = volumes.litres(2)/wmv*100 ; % WM *without* lesion !!!
res.lesionV = struct(...
    'volumes', volumes.litres, ...
    'tiv', tiv, ...
    'wmv', wmv, ...
    'GM_per_tiv', GM_per_tiv, ...
    'WM_per_wmv', WM_per_wmv, ...
    'lesionV_per_tiv', lesionV_per_tiv, ...
    'lesionV_per_wmv', lesionV_per_wmv);

%% 4. extraction of MPM values for the GM/WM/lesion
% and discarding voxels that were thresholded, just in case
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

vMPM = cell(3,nMPM); % #tissue classes x #MPM
for ii=1:nMPM
    v_mpm = spm_read_vols(Vmpm(ii));
    v_mpm = v_mpm(:);
    if nMsk
        v_msk = ~spm_read_vols(Vmsk(ii));
    end
    for jj=1:3 % tc
        tmp = vt_tc123(:,:,:,jj);
        if nMsk
            tmp = tmp & v_msk;
        end
        vMPM{jj,ii} = v_mpm(tmp(:));
    end
end
res.vMPM = vMPM; %#ok<*STRNU>

%% 5. some stats from the MPM values
% Find min-max, mean, median, std, skewness ,kurtosis for each MPM
mM = zeros(nMPM,2); mM(:,1) = Inf; mM(:,2) = -Inf;
meanVal = zeros(3,nMPM); medVal = zeros(3,nMPM); stdVal = zeros(3,nMPM);
skewVal = zeros(3,nMPM); kurtVal = zeros(3,nMPM);
for ii=1:3 % tc
    for jj=1:nMPM % mpm
        % min/max total
        tmp_m = min(res.vMPM{ii,jj});
        if tmp_m<mM(jj,1), mM(jj,1) = tmp_m; end
        tmp_M = max(res.vMPM{ii,jj});
        if tmp_M>mM(jj,2), mM(jj,2) = tmp_M; end
        % mean/meadian/std/skewness/kurtosis
        meanVal(ii,jj) = mean(res.vMPM{ii,jj});
        medVal(ii,jj) = median(res.vMPM{ii,jj});
        stdVal(ii,jj) = std(res.vMPM{ii,jj});
        skewVal(ii,jj) = skewness(res.vMPM{ii,jj});
        kurtVal(ii,jj) = kurtosis(res.vMPM{ii,jj})-3;
    end
end
res.vMPMstats = struct('mM',mM,'meanVal',meanVal,'medVal',medVal,...
    'stdVal',stdVal,'skewVal',skewVal,'kurtVal',kurtVal);

%% 6. save things and pass out fn_out
fn_out.fn_ExParam = fullfile(pth_out,['ExP_',fn]);
save(fn_out.fn_ExParam,'res');


end

%% SUBFUNCTION

function matlabbatch = create_mask(fn_in,fn_out,pth,thr,smK)
% Batch to create a mask from a few segmented tissue classes.
% It works in 3 steps:
% 1. sum up the input images
% 2. smooth the sum
% 3. threshold  the result

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



%% VARIOUS BITS OF CODE FOR DISPLAY
% figure, hist(v_tmsk(v_tmsk>0)),
% min(v_tmsk), max(v_tmsk)
% figure, hist(v_c3(v_c3>0)),
% min(v_c3), max(v_c3)

% figure,
% hist(vMPM{1})


