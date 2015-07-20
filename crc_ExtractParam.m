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
fn_Msk  = char(job.imgMsk); 
fn_MPM  = char(job.imgMPM); 
[pth,fn] = spm_fileparts(fn_MPM(1,:));
if isempty(job.outdir)
    pth_out = pth;
else
    pth_out = job.outdir{1};
end
opt = job.opt;

%% 1. build ICV from the sum of GM, WM, CSF and lesion
Vtc = spm_vol(fn_cImg);
v_tc = spm_read_vols(Vtc(1:4));
v_icv = sum(v_tc,4);
% size(v_tc), size(v_icv)
v_icv = v_icv > opt.thrICV;

Vicv = Vtc(1);
Vicv.fname = spm_file(Vicv.fname,'prefix','icv_');
Vicv = spm_create_vol(Vicv);
Vicv = spm_write_vol(Vicv,v_icv);
res.Picv = Vicv.fname;

%% 2. match between mask and segmented tissue
% In native space
Vtmsk = spm_vol(fn_Msk);
Vc3 = spm_vol(fn_cImg(3,:));
v_tmsk = spm_read_vols(Vtmsk);
v_tmsk = v_tmsk(:) >.5; % thresholding at .5, as it shoul dbe a binary img
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
matlabbatch{1}.spm.util.tvol.outf = 'test_volume';
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
nMPM = size( fn_MPM,1);
Vmpm = spm_vol(fn_MPM);
Vtc = spm_vol(fn_cImg);
v_tc123 = spm_read_vols(Vtc(1:3));
vt_tc123 = v_tc123 > opt.thrTC ;

vMPM = cell(3,1);
for ii=1:nMPM
    v_mpm = spm_read_vols(Vmpm(ii));
    v_mpm = v_mpm(:);
    for jj=1:3
        tmp = vt_tc123(:,:,:,jj);
        vMPM{jj}(:,ii) = v_mpm(tmp(:));
    end
end
res.vMPM = vMPM; %#ok<*STRNU>

fn_out.fn_ExParam = fullfile(pth_out,['ExP_',fn]);
save(fn_out.fn_ExParam,'res');


end

%% VARIOUS BITS OF CODE FOR DISPLAY 
% figure, hist(v_tmsk(v_tmsk>0)),
% min(v_tmsk), max(v_tmsk)
% figure, hist(v_c3(v_c3>0)),
% min(v_c3), max(v_c3)

% figure,
% hist(vMPM{1})


