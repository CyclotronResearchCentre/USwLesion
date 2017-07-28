function fn_msk = crc_fix_msk(fn_msk,opt)
% Fixing a mask image, e.g. an intracranial volume (ICV) mask, by 
% 1/ filling small holes, based on #voxel threshold (1000mm^3 by def.)
% 2/ removing blobs outside brain, i.e. the biggest blob
% 3/ if some inverse warping (MNI->subj) is provided, adjust ICV mask with
%    the SPM one. Simply adding both and checking overlap.
%
% Note that if the mask is fixed it is *overwritten* on disk!
% 
% INPUT:
% - fn_msk  : filename to (ICV) mask volume on disk
% - opt     : option structure
%   .sz_thr  : max size (mm^3) of small blobs to be filled up (1000, def.)
%   .fn_iwarp: filename to inverse warping to use SPM-ICV mask
% 
% OUTPUT:
% - fn_msk  : same filename since the mask is overwritten.
% 
% TO DO:
% more (advanced) procedures could be added, along with some options
% (parameters or switches).
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

%% Deal with input
any_fix = false; % Flag to remember if something was fixed, start as 'no'

% Check options
opt_o = struct(...
    'sz_thr', 1000, ... % maximum size (mm^3) of small blobs to be filled up
    'fn_iwarp', []);
opt = crc_check_flag(opt_o,opt);

% Load in mask
V_msk = spm_vol(fn_msk);
v_msk = spm_read_vols(V_msk);

% Express sz_thr in voxels -> sz_thr_vx
vx_vol = abs(det(V_msk.mat(1:3,1:3))); % in mm^3
sz_thr_vx = ceil(opt.sz_thr/vx_vol);

%% Deal with small holes
[L,num] = spm_bwlabel(double(~v_msk),18);
for ii=1:num
    n_vx = sum(L(:)==ii);
    if n_vx<=sz_thr_vx % if not too big, fix it
        v_msk(L(:)==ii) = 1;
        any_fix = true;
    end
end

%% Remove big blobs outside the biggest one
[L,num] = spm_bwlabel(v_msk);
for ii=1:num
    n_vx = sum(L(:)==ii);
end

if numel(n_vx)>1
    [~,s_ind] = sort(n_vx);
    for ii = s_ind(1:end-1) % clear all but the biggest
        v_msk(L==ii) = 0;
    end
    any_fix = true;
end

%% Write out if something was fixed
if any_fix % Need to save something?
    spm_write_vol(V_msk,v_msk);
end

%% Use inverse warped SPM-ICV to fix ICV mask, if requested
if ~isempty(opt.fn_iwarp)
    % Inverse warp SPM's ICV mask and merge this with the incoming ICV mask
    % + return some matching value (Jaccard index), for informative check.
    % https://en.wikipedia.org/wiki/Jaccard_index
    [matlabbatch,fn_wicvSPM] = create_MB(opt.fn_iwarp,fn_msk);
    spm_jobman('run', matlabbatch);
    % Smoothing a bit the wSPM-ICV to enlarge volume -> play safe!
    V_wicvSPM = spm_vol(fn_wicvSPM);
    Vo = V_wicvSPM; Vo.pinfo(1) = 1; spm_imcalc(V_wicvSPM, Vo, 'i1*255');
    fn_swicvSPM = spm_file(fn_wicvSPM,'prefix','s');
    spm_smooth(fn_wicvSPM,fn_swicvSPM,ones(1,3)*4)
    Vo = V_msk; Vo.pinfo(1) = 1;
    fl_imcalc.interp = 1;
    % Summing both, subj & SPM, ICV such that 1-> subj, 2->SPM, 3->subj&SPM
    Vo = spm_imcalc(char(fn_msk,fn_swicvSPM), Vo, '(i1>0) + 2 * (i2>0)',fl_imcalc);
    v_icv = spm_read_vols(Vo);
    n_123 = [sum(v_icv(:)==1) sum(v_icv(:)==2) sum(v_icv(:)==3)];
    Jind = n_123(3)/sum(n_123);
    fprintf('\nJaccard index for subj-ICV & SPM-ICV : %1.4f\n',Jind)
    % Bring back to binary
    spm_imcalc(Vo, Vo, 'i1>0',fl_imcalc);
    if n_123(2)>0
        % Some voxels added from SPM-ICV
        any_fix = true;
    end
end

%% Write out if something was fixed
if any_fix % Need to save something?
    fprintf('\nUpdate file : %s\n', fn_msk)
end

end

% =======================================================================
%% SUBFUNCTIONS
% =======================================================================

function [matlabbatch,fn_wicvSPM] = create_MB(fn_iwarp,fn_ICV)
% Building matlabbatch for the normalize-write operation, muche easier than
% building the deformation and applying it manually.
% Bring SPM-ICV into subject space -> need to properly define the latter!

prec_round = 1e6;

V_icv = spm_vol(fn_ICV);
vx_sz = round(sqrt(sum(V_icv.mat(1:3,1:3).^2))*prec_round)/prec_round;
p_min = V_icv.mat*ones(4,1);
p_max = V_icv.mat*[V_icv.dim 1]';
img_bb = [-abs(p_min(1:3)') ; abs(p_max(1:3)')];

% Bring in SPM-ICV into subject space
fn_icvSPM = fullfile(spm('dir'),'tpm','mask_ICV.nii');
fn_icvSPM_loc = fullfile(spm_file(fn_iwarp,'path'),'icv_SPM.nii');
copyfile(fn_icvSPM,fn_icvSPM_loc)

matlabbatch{1}.spm.spatial.normalise.write.subj.def(1) = {spm_file(fn_iwarp,'number','')};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {fn_icvSPM_loc};
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = img_bb;
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = vx_sz;
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;

fn_wicvSPM = spm_file(fn_icvSPM_loc,'prefix','w');

end

