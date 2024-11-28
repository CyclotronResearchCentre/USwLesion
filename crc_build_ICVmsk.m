function fn_out = crc_build_ICVmsk(fn_in,opt)
% Function that builds an ICV mask from a set of segmented tissue classes,
% these should be GM, WM & CSF at least.
% This works by 
% - adding the tissue class images together -> all voxels get a p-ICV
% - smoothing (4 mm FHWM) -> cleans up and enlarges things a bit
% - thresholding at .3 to keep those really in the ICV
% - fixing the obtained mask in different ways (see crc_fix_ICV.m) if
%   either a "mask_sz_thr" or "fn_iwarp" is provided.
% 
% Other operations include creating a (smoothed) normalized ICV-mask,
% depending on the options selected.
% 
% INPUT
% - fn_in   : filename of tissue classes to add together (char array)
% - opt     : option structure
%       fn_ref      : reference filename -> ['icv_',fn_ref]
%       smoK        : kernel size (mm) for clean up smoothing [4, by def.]
%       fn_warp     : forward warping -> create warped ICV
%       msk_sz_thr  : max size (mm^3) of small blobs to be filled up [1000, def.]
%       fn_iwarp    : inverse warping -> fix ICV mask with inv-warped SPM-ICV
%       wsmoK       : kernel size (mm) of warped ICV smoothing [4, by def.]
% 
% OUTPUT
% - fn_out
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

%% Deal with input
if nargin<2, opt = struct; end
opt_o = struct( ...
    'fn_ref', [], ...
    'smoK', 4, ...
    'fn_warp', [], ...
    'fn_iwarp', [], ...
    'msk_sz_thr', 1000, ...
    'wsmoK', 4);
opt = crc_check_flag(opt_o,opt);
% Ensures smoothing kernels are 1x3 vectors
if numel(opt.smoK)~=3, opt.smoK = ones(1,3)*opt.smoK(1); end
if numel(opt.wsmoK)~=3, opt.wsmoK = ones(1,3)*opt.wsmoK(1); end

if ~isempty(opt.fn_ref)
    fn_ref = opt.fn_ref;
else
    fn_ref = fn_in(1,:);
end

%% Create ICV mask: sum, smooth, threshold
% 1/ sum
V_in = spm_vol(fn_in);
V_o1 = V_in(1);
fn_tmp = spm_file(V_o1.fname,'filename','tmp.nii');
V_o1.fname = fn_tmp;
V_o1.dt(1) = 2;
imc_fl = struct(...
    'dmtx', 1, ...
    'interp', 1, ...
    'dtype', 2, ...
    'descrip','Sum of TCs');
V_o1 = spm_imcalc(V_in, V_o1, 'sum(X)' ,imc_fl);
% 2/ smooth
fn_stmp = spm_file(fn_tmp,'prefix','s');
spm_smooth(V_o1.fname,fn_stmp,opt.smoK);
% 3/ threshold
fn_ICV = spm_file(fn_ref,'prefix','icv_','number','');
V_o2 = V_o1;
V_o2.fname = fn_ICV;
imc_fl.descrip = 'ICV mask';
imc_fl.dmtx = 0;
V_o2 = spm_imcalc(fn_stmp, V_o2, 'i1>.3' ,imc_fl); %#ok<*NASGU>

% Delete tmp files
delete(fn_tmp), delete(fn_stmp)

% Collect output
fn_out = fn_ICV;

%% Fixing ICV
if (~isempty(opt.msk_sz_thr) && opt.msk_sz_thr>0) || ~isempty(opt.fn_iwarp)
    % proceed with fixing ICV or not, if msk_sz_thr exists and is >0, or if
    % inverse warp is provided
    if ~isempty(opt.msk_sz_thr)
        opt_fx_mask.sz_thr = opt.msk_sz_thr;
    else
        opt_fx_mask.sz_thr = 0;
    end
    if ~isempty(opt.fn_iwarp)
        opt_fx_mask.fn_iwarp = opt.fn_iwarp;
    end
    % do it
    crc_fix_ICV(fn_ICV,opt_fx_mask); % overwriting the file fn_ICV
end

%% Warping and smoothing, if requested
if ~isempty(opt.fn_warp)
    % Apply warping
    matlabbatch = crt_batch_normalize_write(fn_ICV,spm_file(opt.fn_warp,'number',''));
    spm_jobman('run', matlabbatch);
    fn_wICV = spm_file(fn_ICV,'prefix','w');
    fn_out = char(fn_out,fn_wICV);
    if ~isempty(opt.wsmoK) && all(opt.wsmoK>0)
        % Smooth a bit
        fn_swICV = spm_file(fn_ICV,'prefix','sw');
        spm_smooth(fn_wICV,fn_swICV,opt.wsmoK,2);
        fn_out = char(fn_out,fn_swICV);
    end
end

%% Output filename(s) collected on the way -> fn_out
fprintf('\nCreated the following ICV files:\n')
for ii=1:size(fn_out,1)
    fprintf('\t%s\n',fn_out(ii,:));
end

end

% =======================================================================
%% SUBFUNCTIONS
% =======================================================================

function matlabbatch = crt_batch_normalize_write(fn_ICV,fn_warp)
% Building matlabbatch for the normalize-write operation, muche easier than
% building the deformation and applying it manually.
matlabbatch{1}.spm.spatial.normalise.write.subj.def(1) = {fn_warp};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {fn_ICV};
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72 ; 90 90 108];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1.5 1.5 1.5];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;
end


