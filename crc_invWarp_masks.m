function fn_out = crc_invWarp_masks(fn_msk, fn_iwarp, fn_ref)
% Function to inverse-warp, i.e. bring from tempalte (MNI) space into
% subject space, some mask images (or anyother image infact).
%
% INPUT
% - fn_msk   : filename(s) of (mask) image(s) to inv-warp (char array)
% - fn_iwarp : filename to inverse-warr transformation (iy_*.nii image)
% - fn_ref   : filename of reference image in subject space
%
% OUTPUT
% - fn_out   : filename of inv-wapred images
% 
% NOTE
% The created images are placed in the same folder at the reference image
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium


%% Create the matlabbatch job
prec_round = 1e6;
V_ref = spm_vol(fn_ref);
% voxel size rounded to some precision
vx_sz = round(sqrt(sum(V_ref.mat(1:3,1:3).^2))*prec_round)/prec_round;
% defining BB
p_min = -V_ref.mat\[0 0 0 1]' ;
p_max = (V_ref.dim.*vx_sz)' + p_min(1:3) ;
img_bb = [-abs(p_min(1:3)') ; abs(p_max(1:3)')];

matlabbatch{1}.spm.spatial.normalise.write.subj.def(1) = ...
    {spm_file(fn_iwarp,'number','')};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = ...
    cellstr(fn_msk);
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = img_bb;
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = vx_sz;
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'iw';

%% Run matlabbatch job
spm_jobman('run', matlabbatch);
fn_out = spm_file(fn_msk,'prefix','iw');

%% Move files in destination folder & update filenames
pth_dest = spm_file(fn_ref,'path');
n_fn = size(fn_msk,1);
% Difficulty figure if they're from a 4D volume or N 3D volumes
fn_out_noNb = spm_file(fn_out,'number','');
l_match = strcmp(fn_out_noNb(1,:),cellstr(fn_out_noNb));
if all(l_match)
    % Single 4D file case, move only one
    movefile(fn_out_noNb(1,:),pth_dest)
    fn_out = spm_file(fn_out,'path',pth_dest);
else
    % Multiple 3D files -> move 1 by 1
    for ii=1:n_fn
        movefile(fn_out(ii,:),pth_dest)
    end
    fn_out = spm_file(fn_out,'path',pth_dest);
end


end

% % test
% fn_msk = spm_select(Inf,'image');
% fn_iwarp = spm_select(1,'nii','Inverse warp',{},pwd,'^iy_.*\.nii');
% fn_ref = spm_select(1,'nii','Reference image');
% fn_out = crc_invWarp_masks(fn_msk, fn_iwarp, fn_ref)


