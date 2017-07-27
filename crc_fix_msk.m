function fn_msk = crc_fix_msk(fn_msk,opt)
% Fixing a mask image, e.g. an intracranial volume (ICV) mask, by 
% 1/ filling small holes, based on #voxel threshold (1000mm^3 by def.)
% 2/ removing blobs outside brain, i.e. the biggest blob
%
% Note that if the mask is fixed it is *overwritten* on disk!
% 
% INPUT:
% - fn_msk  : filename to (ICV) mask volume on disk
% - opt     : option structure
%   .sz_thr : maximum size (mm^3) of small blobs that can be filled up
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
    'sz_thr', 1000); % maximum size (mm^3) of small blobs to be filled up
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

%% Write out if somethign was fixed
if any_fix % Need to save something
    %     V_msk.private.dat = v_msk;
    spm_write_vol(V_msk,v_msk);
end

end