function fn_tmsk = crc_fix_LesMsk(fn_msk,opt)
% Fixing a lesion mask image by 
% 1/ removing too small blobs
% 2/ optionally removing blobs not meeting an intensity criteria in a
%    complementary image.
% 
% 
% INPUT:
% - fn_msk   : filename to lesion mask image
% - opt      : option structure
%   .minVol  : min volume (mm^3) of blobs to keep [de. 8]
%   .fn_oth  : filename to complementary image [def. empty]
% 
% OUTPUT:
% - fn_tmsk  : file name of 'thresholded' mask image
% 
% NOTES:
% 1/ The 1st fix relies simply on individual lesion volume, by removing 
%    bits that would be too small to really matter according to medical 
%    criteria. 
%    For multiple sclerosis (MS), the usual clinical criteria is
%       "Lesions will ordinarily be larger than 3 mm in cross section"
%    A cubic volume of 2x2x2 mm^3 has a diagonal of sqrt(12)~3.4mm and a 
%    volume of 8 mm^3, hence the default value 'minVol = 8'.
% 2/ The optional 2nd fix was devised for the case of lesion defined from
%    hyper-intenisties in a structural image, e.g. FLAIR image in the case 
%    of MS patient.
%    The aim is thus to improve the original mask by removing bits that are 
%    not really lesion hyper-intensities in the FLAIR image. The intuition 
%    is the following: 
%    a) most of the mask is correct and covers hyper-intense blobs but some 
%       don't and should be eliminated. 
%    b) The lesion blobs will be mostly made of high-intensity voxels but 
%       may also include "a few" low intensity ones.
%    c) Therefore the not-lesion blobs should have a mean/median intensity 
%       lower than the majority of the mask voxels...
%_______________________________________________________________________
% Copyright (C) 2018 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

%% Fixed parameters
prct_threshold = 95; % percentile threshold
use_mean_or_med = 1; % use mean (1) of all lesion voxel intensities or median (2)

%% 0. Check input & Load things
if nargin<2
    opt = struct;
end
opt_def = struct(...
    'minVol', crc_USwL_get_defaults('ImgFix.minVol'), ... % minimum volume of individual lesion (mm^3)
    'fn_oth', []);
opt = crc_check_flag(opt_def,opt);

Vmsk = spm_vol(fn_msk);
vMsk = spm_read_vols(Vmsk);

% Ensures values are [0 1], in case scaling was wrong, e.g. [0 255], or
% there are some tiny negative values, e.g. if mask was resampled
if max(vMsk(:))>1 || min(vMsk(:))<0
    fprintf('WARNING: some bad values in the lesion mask!\n')
    fprintf('\tValue range : [%1.4f %1.4f] -> Setting it to [0 1]\n', ...
        min(vMsk(:)), max(vMsk(:)))
    vMsk(vMsk>=.5) = 1; % values >=.5, set to 1
    vMsk(vMsk<.5) = 0; % anything <.5 set to zero.
end
vMsk = round(vMsk);
% VX coordinates
% XYZvx = round(Vmsk.mat\[XYZmm ; ones(1,size(XYZmm,2))]);
% XYZvx(4,:) = [];

% minVol mm3 -> vx
vx_vol = abs(det(Vmsk.mat(1:3,1:3)));
minNrVx = floor(opt.minVol/vx_vol); 
% Using 'floor' to round off minimum number of voxels

%% 1. removing too small cluster
% Find clusters
[L,num] = spm_bwlabel(vMsk,18);

% Get number of voxels per cluster
N = zeros(1,num);
for ii=1:num
    N(ii) = sum(L(:)==ii);
end

% clusters to remove
l_clust2rem = find(N<minNrVx);
% Create fixed volume w/o clust2rem
tvMsk = zeros(size(vMsk));
for ii=1:num
    if ~any(ii==l_clust2rem)
        tvMsk(L(:)==ii) = 1; % keep cluster
    else
        L(L(:)==ii) = 0; % remove from clusters
    end
end

%% 2. removing clusters based on intensities in adjoining image
if ~isempty(opt.fn_oth)
    % deal with other image
    Voth = spm_vol(opt.fn_oth);
    % check is has same size as lesion mask!
    if all(Vmsk.dim == Voth.dim)
        % proceed with remaining clusters
        clust_list = 1:num; clust_list(l_clust2rem) = [];
        vOth = spm_read_vols(Voth);
        switch use_mean_or_med
            case 1 % use lesion(s) global mean 
                v_threshold = mean(vOth(~~tvMsk(:)));
            case 2 % use lesion(s) global median
                v_threshold = median(vOth(~~tvMsk(:)));
        end
        % load up values per cluster
        n_clust_list = numel(clust_list);
        v_les_fla_clust = cell(1,n_clust_list);
        prct_les_fla_clust = zeros(1,n_clust_list);
        for ii=1:n_clust_list
            v_les_fla_clust{ii} = vOth(L(:)==clust_list(ii));
            prct_les_fla_clust(ii) = prctile(v_les_fla_clust{ii},prct_threshold);
        end
        % actually it's not necessary to save cluster intensities but could
        % be for other checks... so keep it just in case.
        l_clust2rem = clust_list(prct_les_fla_clust<v_threshold);
        % remove these ones
        for ii=1:num
            if ~any(ii==l_clust2rem)
                tvMsk(L(:)==ii) = 1; % keep cluster
            else
                L(L(:)==ii) = 0; % remove from clusters
            end
        end
    else
        beep;
        fprintf(['\nOther image size [%d %d %d] does not match that of '...
            'the mask [%d %d %d].\nNot proceeding any further!\n'], ...
            [Voth.dim Vmsk.dim])
    end
end

%% Save image & return filename

Vtmsk = Vmsk;
Vtmsk.fname = spm_file(Vmsk.fname,'prefix','t');
Vtmsk = spm_create_vol(Vtmsk);
Vtmsk = spm_write_vol(Vtmsk,tvMsk); %#ok<*NASGU>
% output filename
fn_tmsk = Vtmsk.fname;

end
