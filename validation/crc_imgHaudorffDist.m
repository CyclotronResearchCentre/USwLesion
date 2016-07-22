function [mD,D12,D21,nDropped] = crc_imgHaudorffDist(img1,img2,opt)
% Calculate the mean Hausdorff distance between the blobs of 2 binary
% images (img1/2). It works in 2 steps:
% 1. the coordinates of border voxels of each blobs and image are fist 
%   extracted, 
% 2. then the closest distance between pairs of blobs from both images are 
%   calculated (D12/D21) and averaged for each image individually (mD).
% 
% A preprocessing step can also be applied: only blobs that are at least
% partially overlapping across both images are kept for the H-distance
% estimate. In effect non-matched blobs across images tend to lead to 
% spuriously large distances...
%
%
% FORMAT
%   [mD,D12,D21,nDropped] = crc_imgHaudorffDist(img1,img2,opt)
%
% INPUT
%   img1/2  binary 3D array, voxels with 1's are the clusters
%   opt         some options
%       .v2r    voxel-to-realworld transformation
%       .BlobMatchOnly keep only matching blobs or not [false, default]
%
% OUTPUT
%   mD          mean distances between surface in images 1-2 and 2-1
%   D12/D21     distances between surface in images 1-2 and 2-1
%   nDropped    number of non-matching blobs dropped in img1/2 (only if
%               blob matching is requested)
%
% NOTE:
% Assume img1/2 are already binarized and properly loaded. The only extra
% bit needed is the voxel-to-realworld transformation, to account for
% (anisotropic) voxel size.
%__________________________________________________________________________
% Copyright (C) 2016 Cyclotron Research Centre

% Written by Christophe Phillips
% University of Liege, Belgium

%%
% *Check input data*
opt_def = struct('BlobMatchOnly',false,'v2r',eye(4)); % default options

if nargin == 0
    help crc_imgHaudorffDist
    return
    
elseif nargin == 2
    opt = []; % Go with the defaults
    
elseif nargin < 2 || nargin > 3
    error('Two or three inputs are expected - FORMAT: [mD,D12,D21] = crc_imgHaudorffDist(img1,img2,opt)')
end
% Fill the opt structure with defaults
opt = crc_check_flag(opt_def,opt);

%%
% Deal with blob overlap if requested
nDropped = [0 0];
if opt.BlobMatchOnly
    [L2,num2] = spm_bwlabel(double(img2),26);
    if num2>1
        [L1,num1] = spm_bwlabel(double(img1),26);
        if num1>1 % Look for blobs not matching and remove them
            NoMatch1 = false(1,num1); % \_ No match for blobs in img1/2
            NoMatch2 = false(1,num2); % /
            % Find the not-matching blobs
            for ii=1:num1
                NoMatch1(ii) = ~length(intersect(find(img2),find(L1==ii)));
            end
            for ii=1:num2
                NoMatch2(ii) = ~length(intersect(find(img1),find(L2==ii)));
            end
            % Remove not-matching blobs
            for ii=find(NoMatch1)
                img1(L1==ii) = false;
            end
            for ii=find(NoMatch2)
                img2(L2==ii) = false;
            end
            nDropped = [sum(NoMatch1) sum(NoMatch2)];
        end
    end
end

%%
% *get border voxels and calculate H-distance*
[iBx1,iBy1,iBz1] = crc_borderVx(img1);
[iBx2,iBy2,iBz2] = crc_borderVx(img2);

% Get coordinates in mm
Bxyz1_mm = opt.v2r(1:3,1:3)* [iBx1' ; iBy1' ; iBz1'];
Bxyz2_mm = opt.v2r(1:3,1:3)* [iBx2' ; iBy2' ; iBz2'];

% Calculate Hausdorff distance
[mD,D12,D21] = crc_meanHausdorffDist(Bxyz1_mm,Bxyz2_mm); %#ok<*ASGLU>

end

