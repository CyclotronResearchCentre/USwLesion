function [mD,D12,D21] = crc_meanHausdorffDist(xyz1,xyz2)
% Calculate the mean Hausdorff distance between 2 surfaces, based on the
% coordinates of the voxels of these surfaces.
% This implementation only works with 3D volumes but is much faster than
% standard implementations (with for-loops).
%
% FORMAT
%   [mD,D12,D21] = crc_meanHausdorffDist(xyz1,xyz2)
%
% INPUT
%   xyz1/2  coordinates of voxels on borders in both images
% 
% OUTPUT
%   mD      mean distances between surface in images 1-2 and 2-1
%   D12/D21 distances between surface in images 1-2 and 2-1
%__________________________________________________________________________
% Copyright (C) 2016 Cyclotron Research Centre

% Written by Christophe Phillips
% University of Liege, Belgium

% Switch to true to get some stats and plots generated.
fl_disp = false;

% Make sure that the coordinates are in N x 3 format
SZ1 = size(xyz1);
if SZ1(1)<SZ1(2)
    xyz1 = xyz1' ; xyz2 = xyz2';
    SZ1 = size(xyz1); SZ2 = size(xyz2);
else
    SZ2 = size(xyz2); % ??
end

% Initialize the 2 distance vectors
D12 = zeros(SZ1(1),1);
D21 = zeros(SZ2(1),1)+Inf;

% Loop only once for the voxels of the 1st image
for ii=1:SZ1(1)
    Dtmp = sqrt( (xyz2(:,1)-xyz1(ii,1)).^2 + ... % distance from i_th voxel
        (xyz2(:,2)-xyz1(ii,2)).^2 + ...          % in 1st image to all 
        (xyz2(:,3)-xyz1(ii,3)).^2 );             % voxels in 2nd image
    D12(ii) = min(Dtmp); % smallest distance from i_th voxel in 1st image
    D21 = min(D21,Dtmp); % smallest distance for all voxels in 2nd image
end

mD = [mean(D12) mean(D21)];

% provide some stats and evaluation
if fl_disp
    mM = [ min(D12) max(D12) median(D12) ; min(D21) max(D21) median(D12) ]; %#ok<*UNRCH>
    fprintf('\n')
    fprintf('D12 mean/min/max/median: %1.2f / %1.2f / %1.2f / %1.2f',mD(1),mM(1,:))
    fprintf('\n')
    fprintf('D21 mean/min/max/median: %1.2f / %1.2f / %1.2f / %1.2f',mD(2),mM(2,:))
    fprintf('\n')
    
    % Prepare histogram of distances for both images
    Nbins = 10^floor(log10(SZ1(1)/100))/2;
    if mM(1,2)>mM(2,2)
        [H12,bi] = hist(D12,Nbins);
        H21 = hist(D21,bi);
    else
        [H21,bi] = hist(D21,Nbins);
        H12 = hist(D12,bi);
    end
        figure, bar(bi,[H12' H21'])
    legend('D12','D21')
end
end