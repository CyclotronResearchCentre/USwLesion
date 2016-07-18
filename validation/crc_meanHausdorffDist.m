function [mD,D12,D21] = crc_meanHausdorffDist(xyz1,xyz2)
% Calculate the mean Hausdorff distance between 2 surfaces, based on the
% coordinates of the voxels of these surfaces.
% This implementation only works with 3D coorinates but is much faster than
% standard implementations (with for-loops).
% Note that with the crc_borderVx.m function handles 2D images and return
% the coordinates of border pixels as in 3D.
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
%
% NOTE:
% If one of the surfaces is empty, i.e. there was no blob in that image,
% then distances are set to NaN for the empty image and INF for all the
% border voxels of the other.
%__________________________________________________________________________
% Copyright (C) 2016 Cyclotron Research Centre

% Written by Christophe Phillips
% University of Liege, Belgium

% Switch to true to get some stats and plots generated.
fl_disp = false;

% Make sure that the coordinates are in N x 3 format
[xyz1,SZ1,ok1] = check_array(xyz1);
[xyz2,SZ2,ok2] = check_array(xyz2);

if ok1 && ok2 % Coordinates are both ok, proceed
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
else % some border array must be empty
    if ok1 && ~ok2
        D12 = zeros(SZ1(1),1)+Inf;
        D21 = NaN;
    elseif ~ok1 && ok2
        D12 = NaN;
        D21 = zeros(SZ2(1),1)+Inf;
    elseif ~ok1 && ~ok2
        D12 = NaN;
        D21 = NaN;
    end
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
    if ok1 && ok2
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
end

%% SUBFUNCTIONS

% Make sure that the coordinates are in N x 3 format
% If the array is empty then ok=false.
function [xyz,SZ,ok] = check_array(xyz)
SZ = size(xyz);
ok = false;
if any(SZ==3)
    if SZ(1)<SZ(2)
        xyz = xyz' ;
        SZ = size(xyz);
    end
    ok = true;
elseif prod(SZ)
    error('Coordinates should be a 3xN or Nx3 array');
end

end