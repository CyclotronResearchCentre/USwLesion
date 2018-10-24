function varargout = crc_borderVx(vol,verbose)
%
% Function that returns the list of voxels at the 'border' of clusters in a
% 3D binary image. A border voxel is any voxel whose 18 direct neighbours
% are not all surrounded by voxels in that cluster.
%
% FORMAT
%   lBorder = crc_borderVx(vol,verbose)
% or
%   [lBx,lBy,lBz] = crc_borderVx(vol,verbose)
%
% INPUT
%   vol      binary 3D array, voxels with 1's are the clusters
%   verbose  display some numbers about volume and surface [Def., true]
%
% OUTPUT
%   lBorder  list of border voxels
% or
%   [lBx,lBy,lBz] indexes of the border voxels
%
%--------------------------------------------------------------------------
% NOTE:
% - for the voxels in a cluster on the edge of the volume, what's outside 
%   the volume is considered as part of the cluster. In other words, voxels
%   on the surface of the volume are NOT considered at the surface of the
%   cluster.
% - if the image is empty, i.e. no blobs, the output argument(s) will be
%   empty.
%__________________________________________________________________________
% Copyright (C) 2016 Cyclotron Research Centre

% Written by Christophe Phillips
% University of Liège, Belgium

if nargin<2, verbose = true; end

%% Get image sizes, define clique, etc.
vol = vol>0;
DIM = size(vol);
if numel(DIM)==2 % 2D image case
    DIM(3) = 1;
end
lMsk = find(vol(:));
if numel(lMsk)
    [X,Y,Z] = ndgrid(-1:1,-1:1,-1:1);
    clique = [X(:) Y(:) Z(:)];
    sXYZ = sum(abs(clique),2);
    clique((sXYZ==3) | (sXYZ==0),:) = []; % remove centre and 8 outer corners
    nCl = size(clique,1);
    nCl_zero = zeros(nCl,1);
    
    %% Find the border voxles
    lB_tmp = logical(size(lMsk));
    for ii=1:numel(lMsk)
        ii_vx = lMsk(ii);
        [xii,yii,zii] = ind2sub(DIM,ii_vx);
        xyz_neighb = clique + [nCl_zero+xii nCl_zero+yii nCl_zero+zii];
        % check they're inside the volume
        %     l_out = find( ...
        %                     (xyz_neighb(:,1)<1) | (xyz_neighb(:,1)>DIM(1)) | ...)
        %                     (xyz_neighb(:,2)<1) | (xyz_neighb(:,2)>DIM(2)) |...)
        %                     (xyz_neighb(:,3)<1) | (xyz_neighb(:,3)>DIM(3)) );
        % 	xyz_neighb(l_out,:) = [];
        xyz_neighb(...
            (xyz_neighb(:,1)<1) | (xyz_neighb(:,1)>DIM(1)) | ...)
            (xyz_neighb(:,2)<1) | (xyz_neighb(:,2)>DIM(2)) |...)
            (xyz_neighb(:,3)<1) | (xyz_neighb(:,3)>DIM(3)) , :) = [];
        %     l_neighb = sub2ind(DIM,xyz_neighb(:,1),xyz_neighb(:,2),xyz_neighb(:,3));
        %     v_neighb = vol(l_neighb);
        v_neighb = vol(sub2ind(DIM, ...
                        xyz_neighb(:,1),xyz_neighb(:,2),xyz_neighb(:,3)));
        lB_tmp(ii) = ~all(v_neighb);
    end
    lBorder = lMsk(lB_tmp);
else
    lBorder = [];
end
% Number of border voxels:
if verbose
    fprintf('\nThere are %d surface voxels out of %d/%d voxels in clusters/total.\n', ...
        length(lBorder), sum(vol(:)), prod(DIM))
end

%% Deal with output
if nargout==1
    varargout{1} = lBorder;
elseif nargout==3
    [lBx,lBy,lBz] = ind2sub(DIM,lBorder);
    varargout{1} = lBx; varargout{2} = lBy; varargout{3} = lBz;
else
    error('Wrong number of output!')
end

end