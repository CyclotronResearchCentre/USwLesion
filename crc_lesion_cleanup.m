function [fn_out,nDropped] = crc_lesion_cleanup(fn_in,k,prefix)

% This is a silly routine to cleanup anu binary image based on spatial
% extent of clsuters. During segmentation some spurious voxels may end up 
% misclassified. So if we know a clsuter (e.g. lesions) must be at least of
% a given size, we can use spatial extend to clean up the small clusters.
%
% FORMAT 
% fn_out = lesion_crc_cleanup(fn_in,k,prefix)
%
% INPUT 
% fn_in     : the volume name (typically c3 image, the lesion map)
% k         : the spatial extent threshold [Inf, def., i.e largest cluster]
% prefix    : prefix prepended to image name
%
% OUTPUT 
% fn_out    : file name of cleaned image, writen on the disk
% nDropped  : number of clusters that were dropped out and remaining number
%             of clusters.
% 
%
% see also spm_bwlabel
% Cyril Pernet 30 Oct 2015
%_______________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

if nargin<2
    prefix = 'cleaned_';
end
if nargin<1
    k = Inf; % if not specified we keep only the biggest cluster
end
if nargin == 0
    fn_in = spm_select(1,'image','select image to threshold');
end

% get image in
if iscell(fn_in); fn_in = cell2mat(fn_in); end
if ischar(fn_in)
    V_in = spm_vol(fn_in);
    data = spm_read_vols(V_in);
    save_img = true;
else
    data = fn_in;
    save_img = false;
end

[clustered_map,num] = spm_bwlabel(double(data>0),6);
% note we use surface connection only

% check number of elements per clusters
extent_map = zeros(size(data));
clustered_map = clustered_map(:);
nv = histc(clustered_map,0:num);
[~,idxall]=sort(clustered_map,'ascend');
idxall(1:nv(1)) = []; nv(1)=[]; % remove 1st bin ie 0s
ends = cumsum(nv);
inis = ends-nv+1;

% make the new map - if cluster >=k then = 1
if isinf(k) || ischar(k);
    [~,pos] = max(nv);
    idx=idxall(inis(pos):ends(pos));
    extent_map(idx)=1;
    nDropped = [numel(nv)-1 1];
else
    nDropped = 0;
    for i=1:num
        idx=idxall(inis(i):ends(i));
        if nv(i) >= k
            extent_map(idx)=1;
        else
            extent_map(idx)=0;
            nDropped = nDropped + 1;
        end
    end
    nDropped(2) = num-nDropped(1);
end

extent_map = extent_map .* data; % reasign the initial values
if save_img
    fn_out = spm_file(V_in.fname,'prefix',prefix);
    V_out.fname = fn_out;
    V_out.descrip = [V_out.descrip ' thresholded @ k=' num2str(k)];
    V_out.private.descrip = V_out.descrip;
    spm_write_vol(V_out,extent_map);
    fprintf('%s created \n',V_out.fname)
else
    fn_out = extent_map;
end

end
