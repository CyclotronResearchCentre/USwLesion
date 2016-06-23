function Pout = crc_lesion_cleanup(Pin,k,prefix)

% this is a silly routine to cleanup the lesion mask (c3) based on spatial
% extend, during segmentation some voxels may endup misclassified and if
% we know we have lesions of a given size, we can use spatial extend to
% cleanup
%
% FORMAT 
% Pout = lesion_crc_cleanup(Pin,k)
%
% INPUT 
% Pin    : the volume name (typically c3 image, the lesion map)
% k      : the spatial extent threshold [Inf, def., i.e largest cluster]
% prefix : prefix prepended to image name
%
% OUTPUT 
% Pout   : output file name of cleaned image
% clean_P is writen on the disk
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
    Pin = spm_select(1,'image','select image to threshold');
end

% get clusters
if iscell(Pin); Pin = cell2mat(Pin); end    
V = spm_vol(Pin);
data = spm_read_vols(V);
[clustered_map,num] = spm_bwlabel(double(data>0),6); % note we use surface connection only

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
    
else
    for i=1:num
        idx=idxall(inis(i):ends(i));
        if nv(i) >= k
            extent_map(idx)=1;
        else
            extent_map(idx)=0;
        end
    end
end

extent_map = extent_map .* data; % reasign the initial values
Pout = spm_file(V.fname,'prefix',prefix);
V.fname = Pout;
V.descrip = [V.descrip ' thresholded @ k=' num2str(k)];
V.private.descrip = V.descrip;
spm_write_vol(V,extent_map);
fprintf('%s created \n',V.fname)

end
% to make an extent map do
% for i=1:num
%     idx=idxall(inis(i):ends(i));
%     extent_map(idx)=nv(i);
% end
% V.fname = [root filesep 'extent_map_' name ext];
% V.descrip = [V.descrip ' spatial extent map'];
% V.private.descrip = V.descrip;
% spm_write_vol(V,extent_map);
% fprintf('%s created \n',V.fname)
