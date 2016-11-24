function [fn_out,nDropped] = crc_lesion_cleanup(fn_in,opt)

% This is a silly routine to cleanup any binary image based 
% - on eroding/growing operations, and/or
% - on spatial extent of clusters
% 
% During segmentation some spurious voxels may end up misclassified. So if 
% we know a cluster (e.g. lesions) must be at least of a given size, we can
% use spatial extend to clean up the small clusters.
% Similarly it can be useful to erode volumes then grow them to tease out 
% compact-ish clusters as opposed to 'lines'. Indeed partial volume effect
% tends to produce very thin, border line/surface clusters.
%
% FORMAT 
% fn_out = lesion_crc_cleanup(fn_in,k,prefix)
%
% INPUT 
% fn_in     : the volume name (typically c3 image, the lesion map)
% opt       : structure with options
%   .k        the spatial extent threshold [Inf, def., i.e largest cluster]
%   .prefix   prefix prepended to image name [def, 'c']
%   .Neg      number of eroding then growing steps [def, [0 0]]
%   .Nge      number of growing then eroding steps [def, [0 0]]
%   .Nneighb  numbre of neighbours in 3D (6, 18 or 26) [def, 6]
%
% OUTPUT 
% fn_out    : file name of cleaned image, writen on the disk
% nDropped  : number of clusters that were dropped out and remaining number
%             of clusters, after erode/grow.
% 
% NOTE:
% - Eroding-then-growing (EG) or growing-then-eroding (GE)  are of course 
%   mutually exclusive! If both were set up simultaneously, then EG would
%   have priority.
% - the EG/GE operation, if requested, takes place before the cluster size
%   trimming
%
% see also spm_bwlabel
% Cyril Pernet 30 Oct 2015
%_______________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

%% Sort out input
if nargin<2, opt = []; end;
opt_def = struct('k', Inf, 'prefix', 'c', ...
                  'Neg', [], 'Nge', [], 'Nneighb',6);
opt = crc_check_flag(opt_def,opt); % Check and pad filter structure

k = opt.k;
prefix = opt.prefix;
Neg = opt.Neg;
Nge = opt.Nge;
Nneighb = opt.Nneighb;

if nargin == 0
    fn_in = spm_select(1,'image','select image to threshold');
end

if ~isempty(Neg)
    orderEG = 'erode-grow';  % -> erode then grow
    if numel(Neg)~=2
        error('Wrong number of erode-grow steps.');
    end
elseif ~isempty(Nge)
    orderEG = 'grow-erode'; % -> grow then erode
    if numel(Nge)~=2
        error('Wrong number of grow-erode steps.');
    end
else
    orderEG = 0;  % -> no erode/grow operation
end

%% Get image in
if iscell(fn_in); fn_in = cell2mat(fn_in); end
if ischar(fn_in)
    V_in = spm_vol(fn_in);
    data = spm_read_vols(V_in);
    save_img = true;
else
    data = fn_in;
    save_img = false;
end
ddata = double(data>0);

%% Deal with erode/grow
% create neighbourhood
switch Nneighb
    case 26 % 26 neighbours, by corner
        neighb = ones(3,3,3);
    case 18 % 18 neighbours, by edge
        neighb = ones(3,3,3);
        neighb(1,1,1) = 0; neighb(1,1,3) = 0; neighb(1,3,1) = 0; neighb(1,3,3) = 0;
        neighb(3,1,1) = 0; neighb(3,1,3) = 0; neighb(3,3,1) = 0; neighb(3,3,3) = 0; 
    otherwise % 6 neighbours, by surface
        neighb = zeros(3,3,3); 
        neighb(:,2,2) = 1; neighb(2,:,2) = 1; neighb(2,2,:) = 1;
end

switch orderEG
    case 'erode-grow'  % -> erode then grow
        for ii=1:Neg(1) % erode
            ddata = imerode(~~ddata,neighb);
        end
        for ii=1:Neg(2) % grow
            ddata = imdilate(~~ddata,neighb);
        end       
    case 'grow-erode' % -> grow then erode
        for ii=1:Nge(1) % grow
            ddata = imdilate(~~ddata,neighb);
        end       
        for ii=1:Nge(2) % erode
            ddata = imerode(~~ddata,neighb);
        end
    otherwise
        % Do nothing. :-)
end

%% Deal with cluster sizes
[clustered_map,num] = spm_bwlabel(double(ddata),18);
% note we use edge connection

% check number of elements per clusters
extent_map = zeros(size(data));
clustered_map = clustered_map(:);
nv = histc(clustered_map,0:num);
[~,idxall]=sort(clustered_map,'ascend');
idxall(1:nv(1)) = []; nv(1)=[]; % remove 1st bin, i.e. 0s
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

%% Deal with output
extent_map = extent_map .* data; % reasign the initial values
if save_img
    V_out = V_in;
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
