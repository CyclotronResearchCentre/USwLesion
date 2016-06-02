function overlap = image_overlap(img1,img2,opt)

% This function computes the matching between 2 binary 3D images, based on
% different measures:
% - the modified Jaccard index
% - the percentage overlap between the two binary images.
%
% FORMAT:
%   overlap = percent_overlap(img1,img2,opt)
%
% INPUT:
%   - img1 and img2 are two image names or matrices (computes mJ)
%     img1 can also be seen as a source image and img2 as the reference
%     (ground truth) image to computes overlap of source to the reference
%   - opt is a structure with a few processing options
%       .thr is the threshold applied to both images to binarize them,
%            otherwise any non-zero value is considered as 1. By default
%            thr=0, i.e. no thresholding is performed.
%       .mask is a binary image indicating which pixels/voxels
%            have to be taken into account (name or matrix)
%
% OUTPUT:
%   - overlap is a structure with the followign measures
%       .tp:    percentage of img1 roi inside img2 roi this is the true
%               positive rate (sensitivity)
%       .fp: percentage of img1 roi inside img2 0 this is the false
%               positive rate
%       .tn: percentage of img1 0 inside img2 0 this is the true negative
%               rate
%       .fn: percentage of img1 0 inside img2 roi this is the false
%               negative rate (specificity)
%       .mcc: Matthews correlation coefficient  (see Ref here under)
%       .mJ: the modified Jaccard index (see Ref here under)
%
%
% REFERENCES:
% The Jaccard index (1910) is defined as
%   J(A,B) = N(intersect(A,B)) / N(union(A,B))
%   <https://en.wikipedia.org/wiki/Jaccard_index>
%
% The modified Jaccard index follows Maitra (2010)
%   W(A,B) = N(intersect(A,B) / N(A)+N(B)-N(intersect(A,B))
%   <http://www.ncbi.nlm.nih.gov/pubmed/19963068>
%
% The Matthews correlation coefficient is described here
%   <https://en.wikipedia.org/wiki/Matthews_correlation_coefficient>
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% NOTES:
% To execute this function using MRI images, the software SPM needs to be
% in your matlab path - see SPM <http://www.fil.ion.ucl.ac.uk/spm/>
%__________________________________________________________________________
%
% Cyril R Pernet, 6 May 2016
% The University of Edinburgh, NeuroImaging Sciences, UK.
%
% Christophe Phillips, for some added features
% University of Liege, GIGA - Research, Belgium

%%
%  Simple example with synthetic images
if nargin == 0
    display_help_and_example
    return
    
elseif nargin == 2
    % Set opt to [] in order to use the default options
    opt = [];
    
elseif nargin <2 || nargin > 3
    error('Two or three inputs are expected - FORMAT: [mJ,overlap] = percent_overlap(img1,img2,option)')
    
end

%%
% Default values
opt_def = struct('thr',0,'mask',[]);
opt = crc_check_flag(opt_def,opt);

%%
% Check fisrt images in: they need to be of the same size and to have 0s

% check format
if ischar(img1)
    if exist(img1,'file')
        V1 = spm_vol(img1);
        img1 =spm_read_vols(V1);
    else
        error('the file %s doesn''t exist',spm_file(img1,'filename'))
    end
end

if ischar(img2)
    if exist(img2,'file')
        V2 = spm_vol(img2);
        img2 =spm_read_vols(V2);
    else
        error('the file %s doesn''t exist',spm_file(img2,'filename'))
    end
end

% check dimensions
if any(size(img1)~=size(img2))
    error('img1 and img2 are not of the same size')
end

%%
% *Clean-up images using the mask*

mask = opt.mask;  % a mask field will be there, empty or not.
if ~isempty(mask) % load if not empty
    if ischar(mask)
        if exist(mask,'file')
            V3 = spm_vol(mask);
            mask = spm_read_vols(V3);
        else
            error('the file %s doesn''t exist',fileparts(mask))
        end
    end
    
    if any(size(img1)~=size(mask))
        error('the mask is not of the same size as input images')
    end
    
    % invert the mask as logical
    if isnan(sum(mask(:)))
        mask(~isnan(mask)) = 0;
        mask(isnan(mask)) = 1;
    else
        mask = (mask==0);
    end
end

%%
% *Compute the overlap using images as binary vectors*

if ~isempty(mask)
    % remove inverted mask region
    img1 = img1(:); img1(mask(:)) = [];
    img2 = img2(:); img2(mask(:)) = [];
    img1 = img1 > opt.thr;
    img2 = img2 > opt.thr;
else
    img1 = img1(:) > opt.thr;
    img2 = img2(:) > opt.thr;
end

I = sum((img1+img2)==2);
overlap.mJ = I / (sum(img1)+sum(img2)-I);

%%
% *Compute percentage of overalp*

% percentage of img1 roi inside img2 roi
tp = sum(ismember(find(img1),find(img2))) / sum(img2);
overlap.tp = tp*100;

% percentage of img1 0 inside img2 0
tn = sum(ismember(find(img1==0),find(img2==0))) / sum(img2==0);
overlap.tn = tn*100;

% percentage of img1 roi inside img2 0
fp = sum(ismember(find(img1),find(img2==0))) / sum(img2==0);
overlap.fp = fp*100;

% percentage of img1 0 inside img2 roi
fn = sum(ismember(find(img1==0),find(img2))) / sum(img2);
overlap.fn = fn *100;

% Matthews correlation coefficient
overlap.mcc = ((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SUBFUNCTIONS
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function display_help_and_example
% Simple function to display the help and an example of its use.
help image_overlap
disp(' '); disp('a 2D exemple is displayed in the figure')

% this is a simple code check
index1 = 1:6; index2 = 7:12;
figure;
for move = 1:4
    img1 = zeros(12,12); img2 = img1;
    img1(index1,index1) = 1;
    img2(index2,index2) = 1;
    overlap = image_overlap(img1,img2);
    subplot(2,2,move); imagesc(img1+img2); drawnow
    title(['overlap ' num2str(overlap.tp) '% mJ=' num2str(overlap.mJ)])
    index1 = index1+1;
    index2 = index2 -1;
end
end