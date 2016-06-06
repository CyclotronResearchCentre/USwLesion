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
%            otherwise any non-zero value is considered as 1.
%            By default, thr=0, i.e. no thresholding is performed.
%       .mask is a binary image indicating which pixels/voxels
%            have to be taken into account (name or matrix).
%            By default, no masking, i.e. [].
%       .v2r voxel-to-realworld coordinates transformation, for the case
%            where 3D arrays are passed directly.
%            By default, anisotropic of size 1mm3, i.e. eye(4).
%
% OUTPUT:
%   - overlap is a structure with the followign measures
%       .mJ:    the modified Jaccard index (see Ref here under)
%       .cm     confusion matrix [ TP FN ; FP TN ] counts
%       .tp:    percentage of img1 roi inside img2 roi this is the true
%               positive rate (sensitivity)
%       .fp:    percentage of img1 roi inside img2 0 this is the false
%               positive rate
%       .tn:    percentage of img1 0 inside img2 0 this is the true negative
%               rate
%       .fn:    percentage of img1 0 inside img2 roi this is the false
%               negative rate (specificity)
%       .mcc:   Matthews correlation coefficient  (see Ref here under)
%       .CK:    Cohen's Kappa
%       .mHd:   mean Hausdorff distance between the 2 surfaces.
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
% Cohen's Kappa is described here:
%   <https://en.wikipedia.org/wiki/Cohen's_kappa>
%
% Haussdorf distance:
%   <https://en.wikipedia.org/wiki/Hausdorff_measure>
% Here use a modified form where the mean distance between the surfaces is
% calculated, instead of taking the maximum, i.e. the "average symmetric
% surface distance" as in the MS lesion segmentation challenge
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
opt_def = struct('thr',0,'mask',[],'v2r',eye(4));
opt = crc_check_flag(opt_def,opt);

%%
% Check fisrt images in: they need to be of the same size and to have 0s
% vox-to-real world mapping is extracted from the 1st image, or set to the
% value in 'opt', possibly the default value eye(4)

% check format
if ischar(img1)
    if exist(img1,'file')
        V1 = spm_vol(img1);
        img1 =spm_read_vols(V1);
        v2r = V1.mat;
    else
        error('the file %s doesn''t exist',spm_file(img1,'filename'))
    end
else
    v2r = opt.v2r;
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
% *Apply the mask and consider images as binary images & vectors*

img1 = img1 > opt.thr; % \_ binarize
img2 = img2 > opt.thr; % /

if ~isempty(mask)
    % apply mask and vectorize
    img1(~mask(:)) = 0;
    img2(~mask(:)) = 0;
    vimg1 = img1(:); vimg1(mask(:)) = [];
    vimg2 = img2(:); vimg2(mask(:)) = [];
else
    % vectorize
    vimg1 = img1(:);
    vimg2 = img2(:);
end

%%
% *Compute Jaccard coefficient
I = sum((vimg1+vimg2)==2);
overlap.mJ = I / (sum(vimg1)+sum(vimg2)-I);

%%
% *Compute percentage of overalp*

% Get confusion matrix with same convention as in 
%   <https://en.wikipedia.org/wiki/Confusion_matrix>
% that is [ TP FN ; FP TN ]
% cm = [ sum(ismember(find(vimg1),find(vimg2))) ...
%        sum(ismember(find(vimg1==0),find(vimg2))) ; ...
%        sum(ismember(find(vimg1),find(vimg2==0))) ...
%        sum(ismember(find(vimg1==0),find(vimg2==0))) ];
% The following CM calculation is about 10x faster than the previous one.
TP = sum((vimg1+vimg2)==2);
FN = sum((~vimg1+vimg2)==2);
FP = sum((vimg1+~vimg2)==2);
TN = sum((~vimg1+~vimg2)==2);
overlap.cm = [ TP FN ; FP TN ];

% percentage of vimg1 roi inside vimg2 roi
tp = TP / sum(vimg2);
overlap.tp = tp*100;

% percentage of vimg1 0 inside vimg2 0
tn = TN / sum(vimg2==0);
overlap.tn = tn*100;

% percentage of vimg1 roi inside vimg2 0
fp = FP / sum(vimg2==0);
overlap.fp = fp*100;

% percentage of vimg1 0 inside vimg2 roi
fn = FN / sum(vimg2);
overlap.fn = fn *100;

% Matthews correlation coefficient
overlap.mcc = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));

% Cohen's Kappa
Po = (TP+TN)/sum(overlap.cm(:));
Pr = (TP+TN)*(TP+FP)/sum(overlap.cm(:))^2;
overlap.CK = (Po - Pr)/(1 - Pr);

%% Extract the mean Hausdorff distance
% Get border coordinates in vx
[iBx1,iBy1,iBz1] = crc_borderVx(img1);
[iBx2,iBy2,iBz2] = crc_borderVx(img2);

% Get coordinates in mm
Bxyz1_mm = v2r(1:3,1:3)* [iBx1' ; iBy1' ; iBz1'];
Bxyz2_mm = v2r(1:3,1:3)* [iBx2' ; iBy2' ; iBz2'];

[mD,D12,D21] = crc_meanHausdorffDist(Bxyz1_mm,Bxyz2_mm); %#ok<*ASGLU>
overlap.mHd = mean(mD);

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