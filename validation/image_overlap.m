function [mJ,overlap] = image_overlap(img1,img2)

% This function computes the modified Jaccard index and the percentage 
% overlap between two binary 3D images.
%
% FORMAT: overlap = percent_overlap(img1,img2)
% INPUT: img1 and img2 are two image names or matrices (computes mJ)
%        img1 can also be seen as a source image and img2 as the reference
%        (ground truth) image to computes overlap of source to the reference
% OUTPUT: mJ is the modified Jaccard index
%         overlap is a structure with:
%                 overlap.tp: percentage of img1 roi inside img2 roi
%                             this is the true positive rate (sensitivity)
%                 overlap.fp: percentage of img1 roi inside img2 0  
%                             this is the false positive rate
%                 overlap.tn: percentage of img1 0 inside img2 0  
%                             this is the true negative rate
%                 overlap.fn: percentage of img1 0 inside img2 roi
%                             this is the false negative rate (specificity)
%                 overlap.mcc: Matthews correlation coefficient
%          <https://en.wikipedia.org/wiki/Matthews_correlation_coefficient>
%
% To execute this function using MRI images, the software SPM needs to be 
% in your matlab path - see SPM <http://www.fil.ion.ucl.ac.uk/spm/>
%
% Cyril R Pernet, 6 May 2016
% The University of Edinburgh, NeuroImaging Sciences

%%
% The Jaccard index (1910) is defined as
% J(A,B) = N(intersect(A,B)) / N(union(A,B))
% <https://en.wikipedia.org/wiki/Jaccard_index>
%%
% The modified Jaccard index follows Maitra (2010)
% W(A,B) = N(intersect(A,B) / N(A)+N(B)-N(intersect(A,B))
% <http://www.ncbi.nlm.nih.gov/pubmed/19963068>

if nargin == 0
    help image_overlap
    disp(' '); disp('a 2D exemple is displayed in the figure')

    % this is a simple code check
    index1 = 1:6; index2 = 7:12;
    figure; 
    for move = 1:4
        img1 = zeros(12,12); img2 = img1;
        img1(index1,index1) = 1;
        img2(index2,index2) = 1;
        [mJ,overlap] = image_overlap(img1,img2);
        subplot(2,2,move); imagesc(img1+img2); drawnow
        title(['overlap ' num2str(overlap.tp) '% mJ=' num2str(mJ)])
        index1 = index1+1;
        index2 = index2 -1;
    end
    
    return

elseif nargin ~=2
    error('two inputs are expected - FORMAT: overlap = percent_overlap(img1,img2)')
end

%%
% Check fisrt images in: they need to be of the same size and to have 0s

% check format
if ischar(img1)
    if exist(img1,'file')
    V1 = spm_vol(img1); 
    img1 =spm_read_vols(img1);
    else
       error(sprintf('the file %s doesn''t exist',fileparts(img1))) 
    end
end

if ischar(img2)
    if exist(img2,'file')
    V2 = spm_vol(img2); 
    img2 =spm_read_vols(img2);
    else
       error(sprintf('the file %s doesn''t exist',fileparts(img2))) 
    end
end

% check dimensions
if any(size(img1)~=size(img2))
    error('img1 and img2 are not of the same size')
end

%% 
% *Compute the overlap using images as binary vectors*
 
img1 = img1(:) > 0;
img2 = img2(:) > 0;
I = sum((img1+img2)==2);
mJ = I / (sum(img1)+sum(img2)-I);

%% 
% *Compute percentage of overalp*

if nargout == 2
    
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


