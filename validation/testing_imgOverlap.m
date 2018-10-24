%% Small script to test the image overlapping measures
%_______________________________________________________________________
% Copyright (C) 2016 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

%% 0/ Generate 3D-images with prespcified shapes
sz = 100;
% Empty images
img1 = zeros(sz,sz);
img1b = zeros(sz,sz);
img2 = zeros(sz,sz);
img2b = zeros(sz,sz);
img3 = zeros(sz,sz);
img3b = zeros(sz,sz);
[X,Y] = meshgrid(1:sz,1:sz);

% Add elements in img1: 2 disks of radius 20, centered on [25,25] & [75,75]
R = [20 20 10];
Xc = [25 75 25];
Yc = [25 75 75];
Nc = numel(R);
for ii=1:2
    Dc = sqrt((X-Xc(ii)).^2+(Y-Yc(ii)).^2);
    img1(Dc<=R(ii))=1;
end
% img1b is the same as img1 but with one extra small disk, not matching img2
for ii=1:Nc
    Dc = sqrt((X-Xc(ii)).^2+(Y-Yc(ii)).^2);
    img1b(Dc<=R(ii))=1;
end

% Add elements in img2, 4 squares:
% - 1 larger but matching the 1st disk, centered around[25,25], side 41
% - 2 smaller but matching the 2nd disk, side 15, centre [65,75] & [85,75]
% - 1 medium not matching, side 31, centre [75,25]
S = ([41 15 15 31]-1)/2;
Xs = [25 65 85 75];
Ys = [25 75 75 25];
Ns = numel(S);
for ii=1:Ns
    Ds = max(abs(X-Xs(ii)),abs(Y-Ys(ii)));
    img2(Ds<=S(ii))=1;
end
% img2b is like img2 but with only the overlapping parts, i.e. 1st 3 squares
for ii=1:3
    Ds = max(abs(X-Xs(ii)),abs(Y-Ys(ii)));
    img2b(Ds<=S(ii))=1;
end

% Add elements in img3/3b, each has one disk not overlapping with the one
% from the other image.
Dc = sqrt((X-Xc(1)).^2+(Y-Yc(1)).^2);
img3(Dc<=R(1))=1;
Dc = sqrt((X-Xc(2)).^2+(Y-Yc(2)).^2);
img3b(Dc<=R(2))=1;

% Display images 1/1b/2/2b
AllImg = {img1,img1b,img2,img2b,img1+2*img2,img1+2*img2b,img1b+2*img2,img1b+2*img2b};
flag.nRC = [4 2];
flag.labels = {'img1','img1b','img2','img2b','img1+img2','img1+img2b','img1b+img2','img1b+img2b'};
imat(AllImg,flag)

% imat({img1,img2})
% imat({img3,img3b})

%% 1/ Work out how to scale/better use the Hausdorff distance

% Find border voxels
[iBx1,iBy1,iBz1] = crc_borderVx(img1);
[iBx1b,iBy1b,iBz1b] = crc_borderVx(img1b);
[iBx2,iBy2,iBz2] = crc_borderVx(img2);
[iBx2b,iBy2b,iBz2b] = crc_borderVx(img2b);

% Get coordinates in mm
Bxyz1 = [iBx1' ; iBy1' ; iBz1'];
Bxyz1b = [iBx1b' ; iBy1b' ; iBz1b'];
Bxyz2 = [iBx2' ; iBy2' ; iBz2'];
Bxyz2b = [iBx2b' ; iBy2b' ; iBz2b'];

[mD12,D12,D21] = crc_meanHausdorffDist(Bxyz1,Bxyz2); %#ok<*ASGLU>
mHd12 = mean(mD12);
[mD12b,D12b,D2b1] = crc_meanHausdorffDist(Bxyz1,Bxyz2b); %#ok<*ASGLU>
mHd12b = mean(mD12b);
[mD1b2,D1b2,D21b] = crc_meanHausdorffDist(Bxyz1b,Bxyz2); %#ok<*ASGLU>
mHd1b2 = mean(mD1b2);
[mD1b2b,D1b2b,D2bb1] = crc_meanHausdorffDist(Bxyz1b,Bxyz2b); %#ok<*ASGLU>
mHd1b2b = mean(mD1b2b);

fprintf('\n')
fprintf('Mean distance from 1  to 2  : %2.2f and 2  to 1 : %2.2f, average %2.2f\n',mD12,mHd12)
fprintf('Mean distance from 1  to 2b : %2.2f and 2b to 1 : %2.2f, average %2.2f\n',mD12b,mHd12b)
fprintf('Mean distance from 1b to 2  : %2.2f and 2  to 1b: %2.2f, average %2.2f\n',mD1b2,mHd1b2)
fprintf('Mean distance from 1b to 2b : %2.2f and 2b to 1b: %2.2f, average %2.2f\n',mD1b2b,mHd1b2b)

% COMMENTS:
% As expected for the Hausdorff distance:
% - the H-distances 1-to-2 and 1-to-2b are the same, as well as 1b-to-2 and
%   1b-to-2b, 2-to-1 and 2-to-1b, and 2b-to-1 and 2b-to-1b
% - H-distances 2-to-1 and 1b-to-2 are relatively larger because one of the 
%   blob in in the 1st image does not not match any in the 2nd.
% - H-distance 2b-to-1 is of same order as 1-to-2b because all blobs in the 
%   1st image match at least one blob in the 2nd.
% - H-distance 2b-to-1 is larger than 1-to-2b as there are border voxels in 
%   img2b further away from img1 border than the other way round.
% - the average (or mean reciprocal distance) can be difficult to
%   interprete, with the lowest values only when all the blobs in both
%   images match one of the other image!

% Display the histograms
dMax = ceil(max([max(D12) max(D21) max(D1b2) max(D2b1)]));
nBin = 20; edg = (0:nBin)/nBin*dMax;
dHist = zeros(4,nBin+1);
dHist(1,:) = histc(D12,edg); dHist(2,:) = histc(D21,edg);
dHist(3,:) = histc(D1b2,edg); dHist(4,:) = histc(D2b1,edg);
figure
bar(edg,dHist')
legend('D1-2/2b','D2-1/1b','D1b-2/2b','D2b-1/1b')

% All together:
flag.nRC = [5 2];
flag.labels = {'img1','img1b','img2','img2b','img1+img2','img1+img2b','img1b+img2','img1b+img2b'};
imat(AllImg,flag)
figure(gcf)
subplot(5,1,5)
bar(edg,dHist')
legend('D1-2/2b','D2-1/1b','D1b-2/2b','D2b-1/1b')

% Play with images 3/3b
[iBx3,iBy3,iBz3] = crc_borderVx(img3);
[iBx3b,iBy3b,iBz3b] = crc_borderVx(img3b);
Bxyz3 = [iBx3' ; iBy3' ; iBz3'];
Bxyz3b = [iBx3b' ; iBy3b' ; iBz3b'];
[mD33b,D33b,D3b3] = crc_meanHausdorffDist(Bxyz3,Bxyz3b); %#ok<*ASGLU>
mHd33b = mean(mD33b);
figure, plot([D33b,D3b3])

% closest point of 3b to 3:
c = sqrt(R(2)^2/2);
x_cl = Xc(2)-c; y_cl = Yc(2)-c;
% distance from centre of img3 blob to closest point in img3b blob
dC3_cl3b = sqrt((Xc(1)-x_cl)^2+(Yc(1)-y_cl)^2);
fprintf('\n')
fprintf('Distance from centre to other image blob : %2.2f\n',dC3_cl3b);
fprintf('to be compared with mean H-distance : %2.2f\n',mHd33b);

%% 2. Use the image_overlap function
[mJ12,mHd12,overlap12]    = image_overlap(img1,img2);
[mJ12b,mHd12b,overlap12b] = image_overlap(img1,img2b);

fprintf('\n')
fprintf('Modfied Jaccard index, \n\t 1-to-2 : %2.2f and 1-to-2b : %2.2f\n',mJ12,mJ12b);
fprintf('Mean H-distance, \n\t 1-to-2 : %2.2f and 1-to-2b : %2.2f\n',mHd12,mHd12b);
fprintf('Matthews corcoef, \n\t 1-to-2 : %2.2f and 1-to_2b: %2.2f\n',overlap12.voxel.mcc,overlap12b.voxel.mcc);
fprintf('Cohen''s kappa, \n\t 1-to-2 : %2.2f and 1-to_2b: %2.2f\n',overlap12.voxel.CK,overlap12b.voxel.CK);
fprintf('Cluster TP, \n\t 1-to-2 : %2.2f and 1-to_2b: %2.2f\n',overlap12.cluster.tp,overlap12b.cluster.tp);
fprintf('Cluster FP, \n\t 1-to-2 : %2.2f and 1-to_2b: %2.2f\n',overlap12.cluster.fp,overlap12b.cluster.fp);

% COMMENTS:
% As expected the match is better at all level with 1-to-2b than 1-to-2


