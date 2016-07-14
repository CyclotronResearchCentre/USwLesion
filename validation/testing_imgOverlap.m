%% Small script to test the image overlapping measures
%
% 0/ Generate 3D-images with prespcified shapes
sz = 100;
% Empty images
img1 = zeros(sz,sz);
img2 = zeros(sz,sz);
img3 = zeros(sz,sz);
[X,Y] = meshgrid(1:sz,1:sz);

% Add elements in img1: 2 disks of radius 20, centered on [25,25] & [75,75]
R = [20 20];
Xc = [25 75];
Yc = [25 75];
Nc = numel(R);
for ii=1:Nc
    Dc = sqrt((X-Xc(ii)).^2+(Y-Yc(ii)).^2);
    img1(Dc<=R(ii))=1;
end
% Xc = 75; Yc = 75;
% Dc = sqrt((X-Xc).^2+(Y-Yc).^2);
% img1(Dc<=R)=1;

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
% img3 is like img2 but with only the overlapping parts, i.e. 1st 3 squares
for ii=1:3
    Ds = max(abs(X-Xs(ii)),abs(Y-Ys(ii)));
    img3(Ds<=S(ii))=1;
end

% Display images
imat({img1,img2,img3, img1+img2})

% 1/ Work out how to scale/better use the Hausdorff distance

% Find border voxels
[iBx1,iBy1,iBz1] = crc_borderVx(img1);
[iBx2,iBy2,iBz2] = crc_borderVx(img2);
[iBx3,iBy3,iBz3] = crc_borderVx(img3);

% Get coordinates in mm
Bxyz1 = [iBx1' ; iBy1' ; iBz1'];
Bxyz2 = [iBx2' ; iBy2' ; iBz2'];
Bxyz3 = [iBx3' ; iBy3' ; iBz3'];

[mD12,D12,D21] = crc_meanHausdorffDist(Bxyz1,Bxyz2); %#ok<*ASGLU>
mHd12 = mean(mD12);
[mD13,D13,D31] = crc_meanHausdorffDist(Bxyz1,Bxyz3); %#ok<*ASGLU>
mHd13 = mean(mD13);

fprintf('\n')
fprintf('Mean distance from 1 to 2 : %2.2f and 2 to 1: %2.2f\n',mD12)
fprintf('Mean distance from 1 to 3 : %2.2f and 3 to 1: %2.2f\n',mD13)

% As expected :
% - the distance 1to2 and 1to3 are the same.
% - distance 2to1 is relatively larger because one of the blob in img2 does
%   not match any in img1
% - distance 3to1 is of same order as 1to3 because all blobs match at least
%   one blob of img1. It is larger as there are border voxels in img3
%   further away than in the other way round.

% Display the histograms
dMax = ceil(max([max(D12) max(D21) max(D13) max(D31)]));
nBin = 20; edg = (0:nBin)/nBin*dMax;
dHist = zeros(4,nBin+1);
dHist(1,:) = histc(D12,edg); dHist(2,:) = histc(D21,edg);
dHist(3,:) = histc(D13,edg); dHist(4,:) = histc(D31,edg);
figure, bar(edg,dHist')
legend('D12','D21','D13','D31')


