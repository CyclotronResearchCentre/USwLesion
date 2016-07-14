%% Small script to test the image overlapping measures
%
% 0/ Generate two 3D-images with prespcified shapes
sz = 100;
% Empty images
img1 = zeros(sz,sz);
img2 = zeros(sz,sz);
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

% Display images
imat({img1,img2,img1+2*img2})

% 1/ Work out how to scale/better use the Hausdorff distance
[iBx1,iBy1,iBz1] = crc_borderVx(img1);
[iBx2,iBy2,iBz2] = crc_borderVx(img2);

% Get coordinates in mm
Bxyz1 = [iBx1' ; iBy1' ; iBz1'];
Bxyz2 = [iBx2' ; iBy2' ; iBz2'];

[mD,D12,D21] = crc_meanHausdorffDist(Bxyz1,Bxyz2); %#ok<*ASGLU>
mHd = mean(mD);
