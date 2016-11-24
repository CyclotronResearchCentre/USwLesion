function lesvol = crc_lesion_volumes(img)
% Calculate volumic information about the segmented binarized lesion image:
% total lesional volume, number of lesions and volume of each individual 
% lesion.
% 
% FORMAT
%   lesvol = crc_lesion_volumes(img)
% 
% INPUT
% img    3D image array or filename to image
% 
% OUTPUT
% lesvol structure with lesional volumes information
%   .TotVol     total lesional volume
%   .LesVol     number of lesions 
%   .LesNum     volume of each individual lesion
%
% NOTE:
% If a 3D volume array is passed, then volumes are expressed in voxels. 
% If an image filename is passed, then volumes are expressed in mm3 based
% on the image voxel-to-world mapping.
%__________________________________________________________________________
% Copyright (C) 2016 Cyclotron Research Centre

% Written by Christophe Phillips
% University of Liege, Belgium

if ischar(img)
    % load image array and find voxel volume
    V_in = spm_vol(img);
    img = spm_read_vols(V_in);
    vx_vol = abs(det(V_in.mat));
else
    vx_vol = 1;
end

% Total volume
lesvol.TotVol = sum(img(:)>0) * vx_vol;

% Number of lesions and their individual volume
% connectivity criterion: 6 (surface), 18 (edge) or 26 (corner)
[L,num] = spm_bwlabel(double(img),18);
lesvol.LesNum = num;
LesVol = zeros(1,num);
for ii=1:num
    LesVol(ii) = sum(L(:)==ii) * vx_vol;
end
lesvol.LesVol = sort(LesVol);

end
