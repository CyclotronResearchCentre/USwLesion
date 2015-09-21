%% Script to launch the 3D display of lesion over strucutral image
%__________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by C. Phillips, 2015.
% Cyclotron Research Centre, University of Liege, Belgium

% Select the 2 images
% Ps = spm_select; % structural image
% Pl = spm_select; % lesion image, binary or to be thresholded.
% Parameters

clear

if ispc
    Ps = 'C:\Dropbox\Work\3_data\MSdata_ELommers\Processed\kts7374-0005-00001-000176-00_MT_EPIB1.nii';
    Pl = 'C:\Dropbox\Work\3_data\MSdata_ELommers\Processed\c3kts7374-0005-00001-000176-00_MT_EPIB1.nii';
else
    Ps = '/Users/chrisp/Dropbox/Work/3_data/MSdata_ELommers/Processed/kts7374-0005-00001-000176-00_MT_EPIB1.nii';
    Pl = '/Users/chrisp/Dropbox/Work/3_data/MSdata_ELommers/Processed/c3kts7374-0005-00001-000176-00_MT_EPIB1.nii';
end

% Rotate the camera around the main axis, or not
cam_rotation = true;
% Keep on rotating to return to original view
do2rounds = true;
% Filename of video created
fn_vid = 'disp3D.avi';

%% Define Z-cut (i.e. horizontal cut) video
z_step = 1 ; % steps for slices (in mm)
z_list = -75:z_step:75;
nFrames = numel(z_list);

% Define first plot -> to load data and handles to objects
param = struct(...
    'cut_ax', 'z', ... % cuting perpendicular to axis 'x'/'y'/'z'
    'cut_sl', z_list(1), ...  % cuting at slice (in mm) if >0 (<0) removes the > part
    'cut_dr', 1, ... % 1 (resp. -1) removing things >cut_sl (resp. <cut_sl)
    'int_thr', 10e-3 ... % Intensity threshold for surface/caps display
    );
[param_out] = crc_disp_3Dlesion(Ps,Pl,param);
% Define rotating point-of-view (pov)
if cam_rotation
    param.pov = [-120 sign(param.cut_dr)*30];
    az_list = param.pov(1) + (0:nFrames-1)/nFrames*360;
end
% Add the 2nd part
if do2rounds
    nFrames = nFrames*2;
    z_list = [z_list fliplr(z_list)];
    az_list = [az_list az_list];
end

% Prepare video
vidObj = VideoWriter(fn_vid);
vidObj.FrameRate = 20;
vidObj.Quality = 95;
open(vidObj);

% Create frames and save them
param = param_out;
for ii = 1:nFrames
    param.cut_sl = z_list(ii);
    if cam_rotation,
        param.pov(1) = az_list(ii);
    end
    % Update figure
    fprintf('\n Slice %c = %d', param.cut_ax, param.cut_sl);
    crc_disp_3Dlesion(param);
    % Write each frame to the file.
    currFrame = getframe(param.p_fig);
    writeVideo(vidObj,currFrame);
end
fprintf('\n')
% Close the file.
close(vidObj);

