function [p_handles, dat_out] = crc_disp_3Dlesion(Ps,Pl,param,other_in)

% FUNCTION [p_handles, dat_out] = crc_disp_3Dlesion(Ps,Pl,param,other_in)
%
% Plot a 3D visualisation of lesion/blob over a cut structural image. 
% It is also possible to update an exising plot by passign the data, 
% parameters and some extra input. This is particularly useful to create 
% frames for a movie.
%
% INPUT
% - Ps : filename of structural image (or structural image 3D matrix for
%        plot update)
% - Pl : filename of structural image (not needed for plot update)
% - param = set of display parameters
%   .pov : point of view for camera, if left empty -> default value in code
%   .cut_ax : cuting perpendicular to axis 'x'/'y'/'z'
%   .cut_sl : cuting at slice (in mm)
%   .cut_dr : removing things >cut_sl (resp. <cut_sl) if 1 (resp. -1)
%   .int_thr : intensity threshold for surface/caps display
%   .col_isosurf : colour for structural volume isosurface, def. [1 .6 .6]
%   .col_blobs : colour for blobs, def. 'green'
%   .rotate3D : allow manual 3D rotation in plot, def. false
% - other_in = structure with the following possible fields (optional)
%   .xyz_mm : coordinates in mm, of the voxels in struct image
%   .Vs : handle to structural image
%   .pf : pointer/handle to figure to use
%   .pa : pointer/handle to axis object -> update the plot only!
% 	.pl : pointer/handle to light object
%   .p1/p2/p3 : pointer/handle to surface objects
% 
% OUPUT
% - p_handles = handles to elements of the plot
%   .p1/p2/p3 : handles to the structural isosurface, blob isosurface and
%               structural cut surface resp.
%   .pl : handles to the lights
%   .pf : handle to the figure
%   .pa : handle to the axis
% - dat_out = data loaded/created in the execution
%   .Ys : structural image data matrix
%   .xyz_mm : structural image voxel coordinates (in mm)
%   .Vs : mapped structral volume
%   .Yl : blob/lesion image data matrix
%   .range_mm : axis ranges in mm
%__________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by C. Phillips, 2015.
% Cyclotron Research Centre, University of Liege, Belgium

%% Check parameters
param_def = struct(...
    'pov', [], ... % viewing azimuth and elevation angles -> default view
    'cut_ax', 'x', ... % cuting perpendicular to axis 'x'/'y'/'z'
    'cut_sl', 0, ...  % cuting at slice (in mm)
    'cut_dr', 1, ... % 1 (resp. -1) removing things >cut_sl (resp. <cut_sl)
    'int_thr', 10e-3, ... % Intensity threshold for surface/caps display
    'col_isosurf', [1 .6 .6], ... % Colour for structure isosurface
    'col_blobs', 'green', ... % Colour for blobs
    'rotate3D',false ... % allow manual 3D rotation
    );
if nargin>2
    param = crc_check_flag(param_def,param);
else
    param = param_def;
end

if nargin<4
    other_in = struct([]);
end

%% Check image input & load images
if ischar(Ps)
    % then load things
    Vs = spm_vol(Ps);
    [Ys,xyz_mm] = spm_read_vols(Vs);
else
    % Assume we're passing the values directly
    % and pick coordinates from 'other_in
    Ys = Ps;
    xyz_mm = other_in.xyz_mm;
    Vs = other_in.Vs;
end

if ischar(Pl)
    Vl = spm_vol(Pl);
    Yl = spm_read_vols(Vl);
else
    % Assume we're passing the values directly
    Yl = Pl;
end

% Image sizes
sz_vx = size(Ys); % size in voxels
range_mm = [ ...
    floor([min(xyz_mm(1,:)) min(xyz_mm(2,:)) min(xyz_mm(3,:))]) ; ...
    ceil([max(xyz_mm(1,:)) max(xyz_mm(2,:)) max(xyz_mm(3,:))])];
% image axis in mm
vx_size = sqrt(sum(Vs.mat(:,1:3).^2)); % voxel size (in mm) for isosurface
Ys_surf = Ys;
x_mm = reshape(xyz_mm(1,:),sz_vx);
y_mm = reshape(xyz_mm(2,:),sz_vx);
z_mm = reshape(xyz_mm(3,:),sz_vx);

%% Deal with cut position and orientation
switch param.cut_ax
    case 'x'
        % external isosurface
        if param.cut_dr>0
            Ys_surf(xyz_mm(1,:)>param.cut_sl) = NaN;
            l_angle = [90 0];
        else
            Ys_surf(xyz_mm(1,:)<param.cut_sl) = NaN;
            l_angle = [-90 0];
        end
        % cutting plane -> cap image
        [sY,sZ] = meshgrid(range_mm(1,2):vx_size(2):range_mm(2,2), ...
            range_mm(1,3):vx_size(3):range_mm(2,3) );
        sX = zeros(size(sY)) + param.cut_sl ;
        % pov
        if isempty(param.pov)
            param.pov = [90*sign(param.cut_dr)+30 15];
        end
        
    case 'y'
        % external isosurface
        if param.cut_dr>0
            Ys_surf(xyz_mm(2,:)>param.cut_sl) = NaN;
            l_angle = [180 0];
        else
            Ys_surf(xyz_mm(2,:)<param.cut_sl) = NaN;
            l_angle = [0 0];
        end
        % cutting plane -> cap image
        [sX,sZ] = meshgrid(range_mm(1,1):vx_size(1):range_mm(2,1), ...
            range_mm(1,3):vx_size(3):range_mm(2,3) );
        sY = zeros(size(sX)) + param.cut_sl ;
        % pov
        if isempty(param.pov)
            param.pov = [(param.cut_dr>0)*180+30 15];
        end
        
    case 'z'
        % external isosurface
        if param.cut_dr>0
            Ys_surf(xyz_mm(3,:)>param.cut_sl) = NaN;
            l_angle = [0 90];
        else
            Ys_surf(xyz_mm(3,:)<param.cut_sl) = NaN;
            l_angle = [0 -90];
        end
        % cutting plane -> cap image
        [sX,sY] = meshgrid(range_mm(1,1):vx_size(1):range_mm(2,1), ...
            range_mm(1,2):vx_size(2):range_mm(2,2) );
        sZ = zeros(size(sX)) + param.cut_sl ;
        % pov
        if isempty(param.pov)
            param.pov = [-120 sign(param.cut_dr)*30];
        end
        
    otherwise
        fprintf('\n Wrong axis definition!'); beep;
end
% get cut surface image values
sXYZ_vx = Vs.mat\[sX(:) sY(:) sZ(:) ones(numel(sX),1)]';
cap_img = spm_sample_vol(Vs,sXYZ_vx(1,:)',sXYZ_vx(2,:)',sXYZ_vx(3,:)',1);
cap_img = reshape(cap_img,size(sX));
cap_img(cap_img<param.int_thr) = NaN;

%% Dispaly
if isfield(other_in,'pf')
    figure(other_in.pf)
else
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/3 scrsz(4)/2])
end
if isfield(other_in,'pa')
    update_plot = true;
    %     delete(other_in.pa)
else
    update_plot = false;
end
if param.rotate3D
    rotate3d on
end

% Isosurface for structure
if ~update_plot
    p1 = patch(isosurface(x_mm,y_mm,z_mm,Ys_surf, param.int_thr), ...
        'FaceColor', param.col_isosurf, 'EdgeColor', 'none');
    axis tight;
    daspect(vx_size)
    axis(range_mm(:)')
    axis vis3d
    axis off
    xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
    colormap(gray(100))
    hold on
else
    % Update struct patch
    PV = isosurface(x_mm,y_mm,z_mm,Ys_surf, param.int_thr);
    set(other_in.p1,'Faces',PV.faces,'Vertices',PV.vertices)
end

view(param.pov);

if ~update_plot 
    % Isosurface for the lesion
    p2 = patch(isosurface(x_mm,y_mm,z_mm,Yl, .5), ...
        'FaceColor', param.col_blobs, 'EdgeColor', 'none');
end

if ~update_plot
    % Cut surface  = volume cap
    p3 = surface(sX,sY,sZ,cap_img,'LineStyle','none');
    hold off
else
    % update cut slice
    set(other_in.p3,'XData',sX,'YData',sY,'ZData',sZ, ...
        'CData',cap_img);
end

if ~update_plot
    % Add lights
    pl(1) = lightangle(param.pov(1),param.pov(2));
    pl(2) = lightangle(l_angle(1),l_angle(2));
    lighting gouraud
else
    % Update lights' position
    [lX,lY,lZ] = sph2cart(param.pov(1),param.pov(2),2e3);
    set(other_in.pl(1),'Position',[lX,lY,lZ]);
end

%% Output, only if not updating plot.
if nargout>0 && ~update_plot
    p_handles = struct('p1',p1,'p2',p2,'p3',p3,'pl',pl,'pf',gcf,'pa',gca);
end
if nargout>1 && ~update_plot
    % Collect 3D images values + coordinates
    dat_out = struct('Ys',Ys,'xyz_mm',xyz_mm,'Vs',Vs,'Yl',Yl, ...
        'range_mm',range_mm);
end

end

%%
% NOTE
% isonormals(x_mm,y_mm,z_mm,Ys_surf, p); % Not usable here 'cos grid is not
% necessarily uniform as requested for some 3D interpolation step...

% Display a 3D view of a brain image with lesion volume rendering
% This is based on the 'help isocaps' code from Matlab
%
% load mri
% D = squeeze(D);
% D(:,1:60,:) = [];
% p = patch(isosurface(D, 5), 'FaceColor', 'red', 'EdgeColor', 'none');
% p2 = patch(isocaps(D, 5), 'FaceColor', 'interp', 'EdgeColor', 'none');
% view(3); axis tight;  daspect([1 1 .4])
% colormap(gray(100))
% camlight; lighting gouraud
% isonormals(D, p);


%% BROLL

% pl2 = lightangle(l_angle(1),l_angle(2));
% delete(pl2)
% 
% delete(other_in.pl(2))
% other_in.pl(2) = lightangle(l_angle(1),l_angle(2))
% 
% pp = get(other_in.pl(1),'Position')
% [a,e,r] = cart2sph(pp(1),pp(2),pp(3));
% [aa,ee] = [a,e]./pi.*180