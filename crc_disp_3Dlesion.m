function [param_out] = crc_disp_3Dlesion(varargin)
% FUNCTION param_out = crc_disp_3Dlesion(Ps,Pl,param)
% 
% Plot a 3D visualisation of lesion/blob over a cut structural image.
% It is also possible to update an exising plot by passign a set of 
% parameters. This is particularly useful to create frames for a movie.
%
% INPUT for a new plot
% - Ps : filename of structural image (or structural image 3D matrix for
%        plot update)
% - Pl : filename of lesion mask image (not needed for plot update)
% - param = set of display parameters
%   .pov : point of view for camera, if left empty -> default value in code
%   .cut_ax : cuting perpendicular to axis 'x'/'y'/'z'
%   .cut_sl : cuting at slice (in mm)
%   .cut_dr : removing things >cut_sl (resp. <cut_sl) if 1 (resp. -1)
%   .int_thr : intensity threshold for surface/caps display
%   .col_isosurf : colour for structural volume isosurface, def. [1 .6 .6]
%   .col_blobs : colour for blobs, def. 'green'
%   .rotate3D : allow manual 3D rotation in plot, def. false
%   .p_fig : figure number (Def. create a new figure)
%   .fig_pos : figure position/size (Def. ~a quarter of the screen)
%
% OUPUT 
% - param_out : the same as the input parameters with some extra fields
%   .range_mm : figure range (in mm) along the x/y/z axis
%   .Vs : structure array containing structural image volume information
%   .p_fig : figure number
%   .fig_pos : figure position
%   .handles : structure with the following sub-fields
%       .p1/p2/p3 : pointer/handle to surface objects
%       .pl : pointer/handle to light object
%       .pa : pointer/handle to axis object
%
% or
%
% FUNCTION crc_disp_3Dlesion(param)
%
% Up date an existing plot. The update can consist in different cut view 
% (i.e. axis, slice and/or direction) and/or point of view.
%
% INPUT for plot update
% - param = the output parameter from anew plot, with updated
%       + cutting parameters (axis, slice and direction)
%       + poin of view (.pov field)
%       as needed for the updated plot.
%__________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by C. Phillips, 2015.
% Cyclotron Research Centre, University of Liege, Belgium


%% New plot or plot update?
if nargin>1
    new_plot = true;
else
    new_plot = false;
end

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

if new_plot
    if nargin==3
        param = crc_check_flag(param_def,varargin{3});
    else
        param = param_def;
    end
else
    param = varargin{1};
end

%% Check image input & load images
if new_plot
    Ps = varargin{1};
    Pl = varargin{2};
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
    vx_size = sqrt(sum(Vs.mat(:,1:3).^2)); % voxel size (in mm) for isosurf
    x_mm = reshape(xyz_mm(1,:),sz_vx);
    y_mm = reshape(xyz_mm(2,:),sz_vx);
    z_mm = reshape(xyz_mm(3,:),sz_vx);
    
    % Set Figure
    if isfield(param,'p_fig')
        figure(param.p_fig)
    else
        param.p_fig = figure;
    end
    if ~isfield(param,'fig_pos')
        scrsz = get(0,'ScreenSize');
        param.fig_pos = [scrsz(3)/4 scrsz(4)/4 scrsz(3)/3 scrsz(4)/2];
    end
    set(param.p_fig,'Position',param.fig_pos) 

    % Create & display patches for structural surface
    FVs = isosurface(x_mm, y_mm, z_mm, Ys, param.int_thr);
    vert = FVs.vertices;
    Nvert = size(vert,1);
    cdat = ones(size(FVs.vertices,1),1)*param.col_isosurf;
    p1 =  patch(FVs, 'FaceColor', 'flat', 'FaceVertexCData', cdat, ...
                'EdgeColor', 'none');
    axis tight;
    daspect(vx_size)
    axis(range_mm(:)')
    axis vis3d
    axis off
    xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
    colormap(gray(100))
    hold on
    % and lesion surface 
    FVl = isosurface(x_mm,y_mm,z_mm,Yl, .5);
    p2 = patch(FVl, ...
        'FaceColor', param.col_blobs, 'EdgeColor', 'none');
else
    % Collect the necessary bits!
    p1 = param.handles.p1;
    vert = get(p1,'Vertices');
    Nvert = size(vert,1);
    p2 = param.handles.p2;
    p3 = param.handles.p3;
    range_mm = param.range_mm;
    Vs = param.Vs;
    vx_size = sqrt(sum(Vs.mat(:,1:3).^2)); % voxel size (in mm)
end

%% Define surface cut and slice image
switch param.cut_ax
    case 'x'
        % external isosurface
        if param.cut_dr>0
            lTransp = find(vert(:,1)>param.cut_sl);
            l_angle = [90 0];
        else
            lTransp = find(vert(:,1)<param.cut_sl);
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
            lTransp = find(vert(:,2)>param.cut_sl);
            l_angle = [180 0];
        else
            lTransp = find(vert(:,2)<param.cut_sl);
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
            lTransp = find(vert(:,3)>param.cut_sl);
            l_angle = [0 90];
        else
            lTransp = find(vert(:,3)<param.cut_sl);
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

%% Complete plot: add slice + transparent surface
figure(param.p_fig)
% Make transparent the right bit
AlphaData = ones(Nvert,1); AlphaData(lTransp) = 0;
set(p1,'FaceVertexAlphaData',AlphaData, 'FaceAlpha', 'interp');
view(param.pov);
% Deal with cut surface
if new_plot
    % Cut surface  = volume cap
    p3 = surface(sX, sY, sZ, cap_img, 'LineStyle', 'none');
    hold off
else
    % update cut slice
    set(p3, 'XData', sX, 'YData', sY, 'ZData', sZ, 'CData', cap_img);
end

%% Deal with lighting
if new_plot
    % Add lights
    pl(1) = lightangle(param.pov(1),param.pov(2));
    pl(2) = lightangle(l_angle(1),l_angle(2));
    lighting gouraud
else
    % update lights
    pl = param.handles.pl;
    [lX,lY,lZ] = sph2cart(param.pov(1),param.pov(2),2e3);
    set(pl(1),'Position',[lX,lY,lZ]);
end

%% Output, only if not updating plot.
if nargout>0 && new_plot
    param_out = param;
    param_out.range_mm = range_mm;
    param_out.Vs = Vs;
    param_out.handles = struct('p1',p1,'p2',p2,'p3',p3,'pl',pl,'pa',gca);
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
