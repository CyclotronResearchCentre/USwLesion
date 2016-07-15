function varargout = imat(M,flag)

% FORMAT h = imat(M,flag)
%
% Dislay image of matrices in a new figure with or without colour bar.
% If M is a cell array of 2D matrices, the matrices are plot as 2D images 
% in a subplot table.
% If M is a 3D matrix, then each plane in the matrix is plot as a 2D image.
%
% INPUT:
%   M       : matrix(ces) to display, either as 2/3D array or cell array
%   flag    : display parameters
%       .dcbar  : display colorbar [1, default] or not [0]
%       .spcbar : use specified caxis scaling [cmin cmax] or not [default]
%       .eqax   : use 'image' axis [1, default] or not [0]
%       .contour: add contour lines [1] or not [0, default]
%       .cont_v : contour values
%       .daxis  : display axis [1] or not [0, default]
%       .labels : labels for images, one per cell/plane. If left empty 
%                 [default], use plane/cell number.
%       .tform  : format for label none [default] or set as a set of pairs
%                 of {'Property',PropertyValue'}
%       .nRC    : imposed number of rows and columns. [set auto, default]
%
% OUPTUT:
%   h       : handle to each axes
%
%__________________________________________________________________________
% Copyright (C) 2009 Cyclotron Research Centre

% Written by Christophe Phillips, c.phillips@ulg.ac.be
% Cyclotron Research Centre, University of Liege, Belgium
% $Id:$

%% Check input
if nargin<2, flag = []; end
def_flag = struct('dcbar',1,'spcbar',[],'eqax',1,'contour',0, ...
                    'daxis',0,'labels',[],'tform',[],'nRC',[]);
flag = crc_check_flag(def_flag,flag);

if nargin<1
    disp('Nothing to plot')
    help imat
    return
end

%% Number of images, row, col
if iscell(M)
    nM = length(M(:)) ;
    ncol = ceil(sqrt(nM));
    nrow = ceil(nM/ncol);
    mx_M = 1 ;
elseif length(size(M))==3
    nM = size(M,3) ;
    ncol = ceil(sqrt(nM));
    nrow = ceil(nM/ncol);
    mx_M = 2 ;
else
    nM = 1;
    mx_M = 0 ;
    ncol = 1;
    nrow = 1;
end

% Deal with fixed subplots
if ~isempty(flag.nRC)
%     if ncol*nrow==prod(flag.nRC)
    if prod(flag.nRC)>=nM
        nrow = flag.nRC(1);
        ncol = flag.nRC(2);
    end
end

%% Check if sub-image titles are provided
if mx_M && ~isempty(flag.labels)
    if length(flag.labels)~=nM
        flag.labels = [];        
    end
end

%% Display
figure,
if mx_M
    h = zeros(nM,1);
    for i=1:nM
        h(i) = subplot(nrow,ncol,i);
        
        if mx_M==1
            val2plot = M{i};
            if isempty(flag.labels)
                t_text = ['cell nr: ',num2str(i)];
            else
                t_text = flag.labels{i};
            end
        elseif mx_M==2
            val2plot = M(:,:,i);
            if isempty(flag.labels)
                t_text = ['plane nr: ',num2str(i)];
            else
                t_text = flag.labels{i};
            end
        end
        im_plot(val2plot,flag,t_text)
%         imagesc(val2plot)
%         title(t_text)
%         
%         if flag.eqax, axis image, end
%         if flag.dcbar, colorbar, end;
%         if ~isempty(flag.spcbar) && length(flag.spcbar)==2
%             caxis(flag.spcbar)
%         end
%         if flag.contour
%             hold on
%             if isfield(flag,'cont_v')
%                 contour(val2plot,flag.cont_v,'-k')
%             else
%                 contour(val2plot,'-k')
%             end
%         end
    end
else
    h = axes;
    im_plot(M,flag,[])
%     imagesc(M),
%     if flag.eqax, axis image, end
%     if flag.dcbar, colorbar, end;
%     if ~isempty(flag.spcbar) && length(flag.spcbar)==2
%         caxis(flag.spcbar)
%     end
%     if flag.contour
%         hold on
%         if isfield(flag,'cont_v')
%             contour(M,flag.cont_v,'-k')
%         else
%             contour(M,'-k')
%         end
%     end
end

if nargout>0
    varargout{1} = h;
end

return

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% SUBFUNCTION for individual plot
%%

function im_plot(val2plot,flag,t_text)

if any(size(val2plot)==1)
    val2plot = squeeze(val2plot);
end

imagesc(val2plot)
shading interp
if ~isempty(t_text)
    ht = title(t_text);
    if ~isempty(flag.tform)
        for ii=1:size(flag.tform,1)
            set(ht,flag.tform{ii,1},flag.tform{ii,2})
        end
    end
end

if flag.eqax, axis image, end
if flag.dcbar, colorbar, end;
if ~isempty(flag.spcbar) && length(flag.spcbar)==2
    caxis(flag.spcbar)
end
if ~flag.daxis
    axis off
end
if flag.contour
    hold on
    if isfield(flag,'cont_v')
        contour(val2plot,flag.cont_v,'-k')
    else
        contour(val2plot,'-k')
    end
end

return