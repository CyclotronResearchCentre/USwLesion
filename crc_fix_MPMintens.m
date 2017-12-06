function fn_out = crc_fix_MPMintens(fn_in,opt)
% Function to 'fix' some issues with MPM images intensities.
% Those values should not be negative, not be exactly zero (they're 
% quantitiave physical values), and remain below some upper limit.
% Proposed solution (implemented here:
% - take the abs-value if <0
% - put a cap on values. The cap depends on the image type (picked up based
%   on the filename, which typically ends with A, MT, R1 or R2*).
% - zero's are filed by averaging the vaue of non-zero neighbours.
% 
% INPUT
% - fn_in   : filename of MPM images to fix
% - opt     : option structure
%       strMPM   : filename parts used to pick the image type.
%                  Def. {'_A'    '_MT'    '_R1'    '_R2'}
%       thrMPM   : Max cap for the corresponding image type
%                  Def. [200 5 2000 2]
%       prefix   : prefix added to filename. Def. 'fx_'
%       crt_mask : create a map indicating the voxels that were fixed  
%                  Def. true
%       fix_zeros: fix the zero-holes in the image. Def. true.
%
% OUTPUT
% - fn_out  : filename of fixed MPM images
%
% TO DO:
% Create a matlabbatch module for this function -> USwLesion/Utils
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

% Deal with input
if nargin <2, opt = struct; end
opt_o = struct(...
    'prefix', 'fx_', ...
    'strMPM', {crc_USwL_get_defaults('tMPM.strMPM')}, ...
    'thrMPM', crc_USwL_get_defaults('tMPM.thrMPM'), ...
    'crt_mask', true, ...
    'fix_zeros', true);
opt = crc_check_flag(opt_o,opt);

nMPM = size(fn_in,1);
nSt = numel(opt.strMPM);
fn_tmp = [];
for ii=1:nMPM % Loop over MPM files
    mtch = zeros(nSt,1);
    for jj=1:nSt
        tmp = strfind(spm_file(fn_in(ii,:),'filename'),opt.strMPM{jj});
        if ~isempty(tmp), mtch(jj) = tmp(end); end % pick last index if many
    end
    [~,p_mtch] = max(mtch);
    if p_mtch
        if exist(spm_file(fn_in(ii,:),'prefix','t'),'file')
            % Already fixed
            fn_tmp = char( fn_tmp , ...
                spm_file(fn_in(ii,:),'prefix','t'));
        else
            % Do the job
            fn_tmp = char( fn_tmp , ...
                fix_MPMintens(deblank(fn_in(ii,:)), ...
                opt.thrMPM(p_mtch), opt.prefix, ...
                opt.crt_mask, opt.fix_zeros));
        end
    else
        % or return message
        fprintf('\nCould not fix file : %s',fn_in(ii,:))
        fn_tmp = char( fn_tmp , deblank(fn_in(ii,:)));
    end
end

fn_out = fn_tmp(2:end,:);

end

%% =======================================================================
%% SUBFUNCTIONS
%% =======================================================================

function fn_out = fix_MPMintens(fn_in,thrMPM,prefix,crt_mask,fix_zeros)
% Make sure that MPM intensities are within [0 thrMPM] by capping the
% values. The resulting image is written out with the prefix 't'.
% To visually check problematic values, create an info-image of voxels that
% were "fixed", with a value of
% - 1 if the voxel value was <0,
% - 2 if >thrMPM
% - 3 if was equal to zero but got fixed
% - 4 if was and still is equal to zero

% Some flags
sz_thr = 25e3; % arbitrary maximum size of patch to fix -> no big holes
% sz_thr = Inf; % arbitrary maximum size of patch to fix -> no limit

% Load stuff
V = spm_vol(fn_in);
dd = spm_read_vols(V);
sz_dd = size(dd); dd = dd(:);

% Generation info-image
if crt_mask
    ll_fix = (dd<0) + (dd>thrMPM)*2 + (dd==0)*3;
end

% Fix the image for the extreme values, <0 and >thr
dd = abs(dd); % Take the abs-value...
NaboveThr = sum(dd>thrMPM);
dd(dd>thrMPM) = thrMPM * (1 + randn(NaboveThr,1)*1e-3); % cap to max + small rand-value
dd = reshape(dd,sz_dd);

% Fix the zero-holes
if any(dd(:)==0) && fix_zeros
    get_at_it = true;
    n_zeros = sum(dd(:)==0);
    
    % Enlarge image with zero's around it
    % -> if picking outside original image, it's a zero
    dd0 = zeros(sz_dd+2);
    dd0(2:sz_dd(1)+1,2:sz_dd(2)+1,2:sz_dd(3)+1) = dd;

    while get_at_it
        [L,num] = spm_bwlabel(double(~dd),18);
        n_vx        = histc(L(:),(0:num) + 0.5);
        n_vx        = n_vx(1:end-1);
        [sn_vx,li_cl] = sort(n_vx);
        li_cl(sn_vx>sz_thr) = [];
        sn_vx(sn_vx>sz_thr) = [];
        
        % display some info
        n_holes = numel(sn_vx);
        n_disp = min(n_holes,3);
        fprintf('#Holes = %d, #voxels = %d, largest %d : ',n_holes,sum(sn_vx),n_disp)
        for ii = 1:n_disp, fprintf('%d ',sn_vx(end-ii+1)),end
        fprintf('\n');
        
        if isempty(li_cl)
            get_at_it = false;
        else
            for i_cl = li_cl'
                dd0 = fix_ith_hole(find(L(:) == i_cl),dd0,sz_dd);
            end
        end
        
        % recover fixed image
        dd = dd0(2:sz_dd(1)+1,2:sz_dd(2)+1,2:sz_dd(3)+1);

        n_zeros_ith = sum(dd(:)==0);
        if n_zeros == n_zeros_ith
            get_at_it = false;
        else
            n_zeros = n_zeros_ith;
        end
    end
end

% Complete info-image + save image
if crt_mask
    ll_fix = ll_fix + (dd(:)==0);
    Vf = V;
    Vf.dt(1) = 2; % uint8
    Vf.fname = spm_file(V.fname,'prefix',prefix);
    Vf.descrip = 'fx_vals, 1 if <0, 2 if >thrMPM, 3 if =0 fixed, 4 if =0';
    Vf = spm_create_vol(Vf);
    Vf = spm_write_vol(Vf,reshape(ll_fix,sz_dd)); %#ok<*NASGU>
end


% Save results
Vc = V;
Vc.fname = spm_file(V.fname,'prefix','t');
Vc = spm_create_vol(Vc);
Vc = spm_write_vol(Vc,dd);
fn_out = Vc.fname;

end

%% FIXING the MPM maps for their zero-holes in the image
function dd0 = fix_ith_hole(ind_vx,dd0,sz_dd)
% Input:
% - list of voxels to be fixed in current hole
% - data, extended with zeros around
% - data size, original!

% Minimal number of non-zero neighbours for the averaging
min_nz = 10; % -> the larger the more restrictive the filling.
max_zeros = 18-min_nz;

% Define 18-neighbourhood
neighb18 = [ ...
    1 -1 0  0 0  0 1 -1 1 -1 -1  1 -1  1 0  0  0  0 ; ...
    0  0 1 -1 0  0 1 -1 0  0  1 -1  0  0 1 -1  1 -1 ; ...
    0  0 0  0 1 -1 0  0 1 -1  0  0  1 -1 1  1 -1 -1 ]';

% Numbre of voxels in hole to fill
ni_vx = numel(ind_vx);

% xyz-coordinates of voxels to fix + their 18 neighbours
[ix,iy,iz] = ind2sub(sz_dd,ind_vx);
ix_n0 = bsxfun(@plus,ix,neighb18(:,1)')+1; % \
iy_n0 = bsxfun(@plus,iy,neighb18(:,2)')+1; % |- xyz coord in extended data
iz_n0 = bsxfun(@plus,iz,neighb18(:,3)')+1; % /

% Get coordinates in the extended image
l_neighb_in = sub2ind(sz_dd+2,ix_n0,iy_n0,iz_n0);
ind_vx_d0 = sub2ind(sz_dd+2,ix+1,iy+1,iz+1);

go_ahead = true;
while go_ahead
    % Pick up current values in all neighbours
    val_dd0 = reshape(dd0(l_neighb_in(:)),ni_vx,18);
    % Count number of zeros
    n0_val_dd0 = sum(~val_dd0,2);
    
    % Keep those with fewer zero neighbours
    n_min = min(n0_val_dd0);
    l_min = find(n0_val_dd0==n_min);
    
    % Go on or stop
    if ni_vx==0 || n_min>max_zeros
        go_ahead = false;
    else
        % do the job = fix those with mean of non-zero neighbours
        dd0(ind_vx_d0(l_min)) = sum(val_dd0(l_min,:),2)./(18-n0_val_dd0(l_min));
        
        % remove these fixed voxels from the list
        ind_vx_d0(l_min) = [];
        ni_vx = numel(ind_vx_d0);
        l_neighb_in(l_min,:) = [];
    end
end

end

% % Create map of hole-sizes
% nL = L;
% for ii=1:num
%     nL(L(:)==ii) = n_vx(ii);
% end
% Vn = V; Vn.fname = spm_file(V.fname,'prefix','nZ2_');
% Vn.descrip = 'number of zeros per cluster';
% Vn = spm_create_vol(Vn);
% Vn = spm_write_vol(Vn,nL); %#ok<*NASGU>
