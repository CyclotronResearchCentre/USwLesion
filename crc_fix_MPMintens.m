function fn_out = crc_fix_MPMintens(fn_in,opt)
% Function to 'fix' some issues with MPM images intensities.
% Those values should not be negative (they're quantitiave physical values)
% and remain below some upper limit.
% -> take the abs-value if <0 and put a cap on values.
% The cap depends on the image type, which is picked up based on the
% filename, which typically ends with A, MT, R1 or R2*.
%
% INPUT
% - fn_in   : filename of MPM images to fix
% - opt     : option structure
%       strMPM  : filename parts used to pick the image type
%       thrMPM  : Max cap for the corresponding image type
%       prefix  : prefix added to filename ('fx_' by default)
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
    'thrMPM', crc_USwL_get_defaults('tMPM.thrMPM'));
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
        fn_tmp = char( fn_tmp , ...
            fix_MPMintens(deblank(fn_in(ii,:)), ...
            opt.thrMPM(p_mtch), opt.prefix));
    else
        fprintf('\nCould not fix file : %s',fn_in(ii,:))
        fn_tmp = char( fn_tmp , deblank(fn_in(ii,:)));
    end
end

fn_out = fn_tmp(2:end,:);

end

%% =======================================================================
%% SUBFUNCTIONS
%% =======================================================================

function fn_out = fix_MPMintens(fn_in,thrMPM,prefix)
% Make sure that MPM intensities are within [0 thrMPM] by capping the
% values. The resulting image is written out with the prefix 't'.
% To visually check problematic values, create an info-image of voxels that
% were "fixed", with a value of
% - 1 if the voxel value was <0,
% - 2 if >thrMPM
% - 3 if equal to zero

% Some flags
crt_mask = true;
fix_zeros = true;
sz_thr = 25e3; % arbitrary maximum size of patch to fix

% Load stuff
V = spm_vol(fn_in);
dd = spm_read_vols(V);
sz_dd = size(dd); dd = dd(:);

% Generation info-image
if crt_mask
    ll_fix = (dd<0) + (dd>thrMPM)*2 + (dd==0)*3;
    Vf = V;
    Vf.dt(1) = 2; % uint8
    Vf.fname = spm_file(V.fname,'prefix',prefix);
    Vf.descrip = 'fixed voxels, 1 if <0, 2 if >thrMPM, 3 if =0';
    Vf = spm_create_vol(Vf);
    Vf = spm_write_vol(Vf,reshape(ll_fix,sz_dd)); %#ok<*NASGU>
end

% Fix the image for the extreme values, <0 and >thr
dd = abs(dd); % Take the abs-value...
NaboveThr = sum(dd>thrMPM);
dd(dd>thrMPM) = thrMPM * (1 + randn(NaboveThr,1)*1e-3); % cap to max + small rand-value
dd = reshape(dd,sz_dd);

if any(dd(:)==0) && fix_zeros
    get_at_it = true;
    n_zeros = sum(dd(:)==0);
    while get_at_it
        [L,num] = spm_bwlabel(double(~dd),18);
        n_vx = zeros(num,1);
        for ii=1:num
            n_vx(ii) = sum(L(:)==ii);
        end
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
                dd = fix_ith_hole(i_cl,n_vx(i_cl),L,dd,sz_dd);
            end
        end
        n_zeros_ith = sum(dd(:)==0);
        if n_zeros == n_zeros_ith
            get_at_it = false;
        else
            n_zeros = n_zeros_ith;
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

% Save results
Vc = V;
Vc.fname = spm_file(V.fname,'prefix','t');
Vc = spm_create_vol(Vc);
Vc = spm_write_vol(Vc,dd);
fn_out = Vc.fname;

end

%% FIXING the MPM maps for there holes (zero's) in the ICV volume
function dd = fix_ith_hole(i_cl,ni_vx,L,dd,sz_dd)
% i_cl = index of cluster
% ni_vx = number of voxels in cluster
% L = cluster indexes
% dd = image values
% sz_dd = image size

neighb18 = [ ...
    1 -1 0  0 0  0 1 -1 1 -1 -1  1 -1  1 0  0  0  0 ; ...
    0  0 1 -1 0  0 1 -1 0  0  1 -1  0  0 1 -1  1 -1 ; ...
    0  0 0  0 1 -1 0  0 1 -1  0  0  1 -1 1  1 -1 -1 ]';

ind_vx = find(L(:) == i_cl);
[ix,iy,iz] = ind2sub(sz_dd,ind_vx);

% get list of indexes of neighbours, check they're ok, then fix by mean
l_neighb_ind = cell(1,ni_vx);
for jj=1:ni_vx
    % find xyz index
    neighb_ind_xyz = neighb18 + ones(18,1)*[ix(jj) iy(jj) iz(jj)];
    % remove those outside the image
    neighb_ind_xyz( any(neighb_ind_xyz<1,2) | any(neighb_ind_xyz>ones(18,1)*sz_dd,2),:) = [];
    % get indexes
    l_neighb_ind{jj} = sub2ind(sz_dd,neighb_ind_xyz(:,1),neighb_ind_xyz(:,2),neighb_ind_xyz(:,3));
    
    % remove if one of those voxels from a zero-cluster
    l_neighb_ind{jj}(dd(l_neighb_ind{jj})==0) = [];
    % if more than half of possible neighbours OK -> fix with mean
    if numel(l_neighb_ind{jj})>9
        dd(ind_vx(jj)) = mean(dd(l_neighb_ind{jj}));
    else
%         fprintf('PROBLEM with voxel #%d (from %d).\n',ind_vx(jj),ni_vx);
    end
    
end


end

