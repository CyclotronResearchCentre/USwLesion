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
            fix_MPMintens(deblank(fn_in(ii,:)), opt.thrMPM(p_mtch), opt.prefix));
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
% On top, create an info-image of voxels that were "fixed", with a value
% of 1 if the voxel value was <0, or 2 if >thrMPM.

crt_mask = true;

% Load stuff
V = spm_vol(fn_in);
dd = spm_read_vols(V);
sz_dd = size(dd); dd = dd(:);

% Generation info-image
if crt_mask
    ll_fix = (dd<0) + (dd>thrMPM)*2;
    Vf = V;
    Vf.dt(1) = 2; % uint8
    Vf.fname = spm_file(V.fname,'prefix',prefix);
    Vf.descrip = 'fixed voxels, 1 if <0 and 2 if >thrMPM';
    Vf = spm_create_vol(Vf);
    Vf = spm_write_vol(Vf,reshape(ll_fix,sz_dd)); %#ok<*NASGU>
end

% Fix the image for the extreme values, <0 and >thr
dd = abs(dd); % Take the abs-value...
NaboveThr = sum(dd>thrMPM);
dd(dd>thrMPM) = thrMPM * (1 + randn(NaboveThr,1)*1e-3); % cap to max + small rand-value
dd = reshape(dd,sz_dd);

% if any(dd(:)==0)
%     % Fix for the small zero patches in some images:
%     % Principle
%     % do not worry about superlarge chunk of zero as this is outside the head.
%     % replace the zero's by the mean value of at least 9 non-zeros neighbours
%     % or
%     % start with the smaller cluster then iterate till everything is filled up.
%     sz_thr = 1000; % arbitrary maximum size of patch to fix
%     [L,num] = spm_bwlabel(double(~dd),18);
%     any_fix = false;
%     n_vx = zeros(num,1);
%     for ii=1:num
%         n_vx(ii) = sum(L(:)==ii);
%     end
%     
%     [sn_vx,li_cl] = sort(n_vx);
%     li_cl(sn_vx>sz_thr) = [];
%     sn_vx(sn_vx>sz_thr) = [];
%     
%     if ~isempty(li_cl)
%         for i_cl = li_cl'
%             dd = fix_ith_hole(i_cl,L,dd);
%         end
%     end
% end

% Save results
Vc = V;
Vc.fname = spm_file(V.fname,'prefix','t');
Vc = spm_create_vol(Vc);
Vc = spm_write_vol(Vc,dd);
fn_out = Vc.fname;

end

% %% FIXING the MPM maps for there holes (zero's) in the ICV volume
% function dd = fix_ith_hole(i_cl,L,dd)
% % use the image erode/grow to work from the outside in, propagating a "mean
% % value" of the outer values.
% 
% if numel(i_cl)~=1
%     aa=1;
% end
% 
% nvx = sum(L(:)==i_cl);
% if nvx==1
%     % easy case of isolated voxel
% else
%     % deal with other cases
% end
% 
% 
% end
