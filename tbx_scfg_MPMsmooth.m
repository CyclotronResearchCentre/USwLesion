function MPMsmooth = tbx_scfg_MPMsmooth
% MATLABBATCH sub-configuration file.
% Applying tissue spcific smoothing in order to limit partial volume 
% effect. This is specifically useful for the quantitative MPM images.
%_______________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

% ---------------------------------------------------------------------
% wMPM Warped MPM images
% ---------------------------------------------------------------------
wMPM         = cfg_files;
wMPM.tag     = 'wMPM';
wMPM.name    = 'Warped structural quantitative images';
wMPM.help    = {'Select the warped structural quantitative images.'};
wMPM.filter = 'image';
wMPM.ufilter = '^ws.*';
wMPM.num     = [1 Inf];
% wMPM.val       = {''};

% ---------------------------------------------------------------------
% mwcImg Warped tissue class images
% ---------------------------------------------------------------------
mwcImg         = cfg_files;
mwcImg.tag     = 'mwcImg';
mwcImg.name    = 'Modulated warped tissue class images';
mwcImg.help    = {'Select the modulate warped tissue classes.'};
mwcImg.filter = 'image';
% mwcImg.ufilter = '^smwc.*';
mwcImg.ufilter = '^mwc.*';
mwcImg.num     = [2 3];

% ---------------------------------------------------------------------
% tpm_l Subject's TPM with lesion class
% ---------------------------------------------------------------------
tpm_l         = cfg_files;
tpm_l.tag     = 'tpm_l';
tpm_l.name    = 'Subject''s TPM with lesion class';
tpm_l.help    = {'Subject''s TPM with lesion class.'};
tpm_l.filter = 'image';
tpm_l.ufilter = 'TPM.*';
tpm_l.num     = [1 1];

%--------------------------------------------------------------------------
% fwhm FWHM of smoothing kernel
%--------------------------------------------------------------------------
fwhm         = cfg_entry;
fwhm.tag     = 'fwhm';
fwhm.name    = 'FWHM of smoothing kernel';
fwhm.help    = {'FWHM of smoothing kernel to be applied'};
fwhm.strtype = 'r';
fwhm.num     = [1 1];
fwhm.val     = {8};

%% EXEC function
%----------------------------------------------------------------------
% MPMsmooth Partial volume smoothing, accounting for specific tissue 
%           classes "à la Bogdan", of the quantitative MPM images.
% ---------------------------------------------------------------------
MPMsmooth        = cfg_exbranch;
MPMsmooth.tag    = 'MPMsmooth';
MPMsmooth.name   = 'Partial volume smoothing';
MPMsmooth.val    = {wMPM mwcImg tpm_l fwhm};
% MPMsmooth.val    = {wMPM tpm_l fwhm tc_ind};
MPMsmooth.help   = {['Applying tissue spcific smoothing in order to',...
    'limit partial volume effect. This is specifically useful for ',...
    'the quantitative MPM images.']};
MPMsmooth.prog   = @tbx_run_MPMsmooth;
MPMsmooth.vout   = @vout_MPMsmoothn;

end

%% OUTPUT function
%_______________________________________________________________________
function dep = vout_MPMsmoothn(job) %#ok<*INUSD>

img_c = 0;
for ii=1:numel(job.wMPM)
    for jj=1:numel(job.mwcImg)
        img_c = img_c+1;
        cdep(img_c) = cfg_dep; %#ok<*AGROW>
        cdep(img_c).sname  = sprintf('fin%d_c%d image',[ii jj]);
        cdep(img_c).src_output =  ...
            substruct('.','fn','{}',{ii},'()',{jj});
%         cdep(img_c).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
        cdep(img_c).tgt_spec = cfg_findspec({{'filter','nifti'}});
    end
end

dep = cdep;
end
%_______________________________________________________________________

%% RUN function
%_______________________________________________________________________
function fn_out = tbx_run_MPMsmooth(job)
% Applying tissue specific smoothing in order to limit partial volume
% effect. This is specifically useful for the quantitative MPM images.

fn_wMPM = char(job.wMPM);
fn_mwTC = spm_file(char(job.mwcImg),'number','');

% Find list of tissue classes, tc_ind
tc_ind = [];
for ii = 1:size(fn_mwTC,1)
    p = strfind(deblank(fn_mwTC(ii,:)),'wc');
    tc_ind = [tc_ind str2double(deblank(fn_mwTC(ii,p+2)))]; %#ok<*AGROW>
end

% Collect the TPM_l
tpm = spm_file(char(job.tpm_l),'number','');
fn_lTPM = char;
for ii = tc_ind
    fn_lTPM = char( fn_lTPM , [tpm,',',num2str(ii)]);
end
fn_lTPM(1,:) = [];

opt_process = struct( ...
    'fwhm', job.fwhm, ...
    'tpm', fn_lTPM) ; % , ...

% Do it!
% fn_finMPM = crc_unifseg_MPMprocess(fn_wMPM, fn_mwTC, opt_process);
fn_out.fn = crc_uswl_MPMsmooth(fn_wMPM, fn_mwTC, opt_process);

end

