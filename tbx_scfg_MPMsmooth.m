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
% wcImg Warped tissue class images
% ---------------------------------------------------------------------
wcImg         = cfg_files;
wcImg.tag     = 'wcImg';
wcImg.name    = 'Warped tissue class images';
wcImg.help    = {'Select the warped tissue classes.'};
wcImg.filter = 'image';
% wcImg.ufilter = '^smwc.*';
wcImg.ufilter = '^wc.*';
wcImg.num     = [2 3];

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
MPMsmooth.val    = {wMPM wcImg tpm_l fwhm};
% MPMsmooth.val    = {wMPM tpm_l fwhm tc_ind};
MPMsmooth.help   = {['Applying tissue spcific smoothing in order to',...
    'limit partial volume effect. This is specifically useful for ',...
    'the quantitative MPM images.']};
MPMsmooth.prog   = @crc_MPMsmooth;
MPMsmooth.vout   = @vout_MPMsmoothn;

end

%% OUTPUT function
%_______________________________________________________________________
function dep = vout_MPMsmoothn(job) %#ok<*INUSD>

img_c = 0;
for ii=1:numel(job.wMPM)
    for jj=1:numel(job.wcImg)
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


