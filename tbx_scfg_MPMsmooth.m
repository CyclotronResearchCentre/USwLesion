function MPMsmooth = tbx_scfg_MPMsmooth

% ---------------------------------------------------------------------
% wMPM Warped MPM images
% ---------------------------------------------------------------------
wMPM         = cfg_files;
wMPM.tag     = 'wMPM';
wMPM.name    = 'Warped structural quantitative images';
wMPM.help    = {'Select the warped structural quantitative images.'};
wMPM.filter = 'image';
wMPM.ufilter = '^w.*';
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
wcImg.ufilter = '^smwc.*';
wcImg.num     = [2 3];

% ---------------------------------------------------------------------
% tpm_l Subject's TPM with lesion class
% ---------------------------------------------------------------------
tpm_l         = cfg_files;
tpm_l.tag     = 'tpm_l';
tpm_l.name    = 'Warped tissue class images';
tpm_l.help    = {'Subject''s TPM with lesion class.'};
tpm_l.filter = 'image';
tpm_l.ufilter = 'TPM.*';
tpm_l.num     = [1 1];

%--------------------------------------------------------------------------
% smooth Smoothing kernel
%--------------------------------------------------------------------------
smooth         = cfg_entry;
smooth.tag     = 'smooth';
smooth.name    = 'Smoothing kernel';
smooth.help    = {'Size of smoothing kernel to be applied'};
smooth.strtype = 'r';
smooth.num     = [1 1];
smooth.val     = {8};


%% EXEC function
%----------------------------------------------------------------------
% MPMsmooth Partial volume smoothing, accounting for specific tissue 
%           classes "à la Bogdan", of the quantitative MPM images.
% ---------------------------------------------------------------------
MPMsmooth        = cfg_exbranch;
MPMsmooth.tag    = 'MPMsmooth';
MPMsmooth.name   = 'Partial volume smoothing';
MPMsmooth.val    = {wMPM wcImg tpm_l smooth};
MPMsmooth.help   = {['Applying tissue spcific smoothing in order to',...
    'limit partial volume effect. This is specifically useful for ',...
    'the quantitative MPM images.']};
MPMsmooth.prog   = @crc_MPMsmooth;
MPMsmooth.vout   = @vout_MPMsmoothn;

end

%% OUTPUT function
%_______________________________________________________________________
function dep = vout_MPMsmoothn(job) %#ok<*INUSD>
dep            = cfg_dep;
dep.sname      = 'Velocity Kernel';
dep.src_output = substruct('.','fname','()',{':'});
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});
end%_______________________________________________________________________


