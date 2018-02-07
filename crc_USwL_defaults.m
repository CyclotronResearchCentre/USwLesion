function crc_USwL_defaults
% Set the defaults which are used by the USwL toolbox
%__________________________________________________________________________
%
% If you want to customise some defaults for your installation, do not
% modify this file directly, but create a file named crc_my_USwL_defaults.m
% instead, accessible from MATLAB search path; e.g., it can be saved in
% MATLAB Startup Folder: userhome/Documents/MATLAB
% or with your toolbox folder: userhome/Documents/spm/toolbox/USwLesion
%
% Example: create the following file to change the image file extension:
% ------ file C:\Workspm12\toolbox\USwLesion\crc_USwL_my_defaults.m -------
% global uswl_def
% uswl_def.segment.mrf = 0;
%--------------------------------------------------------------------------
%
% crc_USwL_defaults should not be called directly in any script or function
% To get/set the defaults, use crc_uswl_get_defaults.
%
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging
% Copyright (C) 2015 Cyclotron Research Centre

% Volkmar Glauche
% Then modified for use with the USwLesion toolbox by Christophe Phillips
% Cyclotron Research Centre, University of Liege, Belgium

global uswl_def

% Parameters for the segmentation with lesion 
%==========================================================================
uswl_def.segment.imgTpm     = {fullfile(spm('dir'),'tpm','TPM.nii')};
uswl_def.segment.tpm4lesion = 1; % TPM(s) affected by the lesion
%                                   (0, GM; 1, WM; 2, GM+WM; 3, GM+WM+CSF) 
uswl_def.segment.biasreg    = 1e-05; % almost nothing
uswl_def.segment.biasfwhm   = Inf;   % no bias correction
uswl_def.segment.biaswr     = [0 0]; % Not saving bias corrected/field images
uswl_def.segment.NbGaussian = [2 2 3 2 3 4 2]; % GM/WM/lesion/CSF/skull/scalp/air
uswl_def.segment.thrMPM     = 1;
uswl_def.segment.ICVmsk     = 1;
uswl_def.segment.thrLesion  = 0;
uswl_def.segment.mrf        = 2;
uswl_def.segment.cleanup    = 0;
uswl_def.segment.scDefReg   = 1; % scaling of warping regularisation,
% A value <1 leads to more freedom for the warping.
% This is useful for specific subjects where an extra bit of flexibility is
% useful for this preliminary CFM-segmentation step.
% Example:
% Some multiple sclerosis (MS) patients have large lesion along the 
% ventricles, which themselves are enlarged. With the standard CFM and 
% warping regularization, then the true shape of the ventricles is not
% properly captured and the lesion probability overlaps alrgely with prior
% CSF map, leading to erroneous USwL later on.

% Parameters for the segmentation of masked anatomical reference (to build
% the updated TPM)
%==========================================================================
uswl_def.msksegm.imgTpm     = {fullfile(spm('dir'),'tpm','TPM.nii')};
uswl_def.msksegm.biasreg    = 1e-05; % almost nothing, assuming we use MPMs
uswl_def.msksegm.biasfwhm   = Inf;   % no bias correction
uswl_def.msksegm.biaswr     = [0 0]; % Not saving bias corrected/field images
uswl_def.msksegm.NbGaussian = [2 2 2 3 4 4]; % GM/WM/CSF/skull/scalp/air
uswl_def.msksegm.mrf        = 2;
uswl_def.msksegm.cleanup    = 0;
uswl_def.msksegm.native     = [[1 0];[1 0];[1 0];[0 0];[0 0];[0 0]]; 
uswl_def.msksegm.scDefReg   = 1; % scaling of warping regularisation,
% A value <1 leads to more freedom for the warping.
% This is useful for specific subjects where an extra bit of flexibility is
% useful for this preliminary CFM-segmentation step.
% Example:
% Some multiple sclerosis (MS) patients have large lesion along the 
% ventricles, which themselves are enlarged. With the standard CFM and 
% warping regularization, then the true shape of the ventricles is not
% properly captured and the lesion probability overlaps alrgely with prior
% CSF map, leading to erroneous USwL later on.

% Processing parameters for the creation of the updated TPM
%==========================================================================
uswl_def.uTPM.minVol        = 8; % volume of lesion patch must be > minVol (in mm^3)
uswl_def.uTPM.nDilate       = 2; % # of dilation step
uswl_def.uTPM.smoKern       = 2; % smoothing (in mm) of the warped lesion mask
uswl_def.uTPM.tpm_ratio     = 100; % ratio of lesion/tpm
uswl_def.uTPM.min_tpm       = 1e-6; % minimum value of tpm overall, as in standard spm's TPM
uswl_def.uTPM.min_tpm_icv   = 1e-3; % minimum value of tpm in intracranial volume
uswl_def.uTPM.b_write       = [0 0]; % not writing bias corrected images, as in standard spm's TPM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Prevent users from making direct calls to spm_defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent runOnce
try %#ok<*TRYNC>
    if ~isdeployed && isempty(runOnce)
        d = dbstack;
        if isempty(intersect({'spm','crc_USwL_get_defaults'},{d.name}))
            fprintf(['Direct calls to crc_USwL_defauts are deprecated.\n' ...
                'Please use crc_USwL_get_defaults instead.\n']);
            runOnce = 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Execute user-specified defaults files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def = uswl_def;
user_defaults = 'crc_USwL_my_defaults.m';
if exist(user_defaults,'file')
    if isdeployed && exist(fullfile(spm('Dir'),'toolbox','USwLesion',user_defaults),'file')
        user_defaults_file = cellstr(fullfile(spm('Dir'),'toolbox','USwLesion',user_defaults));
    else
        user_defaults_file = which(user_defaults,'-ALL');
    end
    for i=1:numel(user_defaults_file)
        try
            spm('run', user_defaults_file{i});
        catch %#ok<*CTCH>
            lr = lasterror; %#ok<*LERR>
            warning(lr.message);
        end
    end
    if spm_check_version('matlab','8.0') >= 0, my_isequaln = @isequaln;
    else my_isequaln = @isequalwithequalnans; end
    if ~my_isequaln(def,uswl_def)
        fprintf('Defaults settings have been modified by file(s):\n');
        for i=1:numel(user_defaults_file)
            fprintf('  %s\n',user_defaults_file{i});
        end
        fn0 = fieldnames(def);
        mf = fn0(~cellfun(@(x) my_isequaln(def.(x),uswl_def.(x)),fn0));
        if ~isempty(mf)
            fprintf('Modified fields: ');
            for i=1:numel(mf)
                fprintf('%s ',mf{i});
            end
            fprintf('\n');
        end
    end
end

%% NOTICE
% This file is largely inspired from routines available the SPM toolbox:
% http://www.fil.ion.ucl.ac.uk/spm/