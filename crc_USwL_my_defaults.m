function crc_USwL_my_defaults
% These are the CRC specific defaults for the USwLesion toolbox
% In particular, this is the 'best' defaults for the segmentation of MS
% lesions in quantitative MR images
% Reason for those changes:
% - need 3 Gaussians for the GM intensity modeling because of qMRI more
%   detailed (in term of intensity variations) images
% - use a qMRI based TPM set, including some subcortical GM areas
% - no need to clean up as a ICV mask is created and applied 
%_______________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

global uswl_def

% Get path to toolbox
Tbx_pth = spm_file(which('tbx_scfg_USwL.m'),'path');


% Parameters for the segmentation with lesion 
%==========================================================================
uswl_def.segment.NbGaussian = [2  2  3      2   3     4     2   2]; 
                             % GM/WM/lesion/CSF/skull/scalp/air/BG
uswl_def.segment.mrf        = 2;
uswl_def.segment.cleanup    = 0;
uswl_def.segment.imgTpm     = {fullfile(Tbx_pth,'eTPM','eTPM_wBG.nii')};
uswl_def.segment.thrMPM     = 1;
uswl_def.segment.ICVmsk     = 1;

uswl_def.segment.biasreg    = 1e-03; % almost nothing
uswl_def.segment.biasfwhm   = 60; % some bias correction
uswl_def.segment.biaswr     = [1 1]; % Saving bias corrected/field images

% Parameters for the segmentation of masked anatomical reference (to build
% the updated TPM)
%==========================================================================
uswl_def.msksegm.imgTpm     = {fullfile(spm('dir'),'tpm','eTPM.nii')};
uswl_def.msksegm.biasreg    = 1e-03; % small regularisation
uswl_def.msksegm.biasfwhm   = 60; % some bias correction
uswl_def.msksegm.biaswr     = [0 0]; % Not saving bias corrected/field images

return

%% NOTICE
% This file is largely inspired from routines available the SPM toolbox:
% http://www.fil.ion.ucl.ac.uk/spm/
