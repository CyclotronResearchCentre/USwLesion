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
uswl_def.segment.NbGaussian = [3 2 2 2 3 4 2];
uswl_def.segment.mrf        = 1;
uswl_def.segment.cleanup    = 0;
uswl_def.segment.imgTpm     = {fullfile(spm('dir'),'tpm','unwTPM_sl2.nii')};

return

%% NOTICE
% This file is largely inspired from routines available the SPM toolbox:
% http://www.fil.ion.ucl.ac.uk/spm/
