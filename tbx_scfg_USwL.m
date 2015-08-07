function USwL = tbx_scfg_USwL
% MATLABBATCH sub-configuration file.
% Doing the Unified Segmentation with updated TPMs, based on a preliminary
% lesion mask. The latter has to provides by the user!
%_______________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium


%% Input definitions
%_______________________________________________________________________

% ---------------------------------------------------------------------
% imgMsk Mask image
% ---------------------------------------------------------------------
imgMsk         = cfg_files;
imgMsk.tag     = 'imgMsk';
imgMsk.name    = 'Mask image';
imgMsk.help    = {'Select the "lesiosn mask" image. '
    'This should be a binary 1/0 image!' };
imgMsk.filter = 'image';
imgMsk.ufilter = '.*';
imgMsk.num     = [1 1];

% ---------------------------------------------------------------------
% imgRef Structural reference image
% ---------------------------------------------------------------------
imgRef         = cfg_files;
imgRef.tag     = 'imgRef';
imgRef.name    = 'Structural reference image';
imgRef.help    = {'Select the structural reference image. '
    'This image is used 1st to map the lesion mask into MNI and generate '
    'the subjec specific TPM.'};
imgRef.filter = 'image';
imgRef.ufilter = '.*';
imgRef.num     = [1 1];

% ---------------------------------------------------------------------
% imgMPM Structural quantitative images
% ---------------------------------------------------------------------
imgMPM         = cfg_files;
imgMPM.tag     = 'imgMPM';
imgMPM.name    = 'Structural quantitative images';
imgMPM.help    = {'Select the structural quantitative (MPM-VBQ) images .'
    'These will be segmented and warped into MNI space.'};
imgMPM.filter = 'image';
imgMPM.ufilter = '.*';
imgMPM.num     = [0 Inf];
imgMPM.val       = {''};

% ---------------------------------------------------------------------
% imgOth Other structural images
% ---------------------------------------------------------------------
imgOth         = cfg_files;
imgOth.tag     = 'imgOth';
imgOth.name    = 'Other structural images';
imgOth.help    = {'Select the other structural (e.g. FLAIR) images .'
    'These will be segmented and warped into MNI space.'};
imgOth.filter = 'image';
imgOth.ufilter = '.*';
imgOth.num     = [0 Inf];
imgOth.val       = {''};

% ---------------------------------------------------------------------
% img4US Images to use for the segmentation
% ---------------------------------------------------------------------
img4US         = cfg_menu;
img4US.tag     = 'img4US';
img4US.name    = 'Images to use for the segmentation';
img4US.help    = {'Choose which image(s) are used for the segmentation '
    'and estimation of the warping into MNI space.'};
img4US.labels = {
    'Structural reference only'
    'all MPMs [DEF]'
    'all MPMs + others'
    }';
img4US.values = {0 1 2};
img4US.val    = {1};

% ---------------------------------------------------------------------
% tpm4lesion Tissue probability map(s) affected by the lesion
% ---------------------------------------------------------------------
tpm4lesion         = cfg_menu;
tpm4lesion.tag     = 'tpm4lesion';
tpm4lesion.name    = 'Tissue probability map(s) affected by the lesion';
tpm4lesion.help    = {'Choose which tissue class(es) is(are) modified by the lesion'};
tpm4lesion.labels = {
    'GM only'
    'WM only'
    'WM+GM'
    'WM+GM+CSF'
    }';
tpm4lesion.values = {0 1 2 3};
tpm4lesion.val    = {1};

% ---------------------------------------------------------------------
% imgTpm TPM images
% ---------------------------------------------------------------------
imgTpm         = cfg_files;
imgTpm.tag     = 'imgTpm';
imgTpm.name    = 'Tissue probability maps';
imgTpm.help    = {'Select the TPM images. ' };
imgTpm.filter  = 'image';
imgTpm.ufilter = '.*';
imgTpm.num     = [1 1];
imgTpm.val{1}  = {fullfile(spm('dir'),'tpm','unwTPM_sl2.nii')};

% ---------------------------------------------------------------------
% thrMPM Threshold outlier values from the MPM
% ---------------------------------------------------------------------
thrMPM         = cfg_menu;
thrMPM.tag     = 'thrMPM';
thrMPM.name    = 'Thresholding the MPMs';
thrMPM.help    = {['Apply a threshold on the MPM''s to remove outlier ',...
    'values from the images.'],...
    ['Take the absolute value of voxels < 0. Voxels > thr are set to thr ',...
    '+ small random number. The ''thr'' value is defined for each ',...
    'MPM (A/MT/R1/R2) seperately.'],...
    'The modified MPM images are prefixed with ''t''.',...
    'A mask image is created to keep track of those "fixed" voxels".'};
thrMPM.labels = {
    'No'
    'Yes'
    }';
thrMPM.values = {0 1};
thrMPM.val    = {1};

% ---------------------------------------------------------------------
% ICVmsk Create ICV-mask and mask the MPMs
% ---------------------------------------------------------------------
ICVmsk         = cfg_menu;
ICVmsk.tag     = 'ICVmsk';
ICVmsk.name    = 'Create ICV-mask and mask the MPMs';
ICVmsk.help    = {'An ICV mask can be created from the reference structural '
    'image and applied onto the MPMs. This cleans up the images quite a bit and'
    'is equivalent to "skull stripping". This helps, in some cases, the '
    'multi-channel segmentation of the MPMs.'};
ICVmsk.labels = {
    'No'
    'Yes'
    }';
ICVmsk.values = {0 1};
ICVmsk.val    = {1};

% ---------------------------------------------------------------------
% options Options
% ---------------------------------------------------------------------
options         = cfg_branch;
options.tag     = 'options';
options.name    = 'Options';
options.val     = {imgTpm img4US tpm4lesion thrMPM ICVmsk};
options.help    = {'Some processing options.'};
%_______________________________________________________________________

%% EXEC function
% ---------------------------------------------------------------------
% USwL Unified segmentation with lesion mas
% ---------------------------------------------------------------------
USwL         = cfg_exbranch;
USwL.tag     = 'uswl';
USwL.name    = 'US with lesion';
USwL.val     = {imgMsk imgRef imgMPM imgOth options};
USwL.help    = {['Unified segmentation for images with lesions when an ',...
    'approximate mask is also provided. This mask is turned into a ',...
    '"tissue probability map" and added to SPM''s usual TPMs to form ',...
    'a subject-specific set of TPMs. Then "Unified Segmentation ',...
    'is applied.'],...
    '',...
    'All images are assumed to be already coregistered!'};
USwL.prog = @crc_USwL;
USwL.vout = @vout_USwL;

end

%% OUTPUT functions
%_______________________________________________________________________
function dep = vout_USwL(job) %#ok<*INUSD>

% TPM_lesion
cdep = cfg_dep;
cdep.sname      = 'Subject''s TPM with Lesion';
cdep.src_output = substruct('.','TPMl');
% cdep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
cdep.tgt_spec   = cfg_findspec({{'filter','nifti'}});

% ICV mask
cdep(end+1)          = cfg_dep; %#ok<*AGROW>
cdep(end).sname      = 'Subject''s ICV mask, native space';
cdep(end).src_output = substruct('.','ICVmsk');
cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});

% Thresholded MPM + Mask
if job.options.thrMPM
    for ii=1:numel(job.imgMPM)
        cdep(end+1)          = cfg_dep; %#ok<*AGROW>
        cdep(end).sname      = sprintf('Thresholded MPM image #%d',ii);
        cdep(end).src_output = substruct('.',sprintf('thrMPM%d',ii));
        cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});
    end
    for ii=1:numel(job.imgMPM)
        cdep(end+1)          = cfg_dep; %#ok<*AGROW>
%         cdep(end).sname      = sprintf('Mask thresholded MPM #%d',ii);
        cdep(end).sname      = sprintf('Threshold-mask MPM #%d',ii);
        cdep(end).src_output = substruct('.',sprintf('thrMPMmsk%d',ii));
        cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});
    end
end

% ICV-masked MPM
if job.options.ICVmsk
    for ii=1:numel(job.imgMPM)
        cdep(end+1)          = cfg_dep; %#ok<*AGROW>
        cdep(end).sname      = sprintf('ICV-masked MPM #%d',ii);
        cdep(end).src_output = substruct('.',sprintf('kMPM%d',ii));
        cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});
    end
end

% Warped MPM
if ~isempty(job.imgMPM)
    for ii=1:numel(job.imgMPM)
        cdep(end+1)          = cfg_dep; %#ok<*AGROW>
        cdep(end).sname      = sprintf('Warped MPM image #%d',ii);
        cdep(end).src_output = substruct('.',['wMPM',num2str(ii)]);
        cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});
    end
end
% Warped Other
if ~isempty(job.imgOth)
    for ii=1:numel(job.imgOth)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('Warped Other image #%d',ii);
        cdep(end).src_output = substruct('.',['wOth',num2str(ii)]);
        cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});
    end
end

% Segmented images in subject space
for ii=1:4
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = sprintf('c%d image',ii);
    cdep(end).src_output = substruct('.','segmImg','.',['c',num2str(ii)]);
    cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end

% Warped segmented images
for ii=1:4
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = sprintf('wc%d image',ii);
    cdep(end).src_output = substruct('.','segmImg','.',['wc',num2str(ii)]);
    cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end

% Modulated warped segmented images
for ii=1:4
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = sprintf('mwc%d image',ii);
    cdep(end).src_output = substruct('.','segmImg','.',['mwc',num2str(ii)]);
    cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end

dep = cdep;

end
%_______________________________________________________________________
