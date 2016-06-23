function USwL = tbx_scfg_USwL
% MATLABBATCH sub-configuration file.
% Doing the Unified Segmentation with updated TPMs, based on a preliminary
% lesion mask. The latter has to provides by the user!
%_______________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium
% Cyril Pernet fixed typos and updated help
% Edinburgh Imaging, The University of Edinburgh


%% Input definitions
%_______________________________________________________________________

% ---------------------------------------------------------------------
% imgMsk Mask image
% ---------------------------------------------------------------------
imgMsk         = cfg_files;
imgMsk.tag     = 'imgMsk';
imgMsk.name    = 'Mask image';
imgMsk.help    = {'Select the "lesion mask" image. This should be a binary 1/0 image '...
    'the thresholded mask t_Mask based on a min size of 8 voxel and the'...
    ' slighlty dilated mask dt_Msk will be saved on the drive'};
imgMsk.filter = 'image';
imgMsk.ufilter = '.*';
imgMsk.num     = [1 1];

% ---------------------------------------------------------------------
% imgRef Structural reference image
% ---------------------------------------------------------------------
imgRef         = cfg_files;
imgRef.tag     = 'imgRef';
imgRef.name    = 'Structural reference image';
imgRef.help    = {'Select 1 reference structural image.', ...
    ['This image is used 1st to map the lesion mask into MNI and generate ' ...
    'the subjec specific TPM. This would typically be a T1 or MT weighted image. ' ...
    'The masked image k_img and the brain mask icv_k_ing are saved on the drive']};
imgRef.filter = 'image';
imgRef.ufilter = '.*';
imgRef.num     = [1 1];

% ---------------------------------------------------------------------
% imgMPM Structural quantitative images
% ---------------------------------------------------------------------
imgMPM         = cfg_files;
imgMPM.tag     = 'imgMPM';
imgMPM.name    = 'Structural quantitative images';
imgMPM.help    = {'Select the structural quantitative (MPM-VBQ) images .', ...
    'These will be segmented, including the lesioned area, and warped into MNI space.'};
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
imgOth.help    = {'Select the other structural images.', ...
    ['For example the FLAIR image that was used to create the "lesion mask". ' ...
    'These will be segmented and warped into MNI space.']};
imgOth.filter = 'image';
imgOth.ufilter = '.*';
imgOth.num     = [0 Inf];
imgOth.val       = {''};

%% OPTIONS for segmentation with lesion
%_______________________________________________________________________

% ---------------------------------------------------------------------
% imgTpm TPM images
% ---------------------------------------------------------------------
imgTpm         = cfg_files;
imgTpm.tag     = 'imgTpm';
imgTpm.name    = 'Tissue probability maps';
imgTpm.help    = {'Select the TPM images. ', ...
    ['A TPM_les file will be created from this file with lesions as', ...
     ' a new tissue class in 3rd position!']};
imgTpm.filter  = 'image';
imgTpm.ufilter = '.*';
imgTpm.num     = [1 1];
imgTpm.def     = @(val)crc_USwL_get_defaults('segment.imgTpm', val{:});

% ---------------------------------------------------------------------
% img4US Images to use for the segmentation
% ---------------------------------------------------------------------
img4US         = cfg_menu;
img4US.tag     = 'img4US';
img4US.name    = 'Images to use for the segmentation';
img4US.help    = {'Choose which image(s) are used for the segmentation '...
    'and estimation of the warping into MNI space.'};
img4US.labels = {
    'Structural reference only'
    'all MPMs [DEF]'
    'all MPMs + others'
    }';
img4US.values = {0 1 2};
img4US.def     = @(val)crc_USwL_get_defaults('segment.img4US', val{:});

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
tpm4lesion.def     = @(val)crc_USwL_get_defaults('segment.tpm4lesion', val{:});

%--------------------------------------------------------------------------
% biasreg Bias regularisation
%--------------------------------------------------------------------------
biasreg         = cfg_menu;
biasreg.tag     = 'biasreg';
biasreg.name    = 'Bias regularisation';
biasreg.help    = {
    ['MR images are usually corrupted by a smooth, spatially varying ', ...
    'artifact that modulates the intensity of the image (bias). These ', ...
    'artifacts, although not usually a problem for visual inspection, ', ...
    'can impede automated processing of the images.']
    ''
    ['An important issue relates to the distinction between intensity ', ...
    'variations that arise because of bias artifact due to the physics ', ...
    'of MR scanning, and those that arise due to different tissue ', ...
    'properties.  The objective is to model the latter by different ', ...
    'tissue classes, while modelling the former with a bias field. We ', ...
    'know a priori that intensity variations due to MR physics tend to ', ...
    'be spatially smooth, whereas those due to different tissue types ', ...
    'tend to contain more high frequency information. A more accurate ', ...
    'estimate of a bias field can be obtained by including prior ', ...
    'knowledge about the distribution of the fields likely to be ', ...
    'encountered by the correction algorithm. For example, if it is ', ...
    'known that there is little or no intensity non-uniformity, then it ', ...
    'would be wise to penalise large values for the intensity ', ...
    'non-uniformity parameters. This regularisation can be placed within ', ...
    'a Bayesian context, whereby the penalty incurred is the negative ', ...
    'logarithm of a prior probability for any particular pattern of non-uniformity.']
    'Knowing what works best should be a matter of empirical exploration.  For example, if your data has very little intensity non-uniformity artifact, then the bias regularisation should be increased.  This effectively tells the algorithm that there is very little bias in your data, so it does not try to model it.'
    }';
biasreg.labels = {
    'no regularisation (0)'
    'extremely light regularisation (0.00001)'
    'very light regularisation (0.0001)'
    'light regularisation (0.001)'
    'medium regularisation (0.01)'
    'heavy regularisation (0.1)'
    'very heavy regularisation (1)'
    'extremely heavy regularisation (10)'
    }';
biasreg.values = {
    0
    1e-05
    0.0001
    0.001
    0.01
    0.1
    1
    10
    }';
biasreg.def     = @(val)crc_USwL_get_defaults('segment.biasreg', val{:});

%--------------------------------------------------------------------------
% biasfwhm Bias FWHM
%--------------------------------------------------------------------------
biasfwhm         = cfg_menu;
biasfwhm.tag     = 'biasfwhm';
biasfwhm.name    = 'Bias FWHM';
biasfwhm.help    = {'FWHM of Gaussian smoothness of bias. If your ', ...
    'intensity non-uniformity is very smooth, then choose a large FWHM. ', ...
    'This will prevent the algorithm from trying to model out intensity ', ...
    'variation due to different tissue types. The model for intensity ', ...
    'non-uniformity is one of i.i.d. Gaussian noise that has been ', ...
    'smoothed by some amount, before taking the exponential. Note also ', ...
    'that smoother bias fields need fewer parameters to describe them. ', ...
    'This means that the algorithm is faster for smoother intensity non-uniformities.'};
biasfwhm.labels = {
    '30mm cutoff'
    '40mm cutoff'
    '50mm cutoff'
    '60mm cutoff'
    '70mm cutoff'
    '80mm cutoff'
    '90mm cutoff'
    '100mm cutoff'
    '110mm cutoff'
    '120mm cutoff'
    '130mm cutoff'
    '140mm cutoff'
    '150mm cutoff'
    'No correction'
    }';
biasfwhm.values = {
    30
    40
    50
    60
    70
    80
    90
    100
    110
    120
    130
    140
    150
    Inf
    }';
biasfwhm.def     = @(val)crc_USwL_get_defaults('segment.biasfwhm', val{:});

%--------------------------------------------------------------------------
% biaswr Save Bias Corrected
%--------------------------------------------------------------------------
biaswr         = cfg_menu;
biaswr.tag     = 'biaswr';
biaswr.name    = 'Save Bias Corrected';
biaswr.help    = {'This is the option to save a bias corrected version of your images from this channel, or/and the estimated bias field. MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images.  The bias corrected version should have more uniform intensities within the different types of tissues.'};
biaswr.labels = {
                'Save Nothing'
                'Save Bias Corrected'
                'Save Bias Field'
                'Save Field and Corrected'
                }';
biaswr.values = {
                [0 0]
                [0 1]
                [1 0]
                [1 1]
                }';
biaswr.def     = @(val)crc_USwL_get_defaults('segment.biaswr', val{:});

% ---------------------------------------------------------------------
% Number of Gaussians per tissue class used to model the lesion
% ---------------------------------------------------------------------
NbGaussian         = cfg_entry;
NbGaussian.tag     = 'NbGaussian';
NbGaussian.name    = 'Number of Gaussians per tissue class';
NbGaussian.help    = {'Set the number of Gaussians per tissue class to model the intensity histograms.' ,...
    ['The 7 numbers corresponds to the 7 tissue classess in the following order : ', ...
    'GM, WM, Lesion, CSF, Skull, Soft tissues, and Air.'], ...
    ['When using MPMs on "older" subject, it is useful to separately model ',...
    'the pallidum with an 8th tissue class, thus 8 number of Gausians.']};
NbGaussian.strtype = 'n';
NbGaussian.num     = [1 Inf];
NbGaussian.def     = @(val)crc_USwL_get_defaults('segment.NbGaussian', val{:});

% ---------------------------------------------------------------------
% thrMPM Threshold outlier values from the MPM
% ---------------------------------------------------------------------
thrMPM         = cfg_menu;
thrMPM.tag     = 'thrMPM';
thrMPM.name    = 'Thresholding the MPMs';
thrMPM.help    = {['Apply a threshold on the MPM''s to remove outlier ',...
    'values from the images before the segmentation itself.'],...
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
thrMPM.def     = @(val)crc_USwL_get_defaults('segment.thrMPM', val{:});

% ---------------------------------------------------------------------
% ICVmsk Create ICV-mask and mask the MPMs
% ---------------------------------------------------------------------
ICVmsk         = cfg_menu;
ICVmsk.tag     = 'ICVmsk';
ICVmsk.name    = 'Create ICV-mask and mask the MPMs & Other images';
ICVmsk.help    = {['An ICV mask can be created from the reference structural ', ...
    'image and applied onto the MPMs/other images before the segmentation itself. ',...
    'This cleans up the images quite a bit and ', ...
    'is equivalent to "skull stripping". This helps, in some cases, the ', ...
    'multi-channel segmentation of the MPMs.']
    ['Note that the TPMs are also masked so that the images to segment ',...
    'and TPMs match together.']
    'The masked MPM/other images are prefixed with ''k''.'};
ICVmsk.labels = {
    'No'
    'Yes'
    }';
ICVmsk.values = {0 1};
ICVmsk.def     = @(val)crc_USwL_get_defaults('segment.ICVmsk', val{:});

% ---------------------------------------------------------------------
% thrLesion Threshold lesion mask based on spatial extent
% ---------------------------------------------------------------------
thrLesion         = cfg_entry;
thrLesion.tag     = 'thrLesion';
thrLesion.name    = 'Thresholding the lesion mask';
thrLesion.help    = {'Apply a spatial threshold on the lesion tissue class c3: ',...
    '- 0 means no threshold;' ,... 
    '- k indicates that only clusters ofs ize larger or equal to k are kept;', ...
    '- ''Inf'' keeps only the largest cluster.'};
thrLesion.num    = [1 1];
thrLesion.def     = @(val)crc_USwL_get_defaults('segment.thrLesion', val{:});

%--------------------------------------------------------------------------
% mrf MRF Parameter
%--------------------------------------------------------------------------
mrf         = cfg_entry;
mrf.tag     = 'mrf';
mrf.name    = 'MRF Parameter';
mrf.help    = {'When tissue class images are written out, a few iterations of a simple Markov Random Field (MRF) cleanup procedure are run.  This parameter controls the strength of the MRF. Setting the value to zero will disable the cleanup.'};
mrf.strtype = 'r';
mrf.num     = [1 1];
mrf.def     = @(val)crc_USwL_get_defaults('segment.mrf', val{:});

%--------------------------------------------------------------------------
% cleanup Clean up any partitions
%--------------------------------------------------------------------------
cleanup         = cfg_menu;
cleanup.tag     = 'cleanup';
cleanup.name    = 'Clean Up';
cleanup.help    = {
    'This uses a crude routine for extracting the brain from segmented images.  It begins by taking the white matter, and eroding it a couple of times to get rid of any odd voxels.  The algorithm continues on to do conditional dilations for several iterations, where the condition is based upon gray or white matter being present.This identified region is then used to clean up the grey and white matter partitions.  Note that the fluid class will also be cleaned, such that aqueous and vitreous humour in the eyeballs, as well as other assorted fluid regions (except CSF) will be removed.'
    ''
    'If you find pieces of brain being chopped out in your data, then you may wish to disable or tone down the cleanup procedure. Note that the procedure uses a number of assumptions about what each tissue class refers to.  If a different set of tissue priors are used, then this routine should be disabled.'
    }';
cleanup.labels = {
    'Dont do cleanup'
    'Light Clean'
    'Thorough Clean'
    }';
cleanup.values = {0 1 2};
cleanup.def     = @(val)crc_USwL_get_defaults('segment.cleanup', val{:});

% ---------------------------------------------------------------------
% options Options
% ---------------------------------------------------------------------
options         = cfg_branch;
options.tag     = 'options';
options.name    = 'Options';
options.val     = {imgTpm img4US biasreg biasfwhm biaswr NbGaussian tpm4lesion thrMPM ICVmsk mrf cleanup thrLesion};
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
