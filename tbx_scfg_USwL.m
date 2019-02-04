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
    ['This image is only used to map the lesion mask into MNI using a ' ...
    ' "cost function masking approach. This step is necessary to generate ' ...
    'the subjec specific TPM.' , ...
    'This would typically be a T1 or MT weighted image. ' ...
    ]};
imgRef.filter = 'image';
imgRef.ufilter = '.*';
imgRef.num     = [1 1];

% ---------------------------------------------------------------------
% imgStruc Structural images for multi-channel USwL
% ---------------------------------------------------------------------
imgStruc         = cfg_files;
imgStruc.tag     = 'imgStruc';
imgStruc.name    = 'Structural images for US-with-Lesion';
imgStruc.help    = {'Select the structural images .', ...
    ['These will be "multi-channel Unified-Segmented-with-Lesion" and '...
    'warped into MNI space. '] ...
    ['This ''mcUS'' includes a tissue class for the lesioned tissues. ' ...
    'If no image are selected, then the Reference image will be used ' ...
    'for the US-with-Lesion step.']};
imgStruc.filter = 'image';
imgStruc.ufilter = '.*';
imgStruc.num     = [0 Inf];
imgStruc.val       = {''};

% ---------------------------------------------------------------------
% imgOth Other structural images
% ---------------------------------------------------------------------
imgOth         = cfg_files;
imgOth.tag     = 'imgOth';
imgOth.name    = 'Other images';
imgOth.help    = {['Select the other images to be warped along but not ', ...
    'entering the US-with-Lesion.'], ...
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
biaswr.name    = 'Save bias corrected';
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
% bias_yes Apply bias correction
% ---------------------------------------------------------------------
bias_yes         = cfg_branch;
bias_yes.tag     = 'bias_yes';
bias_yes.name    = 'Apply bias correction';
bias_yes.val     = {biasreg biasfwhm biaswr};
bias_yes.help    = {'Bias correction options.'};

% ---------------------------------------------------------------------
% bias_no No bias correction
% ---------------------------------------------------------------------
bias_no         = cfg_const;
bias_no.tag     = 'bias_no';
bias_no.name    = 'No bias correction';
bias_no.val     = {0};
bias_no.help    = {'No bias correction.'};

% ---------------------------------------------------------------------
% bias Should intensity bias correction be applied?
% ---------------------------------------------------------------------
bias        = cfg_choice;
bias.tag    = 'bias';
bias.name   = 'Apply intensity bias correction to structural images in USwL';
bias.values = {bias_no, bias_yes};
bias.val    = {bias_no};
bias.help   = {...
    ['You may wish to correct intensity bias in the structural images ',...
    'during the (multi-channel) US-with-Lesion process.']};

% ---------------------------------------------------------------------
% Number of Gaussians per tissue class used to model the lesion
% ---------------------------------------------------------------------
NbGaussian         = cfg_entry;
NbGaussian.tag     = 'NbGaussian';
NbGaussian.name    = 'Number of Gaussians per tissue class';
NbGaussian.help    = {'Set the number of Gaussians per tissue class to model the intensity histograms.' ,...
    ['For the usual TPM''s as provided by SPM (or the eTPM from hMRI toolbox), ',...
    'the 7 numbers corresponds to the 7 tissue classess in the following order : ', ...
    'GM, WM, Lesion, CSF, Skull, Soft tissues, and Air.'], ...
    ['When using MPMs on "older" subject, it is useful to separately model ',...
    'the pallidum with an 8th tissue class, thus 8 number of Gausians.']};
NbGaussian.strtype = 'w';
NbGaussian.num     = [1 Inf];
NbGaussian.def     = @(val)crc_USwL_get_defaults('segment.NbGaussian', val{:});

% ---------------------------------------------------------------------
% ICVmsk Create ICV-mask and mask the Struct + Other
% ---------------------------------------------------------------------
ICVmsk         = cfg_menu;
ICVmsk.tag     = 'ICVmsk';
ICVmsk.name    = 'Mask the Struct & Other images by created ICV-mask';
ICVmsk.help    = {['An ICV mask is created from the reference structural ', ...
    'image and can be applied onto the Struct images before the segmentation itself. ',...
    'This cleans up the images quite a bit and ', ...
    'is equivalent to "skull stripping". This helps, in some cases, the ', ...
    'multi-channel segmentation of the MPMs.']
    ['After the main USwL step, a final ICV mask is also created and ',...
    'would be applied on all the images and tissue maps.']
    'Note that the (warped) ICV mask images are always created .'
    'The masked Struc/Other images are prefixed with ''k''.'};
ICVmsk.labels = {
    'No'
    'Yes'
    }';
ICVmsk.values = {0 1};
ICVmsk.def     = @(val)crc_USwL_get_defaults('segment.ICVmsk', val{:});

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

% ---------------------------------------------------------------------
% options Options
% ---------------------------------------------------------------------
options         = cfg_branch;
options.tag     = 'options';
options.name    = 'Options';
options.val     = {imgTpm NbGaussian tpm4lesion bias ICVmsk mrf };
options.help    = {'Some processing options.'};
%_______________________________________________________________________

%% EXEC function
% ---------------------------------------------------------------------
% USwL Unified segmentation with lesion mas
% ---------------------------------------------------------------------
USwL         = cfg_exbranch;
USwL.tag     = 'uswl';
USwL.name    = 'US with lesion';
USwL.val     = {imgMsk imgRef imgStruc imgOth options};
USwL.help    = {['Unified segmentation for images with lesions when an ',...
    'approximate mask is also provided. This mask is turned into a ',...
    '"tissue probability map" and added to SPM''s usual TPMs to form ',...
    'a subject-specific set of TPMs. Then "Unified Segmentation ',...
    'is applied.'],...
    '',...
    'All images are assumed to be already coregistered!'};
USwL.prog = @tbx_run_USwL;
USwL.vout = @vout_USwL;

end

%% OUTPUT functions
%_______________________________________________________________________
function dep = vout_USwL(job) %#ok<*INUSD>
% Collecting the output from crc_USwL function:
% - fn_out :
%       ICVmsk      : intra-cranial volume mask, as generated from USwL
%       wICVmsk     : intra-cranial volume mask, as generated from USwL
%       Struc_i     : fixed (masked a/o modulated) i^th structural images, if created
%       Oth_i       : masked i^th other image
%       wStruc_i    : warped (masked/modulated) i^th structural image
%       wOth_i      : warped (masked) i^th other image
%       TPMl        : subject specific TPM with lesion
%       segmImg     : structure with posterior tissue probabilities
%           c(i)    : class #i in subject space (1-4)
%           wc(i)   : class #i in MNI space (1-4)
%           mwc(i)  : modulated class #i in MNI space (1-4)
%           rc(i)   : DARTEL  ready class #i in subject space (1-3)

% ICV mask
cdep            = cfg_dep; %#ok<*AGROW>
cdep.sname      = 'Subject''s ICV mask, native space';
cdep.src_output = substruct('.','ICVmsk');
cdep.tgt_spec   = cfg_findspec({{'filter','nifti'}});

% wICV mask
cdep(end+1)          = cfg_dep; %#ok<*AGROW>
cdep(end).sname      = 'Subject''s ICV mask, template space';
cdep(end).src_output = substruct('.','wICVmsk');
cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});

% ICV-masked Struc & Other
if job.options.ICVmsk || ... % masking or bias corrected
        (isfield(job.options.bias,'bias_yes') && ...
         job.options.bias.bias_yes.biaswr(2))
    for ii=1:max(numel(job.imgStruc),1) % At least one image, the Ref struct.
        cdep(end+1)          = cfg_dep; %#ok<*AGROW>
        cdep(end).sname      = sprintf('Corrected Struc #%d',ii);
        cdep(end).src_output = substruct('.',sprintf('Struc_%d',ii));
        cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});
    end
    if ~isempty(job.imgOth) && ...
            ~( strcmp(job.imgOth,'<UNDEFINED>') || isempty(job.imgOth{1}) )
        for ii=1:numel(job.imgOth)
            cdep(end+1)          = cfg_dep; %#ok<*AGROW>
            cdep(end).sname      = sprintf('Corrected Other #%d',ii);
            cdep(end).src_output = substruct('.',sprintf('Oth_%d',ii));
            cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});
        end
    end
end

% Warped Struc & Other
for ii=1:max(numel(job.imgStruc),1) % At least one image, the Ref struct.
    cdep(end+1)          = cfg_dep; %#ok<*AGROW>
    cdep(end).sname      = sprintf('Warped Struc image #%d',ii);
    cdep(end).src_output = substruct('.',sprintf('wStruc_%d',ii));
    cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end
if ~isempty(job.imgOth) && ...
        ~( strcmp(job.imgOth,'<UNDEFINED>') || isempty(job.imgOth{1}) )
    for ii=1:numel(job.imgOth)
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('Warped Other image #%d',ii);
        cdep(end).src_output = substruct('.',sprintf('wOth_%d',ii));
        cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});
    end
end

% Segmented images in subject space
for ii=1:4
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = sprintf('c%d image',ii);
    cdep(end).src_output = substruct('.','segmImg','.',sprintf('c%d',ii));
    cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end

% Warped segmented images
for ii=1:4
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = sprintf('wc%d image',ii);
    cdep(end).src_output = substruct('.','segmImg','.',sprintf('wc%d',ii));
    cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end

% Modulated warped segmented images
for ii=1:4
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = sprintf('mwc%d image',ii);
    cdep(end).src_output = substruct('.','segmImg','.',sprintf('mwc%d',ii));
    cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end

% Dartel-ready segmented images in subject space
for ii=1:3
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = sprintf('c%d image',ii);
    cdep(end).src_output = substruct('.','segmImg','.',sprintf('rc%d',ii));
    cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});
end

% TPM_lesion
cdep(end+1) = cfg_dep;
cdep(end).sname      = 'Subject''s TPM with Lesion';
cdep(end).src_output = substruct('.','TPMl');
cdep(end).tgt_spec   = cfg_findspec({{'filter','nifti'}});


dep = cdep;

end
%_______________________________________________________________________
