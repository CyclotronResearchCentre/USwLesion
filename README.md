# USwLesion
Unified Segmentation with lesions in the brain

The aim is to extend the "unified segmentation" (US, Ashburner et al. 2005) to brain images with lesional tissue. This was originally developed to process multiple sclerosis MR images. We are using the standard structural MRI but also [quantitative MR images](http://www.fil.ion.ucl.ac.uk/Research/physics_info/QuantMRI_VBM.html), aka. multi-parametric maps or MPM. Becasue we are dealing with VBQ/MPM data we also include the specific smoothing proposed by [Draganski et al, 2011](http://www.ncbi.nlm.nih.gov/pubmed/21277375)

This development should lead to an SPM12 comaptible toolbox with a matlabbatch interface.

Here is how the code is organized:
- the matlabbatch configuration files are all the 'tbx_cfg_*' and 'tbx_scfg_*' files
- the processing files are the 'crc_*' files

## Matlabbatch:
The processing is split into 3 processing modules, called from the main [tbx_cfg_USwithLesion.m](https://github.com/CyclotronResearchCentre/USwLesion/blob/master/tbx_cfg_USwithLesion.m) configuration file:
- [tbx_scfg_USwL.m](https://github.com/CyclotronResearchCentre/USwLesion/blob/master/tbx_scfg_USwL.m) deals with the segmentation of images, accounting for lesions. The lesioned area are indicated through an (approximate) lesion-mask image.
- [tbx_scfg_MPMsmooth.m](https://github.com/CyclotronResearchCentre/USwLesion/blob/master/tbx_scfg_MPMsmooth.m) deals with the smoothing of quantitative MR images/MPM. This follows the paper by [Draganski et al, 2011](http://www.ncbi.nlm.nih.gov/pubmed/21277375)
- [tbx_scfg_ParEx.m](https://github.com/CyclotronResearchCentre/USwLesion/blob/master/tbx_scfg_ParEx.m) deals with the extraction of lesion/healthy tissues parameters such as (relative) volumes and voxel intensities.

## Processing:
Like the matlabbatch, the code is organised in 3 main functions:
- [crc_USwL.m](https://github.com/CyclotronResearchCentre/USwLesion/blob/master/crc_USwL.m) does the segmentation itself. The lesion mask is included in the standard TPM following John Ashburner's advice on how to define TPM: no TPM is ever zero in the image, the TPM's should sum up to 1. My personnel touch is to down scale healthy tissu TPM and give the complement to the lesion TPM. The lesion can be limited to WM, GM, GM+WM, GM+WM+CSF, depending on the the type of lesion, for example with MS then you should select 'WM only'. This ensures that the original healthy tissue TPM as well the lesion TPM are accounted for together (see below for more explanations). The function operates in several steps:
  - **STEP 1**, clean up the mask image by
    - removing small lesion patches from mask (set as a minimum number of voxels per lesion region)
    - then growing the lesion region (set as a numbre of dilation step)
  - **STEP 2**, apply the mask on the reference structural images in order to keep only healthy tissues visible
  - **STEP 3**, segment the masked structure, normalize the cleaned up mask (the one before dilation!) and smooth it. The result will be used as a new tissue probability map for the lesion.
  - **STEP 4**, update the TPMs to include a 7th tissue class, the lesion tissue (Note that the lesion is inserted in *3rd position*, between WM and CSF!). This is done by
    1. loading the original TPM's into memory (simple & easy but memory hungry...)
    2. scaling up the lesion-TPM and adjusting the healthy TPM in the intracranial volume
    3. ensuring a minium value for these 2 TPM's all over the image
    4. updating the TPM's by concatenating the lesion-TPM, replacing the healthy TPM with the new one(s), and adjusting overall (the summ of all TPM's at each voxel should be 1 of course!)
    5. saving the TPM's into a subject-specific TPM file, with the lesion-TPM in 3rd position (between the WM and CSF).
  - **STEP 5**, do the segmentation with the new subject-specific TPM. Dependign on the data at hand this can be a multichannel segmentation (use the 'structural reference only', 'all MPMs', or 'all MPMs + others').
  - **STEP 6**, apply the deformation onto the MPMs and create the warped MPMs
- [crc_MPMsmooth.m](https://github.com/CyclotronResearchCentre/USwLesion/blob/master/crc_MPMsmooth.m) operates the specific smoothing for MPM/VBQ data. See [Draganski et al, 2011](http://www.ncbi.nlm.nih.gov/pubmed/21277375).
- [crc_ExtractParam.m](https://github.com/CyclotronResearchCentre/USwLesion/blob/master/crc_ExtractParam.m) extracts some parameters from the (smoothed) segmented MPM. For the moment only these parameters are extracted:
  1. total intracranial volume (tICV), to be used asreference volume
  2. the match between mask and segmented lesion volume, using the [Dice coefficient](http://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient)
  3. percentage volume of lesion in tICV and in WM volume
  4. MPM voxel values in GM, WM and lesion (for each map!)

New features to be added at some point...

## Rationale for the TPMs updating:
tba.
