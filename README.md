# USwLesion
Unified Segmentation with lesions in the brain

The aim is to extend the "unified segmentation" (US, Ashburner et al. 2005) to brain images with lesional tissue. This was originally developed to process multiple sclerosis MR images. We are using the standard structural MRI but also [quantitative MR images](http://www.fil.ion.ucl.ac.uk/Research/physics_info/QuantMRI_VBM.html), aka. multi-parametric maps or MPM. There now exist a [`hMRI` toolbox](http://hmri.info) for the generation of these quantitative maps

This development should lead to an SPM12 comaptible toolbox with a matlabbatch interface.

Here is how the code is organized:
- the matlabbatch configuration files are all the 'tbx_cfg_\*' and 'tbx_scfg_\*' files
- the processing files are the 'crc_\*' files

## Matlabbatch:
There are a few processing processing modules, called from the main [`tbx_cfg_USwithLesion.m`](https://github.com/CyclotronResearchCentre/USwLesion/blob/master/tbx_cfg_USwithLesion.m) configuration file.
The main module deals with the "unified segmentation with lesion" and is batched through [`tbx_scfg_USwL.m`](https://github.com/CyclotronResearchCentre/USwLesion/blob/master/tbx_scfg_USwL.m). This is the extension of US to account for lesion(s), as indicated through an (approximate) lesion-mask image.

A second [`tbx_scfg_MPMsmooth.m`](https://github.com/CyclotronResearchCentre/USwLesion/blob/master/tbx_scfg_MPMsmooth.m) module 
deals with the tissue specific smoothing of quantitative MR images. This follows the paper by [Draganski et al, 2011](http://www.ncbi.nlm.nih.gov/pubmed/21277375). This is similar to the functionality available in the [`hMRI` toolbox](http://hmri.info).
A bunch of utility modules are also available as sub-modules: 
- fixing the lesion mask image, e.g. ensuring the image is binary (with 0/1 values) or removing "too small" blobs ([`tbx_scfg_Utils_FxLesMsk.m`](https://github.com/CyclotronResearchCentre/USwLesion/blob/master/tbx_scfg_Utils_FxLesMsk.m))
- fixing the qMRI/MPM data, e.g. cleaning up voxels with values <0 or too large occuring through the estimation in some noisy/low signal voxels ([`tbx_scfg_Utils_FxMPM.m`](https://github.com/CyclotronResearchCentre/USwLesion/blob/master/tbx_scfg_Utils_FxMPM.m))
- fixing the generated ICV mask image, e.g. erode/grow to get rid of out-of-the-brain smal blobs ([`tbx_scfg_Utils_FxICVmsk.m`](https://github.com/CyclotronResearchCentre/USwLesion/blob/master/tbx_scfg_Utils_FxICVmsk.m))
- extracting of lesion/healthy tissues parameters such as (relative) volumes and voxel intensities ([`tbx_scfg_ParEx.m`](https://github.com/CyclotronResearchCentre/USwLesion/blob/master/tbx_scfg_ParEx.m))

## Processing:
The main function [crc_USwL.m](https://github.com/CyclotronResearchCentre/USwLesion/blob/master/crc_USwL.m) does the segmentation itself. The central idea is to extend the (healthy brain) tissue probability maps (TPMs) with an extra tissue classe for lesion(s). Afterwards the SPM's US algorithm can be applied on the patient's data with his patient specific updated TPMs, i.e. "US with Lesion"
The function operates in the following steps:
  - **STEP 1**, apply the mask on the reference structural images in order to keep only healthy tissues visible
  - **STEP 2**, segment the masked structure, normalize the mask and smooth it. The result will be used as a new tissue probability map for the lesion.
  - **STEP 3**, update the TPMs to include a 7th tissue class, the lesion tissue (Note that the lesion is inserted in *3rd position*, between WM and CSF!). This is done by
    1. loading the original TPM's into memory (simple & easy but memory hungry...)
    2. scaling up the lesion-TPM and adjusting the healthy TPM in the intracranial volume
    3. ensuring a minium value for these 2 TPM's all over the image
    4. updating the TPM's by concatenating the lesion-TPM, replacing the healthy TPM with the new one(s), and adjusting overall (the summ of all TPM's at each voxel should be 1 of course!)
    5. saving the TPM's into a subject-specific TPM file, with the lesion-TPM in 3rd position (between the WM and CSF).
  - **STEP 4**, do the segmentation with the new subject-specific TPM. Dependign on the data at hand this can be a multichannel segmentation.
  - **STEP 5**, apply the deformation onto the selected images, plus a few others not used for the segmentation if requested, and create the warped images.

The parameter exctraction function [crc_ExtractParam_qMRIs.m](https://github.com/CyclotronResearchCentre/USwLesion/blob/master/crc_ExtractParam_qMRIs.m) aims at providing some 'summary statistics' on the patient images and the extracted tissues, including the lesion:
  1. total intracranial volume (tICV), to be used as reference volume
  2. match between mask and segmented lesion volume
  3. volumic information , in absolute or relative values for GM, WM and lesion
  4. qMRI values in the voxels of different tissue types
  5. some summary statistics: min, max , mean, median, std, skewness, kurtosis, p10, & p90 of the values extracted at 4., for each quantitative maps and tissue types.

## Rationale for the TPMs updating:
The lesion mask is included in the standard TPM following John Ashburner's advice on how to define TPMs. At each voxel: 
- no TPM is ever zero in the image, 
- the TPMs should sum up to 1. 
My personnel touch is to down scale healthy tissu TPM and give the complement to the lesion TPM. The lesion can be limited to WM, GM, GM+WM, GM+WM+CSF, depending on the the type of lesion, for example with MS then you should select 'WM only'. This ensures that the original healthy tissue TPM as well the lesion TPM are accounted for together.
