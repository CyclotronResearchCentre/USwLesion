%-----------------------------------------------------------------------
% Job saved on 22-Jul-2015 15:15:38 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (12.1)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
%__________________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by C. Phillips, 2015.
% Cyclotron Research Centre, University of Liege, Belgium

matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'Mask';
matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {'<UNDEFINED>'};
matlabbatch{2}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'MPM';
matlabbatch{2}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {
                                                                     '<UNDEFINED>'
                                                                     '<UNDEFINED>'
                                                                     '<UNDEFINED>'
                                                                     '<UNDEFINED>'
                                                                     }';
matlabbatch{3}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'Flair';
matlabbatch{3}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {'<UNDEFINED>'};
matlabbatch{4}.spm.tools.USwLtools.uswl.imgMsk(1) = cfg_dep('Named File Selector: Mask(1) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
matlabbatch{4}.spm.tools.USwLtools.uswl.imgRef(1) = cfg_dep('Named File Selector: MPM(1) - Files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
matlabbatch{4}.spm.tools.USwLtools.uswl.imgMPM(1) = cfg_dep('Named File Selector: MPM(1) - Files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
matlabbatch{4}.spm.tools.USwLtools.uswl.imgMPM(2) = cfg_dep('Named File Selector: MPM(2) - Files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{2}));
matlabbatch{4}.spm.tools.USwLtools.uswl.imgMPM(3) = cfg_dep('Named File Selector: MPM(3) - Files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{3}));
matlabbatch{4}.spm.tools.USwLtools.uswl.imgMPM(4) = cfg_dep('Named File Selector: MPM(4) - Files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{4}));
matlabbatch{4}.spm.tools.USwLtools.uswl.imgOth(1) = cfg_dep('Named File Selector: Flair(1) - Files', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
matlabbatch{4}.spm.tools.USwLtools.uswl.options.img4US = 1;
matlabbatch{4}.spm.tools.USwLtools.uswl.options.tpm4lesion = 1;
% matlabbatch{4}.spm.tools.USwLtools.uswl.options.imgTpm = {'C:\Dropbox\Work\1_SPM\FIL_spm12\tpm\unwTPM_sl2.nii'};
matlabbatch{4}.spm.tools.USwLtools.uswl.options.imgTpm = {fullfile(spm('dir'),'tpm','unwTPM_sl2.nii')};
matlabbatch{5}.spm.tools.USwLtools.MPMsmooth.wMPM(1) = cfg_dep('Unified segmentation with lesion mask: Warped MPM image #1', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','wMPM1'));
matlabbatch{5}.spm.tools.USwLtools.MPMsmooth.wMPM(2) = cfg_dep('Unified segmentation with lesion mask: Warped MPM image #2', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','wMPM2'));
matlabbatch{5}.spm.tools.USwLtools.MPMsmooth.wMPM(3) = cfg_dep('Unified segmentation with lesion mask: Warped MPM image #3', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','wMPM3'));
matlabbatch{5}.spm.tools.USwLtools.MPMsmooth.wMPM(4) = cfg_dep('Unified segmentation with lesion mask: Warped MPM image #4', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','wMPM4'));
matlabbatch{5}.spm.tools.USwLtools.MPMsmooth.wcImg(1) = cfg_dep('Unified segmentation with lesion mask: wc1 image', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','wc1'));
matlabbatch{5}.spm.tools.USwLtools.MPMsmooth.wcImg(2) = cfg_dep('Unified segmentation with lesion mask: wc2 image', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','wc2'));
matlabbatch{5}.spm.tools.USwLtools.MPMsmooth.wcImg(3) = cfg_dep('Unified segmentation with lesion mask: wc3 image', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','wc3'));
matlabbatch{5}.spm.tools.USwLtools.MPMsmooth.wcImg(4) = cfg_dep('Unified segmentation with lesion mask: wc4 image', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','wc4'));
matlabbatch{5}.spm.tools.USwLtools.MPMsmooth.tpm_l(1) = cfg_dep('Unified segmentation with lesion mask: Subject''s TPM with Lesion', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','TPMl'));
matlabbatch{5}.spm.tools.USwLtools.MPMsmooth.fwhm = 8;
matlabbatch{6}.spm.tools.USwLtools.ParEx.imgMPM(1) = cfg_dep('Named File Selector: MPM(1) - Files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
matlabbatch{6}.spm.tools.USwLtools.ParEx.imgMPM(2) = cfg_dep('Named File Selector: MPM(2) - Files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{2}));
matlabbatch{6}.spm.tools.USwLtools.ParEx.imgMPM(3) = cfg_dep('Named File Selector: MPM(3) - Files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{3}));
matlabbatch{6}.spm.tools.USwLtools.ParEx.imgMPM(4) = cfg_dep('Named File Selector: MPM(4) - Files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{4}));
matlabbatch{6}.spm.tools.USwLtools.ParEx.cImg(1) = cfg_dep('Unified segmentation with lesion mask: c1 image', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c1'));
matlabbatch{6}.spm.tools.USwLtools.ParEx.cImg(2) = cfg_dep('Unified segmentation with lesion mask: c2 image', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c2'));
matlabbatch{6}.spm.tools.USwLtools.ParEx.cImg(3) = cfg_dep('Unified segmentation with lesion mask: c3 image', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c3'));
matlabbatch{6}.spm.tools.USwLtools.ParEx.cImg(4) = cfg_dep('Unified segmentation with lesion mask: c4 image', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','segmImg', '.','c4'));
matlabbatch{6}.spm.tools.USwLtools.ParEx.imgMsk(1) = cfg_dep('Named File Selector: Mask(1) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
matlabbatch{6}.spm.tools.USwLtools.ParEx.outdir = {''};
matlabbatch{6}.spm.tools.USwLtools.ParEx.opt.thrICV = 0.1;
matlabbatch{6}.spm.tools.USwLtools.ParEx.opt.thrLesion = 0.666666666666667;
matlabbatch{6}.spm.tools.USwLtools.ParEx.opt.thrTC = 0.666666666666667;
