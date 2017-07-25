%-----------------------------------------------------------------------
% Job saved on 08-Jun-2016 17:49:36 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (12.2)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.util.imcalc.input = {'C:\Dropbox\Work\1_SPM\FIL_spm12\tpm\labels_Neuromorphometrics.nii,1'};
matlabbatch{1}.spm.util.imcalc.output = 'b_pallidum.nii';
matlabbatch{1}.spm.util.imcalc.outdir = {'C:\Dropbox\Work\1_SPM\FIL_spm12\toolbox\USwLesion\Script_and_batches'};
matlabbatch{1}.spm.util.imcalc.expression = '(i1==55) + (i1==56)';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 2;
matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('Image Calculator: ImCalc Computed Image: b_pallidum.nii', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{2}.spm.spatial.smooth.fwhm = [4 4 4];
matlabbatch{2}.spm.spatial.smooth.dtype = 16;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's';
