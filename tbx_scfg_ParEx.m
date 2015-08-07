function ParEx = tbx_scfg_ParEx
% MATLABBATCH sub-configuration file.
% Extracting parameters from the segmented images
%__________________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% ---------------------------------------------------------------------
% imgMPM Structural quantitative images
% ---------------------------------------------------------------------
imgMPM         = cfg_files;
imgMPM.tag     = 'imgMPM';
imgMPM.name    = 'Structural quantitative images';
imgMPM.help    = {'Select the structural quantitative (MPM-VBQ) images.'};
imgMPM.filter  = 'image';
imgMPM.ufilter = '.*';
imgMPM.num     = [1 Inf];

% ---------------------------------------------------------------------
% imgMPMmsk Mask of excluded voxels from structural quantitative images
% ---------------------------------------------------------------------
imgMPMmsk         = cfg_files;
imgMPMmsk.tag     = 'imgMPMmsk';
imgMPMmsk.name    = 'Mask of excluded voxels from MPMs';
imgMPMmsk.help    = {'Select the "msk_MPM" images.'};
imgMPMmsk.filter  = 'image';
imgMPMmsk.ufilter = '.*';
imgMPMmsk.num     = [0 Inf];
imgMPMmsk.val       = {''};

% ---------------------------------------------------------------------
% cImg Tissue class images
% ---------------------------------------------------------------------
cImg         = cfg_files;
cImg.tag     = 'cImg';
cImg.name    = 'Tissue class images';
cImg.help    = {'Select the tissue classes of GM/WM/CSF/Lesion'};
cImg.filter  = 'image';
cImg.ufilter = '^c.*';
cImg.num     = [3 4];

% ---------------------------------------------------------------------
% imgMsk Mask image
% ---------------------------------------------------------------------
imgMsk         = cfg_files;
imgMsk.tag     = 'imgMsk';
imgMsk.name    = 'Mask image';
imgMsk.help    = {'Select the "lesiosn mask" image.'};
imgMsk.filter = 'image';
imgMsk.ufilter = '.*';
imgMsk.num     = [0 1];

%--------------------------------------------------------------------------
% outdir Output Directory
%--------------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output Directory';
outdir.val{1}  = {''};
outdir.help    = {'File produced will be written into this output directory. If no directory is given, file will be written to directory of MPM images.'};
outdir.filter  = 'dir';
outdir.ufilter = '.*';
outdir.num     = [0 1];

%--------------------------------------------------------------------------
% thrICV Threshold for ICV definition
%--------------------------------------------------------------------------
thrICV         = cfg_entry;
thrICV.tag     = 'thrICV';
thrICV.name    = 'Threshold for ICV definition';
thrICV.help    = {'Threshold for ICV definition'};
thrICV.strtype = 'r';
thrICV.num     = [1 1];
thrICV.val     = {.5};

%--------------------------------------------------------------------------
% thrLesion Threshold for Lesion volume
%--------------------------------------------------------------------------
thrLesion         = cfg_entry;
thrLesion.tag     = 'thrLesion';
thrLesion.name    = 'Threshold for Lesion volume';
thrLesion.help    = {'Threshold to estimate the lesion volume'};
thrLesion.strtype = 'r';
thrLesion.num     = [1 1];
thrLesion.val     = {.8};

%--------------------------------------------------------------------------
% thrTC Threshold for tissue classes
%--------------------------------------------------------------------------
thrTC         = cfg_entry;
thrTC.tag     = 'thrTC';
thrTC.name    = 'Threshold for tissue classes';
thrTC.help    = {'Threshold for tissue classes, when extracting the MPM values'};
thrTC.strtype = 'r';
thrTC.num     = [1 1];
thrTC.val     = {.8};

% ---------------------------------------------------------------------
% opt Options
% ---------------------------------------------------------------------
opt         = cfg_branch;
opt.tag     = 'opt';
opt.name    = 'Options';
opt.val     = {thrICV thrLesion thrTC};
opt.help    = {'Defining some thresholds for the parameters/values extraction'};
%_______________________________________________________________________


%% EXEC function
% ---------------------------------------------------------------------
% ParEx Unified segmentation with lesion mas
% ---------------------------------------------------------------------
ParEx         = cfg_exbranch;
ParEx.tag     = 'ParEx';
ParEx.name    = 'Parameter extraction for the GM/WM/lesion';
ParEx.val     = {imgMPM imgMPMmsk cImg imgMsk outdir opt}; 
ParEx.help    = {['Extracting some parameters from the MPMs over the ',...
    'GM/WM/lesion tissue classes.'],...
    'See the processing function itself for the details of what''s computed.'};
ParEx.prog    = @crc_ExtractParam;
ParEx.vout    = @vout_ExtractParam;
ParEx.check   = @imgMPMmsk_check;

end

%% OUTPUT function
%_______________________________________________________________________
function dep = vout_ExtractParam(job) %#ok<*INUSD>
dep            = cfg_dep;
dep.sname      = 'Extracted Parameters Mat_file';
dep.src_output = substruct('.','fn_ExParam');
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});

end
%_______________________________________________________________________

%% CHECK function
%_______________________________________________________________________
function t = imgMPMmsk_check(job)
t = {};
if ~isempty(job.imgMPM) && ~isempty(job.imgMPMmsk)
    % Check that numbers do match.
    Nimg = numel(job.imgMPM);
    Nmsk = numel(job.imgMPMmsk);
    if Nimg~=Nmsk
        t = {sprintf('Num MPM images (%d) ~= Num Msk images (%d)', Nimg,Nmsk)};
    end
end
% for i=1:numel(sess.regress)
%     if numel(sess.scans) ~= numel(sess.regress(i).val)
%         t = {t{:}, sprintf('Num scans (%d) ~= Num regress[%d] (%d).',numel(sess.scans),i,numel(sess.regress(i).val))};
%     end
% end

end