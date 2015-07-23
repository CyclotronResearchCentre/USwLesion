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
ParEx.val     = {imgMPM cImg imgMsk outdir opt}; % CP: need to add the input
ParEx.help    = {'Extracting some parameters from the MPMs over the GM/WM/lesion tissue classes'};
ParEx.prog = @crc_ExtractParam;
ParEx.vout = @vout_ExtractParam;

%% OUTPUT function
%_______________________________________________________________________
function dep = vout_ExtractParam(job) %#ok<*INUSD>
dep            = cfg_dep;
dep.sname      = 'Extracted Parameters Mat_file';
dep.src_output = substruct('.','fn_ExParam');
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});
%_______________________________________________________________________
