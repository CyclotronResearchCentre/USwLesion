function FxMPM = tbx_scfg_Utils_FxMPM
% MATLABBATCH sub-configuration file.
% Small tool to fix a lesion mask
%_______________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

% Get some defaults
v_thrMPM = crc_USwL_get_defaults('ImgFix.thrMPM');
n_thrMPM = crc_USwL_get_defaults('ImgFix.strMPM');
if numel(v_thrMPM)~=numel(n_thrMPM)
    error('USwL:fixMPM','Mismatch in qMRI/MPM fixing defaults.');
end

% ---------------------------------------------------------------------
% fnMPM Filename of MPM/qMR images to be fixed
% ---------------------------------------------------------------------
fnMPM         = cfg_files;
fnMPM.tag     = 'fnMPM';
fnMPM.name    = 'MPM/qMR images to be fixed';
fnMPM.help    = {'Select the quantitative MRIs (aka. Multi-Parametric Maps, MPMs) to be fixed.', ...
    'If multiple images are selected, they are processed one by one.'};
fnMPM.filter = 'image';
fnMPM.num     = [1 Inf];

%--------------------------------------------------------------------------
% strMPM  filename suffix used to pick image types (_A, _MT, _R1 & _R2s)
%--------------------------------------------------------------------------
strMPM         = cfg_entry;
strMPM.tag     = 'strMPM';
strMPM.name    = 'Filename suffix for the qMRI/MPM';
strMPM.help    = {'Enter the filename suffix for the corresponding qMRI/MPM.'};
strMPM.strtype = 's';
strMPM.num     = [1 Inf];
strMPM.val     = {char(n_thrMPM)};

%--------------------------------------------------------------------------
% thrMPM  thresholds for the correspondign qMRIs (A, MT, R1 & R2s)
%--------------------------------------------------------------------------
thrMPM         = cfg_entry;
thrMPM.tag     = 'thrMPM';
thrMPM.name    = 'Threshold values for the qMRI/MPM';
thrMPM.help    = {'Enter the maximum value for the corresponding qMRI/MPM.'};
thrMPM.strtype = 'r';
thrMPM.num     = [1 Inf];
thrMPM.val     = {v_thrMPM'};

% ---------------------------------------------------------------------
% fxZeros Fix zeros in the images, or not
% ---------------------------------------------------------------------
fxZeros        = cfg_menu;
fxZeros.tag    = 'fxZeros';
fxZeros.name   = 'Fix zeros in the images';
fxZeros.help   = {
    'Choose wheter zero vlaues in the qMRI/MPMs should be fixed by interpolating from its neighbours.'
    }';
fxZeros.labels = {
    'Yes'
    'No'}';
fxZeros.values = {1 2};
fxZeros.val    = {1};

%--------------------------------------------------------------------------
% fn_prefix  filename suffix used to pick image types (_A, _MT, _R1 & _R2s)
%--------------------------------------------------------------------------
fn_prefix         = cfg_entry;
fn_prefix.tag     = 'fn_prefix';
fn_prefix.name    = 'Filename prefix for the mask';
fn_prefix.help    = {'Enter the prefix for the mask of fixed voxels.'};
fn_prefix.strtype = 's';
fn_prefix.num     = [1 Inf];
fn_prefix.val     = {'fxM_'};

% ---------------------------------------------------------------------
% crFxMsk_yes Create mask of fixed voxels.
% ---------------------------------------------------------------------
crFxMsk_yes         = cfg_branch;
crFxMsk_yes.tag     = 'crFxMsk_yes';
crFxMsk_yes.name    = 'Yes';
crFxMsk_yes.val     = {fn_prefix};
crFxMsk_yes.help    = {'Create mask of fixed voxels.'};

% ---------------------------------------------------------------------
% crFxMsk_no No mask of fixed voxels created.
% ---------------------------------------------------------------------
crFxMsk_no         = cfg_const;
crFxMsk_no.tag     = 'crFxMsk_no';
crFxMsk_no.name    = 'No';
crFxMsk_no.val     = {0};
crFxMsk_no.help    = {'No mask of fixed voxels created.'};

% ---------------------------------------------------------------------
% crFxMsk Should a mask showing the fixed voxels be creates?
% ---------------------------------------------------------------------
crFxMsk        = cfg_choice;
crFxMsk.tag    = 'crFxMsk';
crFxMsk.name   = 'Should a mask showing the fixed voxels be created?';
crFxMsk.values = {crFxMsk_no, crFxMsk_yes};
crFxMsk.val    = {crFxMsk_no};
crFxMsk.help   = {...
    ['You may wish to keep track of which voxels were fixed in the qMRI/MPMs ',...
    ', if so select ''yes'' and pick a prefix for the created image.']};

% ---------------------------------------------------------------------
% options Options
% ---------------------------------------------------------------------
options         = cfg_branch;
options.tag     = 'options';
options.name    = 'Options';
options.val     = {strMPM thrMPM fxZeros crFxMsk};
options.help    = {'Some processing options.'};

%% EXEC function
%----------------------------------------------------------------------
% FxLesMsk Fixing some (lesion) mask image(s)
% ---------------------------------------------------------------------
FxMPM        = cfg_exbranch;
FxMPM.tag    = 'FxMPM';
FxMPM.name   = 'Fixing qMRI/MPMs';
FxMPM.val    = {fnMPM options};
FxMPM.help   = {'Fixing qMRI/MPM', ...
    '1/ capping the values [0 max], depending on image type (based on their filename suffix)', ...
    '2/ if need, fill in the "zeros holes" in the image by interpolating from the neighbours', ...
    'A mask indicating which voxels were fixed can be output too.'};
FxMPM.prog   = @run_FxMPM;
FxMPM.vout   = @vout_FxMPM;
FxMPM.check   = @check_suffix_thresh;

end

%_______________________________________________________________________
%% OUTPUT function
%_______________________________________________________________________
function dep = vout_FxMPM(job) %#ok<*INUSD>

dep(1) = cfg_dep;
dep(1).sname = 'Fixed qMRI/MPMs';
dep(1).src_output = substruct('.','files');
dep(1).tgt_spec   = cfg_findspec({{'filter','image'}});

end
%_______________________________________________________________________
%% RUN function
%_______________________________________________________________________
function out = run_FxMPM(job)
% Need to reorganize the options mainly

% fn_out = crc_fix_MPMintens(fn_in,opt)
%
% INPUT
% - fn_in   : filename of MPM images to fix
% - opt     : option structure
%       strMPM   : filename parts used to pick the image type.
%       thrMPM   : Max cap for the corresponding image type
%       prefix   : prefix added to filename.
%       crt_mask : create a map indicating the voxels that were fixed
%                  [def. false]
%       fix_zeros: fix the zero-holes in the image.
%
% OUTPUT
% - fn_out  : filename of fixed MPM images
opt.strMPM = cellstr(job.options.strMPM);
opt.thrMPM = job.options.thrMPM;
opt.fix_zeros = job.options.fxZeros;

if isfield(job.options.crFxMsk,'crFxMsk_yes')
    opt.crt_mask = true;
    opt.prefix = job.options.crFxMsk.crFxMsk_yes.fn_prefix;
else
    opt.crt_mask = false;
end
fn_in = spm_file(char(job.fnMPM),'number','');
out.files = cellstr(crc_fix_MPMintens(fn_in,opt));

end

%_______________________________________________________________________
%% CHECK function
%_______________________________________________________________________
function t = check_suffix_thresh(job)
% Checking that the data are consistent.
t   = {};

nStr = size(job.options.strMPM,1);
nThr = numel(job.options.thrMPM);

if nStr~=nThr
    t{1} = 'Number of suffixes should match that of thresholds!';
    warndlg(t,'Mismatch in qMRI/MPM thesholds.');
    return
end

end
