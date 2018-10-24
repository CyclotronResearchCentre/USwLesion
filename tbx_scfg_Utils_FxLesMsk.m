function FxLesMsk = tbx_scfg_Utils_FxLesMsk
% MATLABBATCH sub-configuration file.
% Small tool to fix a lesion mask
%_______________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

% ---------------------------------------------------------------------
% fnMsk Filename of mask to be fixed
% ---------------------------------------------------------------------
fnMsk         = cfg_files;
fnMsk.tag     = 'fnMsk';
fnMsk.name    = 'Lesion mask image to be fixed';
fnMsk.help    = {'Select the lesion mask image(s) to be fixed.', ...
    'If multiple images are selected, they are processed one by one.'};
fnMsk.filter = 'image';
fnMsk.num     = [1 Inf];

% ---------------------------------------------------------------------
% fnOth Filename of other image to be used for its intensities
% ---------------------------------------------------------------------
fnOth         = cfg_files;
fnOth.tag     = 'fnOth';
fnOth.name    = 'Other image to use';
fnOth.help    = {'Select  the other image(s) to use for their intensities.', ...
    'If none are selected, then this step will be skipped in the mask fixing', ...
    'If you select some image(s), you should select the same number of ''Mask'' and ''Other'' images!'};
fnOth.filter = 'image';
fnOth.num     = [0 Inf];
fnOth.val       = {''};

%--------------------------------------------------------------------------
% minVol Minimum volume threshold
%--------------------------------------------------------------------------
minVol         = cfg_entry;
minVol.tag     = 'minVol';
minVol.name    = 'Size threshold (mm^3)';
minVol.help    = {'Minimum size (mm^3) of smallest blobs that are kep.'};
minVol.strtype = 'r';
minVol.num     = [1 1];
minVol.def     = @(val)crc_USwL_get_defaults('ImgFix.minVol', val{:});
% minVol.val     = {8};

% ---------------------------------------------------------------------
% options Options
% ---------------------------------------------------------------------
options         = cfg_branch;
options.tag     = 'options';
options.name    = 'Options';
options.val     = {minVol fnOth};
options.help    = {'Some processing options.'};

%% EXEC function
%----------------------------------------------------------------------
% FxLesMsk Fixing some (lesion) mask image(s)
% ---------------------------------------------------------------------
FxLesMsk        = cfg_exbranch;
FxLesMsk.tag    = 'FxLesMsk';
FxLesMsk.name   = 'Fixing lesion mask image';
FxLesMsk.val    = {fnMsk options};
FxLesMsk.help   = {'Fixing lesion mask image(s) by', ...
    '1/ removing small blobs below some minimum volume', ...
    '2/ if available, removing blobs whos intenisty is too low compared to the lesion mean.', ...
    'For details on 2/ see the main function ''crc_fix_LesMsk.m''. Note also that this is qui experimental'};
FxLesMsk.prog   = @run_FxLesMsk;
FxLesMsk.vout   = @vout_FxLesMsk;
FxLesMsk.check   = @check_file_number;

end

%_______________________________________________________________________
%% OUTPUT function
%_______________________________________________________________________
function dep = vout_FxLesMsk(job) %#ok<*INUSD>

dep(1) = cfg_dep;
dep(1).sname = 'Fixed mask image(s)';
dep(1).src_output = substruct('.','files');
dep(1).tgt_spec   = cfg_findspec({{'filter','image'}});

end
%_______________________________________________________________________
%% RUN function
%_______________________________________________________________________
function out = run_FxLesMsk(job)
% Work out case for multiple masks
nMsk = size(job.fnMsk,1);
fn = cell(nMsk,1);
% Get options
opt_ii.minVol = job.options.minVol;
if ~isempty(job.options.fnOth)
    use_other = true;
else
    use_other = false;
end
for ii=1:nMsk
    if use_other
        opt_ii.fn_oth = job.options.fnOth{ii};
    end
    fn{ii} = crc_fix_LesMsk(job.fnMsk{ii},opt_ii);
end

out.files = cellstr(char(fn));

end

%_______________________________________________________________________
%% CHECK function
%_______________________________________________________________________
function t = check_file_number(job)
% Checking that the data are consistent.
t   = {};

nMsk = numel(job.fnMsk);
nOth = numel(job.options.fnOth);
if nOth>0
    if nMsk~=nOth
        t{1} = 'Number of ''Mask'' images not matching that of ''Other''!';
        warndlg(t,'Other images number');
        return
    end
end

end
