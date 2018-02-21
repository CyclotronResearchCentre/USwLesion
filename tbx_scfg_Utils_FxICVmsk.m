function FxMsk = tbx_scfg_Utils_FxMsk
% MATLABBATCH sub-configuration file.
% Applying tissue spcific smoothing in order to limit partial volume 
% effect. This is specifically useful for the quantitative MPM images.
%_______________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

% ---------------------------------------------------------------------
% fnMsk Filename of mask to be fixed
% ---------------------------------------------------------------------
fnMsk         = cfg_files;
fnMsk.tag     = 'fnMsk';
fnMsk.name    = 'Mask image to be fixed';
fnMsk.help    = {'Mask image, e.g. an intracranial volume (ICV) mask to be fixed.', ...
    'If multiple images are selected, they are processed one by one.'};
fnMsk.filter = 'image';
% fnMsk.ufilter = '.*';
fnMsk.num     = [1 Inf];

%--------------------------------------------------------------------------
% SzThr Size threshold
%--------------------------------------------------------------------------
SzThr         = cfg_entry;
SzThr.tag     = 'SzThr';
SzThr.name    = 'Size threshold of small blob (mm^3)';
SzThr.help    = {'Maximum size (mm^3) of small blobs that can be filled up.'};
SzThr.strtype = 'r';
SzThr.num     = [1 1];
SzThr.val     = {1000};

% ---------------------------------------------------------------------
% options Options
% ---------------------------------------------------------------------
options         = cfg_branch;
options.tag     = 'options';
options.name    = 'Options';
options.val     = {SzThr};
options.help    = {'Some processing options.'};

%% EXEC function
%----------------------------------------------------------------------
% FxMsk Fixing some mask image(s)
% ---------------------------------------------------------------------
FxMsk        = cfg_exbranch;
FxMsk.tag    = 'FxMsk';
FxMsk.name   = 'Fixing mask image';
FxMsk.val    = {fnMsk options};
FxMsk.help   = {'Fixing a mask image, e.g. an intracranial volume (ICV) mask, by', ...
    '1/ filling small holes, based on #voxel threshold (1000mm^3 by def.)', ...
    '2/ removing blobs outside brain, i.e. the biggest blob'};
FxMsk.prog   = @crc_run_FxMsk;
FxMsk.vout   = @vout_MPMsmoothn;

end

%% OUTPUT function
%_______________________________________________________________________
function dep = vout_MPMsmoothn(job) %#ok<*INUSD>

cdep(1) = cfg_dep;
cdep(1).sname = 'Fixed maxk image(s)';
% dep(1).src_output = substruct('{}',{1});
cdep(1).src_output = substruct('.','files');
cdep(1).tgt_spec   = cfg_findspec({{'filter','image'}});

% nMsk = size(job.fnMsk,1);
% for ii=1:nMsk
%     dep(ii) = cfg_dep; %#ok<*AGROW>
%     
% end
% img_c = 0;
% for ii=1:numel(job.wMPM)
%     for jj=1:numel(job.wcImg)
%         img_c = img_c+1;
%         cdep(img_c) = cfg_dep; %#ok<*AGROW>
%         cdep(img_c).sname  = sprintf('fin%d_c%d image',[ii jj]);
%         cdep(img_c).src_output =  ...
%             substruct('.','fn','{}',{ii},'()',{jj});
% %         cdep(img_c).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%         cdep(img_c).tgt_spec = cfg_findspec({{'filter','nifti'}});
%     end
% end

dep = cdep;
end
%_______________________________________________________________________
%% RUN function
%_______________________________________________________________________
function out = crc_run_FxMsk(job)

nMsk = size(job.fnMsk,1);
fn = cell(nMsk,1);
for ii=1:nMsk
    fn{ii} = crc_fix_ICV(job.fnMsk{ii},job.options);
end

% out.files = {cellstr(char(fn))};
out.files = cellstr(char(fn));

end

