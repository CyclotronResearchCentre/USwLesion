function [fn_out,res] = crc_ExtractParam_qMRIs(fn_img,opt)
% Extracting parameters from the segmented images.
%
% Which parameters should extracted?
% Here are a few things done so far:
%   1. tICV -> reference volume
%   2. match between mask and segmented lesion volume
%   3. Extract volumic information (see list here under)
%   4. MPMvalues in the follwoing tissue types: cortical GM, deep GM (no
%      cerebellum nor brainstem), cortical NA, lesion
%      Note: cortical = brain - (central, cerebellum, brainstem), as
%      defined from the brainparts mask image.
%   5. some summary statistics: min, max , mean, median, std, skewness,
%      kurtosis, p10, & p90 of the values extracted at 4., for each MPM and
%      tissue types.
%
% INPUT
% - fn_img : structure with necessary filenames
%       .fn_cImg   : tissue classes (GM, WM, Les, CSF) or (GM, WM, CSF)
%       .fn_qMRI    : quantitative MRIs (typically A, MT, R1, R2s)
%       .fn_seg8   : warping & segmentation parameters
%       .fn_BrPart : brain parts mask for value estraction
%       .fn_MskLes : original lesion mask [optional]
% - opt    : structure with some options
%       .outdir    : output directory
%       .subj_Id   : subject Id for the filename of saved parameters
%       .thrICV    : threshold for the creation of the tICV mask [def .5]
%       .thrLesion : lesion map threshold, for match with lesion mask
%       .thrTC     : tissue class threshold, for value collection
%       .prefix    : prefix added to the filename [def. empty]
%       .fn_load   : filename to existing result structure, useful to avoid
%                    recalculating good bits
%       .exe_ops   : flags to decide which operation to perform. A [1x5]
%                    vectors [def. true(1,5)]
%       .qValSScal : scaling of the qMRI values before calculated summary
%                    statistics, a [1x4] vector [def. ones(1,4)]
%
% OUTPUT
% - fn_out  : filename of output matlab file with result structure
% - res     : result structure (the one saved) with the fields
%      .Picv         : filename of final ICV mask used for tICV calculation
%      .match        : struct for the matching between lesion mask and
%                      segmented lesion
%           .DiceC      : Dice coefficient
%           .JaccardI   : Jaccard index
%           .N_tmsk     : #voxels in lesion mask
%           .N_c3       : #voxels in segmented lesion
%           .N_intersec : #voxels in both segmented and lesion mask
%           .N_union    : #voxels in segmented and/or lesion mask
%      .tissueVol    : structu with volumic information (in liter or %)
%           .volumes    : volumes of 3 or 4 tissue classes
%           .tiv        : total intra cranial volume
%           .gmv        : total GM volume
%           .wmv        : total WM volume (normal appearing and lesion)
%           .lesv       : total lesion volume
%           .BrParenchFrac : Brain Parenchymal Fraction over TIV
%           .GMFrac     : GM fraction over TIV
%           .WMFrac     : WM fraction over TIV
%           .LesFrac    : Lesion fraction over TIV
%           .LesFracWM  : Lesion fraction over WM volume
%      .values_qMRI  : values extracted from the qMRIs in tisue classes and
%                      specific brain parts
%           .val_labels : labels of tissue class and brain part, typically
%                         'GM cortical', 'GM central','WM cortical', and
%                         'Lesion'
%           .val_qMRI   : extracted values 1 cell per label. Each cell
%                         contains an array of size [#voxels x #qMRIs]
%           .fn_qMRI    : filename of qMRI files used
%      .sSstats_qMRI :
%           .sstats_label : labels of summary statistics values
%           .sstats_qMRI  : values of summary statistics 1 cell per value
%                          label. Each cell contains an array of size
%                          [#sstats x #qMRIs]
%           .val_labels   : labels of tissue class and brain part,
%                           typically 'GM cortical', 'GM central','WM
%                           cortical', and 'Lesion'
%           .fn_qMRI      : filename of qMRI files used
%
%
% TO DO:
% Check the values picked up from qMRI in each tissue class & brain part???
%
% REFS:
% tICV : http://www.sciencedirect.com/science/article/pii/S1053811914007769
%__________________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by C. Phillips, 2015.
% Cyclotron Research Centre, University of Liege, Belgium

%% 0. deal with input and options
fn_cImg = fn_img.fn_cImg;
ncImg = size(fn_cImg,1); % 4 -> [GM, WM, Lesion, CSF] or 3 -> [GM, WM, CSF]
fn_qMRI  = fn_img.fn_qMRI;
nqMRI = size(fn_qMRI,1);
fn_seg8 = fn_img.fn_seg8;
fn_iwarp = fn_img.fn_iwarp;

if ~isfield(opt,'qValSScal') || isempty(opt.qValSScal) ...
        || numel(opt.qValSScal)~=4
    opt.qValSScal = ones(1,4); % no scaling by default
end

if ~isfield(opt,'exe_ops') || isempty(opt.exe_ops) || numel(opt.exe_ops)~=5
    opt.exe_ops = true(1,5); % run all by default
end

if isfield(fn_img,'fn_MskLes') && ~isempty(fn_img.fn_MskLes)
    fn_MskLes  = fn_img.fn_MskLes;
else
    fn_MskLes  = '';
    opt.exe_ops(2) = false; % Cannot calculate, whatever is requested
end

[pth,fn] = spm_fileparts(fn_qMRI(1,:));
if ~isfield(opt,'prefix')|| isempty(opt.outdir)
    pth_out = pth;
else
    pth_out = opt.outdir;
end

if ~isfield(opt,'prefix') || isempty(opt.prefix)
    fn_prefix = '';
else
    fn_prefix = [opt.prefix,'_'];
end

if ~isfield(opt,'subj_Id') || isempty(opt.subj_Id)
    opt.subj_Id = fn;
end

if ~isfield(opt,'fn_load') || isempty(opt.fn_load)
    loadR_fl = false;
else
    try
        loadR = load(opt.fn_load);
        loadR_fl = true;
    catch
        loadR_fl = false;
        fprintf('\nTried to load file %s but failed.\n', ...
            spm_file(opt.fn_load,'filename'));
    end
end
if loadR_fl
    res = loadR.res;
end


%% 1. build ICV from the sum of GM, WM, CSF and lesion
if opt.exe_ops(1)
    fn_icv = spm_file(spm_file(fn_cImg(1,:),'filename'),'prefix','tICV_');
    matlabbatch = create_mask(fn_cImg(1:end,:),fn_icv,pth,opt.thrICV);
    spm_jobman('run', matlabbatch);
    res.Picv = fn_icv;
end

%% 2. match between lesion mask and segmented tissue
if opt.exe_ops(2)
    % Dice coefficient: http://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient
    % Jaccard index: http://en.wikipedia.org/wiki/Jaccard_index
    
    % In native space
    Vtmsk = spm_vol(fn_MskLes);
    Vc3 = spm_vol(fn_cImg(3,:));
    v_tmsk = spm_read_vols(Vtmsk);
    v_tmsk = v_tmsk(:) >.5; % thresholding at .5, as it should be a binary img
    v_c3 = spm_read_vols(Vc3);
    v_c3 = v_c3(:) > opt.thrLesion;
    
    N_tmsk = sum(v_tmsk);
    N_c3 = sum(v_c3);
    N_intersec = sum(v_tmsk.*v_c3);
    N_union = sum( (v_tmsk+v_c3)>0 );
    
    % Dice coefficient
    DC = 2*N_intersec / (N_tmsk + N_c3) ;
    % Jaccard index
    JI = N_intersec/N_union ;
    
    res.match = struct( ...
        'DiceC', DC, ...
        'JaccardI', JI, ...
        'N_tmsk', N_tmsk, ...
        'N_c3', N_c3, ...
        'N_intersec', N_intersec, ...
        'N_union', N_union);
end

%% 3. Volumes
% Volumes or fractions to extract:
% - TIV
% - GMV = Cortical and Deep Grey matter
% - WMV = NAWM + lesions
% - Brain Parenchymal Fraction = (GMV+WMV)/TIV
% - Grey Matter Fraction = GMV/TIV
% - White Matter Fraction =  WMV/TIV
% - Lesion Fraction = Lesion volume/TIV (more standard than relative than
%   fraction of WMV)
% - Lesion Fraction WM = Lesion volume/WMV

if opt.exe_ops(3)
    clear matlabbatch
    % fn_seg8 = spm_select('FPList',pth,'^kt.*_seg8\.mat$');
    matlabbatch{1}.spm.util.tvol.matfiles = {fn_seg8};
    matlabbatch{1}.spm.util.tvol.tmax = ncImg;
    matlabbatch{1}.spm.util.tvol.mask = {fullfile(spm('dir'),'tpm','mask_ICV.nii,1')};
    % matlabbatch{1}.spm.util.tvol.outf = fullfile(pth_out,['TV_',fn]);
    spm_jobman('run', matlabbatch);
    tmp = load(fn_seg8);
    volumes = tmp.volumes;
    tiv = sum(volumes.litres);
    gmv = volumes.litres(1);
    if ncImg ==4
        wmv = sum(volumes.litres([2 3])); % WM + lesion volume!
        lesv = volumes.litres(3);
    else
        wmv = volumes.litres(2);
    end
    
    BrParenchFrac = (gmv+wmv)/tiv*100; % expressed in percentage
    GMFrac = gmv/tiv*100; % expressed in percentage
    WMFrac = wmv/tiv*100; % expressed in percentage
    if ncImg ==4
        LesFrac = volumes.litres(3)/tiv*100; % expressed in percentage
        LesFracWM = volumes.litres(3)/wmv*100; % expressed in percentage
    end
    
    tissueVol = struct(...
        'volumes', volumes.litres, ...
        'tiv', tiv, ...
        'gmv', gmv, ...
        'wmv', wmv, ...
        'BrParenchFrac', BrParenchFrac, ...
        'GMFrac', GMFrac, ...
        'WMFrac', WMFrac);
    if ncImg ==4
        tissueVol.lesv = lesv;
        tissueVol.LesFrac = LesFrac;
        tissueVol.LesFracWM = LesFracWM;
    end
    
    % Store in main results structure
    res.tissueVol = tissueVol;
end

%% 4. extraction of MPM values for the GM/WM(/lesion)
% Accounting for the brain parts selected.

if opt.exe_ops(4)
    
    % Some definitions & pre-loading
    val_qMRI = cell(1,ncImg);
    val_labels = char('GMcortical', 'GMcentral','WMcortical');
    if ncImg==4
        val_labels = char(val_labels,'Lesion');
    end
    VqMRI = spm_vol(fn_qMRI);
    tval_c123 = zeros(prod(VqMRI(1).dim),ncImg);
    
    % Brain Parts mask into subject subject space
    %--------------------------------------------
    % copy into temporary subdir
    fn_brP = fullfile(spm_file(mfilename('fullpath'),'path'), ...
        'eTPM','msk_BrainParts.nii');
    dr_tmp = fullfile(pth,'tmp');
    if ~exist(dr_tmp,'dir'), mkdir(dr_tmp); end
    copyfile(fn_brP,dr_tmp);
    nBrP = 6; % Use all six images of 4D volume
    fn_brP_loc = '';
    for ii=1:nBrP
        fn_brP_loc = char(fn_brP_loc, ...
            spm_file(fn_brP,'path',dr_tmp,'number',ii));
    end
    fn_brP_loc(1,:) = [];
    
    % Create MatlabBatch & bring BrainParts into subject space
    clear matlabbatch
    matlabbatch = create_MB(fn_iwarp,fn_cImg(1,:),fn_brP_loc);
    spm_jobman('run', matlabbatch);
    
    fn_wbrP_loc = spm_file(fn_brP_loc,'prefix','w');
    Vmsk_cortex = spm_vol(fn_wbrP_loc(2,:));
    Vmsk_central = spm_vol(fn_wbrP_loc(6,:));
    
    % Deal with
    %   * GM: cortical (BP#2) and central (BP#6) only
    %   * WM: cortical (BP#2)
    %   * Lesion: anywhere (if provided)
    %------------------------------------------------------
    % GM, c1
    Vc1 = spm_vol(fn_cImg(1,:));
    [val_c1,XYZmm] = spm_read_vols(Vc1);
    tval_c1 = val_c1(:)>=opt.thrTC; % Keep voxel with p>= thr
    XYZvx_msk = Vmsk_cortex.mat\[XYZmm ; ones(1,numel(tval_c1))];
    % cortical
    val_BPcortex = spm_sample_vol(Vmsk_cortex,XYZvx_msk(1,:)',XYZvx_msk(2,:)',XYZvx_msk(3,:)',1);
    tval_c123(:,1) = tval_c1 .* val_BPcortex>=.5;
    % central
    val_BPcentral = spm_sample_vol(Vmsk_central,XYZvx_msk(1,:)',XYZvx_msk(2,:)',XYZvx_msk(3,:)',1);
    tval_c123(:,2) = tval_c1 .* val_BPcentral>=.5;
    
    % WM, c2
    Vc2 = spm_vol(fn_cImg(2,:));
    val_c2 = spm_read_vols(Vc2);
    tval_c2 = val_c2(:)>=opt.thrTC;
    tval_c123(:,3) = tval_c2 .* val_BPcortex>=.5;
    
    % Lesion, c3
    if ncImg==4
        Vc3 = spm_vol(fn_cImg(3,:));
        val_c3 = spm_read_vols(Vc3);
        tval_c123(:,4) = val_c3(:)>=opt.thrTC;
    end
    tval_c123 = logical(round(tval_c123));
    
    % Initialize values to NaN -> check after if some values not defined!
    for i_tcBP = 1:ncImg
        val_qMRI{i_tcBP} = zeros(sum(tval_c123(:,i_tcBP)),nqMRI) + NaN;
    end
    
    % Load qMRI values, loop over qMRIs and tissue class per brain part
    for j_qMRI = 1:nqMRI
        val_qMRI_i = spm_read_vols(VqMRI(j_qMRI));
        val_qMRI_i = val_qMRI_i(:);
        for i_tcBP = 1:ncImg
            val_qMRI{i_tcBP}(:,j_qMRI) = val_qMRI_i(tval_c123(:,i_tcBP));
        end
    end
    
    % Stacking things into result strucutre
    values_qMRI = struct(...
        'val_labels', {val_labels}, ...
        'val_qMRI', {val_qMRI}, ...
        'fn_qMRI', fn_qMRI);
    res.values_qMRI = values_qMRI;
end

%% 5. some stats from the MPM values
% Find the following summary statistic for intensities in each qMRI in
% each tissue class per brain part:
% mean, median, std, P10, P90, min, max, skewness ,kurtosis

if opt.exe_ops(5)
    sstats_qMRI = cell(1,4); % Make it 4, even if no lesion (-> NaN)
    val_labels = char('GMcortical', 'GMcentral','WMcortical','Lesion');
    sstats_label = char('mean', 'median', 'std', 'P10', 'P90', 'P9010', ...
        'min', 'max', 'skewness' ,'kurtosis');
    nsstats = size(sstats_label,1);
    
    for i_tcBP = 1:ncImg
        sstats_qMRI{i_tcBP} = zeros(nsstats,nqMRI)+NaN;
        for j_qMRI = 1:nqMRI
            tmp_val = res.values_qMRI.val_qMRI{i_tcBP}(:,j_qMRI);
            l_neg = find(tmp_val<0);
            if ~isempty(l_neg)
                tmp_val(l_neg) = []; % remove negative-values
                fprintf('\n Subject %s, tissue %s, #zeros = %d', ...
                    opt.subj_Id, ...
                    deblank(res.values_qMRI.val_labels(i_tcBP,:)), ...
                    numel(l_neg) );
            end
            l_zer = find(tmp_val==0);
            if ~isempty(l_zer)
                tmp_val(l_zer) = []; % remove zero-values
                fprintf('\n Subject %s, tissue %s, #zeros = %d', ...
                    opt.subj_Id, ...
                    deblank(res.values_qMRI.val_labels(i_tcBP,:)), ...
                    numel(l_zer) );
            end
            if any([l_neg ; l_zer]), fprintf('\n'), end
            % Scaling
            tmp_val = tmp_val*opt.qValSScal(j_qMRI);
            % Get all summary statistics
            sstats_qMRI{i_tcBP}(:,j_qMRI) = ...
                get_sstats(tmp_val,sstats_label);
        end
    end
    if ncImg==3 % Add NaNs for the Lesion bits
        sstats_qMRI{4} = NaN(size(sstats_qMRI{3}));
    end
    
    sSstats_qMRI = struct( ...
        'sstats_qMRI', {sstats_qMRI}, ...
        'sstats_label', {sstats_label}, ...
        'val_labels', {val_labels}, ...
        'fn_qMRI', fn_qMRI);
    
    res.sSstats_qMRI = sSstats_qMRI;
end

%% 6. save things and pass out fn_out
fn_out = fullfile(pth_out,['ExP_',opt.subj_Id]);
fn_out = spm_file(fn_out,'prefix',fn_prefix);
save(fn_out,'res');


end

% =======================================================================
%% SUBFUNCTION
% =======================================================================
function matlabbatch = create_mask(fn_in,fn_out,pth,thr,smK)
% Batch to create a mask from a few segmented tissue classes.
% It works in 3 steps:
% 1. sum up the input images
% 2. smooth the sum
% 3. threshold the result

if nargin<5, smK = 8; end
if nargin<4, thr = .5; end

matlabbatch{1}.spm.util.imcalc.input = cellstr(fn_in);
matlabbatch{1}.spm.util.imcalc.output = 'tmp.nii';
matlabbatch{1}.spm.util.imcalc.outdir = {pth};
matlabbatch{1}.spm.util.imcalc.expression = 'sum(X)';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 2;
matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('Image Calculator: ImCalc Computed Image: output', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{2}.spm.spatial.smooth.fwhm = [1 1 1]*smK;
matlabbatch{2}.spm.spatial.smooth.dtype = 0;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's';
matlabbatch{3}.spm.util.imcalc.input(1) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{3}.spm.util.imcalc.output = fn_out;
matlabbatch{3}.spm.util.imcalc.outdir = {pth};
matlabbatch{3}.spm.util.imcalc.expression = sprintf('i1>%f',thr);
matlabbatch{3}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{3}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{3}.spm.util.imcalc.options.mask = 0;
matlabbatch{3}.spm.util.imcalc.options.interp = 1;
matlabbatch{3}.spm.util.imcalc.options.dtype = 2;
matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('Image Calculator: ImCalc Computed Image: output', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_move.files(2) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{4}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;

end

% =======================================================================
function matlabbatch = create_MB(fn_iwarp,fn_subj,fn_brP)
% Building matlabbatch for the normalize-write operation, muche easier than
% building the deformation and applying it manually.
% Bring SPM-ICV into subject space -> need to properly define the latter!

prec_round = 1e6;

V_subj = spm_vol(fn_subj);
% voxel size rounded to some precision
vx_sz = round(sqrt(sum(V_subj.mat(1:3,1:3).^2))*prec_round)/prec_round;
% defining BB
p_min = -V_subj.mat\[0 0 0 1]' ; p_min = p_min(1:3);
p_max = (V_subj.dim.*vx_sz)' + p_min ;
img_bb = round([-abs(p_min') ; abs(p_max')]);

matlabbatch{1}.spm.spatial.normalise.write.subj.def(1) = {spm_file(fn_iwarp,'number','')};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(fn_brP);
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = img_bb;
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = vx_sz;
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;

end

% =======================================================================
function sstats_v = get_sstats(val,labels)
% extracting the summary statistics from a bunch of voxel values, based on
% the labels of the values requested.
% Dealing with mean, median, std, P10, P90, min, max, skewness ,kurtosis.

nlabels = size(labels,1);
sstats_v = zeros(nlabels,1)+NaN;
for ii=1:nlabels
    switch lower(deblank(labels(ii,:)))
        case 'mean'
            sstats_v(ii) = mean(val);
        case 'median'
            sstats_v(ii) = median(val);
        case 'std'
            sstats_v(ii) = std(val);
        case 'p10'
            sstats_v(ii) = prctile(val,10);
        case 'p90'
            sstats_v(ii) = prctile(val,90);
        case 'p9010'
            sstats_v(ii) = prctile(val,90) - prctile(val,10);
        case 'min'
            sstats_v(ii) = min(val);
        case 'max'
            sstats_v(ii) = max(val);
        case 'skewness'
            sstats_v(ii) = skewness(val);
        case 'kurtosis'
            sstats_v(ii) = kurtosis(val);
        otherwise
            fprintf('\nCould not find operation : ''%s''',deblank(labels(ii,:)))
    end
end

end