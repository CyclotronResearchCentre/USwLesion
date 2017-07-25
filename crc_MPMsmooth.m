function fn_out = crc_MPMsmooth(job)
% Applying tissue specific smoothing in order to limit partial volume
% effect. This is specifically useful for the quantitative MPM images.
%_______________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

fn_wMPM = char(job.wMPM);
fn_mwTC = spm_file(char(job.wcImg),'number','');

% Find list of tissue classes, tc_ind
tc_ind = [];
for ii = 1:size(fn_mwTC,1)
    p = strfind(deblank(fn_mwTC(ii,:)),'wc');
    tc_ind = [tc_ind str2double(deblank(fn_mwTC(ii,p+2)))]; %#ok<*AGROW>
end

% Collect the TPM_l
tpm = spm_file(char(job.tpm_l),'number','');
fn_lTPM = char;
for ii = tc_ind
    fn_lTPM = char( fn_lTPM , [tpm,',',num2str(ii)]);
end
fn_lTPM(1,:) = [];

opt_process = struct( ...
    'fwhm', job.fwhm, ...
    'tpm', fn_lTPM) ; % , ...

% Do it!
% fn_finMPM = crc_unifseg_MPMprocess(fn_wMPM, fn_mwTC, opt_process);
fn_out.fn = crc_unifseg_MPMprocess(fn_wMPM, fn_mwTC, opt_process);

end

%__________________________________________________________________________
%
%% SUBFUNCTION
%__________________________________________________________________________
function fn_finMPM = crc_unifseg_MPMprocess(fn_wMPM, fn_mwTC, opt)
%
% FORMAT fn_finMPM = crc_unifseg_MPMprocess(fn_wMPM, fn_TC, opt)
%
% Applyin the "Bogdan treatment" (*) on MPM images that were processed
% with 'unified segmentation' and NOT Dartel (as is usual).
% This is necessary when dealing with patients with brain lesions as GM and
% WM are affected and cannot be properly aligned across subjects.
%
% Since the number of tissue classes is 2 (GM+WM) or 3 (GM+WM+lesion), one
% should be a bit careful how things are processed.
%
% INPUT:
% - fn_wMPM : filename (char array) of the warped MPM
% - fn_mwTC : filename (char array) of the modulated warped tissue classes
% - opt     : structure with some options
%   .fwhm   : kernel size for smoothing
%   .tpm    : tissue probability maps to use, with full path!
%             And one per tissue class passed!
%   .tc_ind : index of tissue classes to be considered, e.g. [1 2 3] for
%             the 1st three classes (e.g. GM/WM/Lesion for the case of
%             brain image with lesion)
%
% OUTPUT:
% - fn_finMPM : filename of the "smoothed tissue masked MPMs".
%
% (*) smoothing and scaling to account for partial volume effect on a
% series of multi-parametric map, e.g. MT, R1, R2*, based on the tissue
% classes, e.g. GM, WM and lesions.
%__________________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by C. Phillips, 2015.
% Cyclotron Research Centre, University of Liege, Belgium

fn_TPM = opt.tpm;

nMPM = size(fn_wMPM,1);
nTC  = size(fn_mwTC,1);
nTPM = size(fn_TPM,1);
if nTC~=nTPM
    error('VBQ:MPM','Wrong number of tissue classes!')
end

ic_flag = struct(...
    'dtype', 16, ... % keep things in floats
    'interp', -4);   % 4th order spline interpolation

fn_finMPM = cell(nMPM,1);
for ii=1:nMPM
    % ii^th MPM to be treated
    fn = fn_wMPM(ii,:) ;
    
    % Get the weighted MPM -> p-images
    p = cell(nTC,1);
    for jj=1:nTC
        % MPM weighted with its own GM/WM/lesion, and a priori>.05
        tmp = char(fn, fn_mwTC(jj,:), fn_TPM(jj,:));
        p_tmp = spm_imcalc(tmp,spm_file(fn,'prefix',['p',num2str(jj),'_']), ...
            '(i1.*i2).*(i3>0.05)',ic_flag);
        p{jj} = p_tmp.fname;
    end
    
    % Smooth TC -> ssmwc1 images
    m = cell(nTC,1);
    for jj=1:nTC
        m{jj} = spm_file(fn_mwTC(jj,:),'prefix','s');
        spm_smooth(fn_mwTC(jj,:),m{jj},opt.fwhm); % smooth mwc(jj)
    end
    
    % Smooth weighted MPM -> sp-images
    n = cell(nTC,1);
    for jj=1:nTC
        n{jj} = spm_file(p{jj},'prefix','s');
        spm_smooth(p{jj},n{jj},opt.fwhm);
    end
    
    % calculate signal, as in paper + masking smoothed TC>.05
    q = cell(nTC,1);
    for jj=1:nTC
        q{jj} = spm_file(p{jj},'prefix','fin_');
        spm_imcalc(char(n{jj},m{jj},m{jj}), q{jj}, ...
            '(i1./i2).*(i3>0.05)',ic_flag);
    end
    
    %     fn_finMPM{ii} = char(q); % saved as char array
    fn_finMPM{ii} = q; % saved as cell array
    fn_2delete = (char(char(p),char(m),char(n)));
    for jj=1:size(fn_2delete,1)
        delete(deblank(fn_2delete(jj,:)));
    end
end


end

