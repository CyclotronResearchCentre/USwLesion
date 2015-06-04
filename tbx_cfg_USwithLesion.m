function USwLtools = tbx_cfg_USwithLesion
% MATLABBATCH Configuration file for toolbox 'USwL', i.e. 'Unified 
% segmentation with lesion'
%_______________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Christophe Phillips

if ~isdeployed,
    addpath(fullfile(spm('Dir'),'toolbox','USwithLesion'));
end

%----------------------------------------------------------------------
% Setting up main choices
%----------------------------------------------------------------------
USwLtools         = cfg_choice;
USwLtools.tag     = 'USwLtools';
USwLtools.name    = 'US with Lesion Tools';
USwLtools.help    = {
    ['zae',...
    'qsdf']
    ['Tiouop ',...
    'ghhgk']
    }';
USwLtools.values  = {tbx_scfg_USwL tbx_scfg_MPMsmooth tbx_scfg_ParEx};

end