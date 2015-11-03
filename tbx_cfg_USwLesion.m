function USwLtools = tbx_cfg_USwLesion
% MATLABBATCH Configuration file for toolbox 'USwL', i.e. 'Unified 
% segmentation with lesion'.
% 
% More details (and updates) can be found on:
% https://github.com/CyclotronResearchCentre/USwLesion
%_______________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

if ~isdeployed,
    addpath(fullfile(spm('Dir'),'toolbox','USwLesion'));
end

%----------------------------------------------------------------------
% Setting up main choices
%----------------------------------------------------------------------
USwLtools         = cfg_choice;
USwLtools.tag     = 'USwLtools';
USwLtools.name    = 'US with Lesion Tools';
USwLtools.help    = {
    ['Tools for the segmentation of MR images with lesions ',...
    'based on the ''unified segementation'' approach.']
    ['One need to provide an (approximate) mask of the lesioned area(s)',...
    'in order to use this tool!']
    }';
USwLtools.values  = {tbx_scfg_USwL tbx_scfg_MPMsmooth tbx_scfg_ParEx};

end