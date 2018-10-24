function fn_out = tbx_run_USwL(varargin)
% Toolbox USwL job execution function
% takes a harvested job data structure and calls tbx function to perform
% computations on the data.
% 
% INPUT:
% job    - harvested job data structure (see matlabbatch help)
% 
% OOUTPUT:
% fn_out - computation results, a struct variable with the filenames of 
%          generated images.
%_______________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

job = varargin{1}; % To follow SPM tradition...

% Deal with intensity bias correction option
if isfield(job.options.bias,'bias_no')
    % bias correction
    job.options.biasreg = 10;
    job.options.biasfwhm = Inf;
    job.options.biaswr = [0 0];
else
    % with bias correction
    job.options.biasreg = job.options.bias.bias_yes.biasreg;
    job.options.biasfwhm = job.options.bias.bias_yes.biasfwhm;
    job.options.biaswr = job.options.bias.bias_yes.biaswr;
end
job.options = rmfield(job.options,'bias');

%% Collect input -> to fit into previously written code. :-)
fn_in{1} = spm_file(job.imgMsk{1},'number',''); % Mask image
fn_in{2} = spm_file(job.imgRef{1},'number',''); % Structural reference
fn_in{3} = char(spm_file(job.imgStruc,'number','')); % All Structurals
fn_in{4} = char(spm_file(job.imgOth,'number','')); % Other images

options = job.options;

% For back compatibility issue
if ~isfield(options,'cleanup') || ...
        isempty(options.cleanup) || strcmp(options.cleanup,'<UNDEFINED>')
    options.cleanup = 0;
end

% Call processing routine and catch output
fn_out = crc_USwL(fn_in,options);

end