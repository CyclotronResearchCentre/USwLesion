function varargout = crc_USwL_get_defaults(defstr, varargin)
% Get/set the defaults values associated with an identifier
% FORMAT uswl_def = crc_USwL_get_defaults
% Return the global "uswl_def" variable defined in crc_USwL_defaults.m.
%
% FORMAT defval = crc_USwL_get_defaults(defstr)
% Return the defaults value associated with identifier "defstr". 
% Currently, this is a '.' subscript reference into the global  
% "uswl_def" variable defined in crc_USwL_defaults.m.
%
% FORMAT spm_get_defaults(defstr, defval)
% Set the defaults value associated with identifier "defstr" to defval.
% The new defaults value applies immediately to:
% * new modules in batch jobs
% * modules in batch jobs that have not been saved yet
% This value will not be saved for future sessions of SPM/USwLesion.
% To make persistent changes, see help section in crc_USwL_defaults.m.
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging
% Copyright (C) 2015 Cyclotron Research Centre

% Volkmar Glauche
% Then modified for use with the USwLesion toolbox by Christophe Phillips
% Cyclotron Research Centre, University of Liege, Belgium


global uswl_def;

if isempty(uswl_def)
    crc_USwL_defaults;
end

if nargin == 0
    varargout{1} = uswl_def;
    return
end

% construct subscript reference struct from dot delimited tag string
tags = textscan(defstr,'%s', 'delimiter','.');
subs = struct('type','.','subs',tags{1}');

if nargin == 1
    varargout{1} = subsref(uswl_def, subs);
else
    uswl_def = subsasgn(uswl_def, subs, varargin{1});
end

%% NOTICE
% This file is largely inspired from routines available the SPM toolbox:
% http://www.fil.ion.ucl.ac.uk/spm/