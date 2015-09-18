function flags = crc_check_flag(flags_o,flags)

% FORMAT flags = crc_check_flag(flags_o,flags)
%
% Function to automatically check the content of a "flag" structure, using
% a "default flag structure", adding the missing fields and putting in the 
% default value if none was provided.
%
% INPUT:
% flags_o   default or reference structure
% flags     input flag/option structure that need to be filled for missing
%           fields with default values
%
% OUPUT:
% flags     filled flag/option structure
%__________________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by Y. Leclercq & C. Phillips, 2008.
% Revised and updated by C. Phillips, 2015
% Cyclotron Research Centre, University of Liege, Belgium

f_names = fieldnames(flags_o);
% list fields in default structure

Nfields = length(f_names);
for ii=1:Nfields
    % Update the output if
    % - a field is missing
    % - the field is empty when it shouldn't
    if ~isfield(flags,f_names{ii}) || ...
            ( isempty(flags.(f_names{ii})) && ~isempty(flags_o.(f_names{ii})) )
        flags.(f_names{ii}) = flags_o.(f_names{ii});
    end
end

end