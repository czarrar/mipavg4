function check_output_file(filename, fidLOG)
% Logs the existence of given output file
%   check_output_file(filename, fidLOG)
%   
%   Created by Zarrar Shehzad on 2012-09-08.
%

if ~exist(filename, 'file')
    fprintf(fidLOG, 'Specified file ''%s'' does not currently exist: Creating file...\n', filename);
else
    fprintf(fidLOG, 'Specified file ''%s'' already exists: Overwriting old file...\n', filename);
end

end %  function
