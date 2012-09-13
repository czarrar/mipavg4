function recodes = read_recode(RCDfiles, fidLOG)
%   READ_RECODE   Read codes from file allowing one to recode events
%     [RECODEs] = READ_RECODE(RCDFILES)
% 
%   Input:
%       RCDfiles    - cell array of filenames and in each file there should
%                     be two columns, one for the original event codes and 
%                     another for the new event codes
%       fidLOG      - pointer to log file
%
%   Output:
%       recodes     - cell array of length RCDfiles with each element being a
%                     2 column matrix of [orig-code new-code]
%   
%   Created by Zarrar Shehzad on 2012-09-12.

    nRCDfiles   = length(RCDfiles);
    recode      = cell(1, nRCDfiles);
    for i=1:nRCDfiles
        recode{i} = dlmread(RCDfiles{i});
        fprintf(fidLOG,'\nRecode file specified: %s\n', RCDfiles);
        fprintf(fidLOG,'Number of recodes specified = %d\n', size(recode{i},1));
    end
    
end %  function