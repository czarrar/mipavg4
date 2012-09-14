function ALLEEG = map_codes_and_labels(ALLEEG, codes, labels, fidLOG)
%   [ALLEEG] = map_codes_and_labels(ALLEEG, codes, labels, fidLOG)
% 
%   Replaces event type values based on mapping in codes, which results in
%   changing [ALLEEG(i).event.type]. Corresponding labels are also added
%   to [ALLEEG(i).urevent.label] and [ALLEEG(i).event.label].
%
%   Inputs:
%       ALLEEG      - struct array of EEG data
%                     contains EEG.event with trial information
%       codes       - a cell array where each element is a 2-column matrix
%                     1st column = input codes, corresponding to [EEG.event.type]
%                     2nd column = output codes, replacing input values
%       labels      - similar to codes, cell array with each element a 2-column matrix
%                     1st column = input labels
%                     2nd column = output labels
%       fidLOG      - optional pointer to log file
%   
%   Outputs:
%       ALLEEG      - struct array of EEG data
%
%   Created by Zarrar Shehzad on 2012-09-06

if nargin < 4; fidLOG = 1; end

fprintf(fidLOG,'\n=========================\n');
fprintf(fidLOG,'\nMapping Codes and Labels\n');

ALLcodes = check_against_ALLEEG(codes);
ALLlabels = check_against_ALLEEG(labels);
[ALLEEG.eventcodes] = deal([]); [ALLEEG.eventlabels] = deal({});
clear codes labels;
n = length(ALLEEG);

% TODO: Should there be more output to the log file?
% TODO: type code is being overwritten as string instead of integer
for i=1:n
    fprintf(fidLOG, '\nProcessing EEG file #%i\n', i);
    
    EEG = ALLEEG(i);
    codes = ALLcodes{i};
    labels = ALLlabels{i};
    
    if size(codes,1) ~= size(labels,1)
       fprintf(fidLOG, '\nSize mismatch between codes and labels\n');
       error('Program cannot proceed');
    end
    
    for j=1:size(codes, 1)
        % Find input codes in EEGLab event data structure
        search = [EEG.event.type] == codes(j,1); 
        if any(search)
            [EEG.urevent([EEG.event(search).urevent]).label] = deal(labels{j,1});
            [EEG.event(search).type] = deal(codes(j,2));
            [EEG.event(search).label] = deal(labels{j,2});
        else
            fprintf(fidLOG, '\tCould not find code %i (label %s)\n', codes(j,1), labels{j,1});
        end
    end
    
    % Were there any codes found?
    if ~isfield(EEG.event, 'label')
        fprintf(fidLOG, '\nNo usuable event codes found for file/run #%d\n', i);
        error('Program cannot proceed');        
    end
    
    % Save unique codes and corresponding labels
    EEG.eventcodes = unique([EEG.event.type]);
    EEG.eventlabels = cell(1, length(EEG.eventcodes));
    for k=1:length(EEG.eventcodes)
        search = find([EEG.event.type] == EEG.eventcodes(k));
        EEG.eventlabels{k} = EEG.event(search(1)).label;
    end
    
    % Remove unique codes/labels with empty labels
    empty_labels = find(cellfun(@isempty, EEG.eventlabels));
    present_labels = find(~cellfun(@isempty, EEG.eventlabels));
    for k=1:length(empty_labels)
        fprintf(fidLOG, '\tEvent %i did not have a matching code/label\n', ...
                    EEG.eventcodes(empty_labels(k)));
    end
    EEG.eventcodes = EEG.eventcodes(present_labels);
    EEG.eventlabels = EEG.eventlabels(present_labels);
    
    ALLEEG(i) = EEG;
end

fprintf(fidLOG,'\n=========================\n');


function variable = check_against_ALLEEG(variable)
%   [variable] = check_against_ALLEEG(variable)
%
%   Checks that variable length matches that of ALLEEG.
%   Will repeat variable if its length is 1 in order to match its length with ALLEEG.
%
    vname=@(x) inputname(1);
    nVAR = length(variable); nEEG = length(ALLEEG);
    if nVAR ~= nEEG;
        if nVAR == 1
          variable = repmat(variable, 1, nEEG); 
        else
          fprintf(fidLOG, '\nLength of %s must match ALLEEG or be 1\n', vname(variable));
          error('Program cannot continue');
        end
    end
end

end %  function