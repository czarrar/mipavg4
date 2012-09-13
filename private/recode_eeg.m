function ALLEEG = recode_eeg(ALLEEG, ALLrecodes, fidLOG)
%   RECODE_EEG   Changes the event codes based on recode matrix
%     [ALLEEG] = RECODE_EEG(ALLEEG, ALLRECODES, FIDLOG)
% 
%   Inputs:
%       ALLEEG      - struct array of EEGLab data structures
%       ALLrecodes  - cell array of recode matrices, each matrix should be 
%                     2 columns and have as many rows as there are events in
%                     the associated EEG object
%       fidLOG      - pointer to the LOG file
%   
%   Created by Zarrar Shehzad on 2012-09-12.

n = length(ALLEEG);

if n ~= length(ALLrecodes)
    fprintf(fidLOG, '\nMismatch between number of EEG files (%d) and recode files (%d)\n', n, length(ALLrecodes));
    fclose('all'); error('Program cannot proceed');
end

for i=1:n
    EEG = ALLEEG(i); recodes = ALLrecodes(i);
    
    norig = length(EEG.event); nrecodes = size(recodes, 1);
    if norig ~= nrecodes
        fprintf(fidLOG,'\nMismatch between number of original and recodes: %d %d', norig, nrecodes);
        fclose('all'); error('Program cannot proceed');
    end
    
    fprintf(fidLOG,'\nSeq\t Old-EEG\t Old-Recode\t New-Recode\n');
    
    for j=1:norig
        fprintf(fidLOG,'%4d\t %4d\t %4d\t %4d\n', j, EEG.event(j), recodes(j,1), recodes(j,2));
        if EEG.event(j) ~= recodes(j,1)
            fprintf(fidLOG,'\nError in recoding - code mismatch\n');
            fclose('all'); error('Program cannot proceed');
        elseif ischar(EEG.event(j).type)
            EEG.event(j).type = char(recodes(j,2));
        else
            EEG.event(j).type = recodes(j,2);
        end
    end
    
    ALLEEG(i) = EEG;
end

end %  function