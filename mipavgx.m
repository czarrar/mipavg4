function mipavgx(varargin)

%%==============================================================================
%%                                              Initialize file and option flags
%%==============================================================================

nEEGfiles = 0; nSGCfiles = 0; nMONfiles = 0; nLOGfiles = 0; 
nOUTprefixes = 0; nARFfiles = 0; nRCDfiles = 0;

EEGfiles = {}; SGCfiles = {}; RCDfiles = {};
OUTprefix = []; ARFfile = [];

lowpass_value = 0; highpass_value = 0; notchfilter_option = 0;
typeEEG = 0; typeCNT = 0; typeSET = 0; verbose = 0;

output = []; output_opts = {'eegad', 'eeglab', 'fieldtrip', 'econnectome', ...
                            'continuous_eeglab', 'continuous_fieldtrip'};
for i=1:length(output_opts); output.(output_opts{i}) = 0; end

conversion_factor = 0; adunitsPerUvolt = 1;

% old stuff...not sure to keep or not
display_option = 0;

% Arguments for pop_loadXXX functions
load_args = {};


%%==============================================================================
%%                                                          Parse user arguments
%%==============================================================================

for j = 1:length(varargin)
    arg_string = lower(char(varargin{j}));
    if(strfind(arg_string,'.eeg'))
        EEGfiles{end+1} = varargin{j};
        typeEEG = 1;
        nEEGfiles = nEEGfiles + 1;
    elseif(strfind(arg_string,'.cnt'))
        EEGfiles{end+1} = varargin{j};
        typeCNT = 1;
        nEEGfiles = nEEGfiles + 1;
    elseif(strfind(arg_string, '.set'))
        EEGfiles{end+1} = varargin{j};
        typeSET = 1;
        nEEGfiles = nEEGfiles + 1;
    elseif(strfind(arg_string,'.sgc'))
        SGCfiles{end+1} = varargin{j};
        nSGCfiles = nSGCfiles + 1;
    elseif(strfind(arg_string,'.log'))
        LOGfile = varargin{j};
        nLOGfiles = nLOGfiles + 1;
    elseif(strfind(arg_string,'.arf'))
        ARFfile = varargin{j};
        nARFfiles = nARFfiles + 1;
    elseif(strfind(arg_string,'.rcd'))
        RCDfiles{end+1} = varargin{j};
        nRCDfiles = nRCDfiles + 1;
    elseif(strfind(arg_string,'-p'))
        OUTprefix = strtrim(varargin{j}(strfind(arg_string, '-p')+2:end));
        nOUTprefixes = nOUTprefixes + 1;
%     elseif(strfind(arg_string,'-d'))
%         display_option = 1;
%         temp = char(varargin{j});
%         display_channel = fix(str2double(temp(strfind(temp,'-d')+2:end)));
%         if(isempty(display_channel)) display_channel = 1; end
%         nDisplayChannels = length(display_channel);
    elseif(strfind(arg_string,'-f'))
        filter_values = arg_string(strfind(arg_string,'-f')+2:end);
        filter_values = str2double(filter_values);
        if length(filter_values) ~= 2; error('Unrecognized option for -f'); end
        lowpass_value = filter_values(1); highpass_value = filter_values(2);
    elseif(strfind(arg_string,'-n'))
        notchfilter_option = 1;
    elseif(strfind(arg_string,'-o'))
        o = lower(strtrim(varargin{j}(strfind(arg_string,'-o')+2:end)));
        if isempty(intersect(o, output_opts))
            error('Unrecognized argument ''%s'' for -o', o); 
        end
        output.(o) = 1;
    % Legacy options for MIP files
    elseif(strfind(arg_string,'-c'))
        temp = arg_string;
        conversion_factor = str2double(temp(strfind(temp,'-c')+2:end));
        if(isempty(conversion_factor)) conversion_factor = 1; end
        load_args{end+1} = 'conversion';
        load_args{end+1} = conversion_factor;
        clear temp conversion_factor;
    % TODO: do i want anyone to have a montage file?
    elseif(strfind(arg_string,'.mon'))
        load_args{end+1} = 'labels_file';
        load_args{end+1} = varargin{j};
        nMONfiles = nMONfiles + 1;
    elseif(strfind(arg_string,'-v'))
        verbose = 1;
    else
        fprintf('\n Invalid argument: %s\n',arg_string);
    end
end



%%==============================================================================
%%                                    Check that we have all the necessary files
%%==============================================================================

if sum([typeEEG typeCNT typeSET]) > 1; error('Cannot use any combination of .cnt, .eeg, or .set files simultaneously'); end

if(~nEEGfiles); error('No EEG file was specified. Program cannot proceed'); end

if (~nSGCfiles)
    error('No SGC file was specified. Program cannot proceed');
elseif (nSGCfiles > 1 && nSGCfiles ~= nEEGfiles)
    error('Insufficient SGC files specified to match EEG files. Program cannot proceed');
end

if nMONfiles > 1
  error('More than one montage file was specified. Program cannot proceed');
end

if ~nOUTprefixes
    fprintf('\nNo output prefix was specified. Using MipAvg.\n');
    OUTprefix = 'MipAvg';
    nOUTprefixes = 1;
end
if nOUTprefixes == 1
    EEGLABcontinuousfile = [OUTprefix '_continuous.set'];
    FIELDTRIPcontinuousfile = [OUTprefix '_ft_continuous.mat'];
    EEGLABprefix = [OUTprefix '_trials'];
    FIELDTRIPprefix = [OUTprefix '_ft_trials'];
    ECONNECOMEfile = [OUTprefix '_econnectome.mat'];
    AVGfile = [OUTprefix '.avg'];
    HDRfile = [OUTprefix '.hdr'];
else
    error('More than one output prefix was specified. Program cannot proceed');
end

if(~nLOGfiles)
    % Maybe the default should be to print to the stdout?
    fprintf('\nNo LOG file was specified. Creating MipAvg.log in working directory\n');
    fidLOG = fopen('MipAvg.log', 'w');
elseif(nLOGfiles == 1)
    fidLOG = fopen(LOGfile, 'w');
else
    error('More than one LOG file specified. Program cannot proceed');
end

if(~nARFfiles) %% IF no arf file, use variance on channel specified in scg file
    fprintf('\nNo ARF file was specified. No artifact rejection will be applied.\n');
elseif(nARFfiles > 1)
    error('More than one ARF file specified. Program cannot proceed');
end


%%==============================================================================
%%                                                       Initialize the log file
%%==============================================================================

fprintf(fidLOG,'\n MIPAVG3 Program executed on %s\n\n', datestr(now));

if (conversion_factor)
    fprintf(fidLOG,'\n Calibration override specified. Microvolt conversion factor set equal to %5.2f\n', conversion_factor);
end


%%==============================================================================
%%                                               Add necessary functions to path
%%==============================================================================

if(exist('pop_importdata.m','file') ~= 2)
    if(exist('/usr/local/packages/MATLABPackages/eeglab/functions/popfunc','dir') == 7)
        fprintf(fidLOG, '\nLoading required EEGLab functions\n');
        addpath('/usr/local/packages/MATLABPackages/eeglab/functions/guifunc');
        addpath('/usr/local/packages/MATLABPackages/eeglab/functions/popfunc');
        addpath('/usr/local/packages/MATLABPackages/eeglab/functions/adminfunc');
        addpath('/usr/local/packages/MATLABPackages/eeglab/functions/sigprocfunc');
    else
        fprintf(fidLOG, '\nRequired eeglab functions cannot be found - check you paths!\n');
        fclose('all');
        return
    end
end


%%==============================================================================
%%                                                              Read in EEG data
%%==============================================================================

fprintf(fidLOG,'\n Reading in EEG files:\n');
ALLEEG = [];
for i=1:nEEGfiles
    fprintf(fidLOG,'\n\t%s\n', EEGfiles{i});
    if typeEEG
        EEG = pop_loadmip(EEGfiles{i}, load_args{:});
    elseif typeCNT
        EEG = pop_loadcnt(EEGfiles{i}, 'dataformat', 'int32');
    elseif typeSET
        EEG = pop_loadset(EEGfiles{i});
    end
    [ALLEEG EEG ~] = eeg_store(ALLEEG, EEG);
end


%%==============================================================================
%%                                                      Recode events (optional)
%%==============================================================================

if(nRCDfiles)
    fprintf(fidLOG, '\n=================\n');
    fprintf(fidLOG,'Recoding section\n');
    
    ALLrecodes  = read_recode(RCDfiles, fidLOG);
    ALLEEG      = recode_eeg(ALLEEG, ALLrecodes, fidLOG);
    
    fprintf(fidLOG, '\n=================\n\n');
end


%%==============================================================================
%%            Read in SGC file with epoch/baseline limits and trial codes/labels
%%==============================================================================

[epoch baseline first_pt_msec ALLcodes ALLlabels] = read_sgc(SGCfiles, fidLOG);


%%==============================================================================
%%                                Map codes and labels from SGC file to EEG file
%%==============================================================================

ALLEEG = map_codes_and_labels(ALLEEG, ALLcodes, ALLlabels, fidLOG);


%%==============================================================================
%%                      Merge data, fix channel labels, and save continuous data
%%==============================================================================


% Merge data together
if nEEGfiles > 1, EEG = pop_mergeset(ALLEEG, 1:nEEGfiles); 
else EEG = ALLEEG{1}; end
clear ALLEEG;

% Ghetto fix of channel labels for cnt files
old_vals = {'CB1', 'CB2', 'HEO', 'VEO'};
new_vals = {'I1', 'I2', 'HEOG', 'VEOG'};
for i=1:length(old_vals)
    search = ismember({EEG.chanlocs.labels}, old_vals{i});
    if any(search); EEG.chanlocs(search).labels = new_vals{i}; end
end

% Save the EEGLab continuous file (will be removed later if user-specified)
check_output_file(EEGLABcontinuousfile, fidLOG);
pop_saveset(EEG, 'filename', EEGLABcontinuousfile);

% Save the fieldtrip continuous file
if output.continuous_fieldtrip
    check_output_file(FIELDTRIPcontinuousfile, fidLOG);
    data            = eeglab2fieldtrip(EEG, 'preprocessing');
    data.hdr        = ft_read_header(EEGLABcontinuousfile);
    data.cfg.event  = ft_read_event(EEGLABcontinuousfile);
    save(FIELDTRIPcontinuousfile, 'data');
end


%%==============================================================================
%%                                                    Epoching and Preprocessing
%%==============================================================================

% Define trials
cfg = [];
cfg.dataset                 = EEGLABcontinuousfile;
cfg.continuous              = 'true';
cfg.trialfun                = 'mip_trialfun';
cfg.trialdef.eventtype      = 'trigger';
cfg.trialdef.eventvalue     = EEG.eventcodes;
cfg.trialdef.prestim        = -1 * epoch(1);
cfg.trialdef.poststim       = epoch(2);
cfg = ft_definetrial(cfg);

% Preprocessing & Creating Epochs
cfg.demean = 'yes';
cfg.baseline = baseline;
if lowpass_value && highpass_value
    cfg.bpfilter    = 'yes';
    cfg.bpfreq      = [lowpass_value highpass_value];
elseif lowpass_value
    cfg.hpfilter    = 'yes';
    cfg.hpfreq      = lowpass_value;
elseif highpass_value
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = highpass_value;
end
if notchfilter_option
    cfg.dftfilter   = 'yes';
    % TODO: What should line filtering be? [50, 100, 150] or [60, 120, 180]
    cfg.dftfreq     = 60;
end
data = ft_preprocessing(cfg);

% Remove continuous EEGLAB file
if ~output.continuous_eeglab
    fprintf(fidLOG, '\nRemoving continuous set file\n');
    [pathstr name] = fileparts(EEGLABcontinuousfile);
    delete(EEGLABcontinuousfile);
    delete([pathstr filesep name '.fdt']);
end


%%==============================================================================
%%                                                            Artifact Rejection
%%==============================================================================

% Read in artifact rejection settings
ar_settings = read_arf(ARFfile, data.label, epoch, fidLOG);

% Loop through and apply each user-defined artifact rejection criteria
cfg = [];
cfg.continuous = 'no';
for i=1:length(ar_settings)
    switch ar_settings(i).method
        case 'flat'
            cfg.trl                             = data.cfg.trl;
            cfg.artfctdef.clip.channel          = ar_settings(i).channels;
            cfg.artfctdef.clip.prestim          = ar_settings(i).prestim;
            cfg.artfctdef.clip.poststim         = ar_settings(i).poststim;
            cfg.artfctdef.clip.thresh           = ar_settings(i).criteria;
            cfg = run_ft_artifact(cfg, data, 'clip');
            fprintf(fidLOG, '\n%d trials rejected with flat', size(cfg.artfctdef.clip.artifact, 1));
        case 'ppa'
            cfg.trl                             = data.cfg.trl;
            % adjust starting sample based on prestim
            if ar_settings(i).prestim < -1*epoch(1)
                cfg.trl(:,1) = cfg.trl(:,1) + (-1*epoch(1) - ar_settings(i).prestim) * data.fsample;
            end
            % adjust ending sample based on poststim
            if ar_settings(i).poststim < epoch(2)
                cfg.trl(:,2) = cfg.trl(:,2) + (ar_settings(i).poststim - epoch(2)) * data.fsample;
            end
            cfg.artfctdef.threshold.channel     = ar_settings(i).channels;
            cfg.artfctdef.threshold.range       = ar_settings(i).criteria;
            cfg.artfctdef.threshold.bpfilter    = 'no';
            for j=1:length(ar_settings(i).opts)
                if strcmpi(ar_settings(i).opts{j}, 'yes-freq')
                    cfg.artfctdef.threshold.bpfilter = 'yes';
                end
            end
            % TODO: are the other filter settings good?
            cfg = run_ft_artifact(cfg, data, 'threshold');
            fprintf(fidLOG, '\n%d trials rejected with ppa', size(cfg.artfctdef.threshold.artifact, 1));
        case 'zthr'
            cfg.trl                             = data.cfg.trl;
            % adjust starting sample based on prestim
            if ar_settings(i).prestim < -1*epoch(1)
                cfg.trl(:,1) = cfg.trl(:,1) + (-1*epoch(1) - ar_settings(i).prestim) * data.fsample;
            end
            % adjust ending sample based on poststim
            if ar_settings(i).poststim < epoch(2)
                cfg.trl(:,2) = cfg.trl(:,2) + (ar_settings(i).poststim - epoch(2)) * data.fsample;
            end
            cfg.artfctdef.eog.channel           = ar_settings(i).channels;
            cfg.artfctdef.eog.cutoff            = ar_settings(i).criteria; % z-value
            cfg.artfctdef.eog.bpfilter          = 'yes';
            cfg.artfctdef.eog.bpfilttype        = 'fir';
            cfg.artfctdef.eog.hilbert           = 'yes';
            cfg.artfctdef.eog.interactive       = 'no';
            cfg.artfctdef.eog.trlpadding        = -0.1; % Padding added for each trial in 
                                                        %   checking for an artifact.
            cfg.artfctdef.eog.fltpadding        = 0.1;  % Padding added before filtering
                                                        %   and removed after filtering.
            cfg.artfctdef.eog.artpadding        = 0.1;  % Padding added around period of time
                                                        %   determined to be an artifact,
            for j=1:length(ar_settings(i).opts)
                if strcmpi(ar_settings(i).opts{j}, 'no-freq')
                    cfg.artfctdef.eog.bpfilter  = 'no';
                    cfg.artfctdef.eog.bpfilttype= 'fir';
                    cfg.artfctdef.eog.hilbert   = 'no';
                elseif strcmpi(ar_settings(i).opts{j}, 'yes-interactive')
                    cfg.artfctdef.eog.interactive = 'yes';
                end
            end
            % TODO: does there need to be an option for the filt freq range?
            cfg = run_ft_artifact(cfg, data, 'eog');
            fprintf(fidLOG, '\n%d trials rejected with zthr', size(cfg.artfctdef.eog.artifact, 1));
        otherwise
            fprintf(fidLOG, '\nUnrecognized artifact rejection approach: %s\n', ...
                     ar_settings(i).method);
    end
end

% Remove offending trials
cfg.artfctdef.reject    = 'complete';
data = ft_rejectartifact(cfg, data);

% Indices of good trials
cleaninds = arrayfun(@(x) find(x == data.cfg.trlold(:,1)), data.cfg.trl(:,1));


%%==============================================================================
%%                                       Write Output (EEGAD, EEGLab, Fieldtrip)
%%==============================================================================

ntpts   = round(diff(epoch) * data.hdr.Fs) + 1;
nchans  = data.hdr.nChans;
ncodes  = length(EEG.eventcodes);

if output.eegad
    % Create average ERP matrix for EEGAD
    avg     = zeros(nchans+1, ntpts, ncodes);

    % AVG file (EEGAD Output)
    fprintf(fidLOG, 'Creating AVG file: %s\n', AVGfile);
    check_output_file(AVGfile, fidLOG);
    fidAVG = fopen(AVGfile,'w');
end

% Split trials based on unique event code
for i=1:ncodes
    code = EEG.eventcodes(i); label = EEG.eventlabels{i};
    
    %% EEGLAB
    if output.eeglab
        % Select clean trials with unique event code and then save
        EEGtrial = pop_epoch(EEG, {code}, epoch, 'eventindices', cleaninds);
        pop_saveset(EEGtrial, 'filename', [EEGLABprefix '_' label '.set']);
    end
    
    
    %% Fieldtrip
    if output.fieldtrip
        % Select trials with unique event code and then save
        trials = data.cfg.trl(:,4) == code;
        tdata = ft_selectdata(data, 'rpt', trials);
        save([FIELDTRIPprefix '_' label '.mat'], 'tdata');
    end
    
    
    %% Average ERP
    
    % Average single trials
    cfg = [];
    timelock = ft_timelockanalysis(cfg, tdata);
    
    % Baseline correct
    cfg = [];
    cfg.baseline = baseline;
    timelock = ft_timelockbaseline(cfg, timelock);
    
    
    %% Econnectome
    if output.econnectome
        tmp = EEG;
        EEG = fieldtrip2econnectome(timelock, label, data.hdr.Fs);
        save(ECONNECOMEfile, 'EEG');
        EEG = tmp; clear tmp;
    end
    
    
    %% EEGAD
    if output.eegad
        % Retain averaged evoked reponse
        avg(1:nchans,:,i) = timelock.avg;
    
        % Add trigger channel info
        % (basically the event code at the onset)
        % TODO: is this done right?
        onset = -epoch(1) * data.hdr.Fs + 1;
        avg(onset,nchans+1,i) = code;
    
        % Save for EEGAD
        fwrite(fidAVG, avg(:,:,i)', 'float32');
    end
end

%% Log Rejected Trials
fprintf(fidLOG,'\n Code\t Label\t Accepted #\t Total #\t Percent Averaged\n');
for j = 1:length(EEG.eventcodes)
    code = EEG.eventcodes(j); label = EEG.eventlabels{j};
    total = sum(data.cfg.trlold(:,4)==code); accepted = sum(data.cfg.trl(:,4)==code);
    fprintf(fidLOG, ' %d\t %s\t %d\t\t %d\t\t %4.2f\n', code, label, accepted, total, accepted/total);
end

%% HDR File (EEGAD Output)
if output.eegad
    fprintf(fidLOG, 'Creating AVG HDR file: %s\n', HDRfile);
    check_output_file(HDRfile, fidLOG);
    fidHDR = fopen(HDRfile, 'w');

    % TODO: WHAT TO DO ABOUT THE EXPERIMENT NAME AND DATE
    % Line 1: expName is the experiment name
    fprintf(fidHDR,'%s\n', 'default'); 
    % Line 2: expDate is the experiment date string
    fprintf(fidHDR,'%s\n', now);
    % Line 3: nChannels is the number of data channels (electrodes).
    fprintf(fidHDR,'%d\n', nchans+1);
    % Line 4: nPoints is the number of data points collected per channel.
    fprintf(fidHDR,'%d\n', ntpts);
    % Line 5: sampling is the sampling rate of the data points in ms/point.
    fprintf(fidHDR,'%d\n', 1000/data.hdr.Fs);
    % Line 6: uvunits is the microvolt conversion factor in raw data units/microvolt.
    fprintf(fidHDR,'%d\n', adunitsPerUvolt);
    % Line 7: onset is the stimulus onset from the first data point in ms.
    fprintf(fidHDR,'%d\n', -first_pt_msec);
    % Line 8+: Event Code Labels
    for k = 1:length(EEG.eventcodes)
        label = char(EEG.eventlabels{k});
        fprintf(fidHDR,'%s\n',label);
    end
    % Line X+: Channel Names
    for k = 1:nchans, fprintf(fidHDR,'%s\n', data.label{k}); end
    fprintf(fidHDR,'%s\n', 'trigger');
end

fclose('all');

end