function [trl, event] = mip_trialfun(cfg)
% mip_trialfun is a trial function for fieldtrip's ft_definetrial.
%  
% It searches for events of type "trigger" and 
% values supplied by the user (through an SGC file)
% 
% It returns similar output as ft_trialfun_general with 2 exceptions:
% - 1st: it adds a 4th column to the trl matrix, which includes the event value
% - 2nd: it adds a field to event called intrl, which indicates if the event is 
%   included as a trial (i.e. is in the trl matrix)

hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

%% event values from eeglab will be string
if isnumeric(cfg.trialdef.eventvalue)
    cfg.trialdef.eventvalue = num2strcell(cfg.trialdef.eventvalue);    
end

%% search for "trigger" events & for specic event values
inds        = strcmp('trigger', {event.type});
inds(inds)  = cellfun(@(c) ~isempty(intersect(c, cfg.trialdef.eventvalue)), {event(inds).value});
codes       = cellfun(@str2num, {event(inds).value});
samples     = [event(inds).sample];

%% determine the number of samples before and after the trigger
pretrig  = -round(cfg.trialdef.prestim  * hdr.Fs);
posttrig =  round(cfg.trialdef.poststim * hdr.Fs);

%% indicate if event is a trial
[event(~inds).intrl] = deal(0); [event(inds).intrl] = deal(1);

%% begin latency; end latency; offset (from 0); event code
trl = [samples + pretrig; samples + posttrig; repmat(pretrig, 1, sum(inds)); codes]';

%% function
function c = num2strcell(n, format)
% num2strcell Convert vector of numbers to cell array of strings
% function c = num2strcell(n, format)
%
% If format is omitted, we use
% c{i} = sprintf('%d', n(i))
  if nargin < 2, format = '%d'; end

  N = length(n);
  c = cell(1,N);
  for i=1:N
    c{i} = sprintf(format, n(i));
  end
end

end