varargin = {'/Users/zarrar/Projects/mipavg/0411121600_abscreen01_run01.cnt', ...
            '/Users/zarrar/Projects/mipavg/0411121600_abscreen01_run02.cnt', ...
            '/Users/zarrar/Dropbox/Programs/current/zshehzad/mipavgx/tests/sample.sgc', ...
            '/Users/zarrar/Dropbox/Programs/current/zshehzad/mipavgx/tests/sample.arf', ...
            '-p /Users/zarrar/Projects/mipavg/output', ...
            '-o eegad', ...
            '-o eeglab', ...
            '-o fieldtrip'}

mipavgx(varargin{:})

cfg = [];
cfg.continuous = 'true';
cfg.dataset = '/Users/zarrar/Projects/tmp.set';
cfg.trialdef.eventtype = 'trigger';
cfg.trialdef.eventvalue = num2strcell(1);
cfg.trialdef.prestim = -1 * epoch(1)/1000;
cfg.trialdef.poststim = epoch(2)/1000;
cfg = ft_definetrial(cfg)

cfg_test = [];
cfg_test.continuous = 'true';
cfg_test.dataset = EEGfiles{1};
cfg_test.dataformat           = 'ns_cnt32';
cfg_test.headerformat         = 'ns_cnt32';
cfg_test.eventformat          = 'ns_cnt32';
cfg_test.trialdef.eventtype = 'trigger';
cfg_test.trialdef.eventvalue = [1:7, 101:166];
cfg_test.trialdef.prestim = epoch(1)/1000;
cfg_test.trialdef.poststim = epoch(2)/1000;
cfg_test = ft_definetrial(cfg_test)
