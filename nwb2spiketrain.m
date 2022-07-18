% addpath(genpath('/Users/hyeyoung/Documents/CODE/matnwb'))
addpath(genpath('/Users/hyeyoung/Downloads/matnwb-master'))
datapath = '/Users/hyeyoung/Documents/OpenScopeData/';
savepath = '/Users/hyeyoung/Documents/OpenScopeAnalyzed/';

nwb = nwbRead('/Users/hyeyoung/Documents/OpenScopeData/000248/sub-Placeholder/sub-Placeholder.nwb');

unit_ids = nwb.units.id.data.load(); % array of unit ids represented within this 
unit_times_data = nwb.units.spike_times.data.load();
unit_times_idx = nwb.units.spike_times_index.data.load();

spiketimes = cell(size(unit_ids));
last_idx = 0;
for ii = 1:length(unit_ids)
    unit_id = unit_ids(ii);
    
%     assert(unit_trials_idx(i) == unit_times_idx(i), 'Expected unit boundaries to match between trials & spike_times jagged arrays')
    start_idx = last_idx + 1;
    end_idx = unit_times_idx(ii);
    
    spiketimes{ii} = unit_times_data(start_idx:end_idx);
    
    last_idx = end_idx;
end