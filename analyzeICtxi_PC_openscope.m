%{
node strength can be calculated with:
-CRF edge potential
-Pearson correlation
-Okun's population coupling
-mutual information between pairs of neurons

sensory sensitivity can be cauculated with:
-CRF AUC
-SP_ICvsRC
-mutual information between vis stim and neural responses
%}
% IMPORTANT NOT TO DO GENPATH
addpath('/Users/hyeyoung/Documents/DATA/OpenScopeData/matnwb_HSLabDesktop')

datadir = '/Users/hyeyoung/Documents/DATA/OpenScopeData/000248/';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

for ises = 1:Nsessions
    clearvars -except datadir nwbsessions ises
    disp(nwbsessions{ises})
    sesclk = tic;

nwbfiles = cat(1, dir([datadir nwbsessions{ises} filesep '*.nwb']), dir([datadir nwbsessions{ises} filesep '*' filesep '*.nwb']));

% take filename  with shortest length or filename that does not contain probe
[~, fileind] = min(cellfun(@length, {nwbfiles.name}));
nwbspikefile = fullfile([nwbfiles(fileind).folder filesep nwbfiles(fileind).name]);
% nwbspikefile = string(nwbspikefile);
disp(nwbspikefile)
nwb = nwbRead(nwbspikefile); %, 'ignorecache');

pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load([pathpp 'postprocessed_probeC.mat'], 'vis')

%% where are the units from? probe position and within-probe electrode
% position? brain area labels?
% There's a 'peak channel' field associated with each unit. 
% The formula for the id is 1000 *probe_val + the local electrode id. 
% This should be the same as the ids that are present in the id field of the 
% electrode section of the NWB (i.e. the 2304th value is 5383 because it is 
% the 383rd value of the 5th probe.)

% nwb.general_extracellular_ephys
electrode_probeid = nwb.general_extracellular_ephys_electrodes.vectordata.get('probe_id').data.load();
electrode_localid = nwb.general_extracellular_ephys_electrodes.vectordata.get('local_index').data.load();
electrode_id = 1000*electrode_probeid + electrode_localid;
electrode_location = nwb.general_extracellular_ephys_electrodes.vectordata.get('location').data.load();

disp(unique(electrode_location)')
% save([pathpp 'info_electrodes.mat'], 'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')

%% extract spike times
unit_ids = nwb.units.id.data.load(); % array of unit ids represented within this session
unit_peakch = nwb.units.vectordata.get('peak_channel_id').data.load();
unit_times_data = nwb.units.spike_times.data.load();
unit_times_idx = nwb.units.spike_times_index.data.load();
% unit_waveform = nwb.units.waveform_mean.data.load();
unit_wfdur = nwb.units.vectordata.get('waveform_duration').data.load();

Nneurons = length(unit_ids);

% all(ismember(unit_peakch, electrode_id))

spiketimes = cell(Nneurons, 1);
last_idx = 0;
for ii = 1:Nneurons
    unit_id = unit_ids(ii);
    
%     assert(unit_trials_idx(i) == unit_times_idx(i), 'Expected unit boundaries to match between trials & spike_times jagged arrays')
    start_idx = last_idx + 1;
    end_idx = unit_times_idx(ii);
    
    spiketimes{ii} = unit_times_data(start_idx:end_idx);
    
    last_idx = end_idx;
end

Tres = 0.001; % 1ms
stlen = ceil((max(unit_times_data)+1)/Tres); % add 1s buffer/padding after the last spike timing

disp([stlen, Nneurons])

% save([pathpp 'info_units.mat'], 'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur') %'unit_times_data', 

%%
tic
elecid = electrode_id+1;
revmapelecid = NaN(max(elecid),1);
revmapelecid(elecid) = 1:numel(elecid);
electrode_location = cellstr(electrode_location);
neuloc = electrode_location(revmapelecid(unit_peakch+1));

% neuctx = [];
% neuctxind = [];
% for iprobe = 1:6
% tempneuoind = find(floor(unit_peakch/1000)==iprobe-1);
% tempneuloc = electrode_location(revmapelecid(unit_peakch(tempneuoind)+1));
% neuctx = cat(1, neuctx, contains(tempneuloc, 'VIS'));
% neuctxind = cat(1, neuctxind, tempneuoind(contains(tempneuloc, 'VIS')));
% end

neuctx = contains(neuloc, 'VIS');
neuctxind = find(neuctx);

Nneuctx = numel(neuctxind);
neulocctx = neuloc(neuctxind);
% isequal(neulocctx, find(contains(neuloc, 'VIS')))

recarealabels = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'};
visarealabels = zeros(Nneuctx,1);
for a = 1:numel(recarealabels)
    if strcmp(recarealabels{a}, 'VISp')
        neuinarea = contains(neulocctx, 'VISp') & ~contains(neulocctx, 'VISpm');
    else
        neuinarea = contains(neulocctx, recarealabels{a});
    end
    visarealabels(neuinarea) = a;
end

kerwinhalf = 25; kersigma = 12/sqrt(2);
kergauss = normpdf( (-kerwinhalf:kerwinhalf), 0,kersigma);
kergauss = (kergauss/sum(kergauss));

ctxspiketrain = false(Nneuctx, stlen);
for ii = 1:Nneuctx
    ci = neuctxind(ii);
    ctxspiketrain(ii, floor(spiketimes{ci}/Tres)+1) = true;
end

ststartend = [floor(min(unit_times_data)/Tres)+1 floor(max(unit_times_data)/Tres)+1];
ICblockstartend = [vis.ICwcfg1_presentations.start_time(1) vis.ICwcfg1_presentations.stop_time(end)];
ICblockstartend = floor(ICblockstartend/Tres)+1;

ctxspiketrain = ctxspiketrain(:,ICblockstartend(1):ICblockstartend(2));

% ctxstgauss = convn(ctxspiketrain, kergauss, 'valid');
Twin = 1000;
Tend = Twin*floor(size(ctxspiketrain,2)/Twin);
ctxstgauss = squeeze(sum(reshape(ctxspiketrain(:,1:Tend), Nneuctx, 1000,Tend/Twin), 2));

save([pathpp 'spiketimes.mat'], 'spiketimes', 'neuloc', 'neuctx', 'neuctxind', 'neulocctx', ...
    'recarealabels', 'visarealabels', 'ststartend', 'ICblockstartend', '-v7.3')

toc
%% population coupling (Okun)

Tlength = size(ctxstgauss,2);
PC = NaN(Nneuctx,1);
PCr = NaN(Nneuctx,1);
tic
for ci = 1:Nneuctx
    vecexcept = true(Nneuctx,1); vecexcept(ci) = false;
    %PC(ci) = (1/sum(ctxspiketrain(ci,:))) * sum(ctxstgauss(ci,:).*sum( ctxstgauss(vecexcept,:)-mean(ctxstgauss(vecexcept,:),2), 1));
    PC(ci) = (1/std(ctxstgauss(ci,:))) * sum(ctxstgauss(ci,:).*sum( ctxstgauss(vecexcept,:)-mean(ctxstgauss(vecexcept,:),2), 1));
    PCr(ci) = corr(ctxstgauss(ci,:)', mean(ctxstgauss(vecexcept,:),1)');
end
toc
disp('done calculating population coupling coefficients')
PCnorm = PC/(Nneuctx*Tlength);
PCz = (PC - mean(PC))./std(PC);

%% population coupling each area
PC_area = NaN(Nneuctx,numel(recarealabels));
PCr_area = NaN(Nneuctx,numel(recarealabels));
tic
for a = 1:numel(recarealabels)
    neuinarea = visarealabels==a;
for ci = 1:Nneuctx
    vecexcept = neuinarea; vecexcept(ci) = false;
    %PC_area(ci,a) = (1/sum(ctxspiketrain(ci,:))) * sum(ctxstgauss(ci,:).*sum( ctxstgauss(vecexcept,:)-mean(ctxstgauss(vecexcept,:),2), 1));
    PC_area(ci,a) = (1/std(ctxstgauss(ci,:))) * sum(ctxstgauss(ci,:).*sum( ctxstgauss(vecexcept,:)-mean(ctxstgauss(vecexcept,:),2), 1));
    PCr_area(ci,a) = corr(ctxstgauss(ci,:)', mean(ctxstgauss(vecexcept,:),1)');
end
end
toc
disp('done calculating inter-area population coupling coefficients')
PCnorm_area = PC_area/(Nneuctx*Tlength);
PCz_area = (PC_area - mean(PC_area,1))./std(PC_area);

save([pathpp 'PC_openscope.mat'], 'unit_peakch', 'neuloc', 'neuctx', 'neulocctx', ...
    'PC', 'PCr', 'PCnorm', 'PCz', 'PC_area', 'PCr_area', 'PCnorm_area', 'PCz_area')

toc(sesclk)
end
