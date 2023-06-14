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

for ises = 2:Nsessions
    clearvars -except datadir nwbsessions ises
    disp(nwbsessions{ises})
    sesclk = tic;
    if ises==4
        continue
    end

    neuopt = 'ctxL23';
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
elecid = electrode_id+1;
revmapelecid = NaN(max(elecid),1);
revmapelecid(elecid) = 1:numel(elecid);
electrode_location = cellstr(electrode_location);
neuloc = electrode_location(revmapelecid(unit_peakch+1));

load([pathpp 'visresponses.mat'])

switch neuopt
    case 'ctx'
neuctx = contains(neuloc, 'VIS');
    case 'ctrctx'
neuctx = contains(neuloc, 'VIS') & RFCIall.RFindclassic==1 & RFCIall.Pkw_rfclassic<0.05;
    case 'ctxL23'
neuctx = contains(neuloc, 'VIS') & contains(neuloc, '2/3');
end
neuctxind = find(neuctx) ;

%%

% neuctx = [];
% neuctxind = [];
% for iprobe = 1:6
% tempneuoind = find(floor(unit_peakch/1000)==iprobe-1);
% tempneuloc = electrode_location(revmapelecid(unit_peakch(tempneuoind)+1));
% neuctx = cat(1, neuctx, contains(tempneuloc, 'VIS'));
% neuctxind = cat(1, neuctxind, tempneuoind(contains(tempneuloc, 'VIS')));
% end

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

%%
ctxspiketrain = false(Nneuctx, stlen);
ststartend = [floor(min(unit_times_data)/Tres)+1 floor(max(unit_times_data)/Tres)+1];
for ii = 1:Nneuctx
    ci = neuctxind(ii);
    ctxspiketrain(ii, floor(spiketimes{ci}/Tres)+1) = true;
end

ICblockstartend = [vis.ICwcfg1_presentations.start_time(1) vis.ICwcfg1_presentations.stop_time(end)];
ICblockstartend = floor(ICblockstartend/Tres)+1;
ctxspiketrain = ctxspiketrain(:,ICblockstartend(1):ICblockstartend(2));
% ctxspiketrain = ctxspiketrain(:,ststartend(1):ststartend(2));

kersm = ones(1,25);
kersm = (kersm/sum(kersm));

ctxstsm = convn(ctxspiketrain, kersm, 'valid'); % this takes ~2min

% save([pathpp 'spiketimes.mat'], 'spiketimes', 'neuloc', 'neuctx', 'neuctxind', 'neulocctx', ...
%     'recarealabels', 'visarealabels', 'ststartend', 'ICblockstartend', '-v7.3')

%%
Thalfwin = 13;
CCGtli = -Thalfwin:Thalfwin;

tic
ctxCCG = NaN(Nneuctx, Nneuctx, length(CCGtli));
for ci = 1:Nneuctx
    for t = 1:length(CCGtli)
        clearvars tempsource temptargmat
        tinds_source = 1+Thalfwin:size(ctxspiketrain,2)-Thalfwin;
        tinds_target = tinds_source + CCGtli(t);
        
        tempsource = ctxspiketrain(ci,tinds_source);
        temptargmat = ctxspiketrain(ci+1:end,tinds_target);
        tempdot = sum(temptargmat(:,tempsource),2);
        %tempdot = ctxspiketrain(ci+1:end,tinds_target)*ctxspiketrain(ci,tinds_source)' ;
        %tempgeomean = 1000*sqrt(mean(ctxspiketrain(ci+1:end,:),2)*mean(ctxspiketrain(ci,:),2)); % normalization used in ctrctxCCG.mat
        tempgeomean = sqrt(sum(ctxspiketrain(ci+1:end,:),2)*sum(ctxspiketrain(ci,:),2));
        tempvec = tempdot./tempgeomean;
        ctxCCG(ci,ci+1:end,t) = tempvec;
        
        tminus = find(CCGtli == -CCGtli(t));
        ctxCCG(ci+1:end,ci,tminus) = tempvec;
    end
end
toc
disp('calculated ctxCCG')

%{
tic
ctxCCGsm = NaN(Nneuctx, Nneuctx, length(CCGtli));
for ci = 1:Nneuctx
    tic
    for t = 1:length(CCGtli)
        tinds_source = 1+Thalfwin:size(ctxstsm,2)-Thalfwin;
        tinds_target = tinds_source + CCGtli(t);
        
        tempdot = ctxstsm(ci+1:end,tinds_target)*ctxstsm(ci,tinds_source)' ;
        tempgeomean = sqrt(sum(ctxstsm(ci+1:end,:),2)*sum(ctxstsm(ci,:),2));
        tempvec = tempdot./tempgeomean;
        ctxCCGsm(ci,ci+1:Nneuctx,t) = tempvec;
        
        tminus = find(CCGtli == -CCGtli(t));
        ctxCCGsm(ci+1:Nneuctx,ci,tminus) = tempvec;
    end
    toc
end
toc
disp('calculated ctrctxCCGsm')

ctxCCGjc = ctxCCG - ctxCCGsm;
ctxCCGweight = squeeze( sum(ctxCCGjc(:,:,CCGtli>0),3) - sum(ctxCCGjc(:,:,CCGtli<0),3) );
save([pathpp 'ctrctxCCG.mat'], 'CCGtli', 'ctrctxCCG', 'ctrctxCCGsm', 'ctrctxCCGjc', 'ctrctxCCGweight', '-v7.3')
%}

ctxCCGweight = squeeze( sum(ctxCCG(:,:,CCGtli>0),3) - sum(ctxCCG(:,:,CCGtli<0),3) );

switch neuopt
    case 'ctx'
        save([pathpp 'ctxCCG.mat'], 'stlen', 'neuctx', 'neulocctx', 'visarealabels', 'CCGtli', 'ctxCCG', 'ctxCCGweight', '-v7.3')
    case 'ctrctx'
        neuctrctx = neuctx;
        neulocctrctx = neulocctx;
        ctrctxCCG = ctxCCG;
        ctrctxCCGweight = ctxCCGweight;
        save([pathpp 'ctrctxCCG.mat'], 'stlen', 'neuctrctx', 'neulocctrctx', 'visarealabels', 'CCGtli', 'ctrctxCCG', 'ctrctxCCGweight', '-v7.3')
    case 'ctxL23'
        neuctxL23 = neuctx;
        neulocctxL23 = neulocctx;
        ctxL23CCG = ctxCCG;
        ctxL23CCGweight = ctxCCGweight;
        save([pathpp 'ctxL23CCG.mat'], 'stlen', 'neuctxL23', 'neulocctxL23', 'visarealabels', 'CCGtli', 'ctxL23CCG', 'ctxL23CCGweight', '-v7.3')
end


toc(sesclk)
end


%%
% n_t = length(tinds_source);
% p = ceil(log2(n_t));
% NFFT = 2^p;
% 
% ci = 18; cj = 119;
% tinds_source = 1+Thalfwin:size(ctrctxspiketrain,2)-Thalfwin;
% tempsource = double( ctrctxspiketrain(ci,tinds_source)' );
% temptarget = double( ctrctxspiketrain(cj,tinds_source)' );
% 
% tempccg = squeeze(ctxCCG(ci,cj,:));
% tic
% CCG = fftshift(ifft(fft(tempsource, NFFT) .* conj(fft(temptarget, NFFT))));
% toc
