%{
conditions = {
    "0": {
        "duration": 1,
        "name": "fast_pulses",
        "condition": "2 ms pulses at 1 Hz"
    },
    "1": {
        "duration": 1,
        "name": "pulse",
        "condition": "a single 10ms pulse"
    },
    "2": {
        "duration": .2,
        "name": "pulse",
        "condition": "1 second of 5Hz pulse train. Each pulse is 2 ms wide"
    },
    "3": {
        "duration": .1,
        "name": "raised_cosine",
        "condition": "half-period of a cosine wave"
    },
    "4": {
        "duration": .05,
        "name": "5 hz pulse train",
        "condition": "Each pulse is 10 ms wide"
    },
    "5": {
        "duration": .033,
        "name": "40 hz pulse train",
        "condition": "Each pulse is 6 ms wide"
    },
    "6": {
        "duration": .025,
        "name": "fast_pulses",
        "condition": "1 second of 40 Hz pulse train. Each pulse is 2 ms wide"
    },
    "7": {
        "duration": 0.02,
        "name": "pulse",
        "condition": "a single square pulse"
    },
    "8": {
        "duration": 0.0167,
        "name": "pulse",
        "condition": "a single square pulse"
    },
    "9": {
        "duration": .0125,
        "name": "raised_cosine",
        "condition": "half-period of a cosine wave"
    },
    "10": {
        "duration": .01,
        "name": "100 hz pulse train",
        "condition": "1 second of 100 Hz pulse train. Each pulse is 2 ms wide"
    },
    "11": {
        "duration": 1.0,
        "name": "Square Pulse",
        "condition": "1 second square pulse: continuously on for 1s"
    }
%}

addpath(genpath(pwd))


datadir = '/Users/hyeyoung/Documents/DATA/OpenScopeData/000248/';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name}; 
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));

for ises = 1:numel(nwbsessions)
clearvars -except ises nwbsessions datadir
sesclk = tic;

% pathpp = 'D:\OpenScopeData\000248\sub-Placeholder\';
% nwb = nwbRead('D:\OpenScopeData\000248\sub-Placeholder\sub-Placeholder.nwb');

nwbfiles = cat(1, dir([datadir nwbsessions{ises} filesep '*.nwb']), dir([datadir nwbsessions{ises} filesep '*' filesep '*.nwb']));

% take filename  with shortest length or filename that does not contain probe
[~, fileind] = min(cellfun(@length, {nwbfiles.name}));
nwbspikefile = fullfile([nwbfiles(fileind).folder filesep nwbfiles(fileind).name]);
% nwbspikefile = string(nwbspikefile);
disp(nwbspikefile)
nwb = nwbRead(nwbspikefile); %, 'ignorecache');

pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];

%ot = nwb.processing.get('optotagging').nwbdatainterface.get('optotagging').data.load();

optocond = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').vectordata.get('condition').data.load();
optostim = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').vectordata.get('stimulus_name').data.load();
optodur = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').vectordata.get('duration').data.load();
optolevel = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').vectordata.get('level').data.load();

%{
o = cellstr([optostim num2str(optodur)]);
v = unique(o);
c = zeros(size(v));
for ii = 1:numel(v)
c(ii) = nnz(strcmp(o,v(ii)));
fprintf('%s %d\n', v{ii},c(ii))
end
%}

%unique(cellstr(horzcat(optocond, optostim)))
opto = struct();
opto.optocond = cellstr(optocond);
opto.optostim = cellstr(optostim);
opto.optodur = optodur;
opto.optolevel = optolevel;

% conditionlist= {'2 ms pulses at 1 Hz', 'a single 10ms pulse', ...
%     '1 second of 5Hz pulse train. Each pulse is 2 ms wide', ...
%     'half-period of a cosine wave', 'Each pulse is 10 ms wide', ...
%     'Each pulse is 6 ms wide', ...
%     '1 second of 40 Hz pulse train. Each pulse is 2 ms wide', ...
%     'a single square pulse', 'a single square pulse', ...
%     'half-period of a cosine wave', ...
%     '1 second of 100 Hz pulse train. Each pulse is 2 ms wide', ...
%     '1 second square pulse: continuously on for 1s'          };
opto.stimlist = {'fast_pulses', 'pulse', 'pulse', 'raised_cosine', ...
    '5 hz pulse train', '40 hz pulse train', 'fast_pulses', 'pulse', ...
    'pulse', 'raised_cosine', '100 hz pulse train', 'Square Pulse'};
opto.durationlist = [1; 1; 0.2; 0.1; 0.05; 0.033; 0.025; 0.02; 0.0167; .0125; 0.01; 1.0];

opto.optotrials = zeros(size(opto.optocond));
for typi = 1:length(stimlist)
    trialsoi = strcmp(opto.optostim, opto.stimlist{typi}) & round(opto.optodur*10^4)==round(opto.durationlist(typi)*10^4) ;
    opto.optotrials(trialsoi) = typi;
end
[v,c]=uniquecnt(opto.optotrials);
disp([v,c])

opto.optostarttime = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').start_time.data.load();
opto.optostoptime = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').stop_time.data.load();

% % below is the order according to my opto script
% condition_list = [0,1,2,3,4,5,6,7,8,9,10,11];
% waveforms_list = {'1Hz', '5Hz', '10Hz', '20Hz', '30Hz', '40Hz', ...
%     '50Hz', '60Hz', '80Hz', '100Hz', 'square1s', 'cosine'};

% % presumed order
opto.presumedcond = {'square1s', '1Hz', '5Hz', '10Hz', '20Hz', '30Hz', '40Hz', ...
    '50Hz', '60Hz', '80Hz', '100Hz', 'cosine'}; % presumably...

% optotimeseries = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').timeseries.data.load();
%%
unit_ids = nwb.units.id.data.load(); % array of unit ids represented within this session
unit_peakch = nwb.units.vectordata.get('peak_channel_id').data.load();
unit_times_data = nwb.units.spike_times.data.load();
unit_times_idx = nwb.units.spike_times_index.data.load();
% unit_waveform = nwb.units.waveform_mean.data.load();
% unit_wfdur = nwb.units.vectordata.get('waveform_duration').data.load();

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

%%
probes = {'A', 'B', 'C', 'D', 'E', 'F'};
optopsthtli = (-500:1500)';
neucnt = 0;
for iprobe = 1:numel(probes)
neuoind = find(floor(unit_peakch/1000)==iprobe-1);
neucnt = neucnt + numel(neuoind);
fprintf('Probe %s: %d\n', probes{iprobe}, numel(neuoind) )

if numel(neuoind)==0
    continue
end

probespiketrain = false(stlen, numel(neuoind));
ststartend = [floor(min(unit_times_data)/Tres)+1 floor(max(unit_times_data)/Tres)+1];
for ii = 1:numel(neuoind)
    ci = neuoind(ii);
    probespiketrain(floor(spiketimes{ci}/Tres)+1, ii) = true;
end


ppfn = sprintf('%spsth_opto_probe%s.mat', pathpp, probes{iprobe});
if exist(ppfn, 'file')
    load(ppfn)
    continue
end


tic
psthtrialinds = floor(opto.optostarttime'/Tres)+1 + optopsthtli;
optopsth = false(length(optopsthtli), length(opto.optostarttime), numel(neuoind));
for ii = 1:numel(neuoind)
    tempST = probespiketrain(:,ii);
    optopsth(:,:,ii) = tempST(psthtrialinds);
end
clear tempST psthtrialinds
toc % takes 6 min for 2500 units


salttrials = opto.optotrials~=find(strcmp(opto.presumedcond, 'cosine'));
saltbasetli = [-floor(0.5/Tres):-1]';
salttesttli = [0:floor(0.009/Tres)]';
probeunits_saltp = NaN(size(neuoind));
probeunits_saltI = NaN(size(neuoind));
for ii = 1:numel(neuoind)
    spt_baseline = squeeze(optopsth(ismember(optopsthtli, saltbasetli), salttrials, ii))';
    spt_test = squeeze(optopsth(ismember(optopsthtli, salttesttli), salttrials, ii))';
    [p I] = salt(spt_baseline,spt_test,Tres);
    probeunits_saltp(ii) = p;
    probeunits_saltI(ii) = I;
end
% appropriate alpha for salt tests: 0.01

save(ppfn, 'neuoind', 'opto', 'Tres', 'optopsthtli', 'optopsth', ...
    'salttrials', 'saltbasetli', 'salttesttli', 'probeunits_saltp', 'probeunits_saltI', '-v7.3')

%{
% % psth averaged across all units, for each opto condition
% relabreord = [opto.relabel(:)', {' '}];
% figure;
for typi = 1:12
    smhalfwin = 0; smwin = smhalfwin*2+1;
    temppsth = convn(optopsth(:, opto.optotrials==typi, :), ones(smwin,1)/smwin, 'valid');
    tempcond = opto.presumedcond{typi};
    %subplot(3,4,typi)
    
    if mod(typi,4)==1
        figure
    end
    subplot(2,2,mod(typi-1,4)+1)
    plot(psthtli(smhalfwin+1:end-smhalfwin), 1000*squeeze(mean(temppsth, [2 3])))
    xlabel('Time (ms)')
    ylabel('Rate (Hz)')
    title(tempcond)
    optohz = strsplit(tempcond, 'Hz');
    optohz = str2double(optohz{1});
    set(gca, 'Xtick', [1000/optohz * (0:optohz)], 'Xgrid', 'on')
end

% % psth averaged across all units, across trials (does baseline decrease across trials?
typi = 3;
smhalfwin = 5; smwin = smhalfwin*2+1; smker = ones(smwin,1)/smwin;
kerwinhalf = 5; kersigma = 2;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
smker = (kergauss/sum(kergauss));
temppsth = convn(optopsth(:, opto.optotrials==typi, :), smker, 'same');
cm = gray(nnz(opto.optotrials==typi)); cm(:,2) = 0; cm(:,3) = 0;
figure
p = plot(psthtli, 1000*squeeze(mean( temppsth, 3)));%, 'Color', cm)
for ii = 1:numel(p)
set(p(ii), 'Color', cm(ii,:))
end
title(opto.presumedcond{typi})

typi = 3;
smhalfwin = 5; smwin = smhalfwin*2+1; smker = ones(smwin,1)/smwin;
kerwinhalf = 5; kersigma = 2;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
smker = (kergauss/sum(kergauss));
temppsth = convn(optopsth(:, opto.optotrials==typi, :), smker, 'same');
cm = gray(nnz(opto.optotrials==typi)); cm(:,2) = 0; cm(:,3) = 0;
tempsalt = (probeunits_saltp<0.01);
figure
p = plot(psthtli, 1000*squeeze(mean( temppsth(:,:,tempsalt), 3)));%, 'Color', cm)
for ii = 1:numel(p)
set(p(ii), 'Color', cm(ii,:))
end
title(opto.presumedcond{typi})

[sv,si]=sort(probeunits_saltp);
for ii = 1:10
    disp([si(ii) nnz(optopsth(:, opto.optotrials==typi, si(ii)))])
end
exunit = 145; typi = 3;
exunit = si(34); typi = 3;
figure; imagesc(squeeze(optopsth(:, opto.optotrials==typi, exunit))')
colormap(flipud(gray))

% CHECK THAT SALT WORKED REASONABLY WELL
[sv,si]=sort(probeunits_saltp);
optopsthavg = squeeze(mean(optopsth(:,salttrials,si), 2));
figure; imagesc(squeeze(optopsthavg)')
colormap(flipud(gray))
caxis([0 0.1])
%}

end


end
