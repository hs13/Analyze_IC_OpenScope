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


if ismac
    addpath(genpath('/Users/hyeyoung/Documents/CODE/Analyze_OpenScope/matnwb'))
    addpath('/Users/hyeyoung/Documents/CODE/helperfunctions')
    datadir = '/Users/hyeyoung/Documents/DATA/OpenScopeData/000248/';
elseif ispc
    rmpath(genpath('C:\Users\USER\GitHub\Analyze_IC_OpenScope\matnwb'))
    addpath(genpath('S:\OpenScopeData\matnwb_HSLabDesktop'))
    addpath('C:\Users\USER\GitHub\helperfunctions')
    datadir = 'S:\OpenScopeData\000248\';
end

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
    %disp(nwbspikefile)
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
    opto.genotype = nwb.general_subject.genotype;
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
    for typi = 1:length(opto.stimlist)
        trialsoi = strcmp(opto.optostim, opto.stimlist{typi}) & round(opto.optodur*10^4)==round(opto.durationlist(typi)*10^4) ;
        opto.optotrials(trialsoi) = typi;
    end
    [v,c]=uniquecnt(opto.optotrials);
    %disp([v,c])

    opto.optostarttime = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').start_time.data.load();
    opto.optostoptime = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').stop_time.data.load();

    % % below is the order according to my opto script
    % condition_list = [0,1,2,3,4,5,6,7,8,9,10,11];
    % waveforms = [data_1Hz_10ms, data_1Hz, data_5Hz, data_10Hz, data_20Hz, data_30Hz, data_40Hz, data_50Hz, data_60Hz, data_80Hz, data_square1s, data_cosine];

    % % presumed order
    % opto.presumedcond = {'square1s', '1Hz', '5Hz', '10Hz', '20Hz', '30Hz', '40Hz', ...
    %     '50Hz', '60Hz', '80Hz', '100Hz', 'cosine'}; % presumably...

    opto.actualcond = {'1Hz_10ms', '1Hz_2ms', '5Hz_2ms', '10Hz_2ms', '20Hz_2ms', ...
        '30Hz_2ms', '40Hz_2ms', '50Hz_2ms', '60Hz_2ms', '80Hz_2ms', 'square1s', 'cosine'};
    fprintf('%s %s\n', nwbsessions{ises}, opto.genotype)

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
        % if exist(ppfn, 'file')
        %     load(ppfn)
        %     continue
        % end


        % tic
        psthtrialinds = floor(opto.optostarttime'/Tres)+1 + optopsthtli;
        optopsth = false(length(optopsthtli), length(opto.optostarttime), numel(neuoind));
        for ii = 1:numel(neuoind)
            tempST = probespiketrain(:,ii);
            optopsth(:,:,ii) = tempST(psthtrialinds);
        end
        clear tempST psthtrialinds
        % toc % takes ~5sec per probe


        salttrials = opto.optotrials~=find(strcmp(opto.actualcond, 'cosine'));
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

    fprintf('%.2f%% (%d/%d)\n', 100*mean(probeunits_saltp<0.01), nnz(probeunits_saltp<0.01), length(probeunits_saltp))

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
    plot(optopsthtli(smhalfwin+1:end-smhalfwin), 1000*squeeze(mean(temppsth, [2 3])))
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
p = plot(optopsthtli, 1000*squeeze(mean( temppsth, 3)));%, 'Color', cm)
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
p = plot(optopsthtli, 1000*squeeze(mean( temppsth(:,:,tempsalt), 3)));%, 'Color', cm)
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

%% plot each session
datadir = 'S:\OpenScopeData\000248\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));

figure
for ises = 1:numel(nwbsessions)
    clearvars('neuoind', 'opto', 'Tres', 'optopsthtli', 'optopsth', ...
        'salttrials', 'saltbasetli', 'salttesttli', 'probeunits_saltp', 'probeunits_saltI')
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    ppfn = [pathpp 'psth_opto_probeC.mat'];
    load(ppfn)

    %opto.genotype = nwb.general_subject.genotype;
    % opto.actualcond = {'1Hz_10ms', '1Hz_2ms', '5Hz_2ms', '10Hz_2ms', '20Hz_2ms', ...
    %     '30Hz_2ms', '40Hz_2ms', '50Hz_2ms', '60Hz_2ms', '80Hz_2ms', 'square1s', 'cosine'};
    % opto = rmfield(opto, 'presumedcond');

    disp(opto.genotype)
    fprintf('%.2f%% (%d/%d)\n', 100*mean(probeunits_saltp<0.01), nnz(probeunits_saltp<0.01), length(probeunits_saltp))

    for typi = 1:12
        subplot(12,numel(nwbsessions),(typi-1)*numel(nwbsessions)+ises)
        smhalfwin = 0; smwin = smhalfwin*2+1;
        temppsth = convn(optopsth(:, opto.optotrials==typi, :), ones(smwin,1)/smwin, 'valid');
        tempcond = opto.actualcond{typi};

        plot(optopsthtli(smhalfwin+1:end-smhalfwin), 1000*squeeze(mean(temppsth, [2 3])))
        if typi==12
        xlabel('Time (ms)')
        end
        if ises==1
        ylabel('Rate (Hz)')
        end
        if typi==1
        title({nwbsessions{ises}, tempcond}, 'interpreter', 'none', 'FontSize', 7, 'FontWeight','normal')
        else
        title(tempcond, 'interpreter', 'none', 'FontSize', 7, 'FontWeight','normal')
        end
        xlim([-250 1250])
        optohz = strsplit(tempcond, 'Hz');
        optohz = str2double(optohz{1});
        xt = [1000/optohz * (0:optohz)];
        if isnan(optohz)
            xt = [0 1000];
        end
        xtl = cell(size(xt));
        xtl{1} = 0;
        xtl{end} = 1;
        set(gca, 'FontSize', 7, 'Xtick', xt, 'XTickLabel', xtl, 'XTickLabelRotation', 0, 'Xgrid', 'on')
    end
end

%% added 230829: RESUME EDITING
% Siegle et al's single unit filter criteria
% isi_violations < 0.5 & amplitude_cutoff < 0.1 & presence_ratio > 0.9
% https://allensdk.readthedocs.io/en/latest/_static/examples/nb/visual_behavior_neuropixels_quality_metrics.html1


% ISI violations: This metric searches for refractory period violations that indicate a unit contains spikes from multiple neurons. 
% The ISI violations metric represents the relative firing rate of contaminating spikes. 
% It is calculated by counting the number of violations of less than 1.5 ms, 
% dividing by the amount of time for potential violations surrounding each spike, 
% and normalizing by the overall spike rate. It is always positive (or 0), but has no upper bound. See ref. 72 for more details.
% Amplitude cutoff: This metric provides an approximation of a unit%s false negative rate. 
% First, a histogram of spike amplitudes is created, and the height of the histogram at the minimum amplitude is extracted. 
% The percentage of spikes above the equivalent amplitude on the opposite side of the histogram peak is then calculated. 
% If the minimum amplitude is equivalent to the histogram peak, the amplitude cutoff is set to 0.5 
% (indicating a high likelihood that more than 50% of spikes are missing). 
% This metric assumes a symmetrical distribution of amplitudes and no drift, so it will not necessarily reflect the true false negative rate.
% Presence ratio: The session was divided into 100 equal-sized blocks; the presence ratio is defined as the fraction of blocks that include one or more spikes from a particular unit. Units with a low presence ratio are likely to have drifted out of the recording, or could not be tracked by Kilosort2 for the duration of the experiment.


addpath(genpath('S:\OpenScopeData\matnwb_HSLabDesktop'))
nwbsessions = {'sub_1171903433','sub_1172968426','sub_1172969394','sub_1174569641', ...
    'sub_1175512783','sub_1176214862','sub_1177693342','sub_1181314060', ...
    'sub_1181585608','sub_1182593231','sub_1183369803','sub_1186544726', ...
    'sub_1189891322','sub_1194090570'};

for ises = 1:numel(nwbsessions)
    clearvars -except nwbsessions ises
    mousedate = nwbsessions{ises};
    datadir = 'S:/OpenScopeData/000248/';
    pathpp = strcat('S:/OpenScopeData/000248/postprocessed/', mousedate, '/');
    disp(mousedate)

    nwbfiles = cat(1, dir([datadir nwbsessions{ises} filesep '*.nwb']), dir([datadir nwbsessions{ises} filesep '*' filesep '*.nwb']));
    % take filename  with shortest length or filename that does not contain probe
    [~, fileind] = min(cellfun(@length, {nwbfiles.name}));
    nwbspikefile = fullfile([nwbfiles(fileind).folder filesep nwbfiles(fileind).name]);
    nwb = nwbRead(nwbspikefile);

    % load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    % load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur'
    % elecid = electrode_id+1;
    % revmapelecid = NaN(max(elecid),1);
    % revmapelecid(elecid) = 1:numel(elecid);
    % neuallloc = electrode_location(revmapelecid(unit_peakch+1));

    unit_ids = nwb.units.id.data.load(); % array of unit ids represented within this session
    unit_peakch = nwb.units.vectordata.get('peak_channel_id').data.load();
    % unit_waveform = nwb.units.waveform_mean.data.load();
    unit_wfdur = nwb.units.vectordata.get('waveform_duration').data.load();
    unit_isi_violations = nwb.units.vectordata.get('isi_violations').data.load();
    unit_amplitude = nwb.units.vectordata.get('amplitude').data.load();
    unit_amplitude_cutoff = nwb.units.vectordata.get('amplitude_cutoff').data.load();
    unit_presence_ratio = nwb.units.vectordata.get('presence_ratio').data.load();

    save([pathpp 'qc_units.mat'], 'unit_wfdur', 'unit_isi_violations', 'unit_amplitude', 'unit_amplitude_cutoff', 'unit_presence_ratio') %'unit_times_data',

    figure('Position', [100 100 1200 400])
    subplot(1,3,1)
    hold all
    histogram(log10(unit_isi_violations+0.00001))
    yl = ylim;
    plot(log10([0.5 0.5]), yl, 'r--')
    title(sprintf('isi_violations<0.5 : %.0f%%', 100*mean(unit_isi_violations<0.5)), 'interpreter', 'none')
    subplot(1,3,2)
    hold all
    histogram(unit_amplitude_cutoff, 'NumBins', 20)
    yl = ylim;
    plot([0.5 0.5], yl, 'r--')
    title(sprintf('amplitude_cutoff<0.1 : %.0f%%', 100*mean(unit_amplitude_cutoff<0.1)), 'interpreter', 'none')
    subplot(1,3,3)
    hold all
    histogram(unit_presence_ratio)
    yl = ylim;
    plot([0.9 0.9], yl, 'r--')
    title(sprintf('presence_ratio>0.9 : %.0f%%', 100*mean(unit_presence_ratio>0.9)), 'interpreter', 'none')

    disp('unit_isi_violations<0.5 & unit_amplitude_cutoff<0.1 & unit_presence_ratio>0.9')
    disp(mean(unit_isi_violations<0.5 & unit_amplitude_cutoff<0.1 & unit_presence_ratio>0.9)) % criteria in Siegle
    disp('unit_isi_violations<0.5 & unit_amplitude_cutoff<0.5 & unit_presence_ratio>0.9')
    disp(mean(unit_isi_violations<0.5 & unit_amplitude_cutoff<0.5 & unit_presence_ratio>0.9)) % relax amplitude_cutoff
    disp('mean(unit_isi_violations<1 & unit_amplitude_cutoff<0.5 & unit_presence_ratio>0.9)')
    disp(mean(unit_isi_violations<1 & unit_amplitude_cutoff<0.5 & unit_presence_ratio>0.9)) % relax amplitude_cutoff and isi_violoations

end

isi5amp1pres9 = zeros(numel(nwbsessions),1);
isi5amp5pres9 = zeros(numel(nwbsessions),1);
isi1amp5pres9 = zeros(numel(nwbsessions),1);
isi5pres9 = zeros(numel(nwbsessions),1);
for ises = 1:numel(nwbsessions)
    mousedate = nwbsessions{ises};
    pathpp = strcat('S:/OpenScopeData/000248/postprocessed/', mousedate, '/');

    load([pathpp 'qc_units.mat'], 'unit_wfdur', 'unit_isi_violations', 'unit_amplitude', 'unit_amplitude_cutoff', 'unit_presence_ratio') %'unit_times_data',
    isi5amp1pres9(ises) = mean(unit_isi_violations<0.5 & unit_amplitude_cutoff<0.1 & unit_presence_ratio>0.9);
    isi5amp5pres9(ises) = mean(unit_isi_violations<0.5 & unit_amplitude_cutoff<0.5 & unit_presence_ratio>0.9);
    isi1amp5pres9(ises) = mean(unit_isi_violations<1 & unit_amplitude_cutoff<0.5 & unit_presence_ratio>0.9);
    isi5pres9(ises) = mean(unit_isi_violations<0.5 & unit_presence_ratio>0.9);

    % disp(mousedate)
    % disp('unit_isi_violations<0.5 & unit_amplitude_cutoff<0.1 & unit_presence_ratio>0.9')
    % disp(mean(unit_isi_violations<0.5 & unit_amplitude_cutoff<0.1 & unit_presence_ratio>0.9)) % criteria in Siegle
    % disp('unit_isi_violations<0.5 & unit_amplitude_cutoff<0.5 & unit_presence_ratio>0.9')
    % disp(mean(unit_isi_violations<0.5 & unit_amplitude_cutoff<0.5 & unit_presence_ratio>0.9)) % relax amplitude_cutoff
    % disp('mean(unit_isi_violations<1 & unit_amplitude_cutoff<0.5 & unit_presence_ratio>0.9)')
    % disp(mean(unit_isi_violations<1 & unit_amplitude_cutoff<0.5 & unit_presence_ratio>0.9)) % relax amplitude_cutoff and isi_violoations
end

figure; hold all
plot(1:4, [isi5amp1pres9 isi5amp5pres9 isi1amp5pres9 isi5pres9], 'o-')
errorbar(1:4, mean([isi5amp1pres9 isi5amp5pres9 isi1amp5pres9 isi5pres9],1), std([isi5amp1pres9 isi5amp5pres9 isi1amp5pres9 isi5pres9],0,1)/sqrt(numel(nwbsessions)), 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 2)
