addpath(genpath('H:\CODE\helperfunctions'))
addpath(genpath('H:\CODE\Analyze_OpenScope'))

% kerwinhalf = 2; kersigma = 1;
% kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
% kergauss = (kergauss/sum(kergauss));

% A-AM, B-PM, C-V1, D-LM, E-AL, F-RL
probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
visind = [6 5 1 2 4 3];

datadir = 'D:\OpenScopeData\000248\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name}; 
nwbsessions = nwbsessions(contains(nwbsessions, 'sub'));

Nsessions = numel(nwbsessions)-1;
probeneuronsagg = cell(size(probes));
neuctxagg = cell(size(probes));
neulocagg = cell(size(probes));
sesneuoiagg = cell(size(probes));
psthsizeCG_ctxagg = cell(size(probes));
psthsizeIG_ctxagg = cell(size(probes));
psthsizeCG4_ctxagg = cell(size(probes));
psthsizeCG64_ctxagg = cell(size(probes));
psthsizeFG_ctxagg = cell(size(probes));
psthRFCG_ctxagg = cell(size(probes));
psthRFIG_ctxagg = cell(size(probes));

% ICsigagg = cell(size(probes));
% oriparamsagg = cell(size(probes));
RFCIfields = {'Rrfclassic','Rrfinverse','RFindclassic','RFindinverse', ...
    'Pkw_rfclassic','Pkw_rfinverse','pRrfclassic','pRrfinverse','pRFclassic','pRFinverse'};
RFCIagg = struct();
for iprobe = 1:numel(probes)
for f= 1:numel(RFCIfields)
RFCIagg(iprobe).(RFCIfields{f}) = [];
end
end
sizeCIfields = {'Rsizeclassic','Rsizeinverse','sizeindclassic','sizeindinverse','Pkw_sizeclassic','Pkw_sizeinverse'};
sizeCIagg = struct();
for iprobe = 1:numel(probes)
for f= 1:numel(sizeCIfields)
sizeCIagg(iprobe).(sizeCIfields{f}) = [];
end
end

for ises = 1:4%Nsessions
    fprintf('Session %d/%d %s\n', ises, Nsessions, nwbsessions{ises} )
pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur'

whichneuctx = 1; % 1: correct definition, 2: >=230, 3: topmost-250um
elecid = electrode_id+1;
revmapelecid = NaN(max(elecid),1);
revmapelecid(elecid) = 1:numel(elecid);

% psthall = struct();
% ICsigall = struct();
% RFCIall = struct();
% sizeCIall = struct();
% oriparamsall = struct();
for iprobe = 1:numel(probes)
%     tic
    load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}))
    % 'neuoind', 'vis', 'Tres', 'psthtli', 'psth'
    load(sprintf('%svisresponses_probe%s.mat', pathpp, probes{iprobe}))
    % 'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'ICsig', 'RFCI', 'sizeCI', 'oriparams'
    
    if ~isequal(unique(floor(unit_peakch(neuoind)/1000)), iprobe-1)
        error('check neuoind')
    end
    
    % check whether CCF registration is correct
    probelocs = electrode_location(ismember(electrode_id, unit_peakch(neuoind)));
    
    neuloc = electrode_location(revmapelecid(unit_peakch(neuoind)+1));
    if ~isequal(unique(probelocs), unique(neuloc))
        disp(unique(neuloc)')
        error('check neuloc')
    end
    
    switch whichneuctx
        case 1
            neuctx = contains(neuloc, 'VIS');
        case 2
            % assume that >=230 is neocortex
            neuctx = mod(unit_peakch(neuoind), 1000)>=230;
            warning('for now, assuming that electrode_localid>=230 is neocortex (until CCF registration is corrected)')
            % 2 electrodes per depth, vertical spacing is 20um
            % typical electrode span in cortex 120 : (120/2)*20 = ~1200um
        case 3
            % to approximately get layer 2/3 cells, take units within 250um of the topmost unit
            neuctx = mod(unit_peakch(neuoind), 1000)>=max(mod(unit_peakch(neuoind), 1000))-25;
            warning('for now, to approximately get layer 2/3 cells, take units within 250um of the topmost unit')
    end

    fprintf('Probe %s Area %s: %d/%d\n', probes{iprobe}, visareas{iprobe}, nnz(neuctx), numel(neuoind) )
    disp(unique(probelocs)')
    
    probeneuronsagg{iprobe} = cat(1, probeneuronsagg{iprobe}, neuoind);
    neuctxagg{iprobe} = cat(1, neuctxagg{iprobe}, neuctx);
    neulocagg{iprobe} = cat(1, neulocagg{iprobe}, neuloc);
    
    % aggregate cortex neurons
    neuoi = neuctx;
    sesneuoiagg{iprobe} = cat(1, sesneuoiagg{iprobe}, ises*ones(nnz(neuoi),1));
%     psthagg = cat(1, psthagg, psth);
%     ICsigagg = cat(1, ICsigagg, ICsig);
%     RFCIagg = cat(1, RFCIagg, RFCI);
%     sizeCIagg = cat(1, sizeCIagg, sizeCI);
%     oriparamsagg = cat(1, oriparamsagg, oriparams);

% RFCIfields = fieldnames(RFCI);
for f= 1:numel(RFCIfields)
RFCIagg(iprobe).(RFCIfields{f}) = cat(1, RFCIagg(iprobe).(RFCIfields{f}), RFCI.(RFCIfields{f})(neuoi,:) );
end

for f= 1:numel(sizeCIfields)
sizeCIagg(iprobe).(sizeCIfields{f}) = cat(1, sizeCIagg(iprobe).(sizeCIfields{f}), sizeCI.(sizeCIfields{f})(neuoi,:) );
end

% size vector [0, 4, 8, 16, 32, 64 ]
for icond = 1:5
    switch icond
        case 1
            trialsoi = floor(vis.sizeCI_presentations.trialorder/1000)==4;
            %tempcol = [0 0 0];
        case 2
            trialsoi = floor(vis.sizeCI_presentations.trialorder/1000)==14;
            %tempcol = [0 1 1];
        case 3
            trialsoi = floor(vis.sizeCI_presentations.trialorder/1000)==2;
        case 4
            trialsoi = floor(vis.sizeCI_presentations.trialorder/1000)==6;
        case 5
            trialsoi = floor(vis.sizeCI_presentations.trialorder/1000)==11;
    end
    temp = 1000*psth.sizeCI_presentations(:,trialsoi,neuoi);
    % plot(psthtli, smooth(squeeze(mean(temp, [2,3])), 5), '-')
    %     temp = convn(temp, kergauss, 'same');
    %plot(psthtli, squeeze(mean(temp, [2,3])), '-', 'Color', tempcol)
    
    switch icond
        case 1
            psthsizeCG_ctxagg{iprobe} = cat(2, psthsizeCG_ctxagg{iprobe}, squeeze(mean(temp, 2)) );
        case 2
            psthsizeIG_ctxagg{iprobe} = cat(2, psthsizeIG_ctxagg{iprobe}, squeeze(mean(temp, 2)) );
        case 3
            psthsizeCG4_ctxagg{iprobe} = cat(2, psthsizeCG4_ctxagg{iprobe}, squeeze(mean(temp, 2)) );
        case 4
            psthsizeCG64_ctxagg{iprobe} = cat(2, psthsizeCG64_ctxagg{iprobe}, squeeze(mean(temp, 2)) );
        case 5
            psthsizeFG_ctxagg{iprobe} = cat(2, psthsizeFG_ctxagg{iprobe}, squeeze(mean(temp, 2)) );
    end
end
  
for icond = 1:2
        switch icond
            case 1
                trialsoi = vis.RFCI_presentations.trialorder==1001;
            tempcol = [0 0 0];
            case 2
                trialsoi = vis.RFCI_presentations.trialorder==11001;
            tempcol = [0 1 1];
        end
    temp = 1000*psth.RFCI_presentations(:,trialsoi,neuoi);
    % plot(psthtli, smooth(squeeze(mean(temp, [2,3])), 5), '-')
    %     temp = convn(temp, kergauss, 'same');
    %plot(psthtli, squeeze(mean(temp, [2,3])), '-', 'Color', tempcol)
    
    switch icond
        case 1
            psthRFCG_ctxagg{iprobe} = cat(2, psthRFCG_ctxagg{iprobe}, squeeze(mean(temp, 2)) );
        case 2
            psthRFIG_ctxagg{iprobe} = cat(2, psthRFIG_ctxagg{iprobe}, squeeze(mean(temp, 2)) );
    end
end  
% toc
end

end

save('G:\My Drive\DATA\OpenScope\psth_CGIG.mat', 'probes', 'visareas', 'visind', 'nwbsessions', ...
    'probeneuronsagg', 'neuctxagg', 'neulocagg', 'sesneuoiagg', 'RFCIagg', 'sizeCIagg', ...
    'psthtli', 'psthsizeCG_ctxagg', 'psthsizeCG4_ctxagg', 'psthsizeCG64_ctxagg', ...
    'psthsizeIG_ctxagg', 'psthsizeFG_ctxagg', 'psthRFCG_ctxagg', 'psthRFIG_ctxagg', '-v7.3')

%% for grant
%sizeCI: each stim is 0.25s, inter-trial interval is 0.5s, drifting grating

addpath(genpath('H:\CODE\helperfunctions'))

kerwinhalf = 2; kersigma = 1;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));

probes2p = [3 4 5 2];
fs=14;
figure('Position', [100 100 1200 300])
% annotation('textbox', [0.1 0.9 0.9 0.1], 'string', 'Neuropixels: center-CRF neurons', 'edgecolor', 'none', 'fontsize', fs)
for isp= 1:length(probes2p)
    iprobe = probes2p(isp);
subplot(1,4, isp)
neuoi = sesneuoiagg{iprobe}<=4 & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
temppsthCG = convn(psthsizeCG_ctxagg{iprobe}(:,neuoi), kergauss, 'same');
hold all
shadedErrorBar(psthtli/1000, squeeze(mean(temppsthCG,2)), squeeze(std(temppsthCG,0,2)/sqrt(nnz(neuoi))), {'k-', 'LineWidth', 1}, 1)
xlim([-.1 .200])
ylim([0 25])
set(gca, 'FontSize', fs)
    xlabel('Time (s)')
    ylabel('Rate (Hz)')
    title(visareas{iprobe})
%     title(sprintf('%s N=%d', visareas{iprobe}, nnz(neuoi)))
end

figure('Position', [100 500 1800 300])
% annotation('textbox', [0.1 0.9 0.9 0.1], 'string', 'Neuropixels: center-CRF neurons', 'edgecolor', 'none', 'fontsize', fs)
for iprobe= 1:6
subplot(1,6, visind(iprobe))
neuoi = sesneuoiagg{iprobe}<=4 & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
temppsthCG = convn(psthsizeCG_ctxagg{iprobe}(:,neuoi), kergauss, 'same');
hold all
shadedErrorBar(psthtli/1000, squeeze(mean(temppsthCG,2)), squeeze(std(temppsthCG,0,2)/sqrt(nnz(neuoi))), {'k-', 'LineWidth', 1}, 1)
xlim([-.1 .200])
ylim([0 25])
set(gca, 'FontSize', fs)
    xlabel('Time (s)')
    ylabel('Rate (Hz)')
    title(sprintf('%s N=%d', visareas{iprobe}, nnz(neuoi)))
end

%%
fs=14;
figure('Position', [100 100 1800 300])
% annotation('textbox', [0.1 0.9 0.9 0.1], 'string', 'Neuropixels: center-CRF neurons', 'edgecolor', 'none', 'fontsize', fs)
for iprobe= 1:6
subplot(1,6, visind(iprobe))
neuoi = sesneuoiagg{iprobe}<=4 & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
temppsthCG = convn(psthsizeCG_ctxagg{iprobe}(:,neuoi), kergauss, 'same');
temppsthIG = convn(psthsizeIG_ctxagg{iprobe}(:,neuoi), kergauss, 'same');
hold all
shadedErrorBar(psthtli/1000, squeeze(mean(temppsthCG,2)), squeeze(std(temppsthCG,0,2)/sqrt(nnz(neuoi))), {'k-', 'LineWidth', 1}, 1)
shadedErrorBar(psthtli/1000, squeeze(mean(temppsthIG,2)), squeeze(std(temppsthIG,0,2)/sqrt(nnz(neuoi))), {'c-', 'LineWidth', 1}, 1)
xlim([-.1 .200])
ylim([0 25])
set(gca, 'FontSize', fs)
    xlabel('Time (s)')
    ylabel('Rate (Hz)')
    title(sprintf('%s N=%d', visareas{iprobe}, nnz(neuoi)))
end

%%
kerwinhalf = 20; kersigma = 10;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));

fs=14;
figure('Position', [100 100 1800 300])
% annotation('textbox', [0.1 0.9 0.9 0.1], 'string', 'Neuropixels: center-CRF neurons', 'edgecolor', 'none', 'fontsize', fs)
for iprobe= 1:6
subplot(1,6, visind(iprobe))
neuoi = sesneuoiagg{iprobe}<=4;
neuoi = sesneuoiagg{iprobe}<=4 & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
%neuoi = sesneuoiagg{iprobe}<=4 & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).RFindinverse==1;
temppsthCG = convn(psthRFCG_ctxagg{iprobe}(:,neuoi), kergauss, 'same');
temppsthIG = convn(psthRFIG_ctxagg{iprobe}(:,neuoi), kergauss, 'same');
hold all
shadedErrorBar(psthtli/1000, squeeze(mean(temppsthCG,2)), squeeze(std(temppsthCG,0,2)/sqrt(nnz(neuoi))), {'k-', 'LineWidth', 0.5}, 1)
shadedErrorBar(psthtli/1000, squeeze(mean(temppsthIG,2)), squeeze(std(temppsthIG,0,2)/sqrt(nnz(neuoi))), {'c-', 'LineWidth', 0.5}, 1)
xlim([-.5 1])
ylim([0 20])
set(gca, 'FontSize', fs)
    xlabel('Time (s)')
    ylabel('Rate (Hz)')
    title(sprintf('%s N=%d', visareas{iprobe}, nnz(neuoi)))
end

%%

fs=14;
figure('Position', [100 100 1800 300])
% annotation('textbox', [0.1 0.9 0.9 0.1], 'string', 'Neuropixels: center-CRF neurons', 'edgecolor', 'none', 'fontsize', fs)
for iprobe= 1:6
subplot(1,6, visind(iprobe))
neuoi = sesneuoiagg{iprobe}<=4 & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
temppsthCG = convn(psthsizeCG_ctxagg{iprobe}(:,neuoi), kergauss, 'same');
temppsthCG64 = convn(psthsizeCG64_ctxagg{iprobe}(:,neuoi), kergauss, 'same');
hold all
shadedErrorBar(psthtli/1000, squeeze(mean(temppsthCG,2)), squeeze(std(temppsthCG,0,2)/sqrt(nnz(neuoi))), {'k-', 'LineWidth', 1}, 1)
shadedErrorBar(psthtli/1000, squeeze(mean(temppsthCG64,2)), squeeze(std(temppsthCG64,0,2)/sqrt(nnz(neuoi))), {'c-', 'LineWidth', 1}, 1)
xlim([-.1 .2])
xlim([-.5 1])
ylim([0 25])
set(gca, 'FontSize', fs)
    xlabel('Time (s)')
    ylabel('Rate (Hz)')
    title(sprintf('%s N=%d', visareas{iprobe}, nnz(neuoi)))
%     title(sprintf('%s N=%d', visareas{iprobe}, nnz(neuoi)))
end

%%
szvec = [0, 4, 8, 16, 32, 64];

figure('Position', [100 100 1800 300])
% annotation('textbox', [0.1 0.9 0.9 0.1], 'string', 'Neuropixels: center-CRF neurons', 'edgecolor', 'none', 'fontsize', fs)
for iprobe= 1:6
subplot(1,6, visind(iprobe))
hold all
neuoi = sesneuoiagg{iprobe}<=4 & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
errorbar(1:length(szvec), mean(sizeCIagg(iprobe).Rsizeclassic(neuoi,:), 1), std(sizeCIagg(iprobe).Rsizeclassic(neuoi,:), 0,1)/sqrt(nnz(neuoi)), 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 2);
errorbar(1:length(szvec), mean(sizeCIagg(iprobe).Rsizeinverse(neuoi,:), 1), std(sizeCIagg(iprobe).Rsizeinverse(neuoi,:), 0,1)/sqrt(nnz(neuoi)), 'co-', 'MarkerFaceColor', 'c', 'LineWidth', 2);
set(gca, 'FontSize', fs, 'Xtick', 1:length(szvec), 'xticklabel', szvec)
xlim([1 length(szvec)])
    xlabel('Size (deg)')
    ylabel('Rate (Hz)')
    title(sprintf('%s N=%d', visareas{iprobe}, nnz(neuoi)))
end
