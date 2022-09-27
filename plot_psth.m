addpath(genpath('H:\CODE\helperfunctions'))
addpath(genpath('H:\CODE\Analyze_OpenScope'))
datadir = 'D:\OpenScopeData\000248\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name}; 
nwbsessions = nwbsessions(contains(nwbsessions, 'sub'));

ises = 4;

pathpp = [datadir nwbsessions{ises} filesep];
load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur'


%% aggregate all probes
% A-AM, B-PM, C-V1, D-LM, E-AL, F-RL
probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
visind = [6 5 1 2 4 3];

% 1000's identify which probe, mod(unit_peakch,1000) identifies electrode position (lower values are deeper) 
% assume that >=230 is neocortex
figure;
for iprobe = 1:numel(probes)
    subplot(2,3,iprobe); hold all
    histogram(unit_peakch(floor(unit_peakch/1000)==iprobe-1), 'BinWidth', 20)
    yl = ylim;
    plot((iprobe-1)*1000+[230 230], yl, 'r--')
    ylim(yl)
end

probeneurons = cell(size(probes));
% psthall = struct();
% ICsigall = struct();
% RFCIall = struct();
% sizeCIall = struct();
% oriparamsall = struct();
for iprobe = 1:numel(probes)
    tic
    load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}))
    % 'neuoind', 'vis', 'Tres', 'psthtli', 'psth'
    load(sprintf('%svisresponses_probe%s.mat', pathpp, probes{iprobe}))
    % 'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'ICsig', 'RFCI', 'sizeCI', 'oriparams'
    
    if ~isequal(unique(floor(unit_peakch(neuoind)/1000)), iprobe-1)
        error('check neuoind')
    end
    fprintf('Probe %s Area %s: %d\n', probes{iprobe}, visareas{iprobe}, numel(neuoind) )
    
    probeneurons{iprobe} = neuoind;
    if iprobe==1
        psthall = psth;
        ICsigall = ICsig;
        RFCIall = RFCI;
        sizeCIall = sizeCI;
        oriparamsall = oriparams;
    else
        psthall = cat(1, psthall, psth);
        ICsigall = cat(1, ICsigall, ICsig);
        RFCIall = cat(1, RFCIall, RFCI);
        sizeCIall = cat(1, sizeCIall, sizeCI);
        oriparamsall = cat(1, oriparamsall, oriparams);
    end
    toc
end

%% check whether CCF registration is correct
figure;
hold all
for iprobe = 1:6
    electrodeoi = electrode_probeid==iprobe-1;
    histogram(electrode_localid(contains(electrode_location(electrodeoi), 'VIS')))
end
legend(probes)
xlabel('electrode localid')

%% neuctxall
whichneuctx = 1; % 1: correct definition, 2: >=230, 3: topmost-250um
elecid = electrode_id+1;
revmapelecid = NaN(max(elecid),1);
revmapelecid(elecid) = 1:numel(elecid);

neuctxall = cell(size(probes));
neulocall = cell(size(probes));
for iprobe = 1:numel(probes)
    neuoind = probeneurons{iprobe};
    
    fprintf('Probe %s\n', probes{iprobe})
    % check whether CCF registration is correct
    probelocs = electrode_location(ismember(electrode_id, unit_peakch(neuoind)));
    disp(unique(probelocs)')
    
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

    neuctxall{iprobe} = neuctx;
    neulocall{iprobe} = neuloc;
end

figure;
for iprobe = 1:6
    neuoind = probeneurons{iprobe};
    neuctx = neuctxall{iprobe};
    neuloc = neulocall{iprobe};
    subplot(2,3,iprobe)
    histogram2(contains(neuloc, 'VIS'), mod(unit_peakch(neuoind), 1000)>=230, 'DisplayStyle', 'tile'); colorbar
    xlabel('contains(neuloc, ''VIS'')')
    ylabel('mod(unit_peakch(neuoind), 1000)>=230', 'interpreter', 'none')
    title(sprintf('Probe %s Area %s N=%d\n', probes{iprobe}, visareas{iprobe}, numel(neuoind) ))
end

%%
whichprobe =3;

load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{whichprobe}))
% 'neuoind', 'vis', 'Tres', 'psthtli', 'psth'
load(sprintf('%svisresponses_probe%s.mat', pathpp, probes{whichprobe}))
% 'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'ICsig', 'RFCI', 'sizeCI', 'oriparams'

visblocks = fieldnames(vis)';
neuctx = neuctxall{whichprobe};

%% response to full screen gratings across areas (all neurons)
%sizeCI: each stim is 0.25s, inter-trial interval is 0.5s, drifting grating
% baseline period: -100 to 0ms (alternatively -200 to 0ms would work)
[~,probehierarchy]=sort(visind);
figure('Position', [100 100 1000 400])
% annotation('textbox', [0.05 0.9 0.8 0.1], 'string', [nwbsessions{ises} ' All Neurons'], 'FontSize', 14, 'EdgeColor', 'none')
subplot(1,2,1)
hold all
for a = 1:numel(probes)
    iprobe = probehierarchy(a);
%     a = visind(iprobe);
%     whicharea = visareas{iprobe};
    neuctx = neuctxall{iprobe};
    
    neuoi = neuctx;
    % neuoi = neuctx & RFCIall(iprobe).Pkw_rfclassic<0.05;
    
    trialsoi = floor(vis.sizeCI_presentations.trialorder/1000)==11;
    temp = 1000*psthall(iprobe).sizeCI_presentations(:,trialsoi,neuoi);
    % plot(psthtli, smooth(squeeze(mean(temp, [2,3])), 5), '-')
    temp = convn(temp, kergauss, 'same');
    plot(psthtli, squeeze(mean(temp, [2,3])))
end
legend(visareas{probehierarchy})
xlim([-50 200])
% xlim([-500 750])
set(gca, 'XGrid', 'on')
title([nwbsessions{ises} ' All Neurons Full Screen Gratings'])
% title(sprintf('%s %s neurons N=%d/%d\nsizeCI_presentations', whicharea, neudesc, nnz(neuoi), nnz(neuctx)), 'interpreter', 'none')
xlabel('Time (ms)')
ylabel('Rate (Hz)')

subplot(1,2,2)
hold all
for a = 1:numel(probes)
    iprobe = probehierarchy(a);
%     a = visind(iprobe);
%     whicharea = visareas{iprobe};
    neuctx = neuctxall{iprobe};
    
    neuoi = neuctx;
    % neuoi = neuctx & RFCIall(iprobe).Pkw_rfclassic<0.05;
    
    trialsoi = floor(vis.sizeCI_presentations.trialorder/1000)==11;
    temp = 1000*psthall(iprobe).sizeCI_presentations(:,trialsoi,neuoi);
    % plot(psthtli, smooth(squeeze(mean(temp, [2,3])), 5), '-')
    temp = convn(temp, kergauss, 'same');
    tempbase = squeeze(mean(temp(psthtli>=-100 & psthtli<0,:,:), 'all'));
    temppeak = max(squeeze(mean(temp, [2,3])));
    plot(psthtli, (squeeze(mean(temp, [2,3]))-tempbase)/(temppeak-tempbase) )
end
legend(visareas{probehierarchy})
xlim([-50 200])
ylim([-0.2 1])
% xlim([-500 750])
set(gca, 'XGrid', 'on')
xlabel('Time (ms)')
ylabel('Normalized Activity')

%% compare CRF and IRF
Nrfs = 9;

figure('Position', [100 100 1000 600])
for iprobe = 1:numel(probes)
    a = visind(iprobe);
    whicharea = visareas{iprobe};
    neuctx = neuctxall{iprobe};
    
neuoi = neuctx;
% neuoi = neuctx & RFCIall(iprobe).Pkw_rfclassic<0.05;

    subplot(2,3,a)
    hold all
histogram2(RFCIall(iprobe).RFindclassic(neuoi), RFCIall(iprobe).RFindinverse(neuoi), 0.5:Nrfs+.5, 0.5:Nrfs+.5, 'displaystyle', 'tile')
plot([0.5 Nrfs+.5], [0.5 Nrfs+.5], 'w--', 'LineWidth', 0.5)
axis([0.5 Nrfs+.5 0.5 Nrfs+.5])
colorbar
xlabel('CRF')
ylabel('IRF')
title('2D Histogram')
title(sprintf('%s neurons N=%d\n2D-Histogram', whicharea, nnz(neuoi)), 'interpreter', 'none')
end

figure('Position', [1100 100 1000 600])
for iprobe = 1:numel(probes)
    a = visind(iprobe);
    whicharea = visareas{iprobe};
    neuctx = neuctxall{iprobe};
    
neuoi = neuctx;
% neuoi = neuctx & RFCIall(iprobe).Pkw_rfclassic<0.05;

histclassic = histcounts(RFCIall(iprobe).RFindclassic(neuoi), 0.5:Nrfs+.5, 'normalization', 'probability');
histinverse = histcounts(RFCIall(iprobe).RFindinverse(neuoi), 0.5:Nrfs+.5, 'normalization', 'probability');
expectedh2prob = histinverse'*histclassic;
actualh2prob = histcounts2(RFCIall(iprobe).RFindclassic(neuoi), RFCIall(iprobe).RFindinverse(neuoi), 0.5:Nrfs+.5, 0.5:Nrfs+.5, 'normalization', 'probability');

if abs(sum(expectedh2prob, 'all')-1)>=2^-32
    error('check expectedh2prob')
end
if abs(sum(actualh2prob, 'all')-1)>=2^-32
    error('check actualh2prob')
end

    subplot(2,3,a)
    hold all
imagesc(Nrfs^2* (actualh2prob - expectedh2prob) )
plot([0.5 Nrfs+.5], [0.5 Nrfs+.5], 'w--', 'LineWidth', 0.5)
axis([0.5 Nrfs+.5 0.5 Nrfs+.5])
set(gca, 'ydir', 'normal')
colorbar
caxis([-2 2])
xlabel('CRF')
ylabel('IRF')
title(sprintf('%s neurons N=%d\nActual-Expected Probability (X Nrfs^2)', whicharea, nnz(neuoi)), 'interpreter', 'none')

% disp(Nrfs^2*diag(actualh2prob - expectedh2prob)')
end

%% ctrCRF neurons PSTH on ctrCRF trials vs ctrIRF trials (ref. Keller et al)
kerwinhalf = 2; kersigma = 1;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));
% figure; plot(-kerwinhalf:kerwinhalf, kergauss)

whichneu=2;

figure('Position', [100 100 1000 600])
for iprobe = 1:numel(probes)
    a = visind(iprobe);
    whicharea = visareas{iprobe};
    neuctx = neuctxall{iprobe};
    
    switch whichneu
        case 1
            neuoi = neuctx;
            neudesc = 'all';
        case 2
            neuoi = neuctx & RFCIall(iprobe).RFindclassic==1 & RFCIall(iprobe).Pkw_rfclassic<0.05;
            neudesc = 'center-CRF';
    end
    
    subplot(2,3,a)
    hold all
    for icond = 1:2
        switch icond
            case 1
                trialsoi = vis.RFCI_presentations.trialorder(1:4:end)==1001;
            tempcol = [0 0 0];
            case 2
                trialsoi = vis.RFCI_presentations.trialorder(1:4:end)==11001;
            tempcol = [0 1 1];
        end
        temp = 1000*psthall(iprobe).RFCI_presentations(:,trialsoi,neuoi);
        % plot(psthtli, smooth(squeeze(mean(temp, [2,3])), 5), '-')
        temp = convn(temp, kergauss, 'same');
        plot(psthtli, squeeze(mean(temp, [2,3])), '-', 'Color', tempcol)
    end
    if a==6
    legend({'Classic', 'Inverse'})
    end
    xlim([-50 300])
set(gca, 'XGrid', 'on')
title(sprintf('%s %s neurons N=%d/%d\nRFCI_presentations', whicharea, neudesc, nnz(neuoi), nnz(neuctx)), 'interpreter', 'none')
xlabel('Time (ms)')
ylabel('Rate (Hz)')
end

%% Classic/Inverse Grating and IC responses in each block
tt2p=[106 107 110 111 506 511];
ttcol = [0 0.5 0; 0.8 0.4 0; 1 0.6 0.2; 0.2 1 0.2; 0 0 0.5; 0.2 0.2 1];
whichneu=2;
if whichneu == 1
    neudesc = 'all';
elseif whichneu == 2
    neudesc = 'center-CRF';
elseif whichneu >= 3
    iprefori = whichneu-2;
    neudesc = sprintf('%ddeg-preferring', vis.sizeCI_presentations.directions(iprefori));
end

figure
annotation('textbox', [0.1 0.9 0.8 0.1], 'string', [nwbsessions{ises} ' ' neudesc], 'FontSize', 14, 'EdgeColor', 'none')
for iprobe = 1:numel(probes)
    a = visind(iprobe);
    whicharea = visareas{iprobe};
    neuctx = neuctxall{iprobe};
    if whichneu == 1
        neuoi = neuctx;
    elseif whichneu == 2
        neuoi = neuctx & RFCIall(iprobe).RFindclassic==1 & RFCIall(iprobe).Pkw_rfclassic<0.05;
    elseif whichneu >= 3
        neuoi = neuctx & mod(oriparamsall(iprobe).prefiori-1,4)+1==iprefori & oriparamsall(iprobe).Pmww_OP<0.05;
    end
    
subplot(numel(probes),5,5*(a-1)+1)
hold all
for icond = 1:2
    switch icond
        case 1
            trialsoi = floor(vis.sizeCI_presentations.trialorder/1000)==4;
            tempcol = [0 0 0];
        case 2
            trialsoi = floor(vis.sizeCI_presentations.trialorder/1000)==14;
            tempcol = [0 1 1];
    end
    temp = 1000*psthall(iprobe).sizeCI_presentations(:,trialsoi,neuoi);
    % plot(psthtli, smooth(squeeze(mean(temp, [2,3])), 5), '-')
    temp = convn(temp, kergauss, 'same');
    plot(psthtli, squeeze(mean(temp, [2,3])), '-', 'Color', tempcol)
end
% legend({'Classic', 'Inverse'})
xlim([-50 200])
set(gca, 'XGrid', 'on')
title(sprintf('%s %s neurons N=%d/%d\nsizeCI_presentations', whicharea, neudesc, nnz(neuoi), nnz(neuctx)), 'interpreter', 'none')
xlabel('Time (ms)')
ylabel('Rate (Hz)')

for b= 1:4
    temptrialorder = ICtrialtypes( vis.(visblocks{b}).trialorder + 1);
    subplot(numel(probes),5,5*(a-1)+1+b)
    hold all
    for typi =1:numel(tt2p)
        trialsoi = temptrialorder==tt2p(typi);
        temp = 1000*psthall(iprobe).(visblocks{b})(:,trialsoi,neuoi);
        temp = convn(temp, kergauss, 'same');
        plot(psthtli, squeeze(mean(temp, [2,3])), '-', 'Color', ttcol(typi,:));%, 'LineWidth', 1)
    end
%     legend(num2str(tt2p'))
    xlim([-50 200])
    set(gca, 'XGrid', 'on')
    title(visblocks{b}, 'interpreter', 'none')
    xlabel('Time (ms)')
    ylabel('Rate (Hz)')
end
end

%% Classic/Inverse Grating and IC responses in each block
tt2p=[106 107 110 111 506 511];
ttcol = [0 0.5 0; 0.8 0.4 0; 1 0.6 0.2; 0.2 1 0.2; 0 0 0.5; 0.2 0.2 1];
whichneu=7;
switch whichneu
    case 1
        neudesc = 'ICencoder1';
    case 2
        neudesc = 'RCencoder1';
    case 3
        neudesc = 'RCencoder2';
    case 4
        neudesc = 'ICencoder2';
    case 5
        neutit = 'RElfaith1';
    case 6
        neutit = 'RElfaith2';
    case 7
        neudesc = 'ICresp1';
end

figure
annotation('textbox', [0.1 0.9 0.8 0.1], 'string', [nwbsessions{ises} ' ' neudesc], 'FontSize', 14, 'EdgeColor', 'none')
for iprobe = 1:numel(probes)
    a = visind(iprobe);
    whicharea = visareas{iprobe};
    neuctx = neuctxall{iprobe};

for b= 1:4
    neuoi = neuctx & ICsigall(iprobe).(visblocks{b}).(neudesc);
    
    temptrialorder = ICtrialtypes( vis.(visblocks{b}).trialorder + 1);
    subplot(numel(probes),4,4*(a-1)+b)
    hold all
    for typi =1:numel(tt2p)
        trialsoi = temptrialorder==tt2p(typi);
        temp = 1000*psthall(iprobe).(visblocks{b})(:,trialsoi,neuoi);
        temp = convn(temp, kergauss, 'same');
        plot(psthtli, squeeze(mean(temp, [2,3])), '-', 'Color', ttcol(typi,:));%, 'LineWidth', 1)
    end
%     legend(num2str(tt2p'))
    xlim([-50 200])
    set(gca, 'XGrid', 'on')
    title(sprintf('%s %s %s N=%d', whicharea, visblocks{b}, neudesc, nnz(neuoi)), 'interpreter', 'none')
    xlabel('Time (ms)')
    ylabel('Rate (Hz)')
end
end

%%
b=1;
tt2p=[106 107 110 111];
ttcol = [0 0.5 0; 0.8 0.4 0; 1 0.6 0.2; 0.2 1 0.2];
temptrialorder = ICtrialtypes( vis.(visblocks{b}).trialorder + 1);
figure
for isp =1:4
    switch isp
        case 1
            neutit = 'ICencoder1';
        case 2
            neutit = 'RCencoder1';
        case 3
            neutit = 'RCencoder2';
        case 4
            neutit = 'ICencoder2';
    end
    neuoi = ICsig.(visblocks{b}).(neutit);
    subplot(2,2,isp)
    hold all
    for t = 1:numel(tt2p)
        trialsoi= temptrialorder==tt2p(t);
        plot(psthtli, 1000*squeeze(mean(psth.(visblocks{b})(:, trialsoi, neuoi), [2,3])), 'Color', ttcol(t,:), 'LineWidth', 1)
    end
    title(sprintf('%s N=%d', neutit, nnz(neuoi)))
    xlim([-100 500])
    xlabel('Time (ms)')
    ylabel('Rate (Hz)')
end


%%
neucol = [0 0.5 0; 0 0 0.5; 0.2 1 0.2; 0.2 0.2 1];
figure
for b= 1:4
    temptrialorder = ICtrialtypes( vis.(visblocks{b}).trialorder + 1);
    trialsoi= temptrialorder==106;
    
    subplot(2,2,b)
    hold all
    for ii =1:4
        switch ii
            case 1
                neutit = 'ICencoder1';
            case 2
                neutit = 'RElfaith1';
            case 3
                neutit = 'ICencoder2';
            case 4
                neutit = 'RElfaith2';
        end
        neuoi = ICsig.(visblocks{b}).(neutit);
        plot(psthtli, smooth(1000*squeeze(mean(psth.(visblocks{b})(:, trialsoi, neuoi), [2,3])), 5), 'Color', neucol(ii,:), 'LineWidth', 1)
    end
    title(visblocks{b}, nnz(neuoi))
    xlim([-100 500])
    xlabel('Time (ms)')
    ylabel('Rate (Hz)')
end


%%
b=1;
ttcol = [0 0.7 0; 0 0 1];
figure
for isp =1:2
    switch isp
        case 1
            neutit = 'ICencoder1';
            tt2p=[106 506];
        case 2
            neutit = 'ICencoder2';
            tt2p=[111 511];
    end
    for b = 1:4
        temptrialorder = ICtrialtypes( vis.(visblocks{b}).trialorder + 1);
        neuoi = ICsig.(visblocks{b}).(neutit);
        subplot(2,4,4*(isp-1)+b)
        hold all
        for t = 1:numel(tt2p)
            trialsoi= temptrialorder==tt2p(t);
            plot(psthtli, 1000*squeeze(mean(psth.(visblocks{b})(:, trialsoi, neuoi), [2,3])), 'Color', ttcol(t,:), 'LineWidth', 1)
        end
        title(sprintf('%s N=%d', neutit, nnz(neuoi)))
        xlim([-100 500])
        xlabel('Time (ms)')
        ylabel('Rate (Hz)')
    end
end