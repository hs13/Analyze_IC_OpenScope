datadir = 'D:\OpenScopeData\000248\';

probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
    'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
spklatencybins = [1 5 10 25];
% vistrialtypes = struct();
% vistrialrep = struct();
spklatencyagg = struct();
spklatencyprobagg = struct();
spklatencyadjagg = struct();
spklatencyadjprobagg = struct();
for b = 1:numel(visblocks)
    for ibin = 1:numel(spklatencybins)
        whichbin = sprintf('bin%dms',spklatencybins(ibin));
        spklatencyagg.(whichbin).(visblocks{b}) = cell(size(probes));
        spklatencyprobagg.(whichbin).(visblocks{b}) = cell(size(probes));
        spklatencyadjagg.(whichbin).(visblocks{b}) = cell(size(probes));
        spklatencyadjprobagg.(whichbin).(visblocks{b}) = cell(size(probes));
    end
end
Nsessions=4;
for ises = 1:Nsessions
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    for iprobe = 1:numel(probes)
        load(sprintf('%sspikelatency_probe%s.mat', pathpp, probes{iprobe}))
        for b = 1:numel(visblocks)
            for ibin = 1:numel(spklatencybins)
                whichbin = sprintf('bin%dms',spklatencybins(ibin));
                spklatencyagg.(whichbin).(visblocks{b}){iprobe} = cat(1, ...
                    spklatencyagg.(whichbin).(visblocks{b}){iprobe}, spklatency.(whichbin).(visblocks{b}));
                spklatencyprobagg.(whichbin).(visblocks{b}){iprobe} = cat(1, ...
                    spklatencyprobagg.(whichbin).(visblocks{b}){iprobe}, spklatencyprob.(whichbin).(visblocks{b}));
                spklatencyadjagg.(whichbin).(visblocks{b}){iprobe} = cat(1, ...
                    spklatencyadjagg.(whichbin).(visblocks{b}){iprobe}, spklatencyadj.(whichbin).(visblocks{b}));
                spklatencyadjprobagg.(whichbin).(visblocks{b}){iprobe} = cat(1, ...
                    spklatencyadjprobagg.(whichbin).(visblocks{b}){iprobe}, spklatencyadjprob.(whichbin).(visblocks{b}));
            end
        end
    end
end

%% determine spklatencyprob threshold: 0.01
iprobe = 3;
typi = ICtrialtypes==106;
whichvisblock = 'ICwcfg1_presentations';
figure;
for ibin = 1:numel(spklatencybins)
    whichbin = sprintf('bin%dms',spklatencybins(ibin));
    subplot(2,2,ibin)
    histogram(spklatencyadjprobagg.(whichbin).(whichvisblock){iprobe}(neuctxagg{iprobe}==1, typi) )%, 'normalization', 'cdf')
    title(whichbin)
end

%% on IC and REl trials, compare IC-encoder vs inducer-encoder latency
adj = true;
if adj
    tempspklatencyprobagg = spklatencyadjprobagg;
    tempspklatencyagg = spklatencyadjagg;
    tempspklatencydesc = 'adjusted spike lstency';
else
    tempspklatencyprobagg = spklatencyprobagg;
    tempspklatencyagg = spklatencyagg;
    tempspklatencydesc = 'spike lstency';
end

iprobe = 3;
typi = ICtrialtypes==106;
whichvisblock = 'ICwcfg1_presentations';
figure;
annotation('textbox', [0.05 0.91 0.8 0.1], 'string', [tempspklatencydesc ' on IC1 trials: IC1-encoder (red) vs pacman13-encoder (blue)'], 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    whichbin = sprintf('bin%dms',spklatencybins(ibin));
    subplot(2,2,ibin)
    hold all
    probeneuctxind = find(neuctxagg{iprobe}==1);
    for ii = 1:2
        switch ii
            case 1
    neuoind = probeneuctxind(ICsigagg.(whichvisblock)(iprobe).indenc13==1);
            case 2
    neuoind = probeneuctxind(ICsigagg.(whichvisblock)(iprobe).ICencoder1==1);
        end
    neuoind = neuoind(tempspklatencyprobagg.(whichbin).(whichvisblock){iprobe}(neuoind, typi)>0.01);
    histogram(tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind, typi), 'binwidth', 5, 'normalization', 'probability')
    end
    title(whichbin)
end

figure;
annotation('textbox', [0.05 0.91 0.8 0.1], 'string', [tempspklatencydesc ' on IC trials: IC-encoder (red) vs pacman-encoder (blue)'], 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    whichbin = sprintf('bin%dms',spklatencybins(ibin));
    subplot(2,2,ibin)
    hold all
    probeneuctxind = find(neuctxagg{iprobe}==1);
    for ii = 1:2
        switch ii
            case 1
                neuoind1 = probeneuctxind(ICsigagg.(whichvisblock)(iprobe).indenc13==1);
                neuoind2 = probeneuctxind(ICsigagg.(whichvisblock)(iprobe).indenc24==1);
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind1, ICtrialtypes==106), ...
                    tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind2, ICtrialtypes==111));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind1 = probeneuctxind(ICsigagg.(whichvisblock)(iprobe).ICencoder1==1);
                neuoind2 = probeneuctxind(ICsigagg.(whichvisblock)(iprobe).ICencoder2==1);
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind1, ICtrialtypes==106), ...
                    tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind2, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
    histogram(tempspklatency, 'binwidth', 5, 'normalization', 'probability')
    end
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = ranksum(tempspklatency1, tempspklatency2);
    title(sprintf('%s p=%.4f', whichbin, p))
end

figure;
annotation('textbox', [0.05 0.91 0.8 0.1], 'string', [tempspklatencydesc ' pacman-encoder: IC (red) vs REl (blue) trials'], 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    whichbin = sprintf('bin%dms',spklatencybins(ibin));
    subplot(2,2,ibin)
    hold all
    probeneuctxind = find(neuctxagg{iprobe}==1);
    for ii = 1:2
        switch ii
            case 1
                neuoind1 = probeneuctxind(ICsigagg.(whichvisblock)(iprobe).indenc13==1);
                neuoind2 = probeneuctxind(ICsigagg.(whichvisblock)(iprobe).indenc24==1);
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind1, ICtrialtypes==506), ...
                    tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind2, ICtrialtypes==511));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind1 = probeneuctxind(ICsigagg.(whichvisblock)(iprobe).indenc13==1);
                neuoind2 = probeneuctxind(ICsigagg.(whichvisblock)(iprobe).indenc24==1);
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind1, ICtrialtypes==106), ...
                    tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind2, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
    histogram(tempspklatency, 'binwidth', 5, 'normalization', 'probability')
    end
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = ranksum(tempspklatency1, tempspklatency2);
    title(sprintf('%s p=%.4f', whichbin, p))
end

figure;
annotation('textbox', [0.05 0.91 0.8 0.1], 'string', [tempspklatencydesc ' IC-encoder: IC (red) vs REl (blue) trials'], 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    whichbin = sprintf('bin%dms',spklatencybins(ibin));
    subplot(2,2,ibin)
    hold all
    probeneuctxind = find(neuctxagg{iprobe}==1);
    for ii = 1:2
        switch ii
            case 1
                neuoind1 = probeneuctxind(ICsigagg.(whichvisblock)(iprobe).ICencoder1==1);
                neuoind2 = probeneuctxind(ICsigagg.(whichvisblock)(iprobe).ICencoder2==1);
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind1, ICtrialtypes==506), ...
                    tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind2, ICtrialtypes==511));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind1 = probeneuctxind(ICsigagg.(whichvisblock)(iprobe).ICencoder1==1);
                neuoind2 = probeneuctxind(ICsigagg.(whichvisblock)(iprobe).ICencoder2==1);
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind1, ICtrialtypes==106), ...
                    tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind2, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
    histogram(tempspklatency, 'binwidth', 5, 'normalization', 'probability')
    end
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = ranksum(tempspklatency1, tempspklatency2);
    title(sprintf('%s p=%.4f', whichbin, p))
end


%% on IC and REl trials, compare ctrRF vs non-ctrRF latency
adj = false;
if adj
    tempspklatencyprobagg = spklatencyadjprobagg;
    tempspklatencyagg = spklatencyadjagg;
    tempspklatencydesc = 'adjusted spike lstency';
else
    tempspklatencyprobagg = spklatencyprobagg;
    tempspklatencyagg = spklatencyagg;
    tempspklatencydesc = 'spike lstency';
end

iprobe = 3;
typi = ICtrialtypes==106;
whichvisblock = 'ICwcfg1_presentations';
figure;
for ibin = 1:numel(spklatencybins)
    whichbin = sprintf('bin%dms',spklatencybins(ibin));
    subplot(2,2,ibin)
    hold all
    probeneuctxind = find(neuctxagg{iprobe}==1);
    for ii = 1:2
        switch ii
            case 1
    neuoind = probeneuctxind(RFCIagg(iprobe).RFindclassic~=1 & RFCIagg(iprobe).pRFclassic<0.05);
            case 2
    neuoind = probeneuctxind(RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05);
        end
    neuoind = neuoind(tempspklatencyprobagg.(whichbin).(whichvisblock){iprobe}(neuoind, typi)>0.01);
    histogram(tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind, typi), 'binwidth', 5, 'normalization', 'probability')
    end
    title(whichbin)
end

figure;
annotation('textbox', [0.05 0.91 0.8 0.1], 'string', [tempspklatencydesc ' on IC trials: center-RF (red) vs non-center-RF (blue)'], 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    whichbin = sprintf('bin%dms',spklatencybins(ibin));
    subplot(2,2,ibin)
    hold all
    probeneuctxind = find(neuctxagg{iprobe}==1);
    for ii = 1:2
        switch ii
            case 1
    neuoind = probeneuctxind(RFCIagg(iprobe).RFindclassic~=1 & RFCIagg(iprobe).pRFclassic<0.05);
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind, ICtrialtypes==106), ...
                    tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind, ICtrialtypes==111));
                tempspklatency1 = tempspklatency;
            case 2
    neuoind = probeneuctxind(RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05);
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind, ICtrialtypes==106), ...
                    tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
    histogram(tempspklatency, 'binwidth', 5, 'normalization', 'probability')
    end
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = ranksum(tempspklatency1, tempspklatency2);
    title(sprintf('%s p=%.4f', whichbin, p))
end

figure;
annotation('textbox', [0.05 0.91 0.8 0.1], 'string', [tempspklatencydesc ' non-center-RF neurons: IC (red) vs REl (blue) trials'], 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    whichbin = sprintf('bin%dms',spklatencybins(ibin));
    subplot(2,2,ibin)
    hold all
    probeneuctxind = find(neuctxagg{iprobe}==1);
    for ii = 1:2
        switch ii
            case 1
    neuoind = probeneuctxind(RFCIagg(iprobe).RFindclassic~=1 & RFCIagg(iprobe).pRFclassic<0.05);
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind, ICtrialtypes==506), ...
                    tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind, ICtrialtypes==511));
                tempspklatency1 = tempspklatency;
            case 2
    neuoind = probeneuctxind(RFCIagg(iprobe).RFindclassic~=1 & RFCIagg(iprobe).pRFclassic<0.05);
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind, ICtrialtypes==106), ...
                    tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
    histogram(tempspklatency, 'binwidth', 5, 'normalization', 'probability')
    end
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = ranksum(tempspklatency1, tempspklatency2);
    title(sprintf('%s p=%.4f', whichbin, p))
end

figure;
annotation('textbox', [0.05 0.91 0.8 0.1], 'string', [tempspklatencydesc ' center-RF neurons: IC (red) vs REl (blue) trials'], 'edgecolor', 'none')
for ibin = 1:numel(spklatencybins)
    whichbin = sprintf('bin%dms',spklatencybins(ibin));
    subplot(2,2,ibin)
    hold all
    probeneuctxind = find(neuctxagg{iprobe}==1);
    for ii = 1:2
        switch ii
            case 1
    neuoind = probeneuctxind(RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05);
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind, ICtrialtypes==506), ...
                    tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind, ICtrialtypes==511));
                tempspklatency1 = tempspklatency;
            case 2
    neuoind = probeneuctxind(RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05);
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind, ICtrialtypes==106), ...
                    tempspklatencyagg.(whichbin).(whichvisblock){iprobe}(neuoind, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
    histogram(tempspklatency, 'binwidth', 5, 'normalization', 'probability')
    end
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = ranksum(tempspklatency1, tempspklatency2);
    title(sprintf('%s p=%.4f', whichbin, p))
end