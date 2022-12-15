% datadir = 'D:\OpenScopeData\000248v220916\';
% nwbdir = dir(datadir);
% nwbsessions = {nwbdir.name};
% nwbsessions = nwbsessions(contains(nwbsessions, 'sub-'));
% Nsessions= 4;%numel(nwbsessions)-1;
% 
% probes = {'A', 'B', 'C', 'D', 'E', 'F'};
% visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
%     'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};

datadir = 'D:\OpenScopeData\000248v221123\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name}; 
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

probes = {'C', 'D'};
visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations'};

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
        spklatencyagg.(whichbin).(visblocks{b}).(probes{iprobe}) = [];
        spklatencyprobagg.(whichbin).(visblocks{b}).(probes{iprobe}) = [];
        spklatencyadjagg.(whichbin).(visblocks{b}).(probes{iprobe}) = [];
        spklatencyadjprobagg.(whichbin).(visblocks{b}).(probes{iprobe}) = [];
    end
end
for ises = 1:Nsessions
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur'
    
    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);

    for iprobe = 1:numel(probes)
        load(sprintf('%sspikelatency_probe%s.mat', pathpp, probes{iprobe}))
    load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}), 'neuoind')
        
    % check whether CCF registration is correct
    probelocs = electrode_location(ismember(electrode_id, unit_peakch(neuoind)));
    
    neuloc = electrode_location(revmapelecid(unit_peakch(neuoind)+1));
    if ~isequal(unique(probelocs), unique(neuloc))
        disp(unique(neuloc)')
        error('check neuloc')
    end
    
%     whichneuctx = 1; % 1: correct definition, 2: >=230, 3: topmost-250um
%     switch whichneuctx
%         case 1
%             neuctx = contains(neuloc, 'VIS');
%         case 2
%             % assume that >=230 is neocortex
%             neuctx = mod(unit_peakch(neuoind), 1000)>=230;
%             warning('for now, assuming that electrode_localid>=230 is neocortex (until CCF registration is corrected)')
%             % 2 electrodes per depth, vertical spacing is 20um
%             % typical electrode span in cortex 120 : (120/2)*20 = ~1200um
%         case 3
%             % to approximately get layer 2/3 cells, take units within 250um of the topmost unit
%             neuctx = mod(unit_peakch(neuoind), 1000)>=max(mod(unit_peakch(neuoind), 1000))-25;
%             warning('for now, to approximately get layer 2/3 cells, take units within 250um of the topmost unit')
%     end
    
        for b = 1:numel(visblocks)
            for ibin = 1:numel(spklatencybins)
                whichbin = sprintf('bin%dms',spklatencybins(ibin));
                spklatencyagg.(whichbin).(visblocks{b}).(probes{iprobe}) = cat(1, ...
                    spklatencyagg.(whichbin).(visblocks{b}).(probes{iprobe}), ...
                    spklatency.(whichbin).(visblocks{b}) ); %(neuctx,:) );
                spklatencyprobagg.(whichbin).(visblocks{b}).(probes{iprobe}) = cat(1, ...
                    spklatencyprobagg.(whichbin).(visblocks{b}).(probes{iprobe}), ...
                    spklatencyprob.(whichbin).(visblocks{b}) ); %(neuctx,:) );
                spklatencyadjagg.(whichbin).(visblocks{b}).(probes{iprobe}) = cat(1, ...
                    spklatencyadjagg.(whichbin).(visblocks{b}).(probes{iprobe}), ...
                    spklatencyadj.(whichbin).(visblocks{b}) ); %(neuctx,:) );
                spklatencyadjprobagg.(whichbin).(visblocks{b}).(probes{iprobe}) = cat(1, ...
                    spklatencyadjprobagg.(whichbin).(visblocks{b}).(probes{iprobe}), ...
                    spklatencyadjprob.(whichbin).(visblocks{b}) ); %(neuctx,:) );
            end
        end
    end
end

%% determine spklatencyprob threshold: 0.01
iprobe = find(strcmp(probes, 'C'));

typi = ICtrialtypes==106;
whichvisblock = 'ICwcfg1_presentations';
figure;
for ibin = 1:numel(spklatencybins)
    whichbin = sprintf('bin%dms',spklatencybins(ibin));
    subplot(2,2,ibin)
    histogram(spklatencyadjprobagg.(whichbin).(whichvisblock).(probes{iprobe})(neuctxagg{iprobe}==1, typi) )%, 'normalization', 'cdf')
    title(whichbin)
end

%% on IC and REl trials, compare IC-encoder vs inducer-encoder latency
bw=10;
adj = false;
if adj
    tempspklatencyprobagg = spklatencyadjprobagg;
    tempspklatencyagg = spklatencyadjagg;
    tempspklatencydesc = 'adjusted spike latency';
else
    tempspklatencyprobagg = spklatencyprobagg;
    tempspklatencyagg = spklatencyagg;
    tempspklatencydesc = 'spike latency';
end

iprobe = find(strcmp(probes, 'C'));
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
    %neuoind = neuoind(tempspklatencyprobagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind, typi)>0.01);
    histogram(tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind, typi), 'binwidth', bw, 'normalization', 'probability')
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
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind1, ICtrialtypes==106), ...
                    tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind2, ICtrialtypes==111));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind1 = probeneuctxind(ICsigagg.(whichvisblock)(iprobe).ICencoder1==1);
                neuoind2 = probeneuctxind(ICsigagg.(whichvisblock)(iprobe).ICencoder2==1);
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind1, ICtrialtypes==106), ...
                    tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind2, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
    histogram(tempspklatency, 'binwidth', bw, 'normalization', 'probability')
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
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind1, ICtrialtypes==506), ...
                    tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind2, ICtrialtypes==511));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind1 = probeneuctxind(ICsigagg.(whichvisblock)(iprobe).indenc13==1);
                neuoind2 = probeneuctxind(ICsigagg.(whichvisblock)(iprobe).indenc24==1);
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind1, ICtrialtypes==106), ...
                    tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind2, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
    histogram(tempspklatency, 'binwidth', bw, 'normalization', 'probability')
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
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind1, ICtrialtypes==506), ...
                    tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind2, ICtrialtypes==511));
                tempspklatency1 = tempspklatency;
            case 2
                neuoind1 = probeneuctxind(ICsigagg.(whichvisblock)(iprobe).ICencoder1==1);
                neuoind2 = probeneuctxind(ICsigagg.(whichvisblock)(iprobe).ICencoder2==1);
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind1, ICtrialtypes==106), ...
                    tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind2, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
    histogram(tempspklatency, 'binwidth', bw, 'normalization', 'probability')
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

iprobe = find(strcmp(probes, 'C'));
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
%     neuoind = neuoind(tempspklatencyprobagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind, typi)>spklatencyprobthr);
    histogram(tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind, typi), 'binwidth', 5, 'normalization', 'probability')
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
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind, ICtrialtypes==106), ...
                    tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind, ICtrialtypes==111));
                tempspklatency1 = tempspklatency;
            case 2
    neuoind = probeneuctxind(RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05);
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind, ICtrialtypes==106), ...
                    tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind, ICtrialtypes==111));
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
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind, ICtrialtypes==506), ...
                    tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind, ICtrialtypes==511));
                tempspklatency1 = tempspklatency;
            case 2
    neuoind = probeneuctxind(RFCIagg(iprobe).RFindclassic~=1 & RFCIagg(iprobe).pRFclassic<0.05);
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind, ICtrialtypes==106), ...
                    tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind, ICtrialtypes==111));
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
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind, ICtrialtypes==506), ...
                    tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind, ICtrialtypes==511));
                tempspklatency1 = tempspklatency;
            case 2
    neuoind = probeneuctxind(RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05);
                tempspklatency = cat(1, tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind, ICtrialtypes==106), ...
                    tempspklatencyagg.(whichbin).(whichvisblock).(probes{iprobe})(neuoind, ICtrialtypes==111));
                tempspklatency2 = tempspklatency;
        end
    histogram(tempspklatency, 'binwidth', 5, 'normalization', 'probability')
    end
    xlabel('Spike Latency (ms)')
    ylabel('Probability')
    p = ranksum(tempspklatency1, tempspklatency2);
    title(sprintf('%s p=%.4f', whichbin, p))
end