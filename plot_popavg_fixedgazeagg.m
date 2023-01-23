load(['G:\My Drive\DATA\ICexpts_submission22\openscope_popavg_fixedgaze_agg.mat'])
Nsessions = numel(nwbsessions);

ICblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations'};
neuprobe_fixedgazeagg = [];
for b = 1:numel(ICblocks)
    for f= 1:numel(ICsigfields)
        ICsig_fixedgazeagg.(ICblocks{b}).all.(ICsigfields{f}) = [];
    end
end
% RFCIfields = fieldnames(RFCI);
for f= 1:numel(RFCIfields)
    RFCI_fixedgazeagg.all.(RFCIfields{f}) = [];
end
for f= 1:numel(RFCIspinfields)
    RFCIspin_fixedgazeagg.all.(RFCIspinfields{f}) = [];
end
% size vector [0, 4, 8, 16, 32, 64 ]
for f= 1:numel(sizeCIfields)
    sizeCI_fixedgazeagg.all.(sizeCIfields{f}) = [];
end
for f= 1:numel(oriparamsfields)
    oriparams_fixedgazeagg.all.(oriparamsfields{f}) = [];
end
for f= 1:numel(ori4paramsfields)
    ori4params_fixedgazeagg.all.(ori4paramsfields{f}) = [];
end

for iprobe = 1:numel(probes)
    neuprobe_fixedgazeagg = cat(1, neuprobe_fixedgazeagg, iprobe*ones(size(neuloc_fixedgazeagg{iprobe})));
    
for b = 1:numel(ICblocks)
    for f= 1:numel(ICsigfields)
        ICsig_fixedgazeagg.(ICblocks{b}).all.(ICsigfields{f}) = cat(1, ICsig_fixedgazeagg.(ICblocks{b}).all.(ICsigfields{f}), ICsig_fixedgazeagg.(ICblocks{b}).(probes{iprobe}).(ICsigfields{f}) );
    end
end

% RFCIfields = fieldnames(RFCI);
for f= 1:numel(RFCIfields)
    RFCI_fixedgazeagg.all.(RFCIfields{f}) = cat(1, RFCI_fixedgazeagg.all.(RFCIfields{f}), RFCI_fixedgazeagg.(probes{iprobe}).(RFCIfields{f}) );
end

for f= 1:numel(RFCIspinfields)
    RFCIspin_fixedgazeagg.all.(RFCIspinfields{f}) = cat(1, RFCIspin_fixedgazeagg.all.(RFCIspinfields{f}), RFCIspin_fixedgazeagg.(probes{iprobe}).(RFCIspinfields{f}) );
end

% size vector [0, 4, 8, 16, 32, 64 ]
for f= 1:numel(sizeCIfields)
    sizeCI_fixedgazeagg.all.(sizeCIfields{f}) = cat(1, sizeCI_fixedgazeagg.all.(sizeCIfields{f}), sizeCI_fixedgazeagg.(probes{iprobe}).(sizeCIfields{f}) );
end

for f= 1:numel(oriparamsfields)
    oriparams_fixedgazeagg.all.(oriparamsfields{f}) = cat(1, oriparams_fixedgazeagg.all.(oriparamsfields{f}), oriparams_fixedgazeagg.(probes{iprobe}).(oriparamsfields{f}) );
end

for f= 1:numel(ori4paramsfields)
    ori4params_fixedgazeagg.all.(ori4paramsfields{f}) = cat(1, ori4params_fixedgazeagg.all.(ori4paramsfields{f}), ori4params_fixedgazeagg.(probes{iprobe}).(ori4paramsfields{f}) );
end
end


%% proportion of neurons in each area
% compare areas
recarealabels = {'LGd', 'LP', 'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam', ...
    'VISp1', 'VISp2/3', 'VISp4', 'VISp5', 'VISp6a', 'VISp6b'};
recarea = {'LGN', 'LP', 'V1', 'LM', 'RL', 'AL', 'PM', 'AM', ...
    'V1L1', 'V1L23', 'V1L4', 'V1L5', 'V1L6a', 'V1L6b'};

% VISp NOT VISpm

% compare layers in V1
neutypes = {'all', 'ICencoder', 'RCencoder', 'inducerencoder', 'indin', 'sigBK', 'sigBKnotsigBI'};
Nneurons_fixedgazeacc = struct();
Nneurons_fixedgazeagg = struct();
for a = 1:numel(recarealabels)
    if strcmp(recarealabels{a}, 'VISp')
        neuinarea = contains(cat(1,neuloc_fixedgazeagg{:}), 'VISp') & ~contains(cat(1,neuloc_fixedgazeagg{:}), 'VISpm');
    else
        neuinarea = contains(cat(1,neuloc_fixedgazeagg{:}), recarealabels{a});
    end
for b= 1:numel(ICblocks)
    for ineu = 1:numel(neutypes)
        switch neutypes{ineu}
            case 'all'
                neuoi = true(size(ICsig_fixedgazeagg.(ICblocks{b}).all.ICencoder));
            case 'ICencoder'
                neuoi = ICsig_fixedgazeagg.(ICblocks{b}).all.ICencoder==1;
            case 'RCencoder'
                neuoi = ICsig_fixedgazeagg.(ICblocks{b}).all.RCencoder==1;
            case 'inducerencoder'
                neuoi = ICsig_fixedgazeagg.(ICblocks{b}).all.inducerencoder==1;
            case 'indin'
                neuoi = ICsig_fixedgazeagg.(ICblocks{b}).all.indin1==1 | ICsig_fixedgazeagg.(ICblocks{b}).all.indin2==1 ...
                    | ICsig_fixedgazeagg.(ICblocks{b}).all.indin3==1 | ICsig_fixedgazeagg.(ICblocks{b}).all.indin4==1;
            case 'sigBK'
                neuoi = ICsig_fixedgazeagg.(ICblocks{b}).all.PkwBK<0.05;
            case 'sigBKnotsigBI'
                neuoi = ICsig_fixedgazeagg.(ICblocks{b}).all.PkwBK<0.05 & ICsig_fixedgazeagg.(ICblocks{b}).all.PkwBI>=0.05;
        end
        neuoi = neuinarea & neuoi;
        Nneurons_fixedgazeacc.(recarea{a}).(ICblocks{b}).(neutypes{ineu}) = nnz(neuoi);
        Nneurons_fixedgazeagg.(recarea{a}).(ICblocks{b}).(neutypes{ineu}) = zeros(Nsessions, 1);
        for ises = 1:Nsessions
            sesneu = cat(1,sesneu_fixedgazeagg{:})==ises;
            Nneurons_fixedgazeagg.(recarea{a}).(ICblocks{b}).(neutypes{ineu})(ises) = nnz(neuoi & sesneu);
        end
    end
end
end

%% compare proportion IC-/RC-encoders across sessions
whichvisblock = 'ICwcfg1_presentations';

pltmeanoragg = 2; % 1 is average proportion of neurons across sessions, 2 is aggregate across sessions

fs = 10;
figure
annotation('textbox', [0.1 0.92 0.9 0.1], 'string', ['Fixed Gaze Trials: ' whichvisblock], 'edgecolor', 'none', 'interpreter', 'none')
for ii = 1:2
    switch ii
        case 1
            neudenom = 'all';
        case 2
            neudenom = 'sigBK';
    end
    for jj = 1:2
    switch jj
        case 1
recareas = {'LGN', 'LP', 'V1', 'LM', 'RL', 'AL', 'PM', 'AM'};
        case 2
recareas = {'V1L1', 'V1L23', 'V1L4', 'V1L5', 'V1L6a', 'V1L6b'};
    end
tempmat = NaN(numel(recareas),2);
for a = 1:numel(recareas)
    switch pltmeanoragg
        case 1
    denom = Nneurons_fixedgazeagg.(recarea{a}).(whichvisblock).(neudenom);
    tempmat(a,1) = nanmean(Nneurons_fixedgazeagg.(recarea{a}).(whichvisblock).ICencoder./denom);
    tempmat(a,2) = nanmean(Nneurons_fixedgazeagg.(recarea{a}).(whichvisblock).RCencoder./denom);
        case 2
    denom = Nneurons_fixedgazeacc.(recarea{a}).(whichvisblock).(neudenom);
    tempmat(a,1) = Nneurons_fixedgazeacc.(recarea{a}).(whichvisblock).ICencoder/denom;
    tempmat(a,2) = Nneurons_fixedgazeacc.(recarea{a}).(whichvisblock).RCencoder/denom;
    end
end
yl = [0 .2];
subplot(2,2,2*(jj-1)+ii)
hold all
b = bar(tempmat);
b(1).FaceColor = [0 0.7 0];
b(1).EdgeColor = 'none';
b(2).FaceColor = [1 0.5 0];
b(2).EdgeColor = 'none';
yl=ylim; yl = [0 1.2*yl(2)];
for a = 1:numel(recareas)
    denom = Nneurons_fixedgazeagg.(recarea{a}).(whichvisblock).(neudenom);
    tempICvec = Nneurons_fixedgazeagg.(recarea{a}).(whichvisblock).ICencoder./denom;
    tempRCvec = Nneurons_fixedgazeagg.(recarea{a}).(whichvisblock).RCencoder./denom;
    plot(a+.1*[-1 1], [tempICvec tempRCvec], '-', 'Color', [0.5 0.5 0.5 0.5])
    try
    p = signrank(tempICvec, tempRCvec);
    if p<0.05
        fontcol = [1 0 0];
    else
        fontcol = [0 0 0];
    end
    t=text(a,yl(2)-.1*range(yl), sprintf('p=%.2f', p), 'Color', fontcol, 'FontSize', fs,'HorizontalAlignment', 'center');
    set(t,'Rotation',45);
    catch
    end
end
ylim(yl)
set(gca, 'FontSize', fs, 'XTick', 1:numel(recareas), 'XTickLabel', recareas, 'XTickLabelRotation', 45)
ylabel('Proportion Neurons')
if ii==2 && jj==2
legend({'IC', 'RC'}, 'Location', 'SouthEast')
end
title(sprintf('out of %s neurons', neudenom))
    end
end


figure
annotation('textbox', [0.1 0.92 0.9 0.1], 'string', ['Fixed Gaze Trials: ' whichvisblock], 'edgecolor', 'none', 'interpreter', 'none')
for ii = 1:2
    switch ii
        case 1
            neudenom = 'all';
        case 2
            neudenom = 'sigBK';
    end
    for jj = 1:2
    switch jj
        case 1
recareas = {'LGN', 'LP', 'V1', 'LM', 'RL', 'AL', 'PM', 'AM'};
        case 2
recareas = {'V1L1', 'V1L23', 'V1L4', 'V1L5', 'V1L6a', 'V1L6b'};
    end
tempmat = NaN(numel(recareas),2);
for a = 1:numel(recareas)
    switch pltmeanoragg
        case 1
    denom = Nneurons_fixedgazeagg.(recarea{a}).(whichvisblock).(neudenom);
    tempmat(a,1) = nanmean(Nneurons_fixedgazeagg.(recarea{a}).(whichvisblock).indin./denom);
    tempmat(a,2) = nanmean(Nneurons_fixedgazeagg.(recarea{a}).(whichvisblock).sigBKnotsigBI./denom);
        case 2
    denom = Nneurons_fixedgazeacc.(recarea{a}).(whichvisblock).(neudenom);
    tempmat(a,1) = Nneurons_fixedgazeacc.(recarea{a}).(whichvisblock).indin/denom;
    tempmat(a,2) = Nneurons_fixedgazeacc.(recarea{a}).(whichvisblock).sigBKnotsigBI/denom;
    end
end
yl = [0 .2];
subplot(2,2,2*(jj-1)+ii)
hold all
b = bar(tempmat);
b(1).FaceColor = [0.5 0.5 0.5];
b(1).EdgeColor = 'none';
b(2).FaceColor = [0 1 1];
b(2).EdgeColor = 'none';
% for a = 1:numel(recareas)
%     denom = Nneurons_fixedgazeagg.(recarea{a}).(whichvisblock).(neudenom);
%     tempvec1 = Nneurons_fixedgazeagg.(recarea{a}).(whichvisblock).indin./denom;
%     tempvec2 = Nneurons_fixedgazeagg.(recarea{a}).(whichvisblock).sigBKnotsigBI./denom;
%     plot(a-.2, tempvec1, 'o', 'MarkerFaceColor', b(1).FaceColor, 'MarkerEdgeColor', 'none')
%     plot(a+.2, tempvec2, 'o', 'MarkerFaceColor', b(2).FaceColor, 'MarkerEdgeColor', 'none')
% end
% ylim(yl)
set(gca, 'FontSize', fs, 'XTick', 1:numel(recareas), 'XTickLabel', recareas, 'XTickLabelRotation', 45)
ylabel('Proportion Neurons')
if ii==2 && jj==2
legend({'Pac-In', 'PkwBK<0.05 & PkwBI>0.05'}, 'Location', 'SouthEast')
end
title(sprintf('out of %s neurons', neudenom))
    end
end

%%
load('D:\OpenScopeData\000248v221123\postprocessed\psthavg_fixedgazeagg.mat')

for b= 1:numel(visblocks)
if ~isequal(vistrialtypes_fixedgazeagg(1).(visblocks{b}), unique(cat(2,vistrialtypes_fixedgazeagg.(visblocks{b}))', 'rows')')
    error([visblocks{b} ' vis trial types NOT the same across sessions'])
end
end

for b = 1:numel(visblocks)
    Ronavg_fixedgazeagg.(visblocks{b}).all = [];
    for iprobe = 1:numel(probes)
        Ronavg_fixedgazeagg.(visblocks{b}).all = cat(2, Ronavg_fixedgazeagg.(visblocks{b}).all, Ronavg_fixedgazeagg.(visblocks{b}).(probes{iprobe}));
    end
end

%% orientation tuning of IC-encoders
% whichprobe = 'all';
% neuinarea = true(size(neuprobeagg));

whichprobe = 'C';
iprobe = find(strcmp(probes, whichprobe));
% neuinarea = strcmp(neuloc_fixedgazeagg{iprobe}, 'VISp2/3');
% neuinarea = contains(neuloc_fixedgazeagg{iprobe}, 'VISp');
neuinarea = neuctx_fixedgazeagg{iprobe}==1;

sigori = ori4params_fixedgazeagg.(whichprobe).Pkw_ori4<0.05;

% 0 cfg0tt111 45 cfg1tt111 90 cfg0tt106 135 cfg1tt106
RonwICavgagg = [Ronavg_fixedgazeagg.ICwcfg0_presentations.(whichprobe)(ICtrialtypes==111,:); ...
    Ronavg_fixedgazeagg.ICwcfg1_presentations.(whichprobe)(ICtrialtypes==111,:); ...
    Ronavg_fixedgazeagg.ICwcfg0_presentations.(whichprobe)(ICtrialtypes==106,:); ...
    Ronavg_fixedgazeagg.ICwcfg1_presentations.(whichprobe)(ICtrialtypes==106,:)]';
RonkICavgagg = [Ronavg_fixedgazeagg.ICkcfg0_presentations.(whichprobe)(ICtrialtypes==111,:); ...
    Ronavg_fixedgazeagg.ICkcfg1_presentations.(whichprobe)(ICtrialtypes==111,:); ...
    Ronavg_fixedgazeagg.ICkcfg0_presentations.(whichprobe)(ICtrialtypes==106,:); ...
    Ronavg_fixedgazeagg.ICkcfg1_presentations.(whichprobe)(ICtrialtypes==106,:)]';
[~, prefwIC] = max(RonwICavgagg, [], 2);
[~, prefkIC] = max(RonkICavgagg, [], 2);

RonwICRCavgagg = [Ronavg_fixedgazeagg.ICwcfg0_presentations.(whichprobe)(ICtrialtypes==111,:); ...
    Ronavg_fixedgazeagg.ICwcfg1_presentations.(whichprobe)(ICtrialtypes==111,:); ...
    Ronavg_fixedgazeagg.ICwcfg0_presentations.(whichprobe)(ICtrialtypes==106,:); ...
    Ronavg_fixedgazeagg.ICwcfg1_presentations.(whichprobe)(ICtrialtypes==106,:); ...
    Ronavg_fixedgazeagg.ICwcfg0_presentations.(whichprobe)(ICtrialtypes==110,:); ...
    Ronavg_fixedgazeagg.ICwcfg1_presentations.(whichprobe)(ICtrialtypes==110,:); ...
    Ronavg_fixedgazeagg.ICwcfg0_presentations.(whichprobe)(ICtrialtypes==107,:); ...
    Ronavg_fixedgazeagg.ICwcfg1_presentations.(whichprobe)(ICtrialtypes==107,:)]';
RonkICRCavgagg = [Ronavg_fixedgazeagg.ICkcfg0_presentations.(whichprobe)(ICtrialtypes==111,:); ...
    Ronavg_fixedgazeagg.ICkcfg1_presentations.(whichprobe)(ICtrialtypes==111,:); ...
    Ronavg_fixedgazeagg.ICkcfg0_presentations.(whichprobe)(ICtrialtypes==106,:); ...
    Ronavg_fixedgazeagg.ICkcfg1_presentations.(whichprobe)(ICtrialtypes==106,:); ...
    Ronavg_fixedgazeagg.ICkcfg0_presentations.(whichprobe)(ICtrialtypes==110,:); ...
    Ronavg_fixedgazeagg.ICkcfg1_presentations.(whichprobe)(ICtrialtypes==110,:); ...
    Ronavg_fixedgazeagg.ICkcfg0_presentations.(whichprobe)(ICtrialtypes==107,:); ...
    Ronavg_fixedgazeagg.ICkcfg1_presentations.(whichprobe)(ICtrialtypes==107,:)]';
[~, prefwICRC] = max(RonwICRCavgagg, [], 2);
[~, prefkICRC] = max(RonkICRCavgagg, [], 2);

Ronwindinagg = cat(3,Ronavg_fixedgazeagg.ICwcfg0_presentations.(whichprobe)(ismember(ICtrialtypes,[1302 1304]),:), ...
    Ronavg_fixedgazeagg.ICwcfg1_presentations.(whichprobe)(ismember(ICtrialtypes,[1302 1304]),:), ...
    Ronavg_fixedgazeagg.ICwcfg0_presentations.(whichprobe)(ismember(ICtrialtypes,[1301 1303]),:), ...
    Ronavg_fixedgazeagg.ICwcfg1_presentations.(whichprobe)(ismember(ICtrialtypes,[1301 1303]),:));
Ronwindinagg = squeeze(mean(Ronwindinagg,1));
Ronkindinagg = cat(3,Ronavg_fixedgazeagg.ICkcfg0_presentations.(whichprobe)(ismember(ICtrialtypes,[1302 1304]),:), ...
    Ronavg_fixedgazeagg.ICkcfg1_presentations.(whichprobe)(ismember(ICtrialtypes,[1302 1304]),:), ...
    Ronavg_fixedgazeagg.ICkcfg0_presentations.(whichprobe)(ismember(ICtrialtypes,[1301 1303]),:), ...
    Ronavg_fixedgazeagg.ICkcfg1_presentations.(whichprobe)(ismember(ICtrialtypes,[1301 1303]),:));
Ronkindinagg = squeeze(mean(Ronkindinagg,1));
[~, prefwindin] = max(Ronwindinagg, [], 2);
[~, prefkindin] = max(Ronkindinagg, [], 2);

wICenc = neuinarea & (ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).ICencoder==1 | ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).ICencoder==1);
kICenc = neuinarea & (ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).ICencoder==1 | ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).ICencoder==1);

wRCenc = neuinarea & (ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).RCencoder==1 | ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).RCencoder==1);
kRCenc = neuinarea & (ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).RCencoder==1 | ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).RCencoder==1);

wemerge = neuinarea & ( (ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).PkwBK<0.05 & ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).PkwBI>0.05) ...
    | (ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).PkwBK<0.05 & ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).PkwBI>0.05) );
kemerge = neuinarea & ( (ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).PkwBK<0.05 & ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).PkwBI>0.05) ...
    | (ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).PkwBK<0.05 & ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).PkwBI>0.05) );

wsigkw = neuinarea & ( (ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).PkwBK<0.05 & ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).PkwBI<0.05) ...
    | (ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).PkwBK<0.05 & ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).PkwBI>0.05) );
ksigkw = neuinarea & ( (ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).PkwBK<0.05 & ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).PkwBI<0.05) ...
    | (ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).PkwBK<0.05 & ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).PkwBI>0.05) );

neuoi = wICenc & sigori;
figure; histogram2(prefkICRC(neuoi), ori4params_fixedgazeagg.(whichprobe).prefiori4(neuoi), 'displaystyle', 'tile')
xlabel('Pref wIC Ori')
ylabel('Pref Grating Ori')

figure('Position', [0 0 2000 600])
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'Fixed Gaze Trials: Significantly Tuned Neurons (rows sum to 1)', 'edgecolor', 'none')
for isp = 1:5
    subplot(2,5,isp)
    hold all
    switch isp
        case 1
            neuoi = neuinarea & sigori;
            neudesc = ['All Neurons in Probe-' whichprobe];
        case 2
            neuoi = wsigkw & sigori;
            neudesc = 'PkwBK<0.05 & PkwBI<0.05';
        case 3
            neuoi = wemerge & sigori;
            neudesc = 'PkwBK<0.05 & PkwBI>0.05';
        case 4
            neuoi = wICenc & sigori;
            neudesc = 'IC-encoder';
        case 5
            neuoi = wRCenc & sigori;
            neudesc = 'RC-encoder';
    end
hc = histcounts2(prefwIC(neuoi), ori4params_fixedgazeagg.(whichprobe).prefiori4(neuoi), 0.5:4.5, 0.5:4.5);
imagesc( (hc./sum(hc,2))' ); colorbar
axis([0.5 4.5 0.5 4.5])
ylabel('Pref Grating Ori')
xlabel('Pref wIC')
title(sprintf('%s N=%d', neudesc, nnz(neuoi)))

    subplot(2,5,isp+5)
    hold all
    switch isp
        case 1
            neuoi = neuinarea & sigori;
            neudesc = ['All Neurons in Probe-' whichprobe];
        case 2
            neuoi = ksigkw & sigori;
            neudesc = 'PkwBK<0.05 & PkwBI<0.05';
        case 3
            neuoi = kemerge & sigori;
            neudesc = 'PkwBK<0.05 & PkwBI>0.05';
        case 4
            neuoi = kICenc & sigori;
            neudesc = 'IC-encoder';
        case 5
            neuoi = kRCenc & sigori;
            neudesc = 'RC-encoder';
    end
hc = histcounts2(prefkIC(neuoi), ori4params_fixedgazeagg.(whichprobe).prefiori4(neuoi), 0.5:4.5, 0.5:4.5);
imagesc( (hc./sum(hc,2))' ); colorbar
axis([0.5 4.5 0.5 4.5])
ylabel('Pref Grating Ori')
xlabel('Pref kIC')
title(sprintf('%s N=%d', neudesc, nnz(neuoi)))
end

figure('Position', [0 0 2000 600])
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'Fixed Gaze Trials: Significantly Tuned Neurons (rows sum to 1)', 'edgecolor', 'none')
for isp = 1:5
    subplot(2,5,isp)
    hold all
    switch isp
        case 1
            neuoi = neuinarea & sigori;
            neudesc = ['All Neurons in Probe-' whichprobe];
        case 2
            neuoi = wsigkw & sigori;
            neudesc = 'PkwBK<0.05 & PkwBI<0.05';
        case 3
            neuoi = wemerge & sigori;
            neudesc = 'PkwBK<0.05 & PkwBI>0.05';
        case 4
            neuoi = wICenc & sigori;
            neudesc = 'IC-encoder';
        case 5
            neuoi = wRCenc & sigori;
            neudesc = 'RC-encoder';
    end
hc = histcounts2(prefwindin(neuoi), ori4params_fixedgazeagg.(whichprobe).prefiori4(neuoi), 0.5:4.5, 0.5:4.5);
imagesc( (hc./sum(hc,2))' ); colorbar
axis([0.5 4.5 0.5 4.5])
ylabel('Pref Grating Ori')
xlabel('Pref wPac-In')
title(sprintf('%s N=%d', neudesc, nnz(neuoi)))

    subplot(2,5,isp+5)
    hold all
    switch isp
        case 1
            neuoi = neuinarea & sigori;
            neudesc = ['All Neurons in Probe-' whichprobe];
        case 2
            neuoi = ksigkw & sigori;
            neudesc = 'PkwBK<0.05 & PkwBI<0.05';
        case 3
            neuoi = kemerge & sigori;
            neudesc = 'PkwBK<0.05 & PkwBI>0.05';
        case 4
            neuoi = kICenc & sigori;
            neudesc = 'IC-encoder';
        case 5
            neuoi = kRCenc & sigori;
            neudesc = 'RC-encoder';
    end
hc = histcounts2(prefkindin(neuoi), ori4params_fixedgazeagg.(whichprobe).prefiori4(neuoi), 0.5:4.5, 0.5:4.5);
imagesc( (hc./sum(hc,2))' ); colorbar
axis([0.5 4.5 0.5 4.5])
ylabel('Pref Grating Ori')
xlabel('Pref kPac-In')
title(sprintf('%s N=%d', neudesc, nnz(neuoi)))
end

figure('Position', [0 0 2000 600])
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'Fixed Gaze Trials: Significantly Tuned Neurons (rows sum to 1)', 'edgecolor', 'none')
for isp = 1:5
    subplot(2,5,isp)
    hold all
    switch isp
        case 1
            neuoi = neuinarea & sigori;
            neudesc = ['All Neurons in Probe-' whichprobe];
        case 2
            neuoi = wsigkw & sigori;
            neudesc = 'PkwBK<0.05 & PkwBI<0.05';
        case 3
            neuoi = wemerge & sigori;
            neudesc = 'PkwBK<0.05 & PkwBI>0.05';
        case 4
            neuoi = wICenc & sigori;
            neudesc = 'IC-encoder';
        case 5
            neuoi = wRCenc & sigori;
            neudesc = 'RC-encoder';
    end
hc = histcounts2(prefwICRC(neuoi), ori4params_fixedgazeagg.(whichprobe).prefiori4(neuoi), 0.5:8.5, 0.5:4.5);
imagesc( (hc./sum(hc,2))' ); colorbar
plot([4.5 4.5], [0.5 4.5], 'w-', 'Linewidth', 2)
axis([0.5 8.5 0.5 4.5])
ylabel('Pref Grating Ori')
xlabel('Pref wIC/RC')
title(sprintf('%s N=%d', neudesc, nnz(neuoi)))

    subplot(2,5,isp+5)
    hold all
    switch isp
        case 1
            neuoi = neuinarea & sigori;
            neudesc = ['All Neurons in Probe-' whichprobe];
        case 2
            neuoi = ksigkw & sigori;
            neudesc = 'PkwBK<0.05 & PkwBI<0.05';
        case 3
            neuoi = kemerge & sigori;
            neudesc = 'PkwBK<0.05 & PkwBI>0.05';
        case 4
            neuoi = kICenc & sigori;
            neudesc = 'IC-encoder';
        case 5
            neuoi = kRCenc & sigori;
            neudesc = 'RC-encoder';
    end
hc = histcounts2(prefkICRC(neuoi), ori4params_fixedgazeagg.(whichprobe).prefiori4(neuoi), 0.5:8.5, 0.5:4.5);
imagesc( (hc./sum(hc,2))' ); colorbar
plot([4.5 4.5], [0.5 4.5], 'w-', 'Linewidth', 2)
axis([0.5 8.5 0.5 4.5])
ylabel('Pref Grating Ori')
xlabel('Pref kIC/RC')
title(sprintf('%s N=%d', neudesc, nnz(neuoi)))
end

% % does OP correlate with SP_ICvsRC for "emergent" cells? -- no, not at all obvious
% neuoi = wemerge & (ori4params_fixedgazeagg.(whichprobe).prefiori4==1 | ori4params_fixedgazeagg.(whichprobe).prefiori4==3);
% corr(ori4params_fixedgazeagg.(whichprobe).OP4678(neuoi), ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).SP_ICvsRC(neuoi), 'type', 'spearman')
% figure; plot(ori4params_fixedgazeagg.(whichprobe).OP4678(neuoi), ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).SP_ICvsRC(neuoi), 'o')
% figure; histogram2(ori4params_fixedgazeagg.(whichprobe).OP4678(neuoi), ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).SP_ICvsRC(neuoi), 'displaystyle', 'tile')
% 
% neuoi = wemerge & (ori4params_fixedgazeagg.(whichprobe).prefiori4==2 | ori4params_fixedgazeagg.(whichprobe).prefiori4==4);
% corr(ori4params_fixedgazeagg.(whichprobe).OP4678(neuoi), ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).SP_ICvsRC(neuoi), 'type', 'spearman')
% figure; plot(ori4params_fixedgazeagg.(whichprobe).OP4678(neuoi), ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).SP_ICvsRC(neuoi), 'o')
% figure; histogram2(ori4params_fixedgazeagg.(whichprobe).OP4678(neuoi), ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).SP_ICvsRC(neuoi), 'displaystyle', 'tile')

%% plot each IC block
kerwinhalf = 5; kersigma = 2;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));

whichprobe = 'C';

iprobe = find(strcmp(probes, whichprobe));
% neuinarea = strcmp(neuloc_fixedgazeagg{iprobe}, 'VISp2/3');
% neuinarea = contains(neuloc_fixedgazeagg{iprobe}, 'VISp');
neuinarea = neuctx_fixedgazeagg{iprobe}==1;
% neuinarea = true(size(neuctx_fixedgazeagg{iprobe}));

tt2p = [106 107 110 111 1105 1109]';
tt2p = [106 107 110 111]';
ttcol = [0 0.5 0;
    0.75 0.5 0;
    1 0.75 0.25;
    0.25 0.75 0.25;
    0 0 0.5;
    0.25 0.25 0.75];
fs=10;

for ineu = 1:4
switch ineu
    case 1
        neudesc = 'All Neurons';
        neuoi = true(size(RFCI_fixedgazeagg.(whichprobe).RFindclassic));
    case 2
        neudesc = 'ctrCRF';
        neuoi = RFCI_fixedgazeagg.(whichprobe).RFindclassic==1 & RFCI_fixedgazeagg.(whichprobe).Pkw_rfclassic<0.05;
    case 3
        neudesc = 'ctrCRFsigexcl';
        neuoi = RFCI_fixedgazeagg.(whichprobe).RFindclassic==1 & RFCI_fixedgazeagg.(whichprobe).RFsigexclclassic==1;
    case 4
        neudesc = 'ctrCRFexclsig';
        neuoi = RFCI_fixedgazeagg.(whichprobe).RFindclassic==1 & RFCI_fixedgazeagg.(whichprobe).RFexclsigclassic==1;
end
neuoi = neuinarea & neuoi; 
figure%('Position', [100 0 1800 1200])
annotation('textbox', [0.1 0.92 0.9 0.1], 'string', ['Fixed Gaze Trials: ' visareas{iprobe} ' ' neudesc ' N=' num2str(nnz(neuoi))], 'edgecolor', 'none', 'fontsize', fs, 'interpreter', 'none')
for b= 1:numel(ICblocks)
    whichvisblock = ICblocks{b};
    subplot(2,2, b)
    hold all
    for ii = 1:numel(tt2p)
        typi = ICtrialtypes==tt2p(ii);
        temppsth = convn(squeeze(psthavg_fixedgazeagg.(whichvisblock).(whichprobe)(:,typi,neuoi)), kergauss, 'same');
        plot(psthtli/1000, squeeze(nanmean(temppsth,2)), '-', 'Color', ttcol(ii,:), 'LineWidth', 1)
        %shadedErrorBar(psthtli/1000, squeeze(mean(temppsth,2)), squeeze(std(temppsth,0,2)/sqrt(nnz(neuoi))), {'-', 'Color', ttcol(ii,:), 'LineWidth', 1}, 1)
    end
    typoi = ismember(ICtrialtypes, 1301:1308);
    temppsth = convn(squeeze( mean(psthavg_fixedgazeagg.(whichvisblock).(whichprobe)(:,typoi,neuoi),2) ), kergauss, 'same');
    plot(psthtli/1000, squeeze(nanmean(temppsth,2)), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
    xlim([-.1 .500])
    %     ylim([0 25])
    set(gca, 'FontSize', fs)
    xlabel('Time (s)')
    ylabel('Rate (Hz)')
    title(whichvisblock, 'interpreter', 'none')
end
end

figure%('Position', [100 0 1800 1200])
annotation('textbox', [0.1 0.92 0.9 0.1], 'string', ['Fixed Gaze Trials: ' visareas{iprobe} ' PkwBK<0.05 & PkwBI>0.05'], 'edgecolor', 'none', 'fontsize', fs, 'interpreter', 'none')
for b= 1:numel(ICblocks)
    whichvisblock = ICblocks{b};
    neuoi = neuinarea & ICsig_fixedgazeagg.(whichvisblock).(whichprobe).PkwBK<0.05 & ICsig_fixedgazeagg.(whichvisblock).(whichprobe).PkwBI>0.05;
    subplot(2,2, b)
    hold all
    for ii = 1:numel(tt2p)
        typi = ICtrialtypes==tt2p(ii);
        temppsth = convn(squeeze(psthavg_fixedgazeagg.(whichvisblock).(whichprobe)(:,typi,neuoi)), kergauss, 'same');
        plot(psthtli/1000, squeeze(mean(temppsth,2)), '-', 'Color', ttcol(ii,:), 'LineWidth', 1)
        %shadedErrorBar(psthtli/1000, squeeze(mean(temppsth,2)), squeeze(std(temppsth,0,2)/sqrt(nnz(neuoi))), {'-', 'Color', ttcol(ii,:), 'LineWidth', 1}, 1)
    end
    typoi = ismember(ICtrialtypes, 1301:1308);
    temppsth = convn(squeeze( mean(psthavg_fixedgazeagg.(whichvisblock).(whichprobe)(:,typoi,neuoi),2) ), kergauss, 'same');
    plot(psthtli/1000, squeeze(nanmean(temppsth,2)), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
    xlim([-.1 .500])
    %     ylim([0 25])
    set(gca, 'FontSize', fs)
    xlabel('Time (s)')
    ylabel('Rate (Hz)')
    title(sprintf('%s N=%d', whichvisblock, nnz(neuoi)), 'interpreter', 'none')
end

%% are there emergent responses? for neurons tuned to each orientation, see whether IC evokes responses above and beyond the sum of ind-in trials
kerwinhalf = 10; kersigma = 5;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));

orivec = 0:45:180-1;
oricol = linspace(1,0,4)' * [1 0 0];

neudesc = 'PkwBK<0.05 & PkwBI>0.05';

figure('Position', [100 0 1600 900])
annotation('textbox', [0.1 0.9 0.9 0.1], 'string', ['Fixed Gaze Trials: ' visareas{iprobe} ' ' neudesc], 'edgecolor', 'none', 'fontsize', fs, 'interpreter', 'none')
for iori= 1:4
    neuori = neuinarea & (ori4params_fixedgazeagg.(whichprobe).prefiori4==iori & ori4params_fixedgazeagg.(whichprobe).Pkw_ori4<0.05);
    for jj = 1:2
    subplot(2,4, (jj-1)*4+iori)
    hold all
    leg = cell(4,1);
    for ii = 1:4
        switch ii
            case 1
                whichvisblock = 'ICwcfg0_presentations';
                typi = ICtrialtypes==111;
                typoi = ismember(ICtrialtypes, [1302 1304]);
            case 2
                whichvisblock = 'ICwcfg1_presentations';
                typi = ICtrialtypes==111;
                typoi = ismember(ICtrialtypes, [1302 1304]);
            case 3
                whichvisblock = 'ICwcfg0_presentations';
                typi = ICtrialtypes==106;
                typoi = ismember(ICtrialtypes, [1301 1303]);
            case 4
                whichvisblock = 'ICwcfg1_presentations';
                typi = ICtrialtypes==106;
                typoi = ismember(ICtrialtypes, [1301 1303]);
        end
        switch neudesc
            case 'IC-encoder'
                neuoi = neuori & ICsig_fixedgazeagg.(whichvisblock).(whichprobe).ICencoder==1;
            case 'PkwBK<0.05 & PkwBI>0.05'
                neuoi = neuori & ICsig_fixedgazeagg.(whichvisblock).(whichprobe).PkwBK<0.05 & ICsig_fixedgazeagg.(whichvisblock).(whichprobe).PkwBI>0.05;
            case 'All Neurons'
                neuoi = neuori;
            case 'ctrCRF'
                neuoi = neuori & RFCI_fixedgazeagg.(whichprobe).RFindclassic==1 & RFCI_fixedgazeagg.(whichprobe).Pkw_rfclassic<0.05;
            case 'ctrCRFsigexcl'
                neuoi = neuori & RFCI_fixedgazeagg.(whichprobe).RFindclassic==1 & RFCI_fixedgazeagg.(whichprobe).RFsigexclclassic==1;
            case 'ctrCRFexclsig'
                neuoi = neuori & RFCI_fixedgazeagg.(whichprobe).RFindclassic==1 & RFCI_fixedgazeagg.(whichprobe).RFexclsigclassic==1;
        end
        if jj==1
        temppsth = convn(squeeze(psthavg_fixedgazeagg.(whichvisblock).(whichprobe)(:,typi,neuoi)), kergauss, 'same');
        plot(psthtli/1000, squeeze(mean(temppsth,2)), '-', 'Color', oricol(ii,:), 'LineWidth', 0.5)
        %shadedErrorBar(psthtli/1000, squeeze(mean(temppsth,2)), squeeze(std(temppsth,0,2)/sqrt(nnz(neuoi))), {'-', 'Color', ttcol(ii,:), 'LineWidth', 1}, 1)
        leg{ii} = sprintf('IC%d', orivec(ii));
        else
        temppsth = convn( psthavg_fixedgazeagg.(whichvisblock).(whichprobe)(:,typoi,neuoi), kergauss, 'same');
        plot(psthtli/1000, squeeze(nanmean(temppsth,[2 3])), '-', 'Color', oricol(ii,:), 'LineWidth', 0.5)
        leg{ii} = sprintf('indin%d', orivec(ii));
        end
    end
    legend(leg)
    xlim([-.1 .500])
    ylim([0 20])
    set(gca, 'FontSize', fs)
    xlabel('Time (s)')
    ylabel('Rate (Hz)')
    title(sprintf('Neurons Tuned to %ddeg N=%d', orivec(iori), nnz(neuoi)), 'interpreter', 'none')    
    end
end

%% are there emergent responses? for neurons tuned to each orientation, see whether IC evokes responses above and beyond the sum of ind-in trials

inNout = false; 

figure
annotation('textbox', [0.1 0.9 0.9 0.1], 'string', ['Fixed Gaze Trials: ' visareas{iprobe} ' Blank-Subtracted Evoked Responses'], 'edgecolor', 'none', 'fontsize', fs, 'interpreter', 'none')
for jj = 1:2
    switch jj
        case 1
            indcol = 'w';
        case 2
            indcol = 'k';
    end
for ii = 1:4
    switch ii
        case 1
            whichvisblock = ['IC' indcol 'cfg0_presentations'];
            typi = ICtrialtypes==111;
            if inNout
                typoi = ismember(ICtrialtypes, [1302 1304 1305 1307]);
            else
                typoi = ismember(ICtrialtypes, [1302 1304]);
            end
        case 2
            whichvisblock = ['IC' indcol 'cfg1_presentations'];
            typi = ICtrialtypes==111;
            if inNout
                typoi = ismember(ICtrialtypes, [1302 1304 1305 1307]);
            else
                typoi = ismember(ICtrialtypes, [1302 1304]);
            end
        case 3
            whichvisblock = ['IC' indcol 'cfg0_presentations'];
            typi = ICtrialtypes==106;
            if inNout
                typoi = ismember(ICtrialtypes, [1301 1303 1306 1308]);
            else
                typoi = ismember(ICtrialtypes, [1301 1303]);
            end
        case 4
            whichvisblock = ['IC' indcol 'cfg1_presentations'];
            typi = ICtrialtypes==106;
            if inNout
                typoi = ismember(ICtrialtypes, [1301 1303 1306 1308]);
            else
                typoi = ismember(ICtrialtypes, [1301 1303]);
            end
    end
    neuoi = neuinarea & ICsig_fixedgazeagg.(whichvisblock).(whichprobe).PkwBK<0.05;

    Rblank = Ronavg_fixedgazeagg.(whichvisblock).(whichprobe)(ICtrialtypes==0,:);
    RIC = Ronavg_fixedgazeagg.(whichvisblock).(whichprobe)(typi,:)-Rblank;
    Rindin = sum(Ronavg_fixedgazeagg.(whichvisblock).(whichprobe)(typoi,:)-Rblank,1);
    subplot(2,4,4*(jj-1)+ii)
    hold all
    scatter(RIC(neuoi), Rindin(neuoi), 'o', 'MarkerFaceAlpha', 0.1, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k')
    xl = xlim; yl = ylim; al = [min([xl yl]) max([xl yl])];
    plot(al, al, 'r--')
    axis([al al])
    xlabel('IC trials')
    if inNout
    ylabel('Sum of Pacmans')
    else
    ylabel('Sum of Pac-Ins')
    end
    title(sprintf('%d-deg %sIC', orivec(ii), indcol))
end
end

figure
annotation('textbox', [0.1 0.9 0.9 0.1], 'string', ['Fixed Gaze Trials: ' visareas{iprobe} ' Blank-Subtracted Evoked Responses'], 'edgecolor', 'none', 'fontsize', fs, 'interpreter', 'none')
for jj = 1:2
    switch jj
        case 1
            indcol = 'w';
        case 2
            indcol = 'k';
    end
for ii = 1:4
    switch ii
        case 1
            whichvisblock = ['IC' indcol 'cfg0_presentations'];
            typi = ICtrialtypes==111;
            if inNout
                typoi = ismember(ICtrialtypes, [1302 1304 1305 1307]);
            else
                typoi = ismember(ICtrialtypes, [1302 1304]);
            end
        case 2
            whichvisblock = ['IC' indcol 'cfg1_presentations'];
            typi = ICtrialtypes==111;
            if inNout
                typoi = ismember(ICtrialtypes, [1302 1304 1305 1307]);
            else
                typoi = ismember(ICtrialtypes, [1302 1304]);
            end
        case 3
            whichvisblock = ['IC' indcol 'cfg0_presentations'];
            typi = ICtrialtypes==106;
            if inNout
                typoi = ismember(ICtrialtypes, [1301 1303 1306 1308]);
            else
                typoi = ismember(ICtrialtypes, [1301 1303]);
            end
        case 4
            whichvisblock = ['IC' indcol 'cfg1_presentations'];
            typi = ICtrialtypes==106;
            if inNout
                typoi = ismember(ICtrialtypes, [1301 1303 1306 1308]);
            else
                typoi = ismember(ICtrialtypes, [1301 1303]);
            end
    end

    Rblank = Ronavg_fixedgazeagg.(whichvisblock).(whichprobe)(ICtrialtypes==0,:);
    RIC = Ronavg_fixedgazeagg.(whichvisblock).(whichprobe)(typi,:)-Rblank;
    Rindin = sum(Ronavg_fixedgazeagg.(whichvisblock).(whichprobe)(typoi,:)-Rblank,1);
    subplot(2,4,4*(jj-1)+ii)
    hold all
    neuoi = neuinarea & ICsig_fixedgazeagg.(whichvisblock).(whichprobe).PkwBK<0.05;
    histogram(RIC(neuoi)-Rindin(neuoi))
    neuoi = neuinarea & ICsig_fixedgazeagg.(whichvisblock).(whichprobe).PkwBK<0.05 & ICsig_fixedgazeagg.(whichvisblock).(whichprobe).PkwBI>0.05;
    %neuoi = neuinarea & ICsig_fixedgazeagg.(whichvisblock).(whichprobe).ICencoder==1;
    histogram(RIC(neuoi)-Rindin(neuoi))
    if inNout
    xlabel('IC - (Sum of Pacmans)')
    else
    xlabel('IC - (Sum of Pac-Ins)')
    end
    ylabel('# Neurons')
    title(sprintf('%d-deg %sIC', orivec(ii), indcol))
end
end

%% compare population evoked response of IC vs RC trials for ctrCRFsigexcl neurons
whichprobe = 'C';
iprobe = find(strcmp(probes, whichprobe));
% neuinarea = strcmp(neuloc_fixedgazeagg{iprobe}, 'VISp2/3');
% neuinarea = contains(neuloc_fixedgazeagg{iprobe}, 'VISp');
neuinarea = neuctx_fixedgazeagg{iprobe}==1;
% neuinarea = true(size(neuctx_fixedgazeagg{iprobe}));

neuoi = neuinarea & RFCI_fixedgazeagg.(whichprobe).RFindclassic==1 & RFCI_fixedgazeagg.(whichprobe).RFsigexclclassic==1;
%neuoi = neuinarea & RFCI_fixedgazeagg.(whichprobe).RFindclassic==1 & RFCI_fixedgazeagg.(whichprobe).RFexclsigclassic==1;

typi = ICtrialtypes;
% psthavg_fixedgazeagg.ICwcfg1_presentations.(whichprobe).
