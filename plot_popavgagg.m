load(['G:\My Drive\DATA\ICexpts_submission22\openscope_popavg_agg.mat'])
Nsessions = numel(nwbsessions);

probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
visind = [6 5 1 2 4 3];
visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
    'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
ICblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations'};

%%
neuprobeagg = [];
for b = 1:numel(ICblocks)
    for f= 1:numel(ICsigfields)
        ICsigagg.(ICblocks{b}).all.(ICsigfields{f}) = [];
    end
end
% RFCIfields = fieldnames(RFCI);
for f= 1:numel(RFCIfields)
    RFCIagg.all.(RFCIfields{f}) = [];
end
for f= 1:numel(RFCIspinfields)
    RFCIspinagg.all.(RFCIspinfields{f}) = [];
end
% size vector [0, 4, 8, 16, 32, 64 ]
for f= 1:numel(sizeCIfields)
    sizeCIagg.all.(sizeCIfields{f}) = [];
end
for f= 1:numel(oriparamsfields)
    oriparamsagg.all.(oriparamsfields{f}) = [];
end
for f= 1:numel(ori4paramsfields)
    ori4paramsagg.all.(ori4paramsfields{f}) = [];
end

for iprobe = 1:numel(probes)
    neuprobeagg = cat(1, neuprobeagg, iprobe*ones(size(neulocagg{iprobe})));
    
for b = 1:numel(ICblocks)
    for f= 1:numel(ICsigfields)
        ICsigagg.(ICblocks{b}).all.(ICsigfields{f}) = cat(1, ICsigagg.(ICblocks{b}).all.(ICsigfields{f}), ICsigagg.(ICblocks{b}).(probes{iprobe}).(ICsigfields{f}) );
    end
end

% RFCIfields = fieldnames(RFCI);
for f= 1:numel(RFCIfields)
    RFCIagg.all.(RFCIfields{f}) = cat(1, RFCIagg.all.(RFCIfields{f}), RFCIagg.(probes{iprobe}).(RFCIfields{f}) );
end

for f= 1:numel(RFCIspinfields)
    RFCIspinagg.all.(RFCIspinfields{f}) = cat(1, RFCIspinagg.all.(RFCIspinfields{f}), RFCIspinagg.(probes{iprobe}).(RFCIspinfields{f}) );
end

% size vector [0, 4, 8, 16, 32, 64 ]
for f= 1:numel(sizeCIfields)
    sizeCIagg.all.(sizeCIfields{f}) = cat(1, sizeCIagg.all.(sizeCIfields{f}), sizeCIagg.(probes{iprobe}).(sizeCIfields{f}) );
end

for f= 1:numel(oriparamsfields)
    oriparamsagg.all.(oriparamsfields{f}) = cat(1, oriparamsagg.all.(oriparamsfields{f}), oriparamsagg.(probes{iprobe}).(oriparamsfields{f}) );
end

for f= 1:numel(ori4paramsfields)
    ori4paramsagg.all.(ori4paramsfields{f}) = cat(1, ori4paramsagg.all.(ori4paramsfields{f}), ori4paramsagg.(probes{iprobe}).(ori4paramsfields{f}) );
end
end


%% proportion of neurons in each area
% compare areas
recarealabels = {'LGd', 'LP', 'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam', ...
    'VISp1', 'VISp2/3', 'VISp4', 'VISp5', 'VISp6a', 'VISp6b'};
recareas = {'LGN', 'LP', 'V1', 'LM', 'RL', 'AL', 'PM', 'AM', ...
    'V1L1', 'V1L23', 'V1L4', 'V1L5', 'V1L6a', 'V1L6b'};

% VISp NOT VISpm

% compare layers in V1
neutypes = {'all', 'ICencoder', 'RCencoder', 'inducerencoder', 'indin', 'sigBK', 'sigBKnotsigBI'};
Nneuronsacc = struct();
Nneuronsagg = struct();
for a = 1:numel(recarealabels)
    if strcmp(recarealabels{a}, 'VISp')
        neuinarea = contains(cat(1,neulocagg{:}), 'VISp') & ~contains(cat(1,neulocagg{:}), 'VISpm');
    else
        neuinarea = contains(cat(1,neulocagg{:}), recarealabels{a});
    end
for b= 1:numel(ICblocks)
    for ineu = 1:numel(neutypes)
        switch neutypes{ineu}
            case 'all'
                neuoi = true(size(ICsigagg.(ICblocks{b}).all.ICencoder));
            case 'ICencoder'
                neuoi = ICsigagg.(ICblocks{b}).all.ICencoder==1;
            case 'RCencoder'
                neuoi = ICsigagg.(ICblocks{b}).all.RCencoder==1;
            case 'inducerencoder'
                neuoi = ICsigagg.(ICblocks{b}).all.inducerencoder==1;
            case 'indin'
                neuoi = ICsigagg.(ICblocks{b}).all.indin1==1 | ICsigagg.(ICblocks{b}).all.indin2==1 ...
                    | ICsigagg.(ICblocks{b}).all.indin3==1 | ICsigagg.(ICblocks{b}).all.indin4==1;
            case 'sigBK'
                neuoi = ICsigagg.(ICblocks{b}).all.PkwBK<0.05;
            case 'sigBKnotsigBI'
                neuoi = ICsigagg.(ICblocks{b}).all.PkwBK<0.05 & ICsigagg.(ICblocks{b}).all.PkwBI>=0.05;
        end
        neuoi = neuinarea & neuoi;
        Nneuronsacc.(recareas{a}).(ICblocks{b}).(neutypes{ineu}) = nnz(neuoi);
        Nneuronsagg.(recareas{a}).(ICblocks{b}).(neutypes{ineu}) = zeros(Nsessions, 1);
        for ises = 1:Nsessions
            sesneu = cat(1,sesneuagg{:})==ises;
            Nneuronsagg.(recareas{a}).(ICblocks{b}).(neutypes{ineu})(ises) = nnz(neuoi & sesneu);
        end
    end
end
end

%% compare proportion IC-/RC-encoders across sessions
whichvisblock = 'ICwcfg1_presentations';

pltmeanoragg = 2; % 1 is average proportion of neurons across sessions, 2 is aggregate across sessions

fs = 10;
figure
annotation('textbox', [0.1 0.92 0.9 0.1], 'string', ['All Trials: ' whichvisblock], 'edgecolor', 'none', 'interpreter', 'none')
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
    denom = Nneuronsagg.(recareas{a}).(whichvisblock).(neudenom);
    tempmat(a,1) = nanmean(Nneuronsagg.(recareas{a}).(whichvisblock).ICencoder./denom);
    tempmat(a,2) = nanmean(Nneuronsagg.(recareas{a}).(whichvisblock).RCencoder./denom);
        case 2
    denom = Nneuronsacc.(recareas{a}).(whichvisblock).(neudenom);
    tempmat(a,1) = Nneuronsacc.(recareas{a}).(whichvisblock).ICencoder/denom;
    tempmat(a,2) = Nneuronsacc.(recareas{a}).(whichvisblock).RCencoder/denom;
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
    denom = Nneuronsagg.(recareas{a}).(whichvisblock).(neudenom);
    tempICvec = Nneuronsagg.(recareas{a}).(whichvisblock).ICencoder./denom;
    tempRCvec = Nneuronsagg.(recareas{a}).(whichvisblock).RCencoder./denom;
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
annotation('textbox', [0.1 0.92 0.9 0.1], 'string', ['All Trials: ' whichvisblock], 'edgecolor', 'none', 'interpreter', 'none')
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
    denom = Nneuronsagg.(recareas{a}).(whichvisblock).(neudenom);
    tempmat(a,1) = nanmean(Nneuronsagg.(recareas{a}).(whichvisblock).indin./denom);
    tempmat(a,2) = nanmean(Nneuronsagg.(recareas{a}).(whichvisblock).sigBKnotsigBI./denom);
        case 2
    denom = Nneuronsacc.(recareas{a}).(whichvisblock).(neudenom);
    tempmat(a,1) = Nneuronsacc.(recareas{a}).(whichvisblock).indin/denom;
    tempmat(a,2) = Nneuronsacc.(recareas{a}).(whichvisblock).sigBKnotsigBI/denom;
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
