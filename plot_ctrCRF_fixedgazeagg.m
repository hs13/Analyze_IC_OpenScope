addpath(genpath('/Users/hyeyoung/Documents/CODE/helperfunctions/'))

if exist('/Users/hyeyoung/Documents/', 'dir')
    drivepath = '/Users/hyeyoung/Documents/';
elseif exist('/Volumes/GoogleDrive-116160770365018316196/My Drive/', 'dir')
    drivepath = '/Volumes/GoogleDrive-116160770365018316196/My Drive/';
elseif exist('G:/My Drive/', 'dir')
    drivepath = 'G:/My Drive/';
else
    error('check drivepath in this computer')
end
load([drivepath, 'DATA/ICexpts_submission22/openscope_popavg_fixedgazeagg.mat'])
load([drivepath, 'DATA/ICexpts_submission22/openscope_psthavg_fixedgazeagg.mat'])

probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
vislocs = {'VISam', 'VISpm', 'VISp', 'VISl', 'VISal', 'VISrl'};
visind = [6 5 1 2 4 3];

visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
    'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
ICblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations'};

Nrfs = size(RFCI_fixedgazeagg.C.Rrfclassic, 2);
Nszs = size(sizeCI_fixedgazeagg.C.Rsizeclassic, 2);
Ndirs = size(oriparams_fixedgazeagg.C.Rori, 2);
Noris = Ndirs/2;

dirvec = vis.sizeCI_presentations.directions;
if length(dirvec)~=Ndirs
    error('check sizeCI_presentations directions')
end
orivec = vis.sizeCI_presentations.directions(1:Noris);
szvec = [0, 4, 8, 16, 32, 64];

kerwinhalf = 5; kersigma = 2;
kerwinhalf = 25; kersigma = 5;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));

fs=18;
%%
aggneuloc_fixedgaze = cat(1,neuloc_fixedgazeagg{:});
Nneuronsall = size(aggneuloc_fixedgaze, 1);

neuctxaggall_fixedgaze = zeros(Nneuronsall,1);
neucnt = 0;
for iprobe = 1:numel(probes)
    neuctxaggall_fixedgaze( neucnt+find(neuctx_fixedgazeagg{iprobe}==1) ) = iprobe;
    neuctxaggall_fixedgaze( neucnt+find(neuctx_fixedgazeagg{iprobe}==0) ) = -iprobe;
    neucnt = neucnt+numel(neuctx_fixedgazeagg{iprobe});
end
if nnz(neuctxaggall_fixedgaze==0) ~= 0
    error('check neuctxaggall')
end

neuvisaggall_fixedgaze = zeros(Nneuronsall,1);
for iprobe = 1:numel(probes)
    if strcmp(vislocs{iprobe}, 'VISp')
        neuoi = contains(aggneuloc_fixedgaze, 'VISp') & ~contains(aggneuloc_fixedgaze, 'VISpm');
        neuvisaggall_fixedgaze(neuoi) = iprobe;
    else
        neuoi = contains(aggneuloc_fixedgaze, vislocs{iprobe});
        neuvisaggall_fixedgaze(neuoi) = iprobe;
    end
end

for ii = 1:numel(visind)
    probeind = find(visind==ii);
    fprintf('%s neurons in the right probe %.2f (%d/%d)\n', vislocs{probeind}, ...
        mean(neuctxaggall_fixedgaze(neuvisaggall_fixedgaze==probeind)==probeind), nnz((neuctxaggall_fixedgaze(neuvisaggall_fixedgaze==probeind)==probeind)), nnz(neuvisaggall_fixedgaze==probeind) )
end

for iprobe = 1:numel(probes)
    fprintf('%s %s\n', probes{iprobe}, vislocs{iprobe})
    fprintf('%s ', unique(aggneuloc_fixedgaze(neuctxaggall_fixedgaze==iprobe & neuvisaggall_fixedgaze~=iprobe))' )
    fprintf('\n')
end

%%
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


for b = 1:numel(visblocks)
    Ronavg_fixedgazeagg.(visblocks{b}).all = [];
    for iprobe = 1:numel(probes)
        Ronavg_fixedgazeagg.(visblocks{b}).all = cat(2, Ronavg_fixedgazeagg.(visblocks{b}).all, Ronavg_fixedgazeagg.(visblocks{b}).(probes{iprobe}));
    end
end

%% ctrCRF neurons

neudesc = 'ctrCRF9sigexcl';
switch neudesc
    case 'ctrCG9'
        neuctr = RFCI_fixedgazeagg.all.RFindclassic==1 & RFCI_fixedgazeagg.all.pRFclassic<0.05;
        neuctrdesc = 'center-RF (Mann-Whitney U)';
    case 'ctrCRF9'
        neuctr = RFCI_fixedgazeagg.all.RFindclassic==1 & RFCI_fixedgazeagg.all.Pkw_rfclassic<0.05;
        neuctrdesc = 'center-RF (Kruskal-Wallis)';
    case 'ctrCRF9sigexcl' % Bonferroni-Holms corrected
        neuctr = RFCI_fixedgazeagg.all.RFindclassic==1 & RFCI_fixedgazeagg.all.RFsigexclclassic==1;
        neuctrdesc = 'center-RF (Bonferroni-Holms)';
    case 'ctrCRF9exclsig' % only one RF has p<0.05
        neuctr = RFCI_fixedgazeagg.all.RFindclassic==1 & RFCI_fixedgazeagg.all.RFexclsigclassic==1;
        neuctrdesc = 'exclusive center-RF';
end

%% preferred IC/REl
RwICagg = [Ronavg_fixedgazeagg.ICwcfg0_presentations.all(ICtrialtypes==111,:); ...
    Ronavg_fixedgazeagg.ICwcfg1_presentations.all(ICtrialtypes==111,:); ...
    Ronavg_fixedgazeagg.ICwcfg0_presentations.all(ICtrialtypes==106,:); ...
    Ronavg_fixedgazeagg.ICwcfg1_presentations.all(ICtrialtypes==106,:)]';
RkICagg = [Ronavg_fixedgazeagg.ICkcfg0_presentations.all(ICtrialtypes==111,:); ...
    Ronavg_fixedgazeagg.ICkcfg1_presentations.all(ICtrialtypes==111,:); ...
    Ronavg_fixedgazeagg.ICkcfg0_presentations.all(ICtrialtypes==106,:); ...
    Ronavg_fixedgazeagg.ICkcfg1_presentations.all(ICtrialtypes==106,:)]';

RwICRCagg = [Ronavg_fixedgazeagg.ICwcfg0_presentations.all(ICtrialtypes==111,:); ...
    Ronavg_fixedgazeagg.ICwcfg1_presentations.all(ICtrialtypes==111,:); ...
    Ronavg_fixedgazeagg.ICwcfg0_presentations.all(ICtrialtypes==106,:); ...
    Ronavg_fixedgazeagg.ICwcfg1_presentations.all(ICtrialtypes==106,:); ...
    Ronavg_fixedgazeagg.ICwcfg0_presentations.all(ICtrialtypes==110,:); ...
    Ronavg_fixedgazeagg.ICwcfg1_presentations.all(ICtrialtypes==110,:); ...
    Ronavg_fixedgazeagg.ICwcfg0_presentations.all(ICtrialtypes==107,:); ...
    Ronavg_fixedgazeagg.ICwcfg1_presentations.all(ICtrialtypes==107,:)]';
RkICRCagg = [Ronavg_fixedgazeagg.ICkcfg0_presentations.all(ICtrialtypes==111,:); ...
    Ronavg_fixedgazeagg.ICkcfg1_presentations.all(ICtrialtypes==111,:); ...
    Ronavg_fixedgazeagg.ICkcfg0_presentations.all(ICtrialtypes==106,:); ...
    Ronavg_fixedgazeagg.ICkcfg1_presentations.all(ICtrialtypes==106,:); ...
    Ronavg_fixedgazeagg.ICkcfg0_presentations.all(ICtrialtypes==110,:); ...
    Ronavg_fixedgazeagg.ICkcfg1_presentations.all(ICtrialtypes==110,:); ...
    Ronavg_fixedgazeagg.ICkcfg0_presentations.all(ICtrialtypes==107,:); ...
    Ronavg_fixedgazeagg.ICkcfg1_presentations.all(ICtrialtypes==107,:)]';

RwRElagg = [Ronavg_fixedgazeagg.ICwcfg0_presentations.all(ICtrialtypes==511,:); ...
    Ronavg_fixedgazeagg.ICwcfg1_presentations.all(ICtrialtypes==511,:); ...
    Ronavg_fixedgazeagg.ICwcfg0_presentations.all(ICtrialtypes==506,:); ...
    Ronavg_fixedgazeagg.ICwcfg1_presentations.all(ICtrialtypes==506,:)]';
RkRElagg = [Ronavg_fixedgazeagg.ICkcfg0_presentations.all(ICtrialtypes==511,:); ...
    Ronavg_fixedgazeagg.ICkcfg1_presentations.all(ICtrialtypes==511,:); ...
    Ronavg_fixedgazeagg.ICkcfg0_presentations.all(ICtrialtypes==506,:); ...
    Ronavg_fixedgazeagg.ICkcfg1_presentations.all(ICtrialtypes==506,:)]';

RwREtagg = [Ronavg_fixedgazeagg.ICwcfg0_presentations.all(ICtrialtypes==1109,:); ...
    Ronavg_fixedgazeagg.ICwcfg1_presentations.all(ICtrialtypes==1109,:); ...
    Ronavg_fixedgazeagg.ICwcfg0_presentations.all(ICtrialtypes==1105,:); ...
    Ronavg_fixedgazeagg.ICwcfg1_presentations.all(ICtrialtypes==1105,:)]';
RkREtagg = [Ronavg_fixedgazeagg.ICkcfg0_presentations.all(ICtrialtypes==1109,:); ...
    Ronavg_fixedgazeagg.ICkcfg1_presentations.all(ICtrialtypes==1109,:); ...
    Ronavg_fixedgazeagg.ICkcfg0_presentations.all(ICtrialtypes==1105,:); ...
    Ronavg_fixedgazeagg.ICkcfg1_presentations.all(ICtrialtypes==1105,:)]';

[~,prefiwICagg] = max(RwICagg,[],2);
[~,prefikICagg] = max(RkICagg,[],2);
[~,prefiwICRCagg] = max(RwICRCagg,[],2);
[~,prefikICRCagg] = max(RkICRCagg,[],2);

[~,prefiwRElagg] = max(RwRElagg,[],2);
[~,prefikRElagg] = max(RkRElagg,[],2);
[~,prefiwREtagg] = max(RwREtagg,[],2);
[~,prefikREtagg] = max(RkREtagg,[],2);

neuinarea = neuvisaggall_fixedgaze==find(strcmp(vislocs, 'VISp'));
figure('Position', [100 100 300 240])
histogram2(prefiwICagg(neuinarea & neuctr), prefiwRElagg(neuinarea & neuctr), 'displaystyle', 'tile', 'showemptybins', 'on')
colorbar
colormap redblue
set(gca,'FontSize', fs, 'XTick', 1:4, 'XTickLabel', 0:45:180-1, 'YTick', 1:4, 'YTickLabel', 0:45:180-1)
xlabel('Pref. IC Orientation')
ylabel('Pref. RE_I Orientation')
% title(sprintf('V1 %s (fixed-gaze, N=%d)', neuctrdesc, nnz(neuinarea & neuctr) ), 'fontweight', 'normal')
fprintf('V1 %s (N=%d)', neuctrdesc, nnz(neuinarea & neuctr) )

figure; annotation('textbox', [0.1 0.9 0.8 0.1], 'string', sprintf('V1 %s (fixed-gaze, N=%d)', neuctrdesc, nnz(neuinarea & neuctr) ), 'fontsize', fs, 'edgecolor', 'none' )
subplot(2,2,1); 
h = histogram2(prefiwICagg(neuinarea & neuctr), prefiwRElagg(neuinarea & neuctr), 'displaystyle', 'tile', 'showemptybins', 'on')
colorbar
set(gca,'FontSize', fs, 'XTick', 1:4, 'XTickLabel', 0:45:180-1, 'YTick', 1:4, 'YTickLabel', 0:45:180-1)
xlabel('Pref wIC Ori')
ylabel('Pref wRE_I Ori')
subplot(2,2,2); histogram2(prefiwICagg(neuinarea & neuctr), prefiwREtagg(neuinarea & neuctr), 'displaystyle', 'tile', 'showemptybins', 'on')
colorbar
set(gca,'FontSize', fs, 'XTick', 1:4, 'XTickLabel', 0:45:180-1, 'YTick', 1:4, 'YTickLabel', 0:45:180-1)
xlabel('Pref wIC Ori')
ylabel('Pref wRE_T Ori')
subplot(2,2,3); histogram2(prefiwICRCagg(neuinarea & neuctr), prefiwREtagg(neuinarea & neuctr), 'displaystyle', 'tile', 'showemptybins', 'on')
colorbar
set(gca,'FontSize', fs, 'XTick', 1:4, 'XTickLabel', 0:45:180-1, 'YTick', 1:4, 'YTickLabel', 0:45:180-1)
xlabel('Pref wIC/RC Ori')
ylabel('Pref wRE_T Ori')
subplot(2,2,4); histogram2(prefiwICagg(neuinarea & neuctr), ori4params_fixedgazeagg.all.prefiori4(neuinarea & neuctr), 'displaystyle', 'tile', 'showemptybins', 'on')
colorbar
set(gca,'FontSize', fs, 'XTick', 1:4, 'XTickLabel', 0:45:180-1, 'YTick', 1:4, 'YTickLabel', 0:45:180-1)
xlabel('Pref wIC Ori')
ylabel('Pref Grating Ori')
colormap(redblue)

for icol = 0:1
    switch icol
        case 0
tempmatch = prefikICagg==prefikRElagg;
indcol = 'k';
        case 1
tempmatch = prefiwICagg==prefiwRElagg;
indcol = 'w';
    end
for ii = 1:numel(probes)
    probeind = find(visind==ii);
    neuoi = neuvisaggall_fixedgaze==probeind & neuctr;
    fprintf('%s pref %sIC-%sREl match %d/%d=%.4f\n', visareas{probeind}, indcol, indcol, ...
        nnz(tempmatch(neuoi)), nnz(neuoi), mean(tempmatch(neuoi)) )
end
fprintf('\n')
end

%% size tuning
figure('Position', [100 100 300 240])
hold all
legs = cell(size(probes));
for ii = 1:numel(probes)
    probeind = find(visind==ii);
    neuoi = neuvisaggall_fixedgaze==probeind & neuctr;
    errorbar(szvec, mean(sizeCI_fixedgazeagg.all.Rsizeclassic(neuoi,:), 1), std(sizeCI_fixedgazeagg.all.Rsizeclassic(neuoi,:),0,1)/sqrt(nnz(neuoi)) )
    legs{ii} = visareas{probeind};
end
legend(legs)
set(gca, 'XTick', szvec, 'XGrid', 'on')

%% psth
ttinds = [506 511; 107 110; 106 111];
ttcol = [0 0 1; 1 0.5 0; 0 0.7 0];
yl = [];
figure('Position', [100 100 1440 240])
hold all
for ii = 1:numel(probes)
    probeind = find(visind==ii);
neuoi = neuvisaggall_fixedgaze==probeind & neuctr;
%subplot(2,3,ii)
subplot(1,6,ii)
hold all
for typi = 1:3
    temppsth_ICwcfg1 = [];
    temppsth_ICwcfg0 = [];
    for iprobe = 1:numel(probes)
        ttoi = ismember(ICtrialtypes, ttinds(typi,:));
    temppsth_ICwcfg1 = cat(3, temppsth_ICwcfg1, psthavg_fixedgazeagg.ICwcfg1_presentations.(probes{iprobe})(:,ttoi,neuoi(neuprobe_fixedgazeagg==iprobe)));
    temppsth_ICwcfg0 = cat(3, temppsth_ICwcfg0, psthavg_fixedgazeagg.ICwcfg0_presentations.(probes{iprobe})(:,ttoi,neuoi(neuprobe_fixedgazeagg==iprobe)));
    end
    temppsth = squeeze(mean(cat(2,temppsth_ICwcfg1, temppsth_ICwcfg0),2));
    if size(temppsth,2) ~= nnz(neuoi)
        error('check temppsth')
    end
    shadedErrorBar(psthtli, mean(convn(temppsth,kergauss,'same'),2), std(convn(temppsth,kergauss,'same'),0,2)/sqrt(nnz(neuoi)), {'color', ttcol(typi,:), 'linewidth', 2}, 1 )
end
xlim([-100 400])
if ~isempty(yl)
    ylim(yl)
end
xlabel('Time (ms)')
ylabel('Firing Rate (Hz)')
title(sprintf('%s center-CRF Neurons N=%d', visareas{probeind}, nnz(neuoi) ))
end

% isequal( mean(convn(temppsth,kergauss,'same'),2), conv(mean(temppsth,2),kergauss,'same') )
% max(abs( mean(convn(temppsth,kergauss,'same'),2)-conv(mean(temppsth,2),kergauss,'same') ))
