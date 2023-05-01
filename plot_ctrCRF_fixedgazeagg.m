% checked version match on 230317
% olderver = load([drivepath, 'DATA/ICexpts_submission22/openscope_popavg_fixedgazeagg.mat'])
% newerver = load('/Users/hyeyoung/Documents/OpenScopeData/000248/postprocessed/openscope_popavg_fixedgazeagg.mat')
% !!! isequaln(olderver, newerver) returns FALSE!!!!
% isequaln(olderver.ICsig_fixedgazeagg, newerver.ICsig_fixedgazeagg) returns true
% !!! isequaln(olderver.RFCI_fixedgazeagg, newerver.RFCI_fixedgazeagg) returns FALSE!!!!
% newerver.RFCI_fixedgazeagg is empty, therefore should stick with the olderver

% olderpsth = load([drivepath, 'DATA/ICexpts_submission22/openscope_psthavg_fixedgazeagg.mat'])
% newerpsth = load('/Users/hyeyoung/Documents/OpenScopeData/000248/postprocessed/openscope_psthavg_fixedgazeagg.mat');
% isequaln(olderpsth, newerpsth) returns true

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
    case 'emergent'
        neuctr = ICsig_fixedgazeagg.ICwcfg1_presentations.all.PkwBK<0.05 & ICsig_fixedgazeagg.ICwcfg1_presentations.all.PkwBI>=0.05;
        neuctrdesc = 'PkwBK<0.05 & PkwBI>=0.05';
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
cb = colorbar;
cb.Label.String = '# Neurons';
cb.Label.Rotation = 270;
cb.Label.FontSize = fs;
cb.Label.Position(1) = 3.8;
colormap redblue
set(gca,'FontSize', fs, 'XGrid', 'off', 'YGrid', 'off', 'YDir', 'reverse', 'XTick', 1:4, 'XTickLabel', 0:45:180-1, 'YTick', 1:4, 'YTickLabel', 0:45:180-1)
xlabel('Pref. I_C Orientation','FontSize', fs)
ylabel('Pref. I_R_E Orientation','FontSize', fs)
% title(sprintf('V1 %s (fixed-gaze, N=%d)', neuctrdesc, nnz(neuinarea & neuctr) ), 'fontweight', 'normal')
fprintf('V1 %s (N=%d)', neuctrdesc, nnz(neuinarea & neuctr) )

figure; annotation('textbox', [0.1 0.9 0.8 0.1], 'string', sprintf('V1 %s (fixed-gaze, N=%d)', neuctrdesc, nnz(neuinarea & neuctr) ), 'fontsize', fs, 'edgecolor', 'none' )
subplot(2,2,1); 
h = histogram2(prefiwICagg(neuinarea & neuctr), prefiwRElagg(neuinarea & neuctr), 'displaystyle', 'tile', 'showemptybins', 'on')
colorbar
set(gca,'FontSize', fs, 'YDir', 'reverse', 'XTick', 1:4, 'XTickLabel', 0:45:180-1, 'YTick', 1:4, 'YTickLabel', 0:45:180-1)
xlabel('Pref wIC Ori')
ylabel('Pref wRE_I Ori')
subplot(2,2,2); histogram2(prefiwICagg(neuinarea & neuctr), prefiwREtagg(neuinarea & neuctr), 'displaystyle', 'tile', 'showemptybins', 'on')
colorbar
set(gca,'FontSize', fs, 'YDir', 'reverse', 'XTick', 1:4, 'XTickLabel', 0:45:180-1, 'YTick', 1:4, 'YTickLabel', 0:45:180-1)
xlabel('Pref wIC Ori')
ylabel('Pref wRE_T Ori')
subplot(2,2,3); histogram2(prefiwICRCagg(neuinarea & neuctr), prefiwREtagg(neuinarea & neuctr), 'displaystyle', 'tile', 'showemptybins', 'on')
colorbar
set(gca,'FontSize', fs, 'YDir', 'reverse', 'XTick', 1:4, 'XTickLabel', 0:45:180-1, 'YTick', 1:4, 'YTickLabel', 0:45:180-1)
xlabel('Pref wIC/RC Ori')
ylabel('Pref wRE_T Ori')
subplot(2,2,4); histogram2(prefiwICagg(neuinarea & neuctr), ori4params_fixedgazeagg.all.prefiori4(neuinarea & neuctr), 'displaystyle', 'tile', 'showemptybins', 'on')
colorbar
set(gca,'FontSize', fs, 'YDir', 'reverse', 'XTick', 1:4, 'XTickLabel', 0:45:180-1, 'YTick', 1:4, 'YTickLabel', 0:45:180-1)
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

%%
neuinarea = neuvisaggall_fixedgaze==find(strcmp(vislocs, 'VISp'));
figure('Position', [100 100 300 240])
histogram2(prefiwICagg(neuinarea & neuctr), ori4params_fixedgazeagg.all.prefiori4(neuinarea & neuctr), 'displaystyle', 'tile', 'showemptybins', 'on')
colorbar
colormap redblue
set(gca,'FontSize', fs, 'XGrid', 'off', 'YGrid', 'off', 'YDir', 'reverse', 'XTick', 1:4, 'XTickLabel', 0:45:180-1, 'YTick', 1:4, 'YTickLabel', 0:45:180-1)
xlabel('Pref. I_C Orientation')
ylabel('Pref. Grating Orientation')
% title(sprintf('V1 %s (fixed-gaze, N=%d)', neuctrdesc, nnz(neuinarea & neuctr) ), 'fontweight', 'normal')
fprintf('V1 %s (N=%d)', neuctrdesc, nnz(neuinarea & neuctr) )

%% size tuning
lw = 2;
figure('Position', [100 100 300 240])
hold all
legs = cell(size(probes));
for ii = 1:numel(probes)
    probeind = find(visind==ii);
    neuoi = neuvisaggall_fixedgaze==probeind & neuctr;
    errorbar(szvec, mean(sizeCI_fixedgazeagg.all.Rsizeclassic(neuoi,:), 1), ...
        std(sizeCI_fixedgazeagg.all.Rsizeclassic(neuoi,:),0,1)/sqrt(nnz(neuoi)), 'linewidth', lw )
    legs{ii} = visareas{probeind};
end
legend(legs)
set(gca, 'XTick', szvec, 'XGrid', 'on')

%% psth
% % REl, RC, IC
% ttinds = {[506 511], [107 110], [106 111]};
% ttcol = [0 0 1; 1 0.5 0; 0 0.7 0];
% REl, inducers, IC
ttinds = {[506 511], [1301:1308], [106 111]};
ttcol = [0 0 1; 0 0 0; 0 0.7 0];
yl = [];
% yl = [0 30];
figure%('Position', [100 100 1440 240])
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'fixed-gaze trials', 'edgecolor', 'none')%, 'fontsize', fs)
hold all
for ii = 1:numel(probes)
    probeind = find(visind==ii);
neuoi = neuvisaggall_fixedgaze==probeind & neuctr;
subplot(2,3,ii)
% subplot(1,6,ii)
hold all
for typi = 1:numel(ttinds)
    temppsth_ICwcfg1 = [];
    temppsth_ICwcfg0 = [];
    for iprobe = 1:numel(probes)
        ttoi = ismember(ICtrialtypes, ttinds{typi});
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

%% psth just V1
% ttinds = {[1305:1308], [1301:1304], [506 511], [106 111]};
% ttcol = [0.5 0.5 0.5; 0 0 0; 0 0 1; 0 0.7 0];
% ttlabs = {'P_o_u_t_1_-_4', 'P_i_n_1_-_4', 'I_R_E_1_-_2', 'I_C_1_-_2'};

% ttinds = {[1301:1308], [506 511], [106 111]};
% ttcol = [0 0 0; 0 0 1; 0 0.7 0];
% ttlabs = {'P_i_n_/_o_u_t_1_-_4', 'I_R_E_1_-_2', 'I_C_1_-_2'};

ttinds = {[506 511], [106 111]};
ttcol = [0 0 1; 0 0.7 0];
ttlabs = {'I_R_E_1_-_2', 'I_C_1_-_2'};


yl = [];

figure('Position', [100 100 300 240])
hold all
for ii = 1%:numel(probes)
    probeind = find(visind==ii);
neuoi = neuvisaggall_fixedgaze==probeind & neuctr;
hold all
for typi = 1:numel(ttinds)
    temppsth_ICwcfg1 = [];
    temppsth_ICwcfg0 = [];
    for iprobe = 1:numel(probes)
        ttoi = ismember(ICtrialtypes, ttinds{typi});
    temppsth_ICwcfg1 = cat(3, temppsth_ICwcfg1, psthavg_fixedgazeagg.ICwcfg1_presentations.(probes{iprobe})(:,ttoi,neuoi(neuprobe_fixedgazeagg==iprobe)));
    temppsth_ICwcfg0 = cat(3, temppsth_ICwcfg0, psthavg_fixedgazeagg.ICwcfg0_presentations.(probes{iprobe})(:,ttoi,neuoi(neuprobe_fixedgazeagg==iprobe)));
    end
    temppsth = squeeze(mean(cat(2,temppsth_ICwcfg1, temppsth_ICwcfg0),2));
    if size(temppsth,2) ~= nnz(neuoi)
        error('check temppsth')
    end
    shadedErrorBar(psthtli/1000, mean(convn(temppsth,kergauss,'same'),2), std(convn(temppsth,kergauss,'same'),0,2)/sqrt(nnz(neuoi)), {'color', ttcol(typi,:), 'linewidth', 2}, 1 )
end
xl = [-100 400]/1000;
xlim(xl)
if ~isempty(yl)
    ylim(yl)
end
yl = ylim;
for typi = 1:numel(ttinds)
    text(xl(1)+0.75*range(xl), yl(2)-(numel(ttinds)-typi)*0.15*range(yl), ttlabs{typi}, 'fontweight', 'bold', 'Color', ttcol(typi,:), 'FontSize', fs, 'VerticalAlignment', 'top')
end
text(xl(1)+0.75*range(xl), yl(1), sprintf('N=%d', nnz(neuoi)), 'fontweight', 'bold', 'FontSize', fs, 'VerticalAlignment', 'bottom');%, 'HorizontalAlignment', 'right')
ylim(yl)
set(gca, 'FontSize', fs)
xlabel('Time (s)','FontSize', fs)
ylabel('Firing Rate (Hz)','FontSize', fs)
% title(sprintf('%s center-CRF Neurons N=%d', visareas{probeind}, nnz(neuoi) ))
end

%% stacked psth just V1 individual cells REl vs IC trials
ttinds = {[106 111], [506 511]};
ttcol = [0 0.7 0; 0 0 1];
%ttlabs = {'I_C_1_-_2', 'I_R_E_1_-_2'};
ttlabs = {'I_C', 'I_R_E'};
cl = 30*[-1 1];

probeind = find(visind==1);
neuoi = neuvisaggall_fixedgaze==probeind & neuctr;
temppsth_ICwcfg1 = [];
temppsth_ICwcfg0 = [];
for iprobe = 1:numel(probes)
    temppsth_ICwcfg1 = cat(3, temppsth_ICwcfg1, psthavg_fixedgazeagg.ICwcfg1_presentations.(probes{iprobe})(:,ICtrialtypes==0,neuoi(neuprobe_fixedgazeagg==iprobe)));
    temppsth_ICwcfg0 = cat(3, temppsth_ICwcfg0, psthavg_fixedgazeagg.ICwcfg0_presentations.(probes{iprobe})(:,ICtrialtypes==0,neuoi(neuprobe_fixedgazeagg==iprobe)));
end
temppsth = squeeze(mean(cat(2,temppsth_ICwcfg1, temppsth_ICwcfg0),2));
Rblank = mean(temppsth(psthtli>0 & psthtli<=400,:),1);

plotcb = true;
figure('Position', [400 400 360 240])
for typi = 1:2
    subplot(1,2,typi)
temppsth_ICwcfg1 = [];
temppsth_ICwcfg0 = [];
for iprobe = 1:numel(probes)
    ttoi = ismember(ICtrialtypes, ttinds{typi});
    temppsth_ICwcfg1 = cat(3, temppsth_ICwcfg1, psthavg_fixedgazeagg.ICwcfg1_presentations.(probes{iprobe})(:,ttoi,neuoi(neuprobe_fixedgazeagg==iprobe)));
    temppsth_ICwcfg0 = cat(3, temppsth_ICwcfg0, psthavg_fixedgazeagg.ICwcfg0_presentations.(probes{iprobe})(:,ttoi,neuoi(neuprobe_fixedgazeagg==iprobe)));
end
temppsth = squeeze(mean(cat(2,temppsth_ICwcfg1, temppsth_ICwcfg0),2));

temppsth = temppsth - Rblank;
Revoked = mean(temppsth(psthtli>0 & psthtli<=400,:),1);
if typi ==1
[sv,si]=sort(Revoked, 'descend');
end

imagesc(psthtli/1000, 1:nnz(neuoi), convn(temppsth(:,si), kergauss, 'same')')
set(gca, 'YDir', 'reverse', 'FontSize', fs)
colormap redblue
xl = [0 400]/1000;
xlim(xl)
caxis(cl)
xlabel('Time (s)')
if typi==1
ylabel('Neurons')
else
    set(gca, 'YTickLabel', {})
end
if plotcb
cb = colorbar;
% cb.Position(1) = 0.75;
cb.Label.String = '\DeltaFiring Rate (Hz)';
cb.Label.FontSize = fs;
cb.Label.Rotation = 270;
disp(cb.Label.Position)
cb.Label.Position(1) = cb.Label.Position(1)+4;
title(sprintf('%s', ttlabs{typi}), 'color', ttcol(typi,:))
else
title(sprintf('%s Trials', ttlabs{typi}), 'color', ttcol(typi,:))
end
end

%% stacked psth just V1 individual cells IC trials
ttinds = {[106 111]};
ttcol = [0 0.7 0];
ttlabs = {'I_C_1_-_2'};

probeind = find(visind==1);
neuoi = neuvisaggall_fixedgaze==probeind & neuctr;
temppsth_ICwcfg1 = [];
temppsth_ICwcfg0 = [];
for iprobe = 1:numel(probes)
    temppsth_ICwcfg1 = cat(3, temppsth_ICwcfg1, psthavg_fixedgazeagg.ICwcfg1_presentations.(probes{iprobe})(:,ICtrialtypes==0,neuoi(neuprobe_fixedgazeagg==iprobe)));
    temppsth_ICwcfg0 = cat(3, temppsth_ICwcfg0, psthavg_fixedgazeagg.ICwcfg0_presentations.(probes{iprobe})(:,ICtrialtypes==0,neuoi(neuprobe_fixedgazeagg==iprobe)));
end
temppsth = squeeze(mean(cat(2,temppsth_ICwcfg1, temppsth_ICwcfg0),2));
Rblank = mean(temppsth(psthtli>0 & psthtli<=400,:),1);

typi = 1;
temppsth_ICwcfg1 = [];
temppsth_ICwcfg0 = [];
for iprobe = 1:numel(probes)
    ttoi = ismember(ICtrialtypes, ttinds{typi});
    temppsth_ICwcfg1 = cat(3, temppsth_ICwcfg1, psthavg_fixedgazeagg.ICwcfg1_presentations.(probes{iprobe})(:,ttoi,neuoi(neuprobe_fixedgazeagg==iprobe)));
    temppsth_ICwcfg0 = cat(3, temppsth_ICwcfg0, psthavg_fixedgazeagg.ICwcfg0_presentations.(probes{iprobe})(:,ttoi,neuoi(neuprobe_fixedgazeagg==iprobe)));
end
temppsth = squeeze(mean(cat(2,temppsth_ICwcfg1, temppsth_ICwcfg0),2));

temppsth = temppsth - Rblank;
Revoked = mean(temppsth(psthtli>0 & psthtli<=400,:),1);
[sv,si]=sort(Revoked, 'descend');

figure('Position', [400 100 330 240])
imagesc(psthtli/1000, 1:nnz(neuoi), convn(temppsth(:,si), kergauss, 'same')')
set(gca, 'YDir', 'reverse', 'FontSize', fs)
colormap redblue
xl = [0 400]/1000;
xlim(xl)
caxis(cl)
xlabel('Time (s)')
ylabel('Neurons')
cb = colorbar;
% cb.Position(1) = 0.75;
cb.Label.String = '\DeltaFiring Rate (Hz)';
cb.Label.FontSize = fs;
cb.Label.Rotation = 270;
cb.Label.Position(1) = 4.1;


yl = [];
figure('Position', [100 100 300 240])
hold all
for ii = 1%:numel(probes)
    probeind = find(visind==ii);
neuoi = neuvisaggall_fixedgaze==probeind & neuctr;
hold all
for typi = 1:numel(ttinds)
    temppsth_ICwcfg1 = [];
    temppsth_ICwcfg0 = [];
    for iprobe = 1:numel(probes)
        ttoi = ismember(ICtrialtypes, ttinds{typi});
    temppsth_ICwcfg1 = cat(3, temppsth_ICwcfg1, psthavg_fixedgazeagg.ICwcfg1_presentations.(probes{iprobe})(:,ttoi,neuoi(neuprobe_fixedgazeagg==iprobe)));
    temppsth_ICwcfg0 = cat(3, temppsth_ICwcfg0, psthavg_fixedgazeagg.ICwcfg0_presentations.(probes{iprobe})(:,ttoi,neuoi(neuprobe_fixedgazeagg==iprobe)));
    end
    temppsth = squeeze(mean(cat(2,temppsth_ICwcfg1, temppsth_ICwcfg0),2));
    if size(temppsth,2) ~= nnz(neuoi)
        error('check temppsth')
    end
    %shadedErrorBar(psthtli/1000, mean(convn(temppsth,kergauss,'same'),2), std(convn(temppsth,kergauss,'same'),0,2)/sqrt(nnz(neuoi)), {'color', ttcol(typi,:), 'linewidth', 2}, 1 )
    plot(psthtli/1000, convn(temppsth,kergauss,'same') )
    plot(psthtli/1000, mean(convn(temppsth,kergauss,'same'),2), 'k-', 'LineWidth', 3)
end
xl = [-100 400]/1000;
xlim(xl)
if ~isempty(yl)
    ylim(yl)
end
yl = ylim;
for typi = 1:numel(ttinds)
    text(xl(1)+0.75*range(xl), yl(2)-(numel(ttinds)-typi)*0.15*range(yl), ttlabs{typi}, 'fontweight', 'bold', 'Color', ttcol(typi,:), 'FontSize', fs, 'VerticalAlignment', 'top')
end
text(xl(1)+0.75*range(xl), yl(1), sprintf('N=%d', nnz(neuoi)), 'fontweight', 'bold', 'FontSize', fs, 'VerticalAlignment', 'bottom');%, 'HorizontalAlignment', 'right')
ylim(yl)
set(gca, 'FontSize', fs)
xlabel('Time (s)','FontSize', fs)
ylabel('Firing Rate (Hz)','FontSize', fs)
% title(sprintf('%s center-CRF Neurons N=%d', visareas{probeind}, nnz(neuoi) ))
end

%% proportion of IC responsive neurons among ctrRF neurons
probeind = find(visind==1);
neuoi = neuvisaggall_fixedgaze==probeind & neuctr;

sigBKw = ICsig_fixedgazeagg.ICwcfg0_presentations.all.PkwBK(:,1)<0.05 | ICsig_fixedgazeagg.ICwcfg1_presentations.all.PkwBK(:,1)<0.05;
mean(sigBKw(neuoi))

ICresp11 = ICsig_fixedgazeagg.ICwcfg1_presentations.all.SP_BK(:,1)>0.5 & ICsig_fixedgazeagg.ICwcfg1_presentations.all.Pmww_BK(:,1)<0.05;
ICresp12 = ICsig_fixedgazeagg.ICwcfg1_presentations.all.SP_BK(:,end)>0.5 & ICsig_fixedgazeagg.ICwcfg1_presentations.all.Pmww_BK(:,end)<0.05;
ICresp01 = ICsig_fixedgazeagg.ICwcfg0_presentations.all.SP_BK(:,1)>0.5 & ICsig_fixedgazeagg.ICwcfg0_presentations.all.Pmww_BK(:,1)<0.05;
ICresp02 = ICsig_fixedgazeagg.ICwcfg0_presentations.all.SP_BK(:,end)>0.5 & ICsig_fixedgazeagg.ICwcfg0_presentations.all.Pmww_BK(:,end)<0.05;

ICresp = ICresp11 | ICresp12 | ICresp01 | ICresp02;

SP_BIC = cat(2, ICsig_fixedgazeagg.ICwcfg0_presentations.all.SP_BK(:,1), ICsig_fixedgazeagg.ICwcfg0_presentations.all.SP_BK(:,end), ...
    ICsig_fixedgazeagg.ICwcfg1_presentations.all.SP_BK(:,1), ICsig_fixedgazeagg.ICwcfg1_presentations.all.SP_BK(:,end));
Pmww_BIC = cat(2, ICsig_fixedgazeagg.ICwcfg0_presentations.all.Pmww_BK(:,1), ICsig_fixedgazeagg.ICwcfg0_presentations.all.Pmww_BK(:,end), ...
    ICsig_fixedgazeagg.ICwcfg1_presentations.all.Pmww_BK(:,1), ICsig_fixedgazeagg.ICwcfg1_presentations.all.Pmww_BK(:,end));
[pmat, si] = sort(Pmww_BIC, 2);
reordind = sub2ind(size(pmat), repmat((1:size(pmat,1))',1,size(pmat,2)), si);
if ~isequaln(Pmww_BIC(reordind), pmat)
    error('check reordind')
end
ICrespbonfholms = any(SP_BIC(reordind)>0.5 & pmat.*[4:-1:1]<0.1, 2); % one-tail
ICrespbonfholms2 = any(SP_BIC(reordind)>0.5 & pmat.*[4:-1:1]<0.05, 2); % two-tailed

disp([mean(ICresp(neuoi)) mean(ICrespbonfholms(neuoi)) mean(ICrespbonfholms2(neuoi))])
fprintf('%.4f is the proportion of IC responsive neurons among ctrRF neurons\n', mean(ICresp(neuoi)))
% 0.7917

M = ICsig_fixedgazeagg.ICkcfg1_presentations.all.sigmcBICREl2;
[Mu,ia,ic] = unique(M, 'rows', 'stable'); % Unique Values By Row, Retaining Original Order
h = accumarray(ic, 1); % Count Occurrences
% maph = h(ic); % Map Occurrences To ‘ic’ Values
disp([Mu, h])

sigmcBICRElw = cat(2, ICsig_fixedgazeagg.ICwcfg0_presentations.all.sigmcBICREl2, ...
    ICsig_fixedgazeagg.ICwcfg1_presentations.all.sigmcBICREl2, ...
    ICsig_fixedgazeagg.ICwcfg0_presentations.all.sigmcBICREl1, ...
    ICsig_fixedgazeagg.ICwcfg1_presentations.all.sigmcBICREl1);
mean(any(sigmcBICRElw(neuoi, 1:2:end)==1, 2))

fprintf('%.4f vs %.4f is the proportion of IC vs REl responsive neurons among %s neurons\n', ...
    mean(any(sigmcBICRElw(neuoi, 1:2:end)==1, 2)), mean(any(sigmcBICRElw(neuoi, 2:2:end)==1, 2)), neudesc )

%% 2d-histogram preferred IC vs REl: IC-responsive center-RF neurons
neuinarea = neuvisaggall_fixedgaze==find(strcmp(vislocs, 'VISp'));
neuoi = neuinarea & neuctr & any(sigmcBICRElw(:, 1:2:end)==1, 2);% & any(sigmcBICRElw(:, 2:2:end)==1, 2);

tempmatch = prefiwICagg==prefiwRElagg;
fprintf('Pref IC-REl match %d/%d=%.4f\n', ...
    nnz(tempmatch(neuoi)), nnz(neuoi), mean(tempmatch(neuoi)) )

figure('Position', [100 100 300 240])
histogram2(prefiwICagg(neuoi), ori4params_fixedgazeagg.all.prefiori4(neuoi), 'displaystyle', 'tile', 'showemptybins', 'on')
colorbar
colormap redblue
set(gca,'FontSize', fs, 'XGrid', 'off', 'YGrid', 'off', 'YDir', 'reverse', 'XTick', 1:4, 'XTickLabel', 0:45:180-1, 'YTick', 1:4, 'YTickLabel', 0:45:180-1)
xlabel('Pref. I_C Orientation')
ylabel('Pref. Grating Orientation')
% title(sprintf('V1 %s (fixed-gaze, N=%d)', neuctrdesc, nnz(neuinarea & neuctr) ), 'fontweight', 'normal')
fprintf('V1 IC-responsive %s (N=%d)\n', neuctrdesc, nnz(neuoi) )


figure('Position', [400 100 300 240])
histogram2(prefiwICagg(neuoi), prefiwRElagg(neuoi), 'displaystyle', 'tile', 'showemptybins', 'on')
cb = colorbar;
cb.Label.String = '# Neurons';
cb.Label.Rotation = 270;
cb.Label.FontSize = fs;
cb.Label.Position(1) = 3.8;
colormap redblue
set(gca,'FontSize', fs, 'XGrid', 'off', 'YGrid', 'off', 'YDir', 'reverse', 'XTick', 1:4, 'XTickLabel', 0:45:180-1, 'YTick', 1:4, 'YTickLabel', 0:45:180-1)
xlabel('Pref. I_C Orientation','FontSize', fs)
ylabel('Pref. I_R_E Orientation','FontSize', fs)
% title(sprintf('V1 %s (fixed-gaze, N=%d)', neuctrdesc, nnz(neuinarea & neuctr) ), 'fontweight', 'normal')
fprintf('V1 IC-responsive %s (N=%d)\n', neuctrdesc, nnz(neuoi) )

%% sum of parts vs whole
ii=1; probeind = find(visind==ii);
neuoi = neuvisaggall_fixedgaze==probeind & neuctr;
suminds = sum(Ronavg_fixedgazeagg.ICwcfg1_presentations.all(ismember(ICtrialtypes, [1301:1308]), neuoi) - Ronavg_fixedgazeagg.ICwcfg1_presentations.all(ICtrialtypes==0, neuoi), 1);
sumICs = sum(Ronavg_fixedgazeagg.ICwcfg1_presentations.all(ismember(ICtrialtypes, [106 111]), neuoi) - Ronavg_fixedgazeagg.ICwcfg1_presentations.all(ICtrialtypes==0, neuoi), 1);

figure; hold all
plot(suminds, sumICs, 'o')
xl = xlim;
plot(xl, xl, 'k--')


%% RF heatmap
neuinarea = neuvisaggall_fixedgaze==find(strcmp(vislocs, 'VISp'));
neuoi = neuinarea & neuctr;

RFCItrialtypes = vistrialtypes_fixedgazeagg(1).RFCI_presentations;
crftrialtypes = RFCItrialtypes<10000;

Rspon = cat(1,sponFRvec_fixedgazeagg{:});%meanFRvec_fixedgazeagg
Rcrf9agg = Ronavg_fixedgazeagg.RFCI_presentations.all(crftrialtypes,:)';
if ~( size(Rcrf9agg,1)==length(neuctr) )
    error('mismatch between neuctr and RFCIagg')
end
% tempR = mean(Rcrf9agg(neuctr,:), 1);
tempR = mean(Rcrf9agg(neuctr,:)-Rspon(neuctr,:), 1);

RFcentersrel9 = [0 0; vis.RFCI_presentations.RFcentersrel];

% figure('Position', [800 300 300 200])
figure('Position', [800 300 360 240])
hold all
h=scatter(RFcentersrel9(:,2), RFcentersrel9(:,1), [], tempR, 'filled', 'MarkerFaceAlpha',.7);
for c = 1:size(RFcentersrel9,1)
    text(RFcentersrel9(c,2), RFcentersrel9(c,1), num2str(c), ...
        'FontSize', fs, 'FontWeight', 'bold', 'Color','c', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
end
caxis([0 8])
cb = colorbar;
cb.Label.String = '\DeltaFiring Rate (Hz)';
cb.Label.Rotation = 270;
cb.Label.FontSize = fs;
cb.Label.Position(1) = 4;
axis off
axis square
set(gca, 'FontSize', fs, 'YDir', 'reverse')
axis([-1.5 1.5 -1.5 1.5])
s=.82;
s=.8;
currentunits = get(gca,'Units');
set(gca, 'Units', 'Points');
axpos = get(gca,'Position');
set(gca, 'Units', currentunits);
markerWidth = s/3*axpos(3); % Calculate Marker width in points
set(h, 'SizeData', markerWidth^2)
% colormap gray
colormap redblue
