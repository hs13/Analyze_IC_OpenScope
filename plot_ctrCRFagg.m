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
load([drivepath, 'DATA/ICexpts_submission22/openscope_popavg_agg.mat'])
load([drivepath, 'DATA/ICexpts_submission22/openscope_psthavgagg.mat'])

probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
vislocs = {'VISam', 'VISpm', 'VISp', 'VISl', 'VISal', 'VISrl'};
visind = [6 5 1 2 4 3];

visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
    'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
ICblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations'};

Nrfs = size(RFCIagg.C.Rrfclassic, 2);
Nszs = size(sizeCIagg.C.Rsizeclassic, 2);
Ndirs = size(oriparamsagg.C.Rori, 2);
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
aggneuloc = cat(1,neulocagg{:});
Nneuronsall = size(aggneuloc, 1);

neuctxaggall = zeros(Nneuronsall,1);
neucnt = 0;
for iprobe = 1:numel(probes)
    neuctxaggall( neucnt+find(neuctxagg{iprobe}==1) ) = iprobe;
    neuctxaggall( neucnt+find(neuctxagg{iprobe}==0) ) = -iprobe;
    neucnt = neucnt+numel(neuctxagg{iprobe});
end
if nnz(neuctxaggall==0) ~= 0
    error('check neuctxaggall')
end

neuvisaggall = zeros(Nneuronsall,1);
for iprobe = 1:numel(probes)
    if strcmp(vislocs{iprobe}, 'VISp')
        neuoi = contains(aggneuloc, 'VISp') & ~contains(aggneuloc, 'VISpm');
        neuvisaggall(neuoi) = iprobe;
    else
        neuoi = contains(aggneuloc, vislocs{iprobe});
        neuvisaggall(neuoi) = iprobe;
    end
end

for ii = 1:numel(visind)
    probeind = find(visind==ii);
    fprintf('%s neurons in the right probe %.2f (%d/%d)\n', vislocs{probeind}, ...
        mean(neuctxaggall(neuvisaggall==probeind)==probeind), nnz((neuctxaggall(neuvisaggall==probeind)==probeind)), nnz(neuvisaggall==probeind) )
end

for iprobe = 1:numel(probes)
    fprintf('%s %s\n', probes{iprobe}, vislocs{iprobe})
    fprintf('%s ', unique(aggneuloc(neuctxaggall==iprobe & neuvisaggall~=iprobe))' )
    fprintf('\n')
end

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


for b = 1:numel(visblocks)
    Ronavgagg.(visblocks{b}).all = [];
    for iprobe = 1:numel(probes)
        Ronavgagg.(visblocks{b}).all = cat(2, Ronavgagg.(visblocks{b}).all, Ronavgagg.(visblocks{b}).(probes{iprobe}));
    end
end

%% ctrCRF neurons

neudesc = 'ctrCRF9sigexcl';
switch neudesc
    case 'ctrCG9'
        neuctr = RFCIagg.all.RFindclassic==1 & RFCIagg.all.pRFclassic<0.05;
        neuctrdesc = 'center-RF (Mann-Whitney U)';
    case 'ctrCRF9'
        neuctr = RFCIagg.all.RFindclassic==1 & RFCIagg.all.Pkw_rfclassic<0.05;
        neuctrdesc = 'center-RF (Kruskal-Wallis)';
    case 'ctrCRF9sigexcl' % Bonferroni-Holms corrected
        neuctr = RFCIagg.all.RFindclassic==1 & RFCIagg.all.RFsigexclclassic==1;
        neuctrdesc = 'center-RF (Bonferroni-Holms)';
    case 'ctrCRF9exclsig' % only one RF has p<0.05
        neuctr = RFCIagg.all.RFindclassic==1 & RFCIagg.all.RFexclsigclassic==1;
        neuctrdesc = 'exclusive center-RF';
end

%% preferred IC/REl

RwICagg = [Ronavgagg.ICwcfg0_presentations.all(ICtrialtypes==111,:); ...
    Ronavgagg.ICwcfg1_presentations.all(ICtrialtypes==111,:); ...
    Ronavgagg.ICwcfg0_presentations.all(ICtrialtypes==106,:); ...
    Ronavgagg.ICwcfg1_presentations.all(ICtrialtypes==106,:)]';
RkICagg = [Ronavgagg.ICkcfg0_presentations.all(ICtrialtypes==111,:); ...
    Ronavgagg.ICkcfg1_presentations.all(ICtrialtypes==111,:); ...
    Ronavgagg.ICkcfg0_presentations.all(ICtrialtypes==106,:); ...
    Ronavgagg.ICkcfg1_presentations.all(ICtrialtypes==106,:)]';

RwICRCagg = [Ronavgagg.ICwcfg0_presentations.all(ICtrialtypes==111,:); ...
    Ronavgagg.ICwcfg1_presentations.all(ICtrialtypes==111,:); ...
    Ronavgagg.ICwcfg0_presentations.all(ICtrialtypes==106,:); ...
    Ronavgagg.ICwcfg1_presentations.all(ICtrialtypes==106,:); ...
    Ronavgagg.ICwcfg0_presentations.all(ICtrialtypes==110,:); ...
    Ronavgagg.ICwcfg1_presentations.all(ICtrialtypes==110,:); ...
    Ronavgagg.ICwcfg0_presentations.all(ICtrialtypes==107,:); ...
    Ronavgagg.ICwcfg1_presentations.all(ICtrialtypes==107,:)]';
RkICRCagg = [Ronavgagg.ICkcfg0_presentations.all(ICtrialtypes==111,:); ...
    Ronavgagg.ICkcfg1_presentations.all(ICtrialtypes==111,:); ...
    Ronavgagg.ICkcfg0_presentations.all(ICtrialtypes==106,:); ...
    Ronavgagg.ICkcfg1_presentations.all(ICtrialtypes==106,:); ...
    Ronavgagg.ICkcfg0_presentations.all(ICtrialtypes==110,:); ...
    Ronavgagg.ICkcfg1_presentations.all(ICtrialtypes==110,:); ...
    Ronavgagg.ICkcfg0_presentations.all(ICtrialtypes==107,:); ...
    Ronavgagg.ICkcfg1_presentations.all(ICtrialtypes==107,:)]';

RwRElagg = [Ronavgagg.ICwcfg0_presentations.all(ICtrialtypes==511,:); ...
    Ronavgagg.ICwcfg1_presentations.all(ICtrialtypes==511,:); ...
    Ronavgagg.ICwcfg0_presentations.all(ICtrialtypes==506,:); ...
    Ronavgagg.ICwcfg1_presentations.all(ICtrialtypes==506,:)]';
RkRElagg = [Ronavgagg.ICkcfg0_presentations.all(ICtrialtypes==511,:); ...
    Ronavgagg.ICkcfg1_presentations.all(ICtrialtypes==511,:); ...
    Ronavgagg.ICkcfg0_presentations.all(ICtrialtypes==506,:); ...
    Ronavgagg.ICkcfg1_presentations.all(ICtrialtypes==506,:)]';

RwREtagg = [Ronavgagg.ICwcfg0_presentations.all(ICtrialtypes==1109,:); ...
    Ronavgagg.ICwcfg1_presentations.all(ICtrialtypes==1109,:); ...
    Ronavgagg.ICwcfg0_presentations.all(ICtrialtypes==1105,:); ...
    Ronavgagg.ICwcfg1_presentations.all(ICtrialtypes==1105,:)]';
RkREtagg = [Ronavgagg.ICkcfg0_presentations.all(ICtrialtypes==1109,:); ...
    Ronavgagg.ICkcfg1_presentations.all(ICtrialtypes==1109,:); ...
    Ronavgagg.ICkcfg0_presentations.all(ICtrialtypes==1105,:); ...
    Ronavgagg.ICkcfg1_presentations.all(ICtrialtypes==1105,:)]';

[~,prefiwICagg] = max(RwICagg,[],2);
[~,prefikICagg] = max(RkICagg,[],2);
[~,prefiwICRCagg] = max(RwICRCagg,[],2);
[~,prefikICRCagg] = max(RkICRCagg,[],2);

[~,prefiwRElagg] = max(RwRElagg,[],2);
[~,prefikRElagg] = max(RkRElagg,[],2);
[~,prefiwREtagg] = max(RwREtagg,[],2);
[~,prefikREtagg] = max(RkREtagg,[],2);

neuinarea = neuvisaggall==find(strcmp(vislocs, 'VISp'));
figure; histogram2(prefiwICagg(neuinarea & neuctr), prefiwRElagg(neuinarea & neuctr), 'displaystyle', 'tile')
figure; histogram2(prefiwICagg(neuinarea & neuctr), prefiwREtagg(neuinarea & neuctr), 'displaystyle', 'tile')

figure; histogram2(prefiwICRCagg(neuinarea & neuctr), prefiwREtagg(neuinarea & neuctr), 'displaystyle', 'tile', 'showemptybins', 'on')
colormap(redblue)

neuinarea = neuvisaggall==find(strcmp(vislocs, 'VISp'));
figure('Position', [100 100 300 240])
histogram2(prefiwICagg(neuinarea & neuctr), prefiwRElagg(neuinarea & neuctr), 'displaystyle', 'tile', 'showemptybins', 'on')
colorbar
colormap redblue
set(gca,'FontSize', fs, 'YDir', 'reverse', 'XTick', 1:4, 'XTickLabel', 0:45:180-1, 'YTick', 1:4, 'YTickLabel', 0:45:180-1)
xlabel('Pref. IC Orientation')
ylabel('Pref. RE_I Orientation')
% title(sprintf('V1 %s (fixed-gaze, N=%d)', neuctrdesc, nnz(neuinarea & neuctr) ), 'fontweight', 'normal')
fprintf('V1 %s (N=%d)', neuctrdesc, nnz(neuinarea & neuctr) )

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
    neuoi = neuvisaggall==probeind & neuctr;
    fprintf('%s pref %sIC-%sREl match %d/%d=%.4f\n', visareas{probeind}, indcol, indcol, ...
        nnz(tempmatch(neuoi)), nnz(neuoi), mean(tempmatch(neuoi)) )
end
fprintf('\n')
end

%% size tuning
figure
hold all
legs = cell(size(probes));
for ii = 1:numel(probes)
    probeind = find(visind==ii);
    neuoi = neuvisaggall==probeind & neuctr;
    errorbar(szvec, mean(sizeCIagg.all.Rsizeclassic(neuoi,:), 1), std(sizeCIagg.all.Rsizeclassic(neuoi,:),0,1)/sqrt(nnz(neuoi)) )
    legs{ii} = visareas{probeind};
end
legend(legs)
set(gca, 'XTick', szvec, 'XGrid', 'on')

%% psth
% fs=10;
% % REl, RC, IC
% ttinds = {[506 511], [107 110], [106 111]};
% ttcol = [0 0 1; 1 0.5 0; 0 0.7 0];
% REl, inducers, IC
ttinds = {[506 511], [1301:1308], [106 111]};
ttcol = [0 0 1; 0 0 0; 0 0.7 0];
yl = [2 18];
yl = [0 30];
figure
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'all trials', 'edgecolor', 'none')%, 'fontsize', fs)
hold all
for ii = 1:numel(probes)
    probeind = find(visind==ii);
neuoi = neuvisaggall==probeind & neuctr;
subplot(2,3,ii)
hold all
for typi = 1:numel(ttinds)
    temppsth_ICwcfg1 = [];
    temppsth_ICwcfg0 = [];
    for iprobe = 1:numel(probes)
        ttoi = ismember(ICtrialtypes, ttinds{typi});
    temppsth_ICwcfg1 = cat(3, temppsth_ICwcfg1, psthavgagg.ICwcfg1_presentations.(probes{iprobe})(:,ttoi,neuoi(neuprobeagg==iprobe)));
    temppsth_ICwcfg0 = cat(3, temppsth_ICwcfg0, psthavgagg.ICwcfg0_presentations.(probes{iprobe})(:,ttoi,neuoi(neuprobeagg==iprobe)));
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
% set(gca, 'fontsize', fs)
xlabel('Time (ms)')
ylabel('Firing Rate (Hz)')
title(sprintf('%s center-CRF Neurons N=%d', visareas{probeind}, nnz(neuoi) ))
end

% isequal( mean(convn(temppsth,kergauss,'same'),2), conv(mean(temppsth,2),kergauss,'same') )
% max(abs( mean(convn(temppsth,kergauss,'same'),2)-conv(mean(temppsth,2),kergauss,'same') ))

%% psth just V1
ttinds = {[1301:1304], [506 511], [106 111]};
ttcol = [0 0 0; 0 0 1; 0 0.7 0];
ttlabs = {'P1-8', 'I_R_E1-2', 'I1-2'};
yl = [];

figure('Position', [100 100 300 240])
hold all
for ii = 1%:numel(probes)
    probeind = find(visind==ii);
neuoi = neuvisaggall==probeind & neuctr;
hold all
for typi = 1:numel(ttinds)
    temppsth_ICwcfg1 = [];
    temppsth_ICwcfg0 = [];
    for iprobe = 1:numel(probes)
        ttoi = ismember(ICtrialtypes, ttinds{typi});
    temppsth_ICwcfg1 = cat(3, temppsth_ICwcfg1, psthavgagg.ICwcfg1_presentations.(probes{iprobe})(:,ttoi,neuoi(neuprobeagg==iprobe)));
    temppsth_ICwcfg0 = cat(3, temppsth_ICwcfg0, psthavgagg.ICwcfg0_presentations.(probes{iprobe})(:,ttoi,neuoi(neuprobeagg==iprobe)));
    end
    temppsth = squeeze(mean(cat(2,temppsth_ICwcfg1, temppsth_ICwcfg0),2));
    if size(temppsth,2) ~= nnz(neuoi)
        error('check temppsth')
    end
    shadedErrorBar(psthtli/1000, mean(convn(temppsth,kergauss,'same'),2), std(convn(temppsth,kergauss,'same'),0,2)/sqrt(nnz(neuoi)), {'color', ttcol(typi,:), 'linewidth', 2}, 1 )
end
xlim([-100 400]/1000)
if ~isempty(yl)
    ylim(yl)
end
yl = ylim;
for typi = 1:numel(ttinds)
    text(250/1000, yl(2)-(numel(ttinds)-typi)*0.15*range(yl), ttlabs{typi}, 'fontweight', 'bold', 'Color', ttcol(typi,:), 'FontSize', fs, 'VerticalAlignment', 'top')
end
ylim(yl)
set(gca, 'FontSize', fs)
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
% title(sprintf('%s center-CRF Neurons N=%d', visareas{probeind}, nnz(neuoi) ))
end

%% stacked psth just V1 individual cells REl vs IC trials
ttinds = {[106 111], [506 511]};
ttcol = [0 0.7 0; 0 0 1];
ttlabs = {'I_C_1_-_2', 'I_R_E_1_-_2'};
cl = 30*[-1 1];

probeind = find(visind==1);
neuoi = neuvisaggall==probeind & neuctr;
temppsth_ICwcfg1 = [];
temppsth_ICwcfg0 = [];
for iprobe = 1:numel(probes)
    temppsth_ICwcfg1 = cat(3, temppsth_ICwcfg1, psthavgagg.ICwcfg1_presentations.(probes{iprobe})(:,ICtrialtypes==0,neuoi(neuprobeagg==iprobe)));
    temppsth_ICwcfg0 = cat(3, temppsth_ICwcfg0, psthavgagg.ICwcfg0_presentations.(probes{iprobe})(:,ICtrialtypes==0,neuoi(neuprobeagg==iprobe)));
end
temppsth = squeeze(mean(cat(2,temppsth_ICwcfg1, temppsth_ICwcfg0),2));
Rblank = mean(temppsth(psthtli>0 & psthtli<=400,:),1);

plotcb = false;
figure('Position', [400 400 360 240])
for typi = 1:2
    subplot(1,2,typi)
temppsth_ICwcfg1 = [];
temppsth_ICwcfg0 = [];
for iprobe = 1:numel(probes)
    ttoi = ismember(ICtrialtypes, ttinds{typi});
    temppsth_ICwcfg1 = cat(3, temppsth_ICwcfg1, psthavgagg.ICwcfg1_presentations.(probes{iprobe})(:,ttoi,neuoi(neuprobeagg==iprobe)));
    temppsth_ICwcfg0 = cat(3, temppsth_ICwcfg0, psthavgagg.ICwcfg0_presentations.(probes{iprobe})(:,ttoi,neuoi(neuprobeagg==iprobe)));
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
cb.Label.Position(1) = 9.5;
title(sprintf('%s', ttlabs{typi}), 'color', ttcol(typi,:))
else
title(sprintf('%s Trials', ttlabs{typi}), 'color', ttcol(typi,:))
end
end

%% proportion of IC responsive neurons among ctrRF neurons
probeind = find(visind==1);
neuoi = neuvisaggall==probeind & neuctr;

sigBKw = ICsigagg.ICwcfg0_presentations.all.PkwBK(:,1)<0.05 | ICsigagg.ICwcfg1_presentations.all.PkwBK(:,1)<0.05;
mean(sigBKw(neuoi))

ICresp11 = ICsigagg.ICwcfg1_presentations.all.SP_BK(:,1)>0.5 & ICsigagg.ICwcfg1_presentations.all.Pmww_BK(:,1)<0.05;
ICresp12 = ICsigagg.ICwcfg1_presentations.all.SP_BK(:,end)>0.5 & ICsigagg.ICwcfg1_presentations.all.Pmww_BK(:,end)<0.05;
ICresp01 = ICsigagg.ICwcfg0_presentations.all.SP_BK(:,1)>0.5 & ICsigagg.ICwcfg0_presentations.all.Pmww_BK(:,1)<0.05;
ICresp02 = ICsigagg.ICwcfg0_presentations.all.SP_BK(:,end)>0.5 & ICsigagg.ICwcfg0_presentations.all.Pmww_BK(:,end)<0.05;

ICresp = ICresp11 | ICresp12 | ICresp01 | ICresp02;

SP_BIC = cat(2, ICsigagg.ICwcfg0_presentations.all.SP_BK(:,1), ICsigagg.ICwcfg0_presentations.all.SP_BK(:,end), ...
    ICsigagg.ICwcfg1_presentations.all.SP_BK(:,1), ICsigagg.ICwcfg1_presentations.all.SP_BK(:,end));
Pmww_BIC = cat(2, ICsigagg.ICwcfg0_presentations.all.Pmww_BK(:,1), ICsigagg.ICwcfg0_presentations.all.Pmww_BK(:,end), ...
    ICsigagg.ICwcfg1_presentations.all.Pmww_BK(:,1), ICsigagg.ICwcfg1_presentations.all.Pmww_BK(:,end));
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

M = ICsigagg.ICkcfg1_presentations.all.sigmcBICREl2;
[Mu,ia,ic] = unique(M, 'rows', 'stable'); % Unique Values By Row, Retaining Original Order
h = accumarray(ic, 1); % Count Occurrences
% maph = h(ic); % Map Occurrences To ‘ic’ Values
disp([Mu, h])

sigmcBICRElw = cat(2, ICsigagg.ICwcfg0_presentations.all.sigmcBICREl2, ...
    ICsigagg.ICwcfg1_presentations.all.sigmcBICREl2, ...
    ICsigagg.ICwcfg0_presentations.all.sigmcBICREl1, ...
    ICsigagg.ICwcfg1_presentations.all.sigmcBICREl1);
mean(any(sigmcBICRElw(neuoi, 1:2:end)==1, 2))

fprintf('%.4f vs %.4f is the proportion of IC vs REl responsive neurons among %s neurons\n', ...
    mean(any(sigmcBICRElw(neuoi, 1:2:end)==1, 2)), mean(any(sigmcBICRElw(neuoi, 2:2:end)==1, 2)), neudesc )

%% 2d-histogram preferred IC vs REl: IC-responsive center-RF neurons
neuinarea = neuvisaggall==find(strcmp(vislocs, 'VISp'));
neuoi = neuinarea & neuctr & any(sigmcBICRElw(:, 1:2:end)==1, 2);% & any(sigmcBICRElw(:, 2:2:end)==1, 2);

figure('Position', [100 100 300 240])
histogram2(prefiwICagg(neuoi), ori4paramsagg.all.prefiori4(neuoi), 'displaystyle', 'tile', 'showemptybins', 'on')
colorbar
colormap redblue
set(gca,'FontSize', fs, 'XGrid', 'off', 'YGrid', 'off', 'YDir', 'reverse', 'XTick', 1:4, 'XTickLabel', 0:45:180-1, 'YTick', 1:4, 'YTickLabel', 0:45:180-1)
xlabel('Pref. I_C Orientation')
ylabel('Pref. Grating Orientation')
% title(sprintf('V1 %s (fixed-gaze, N=%d)', neuctrdesc, nnz(neuinarea & neuctr) ), 'fontweight', 'normal')
fprintf('V1 IC-responsive %s (N=%d/%d)\n', neuctrdesc, nnz(neuoi), nnz(neuinarea & neuctr) )


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
