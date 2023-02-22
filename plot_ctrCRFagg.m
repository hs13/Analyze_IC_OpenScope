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

Nrfs = size(RFCIagg.all.Rrfclassic, 2);
Nszs = size(sizeCIagg.all.Rsizeclassic, 2);
Ndirs = size(oriparamsagg.all.Rori, 2);
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

neudesc = 'ctrCG9';
switch neudesc
    case 'ctrCG9'
        neuctr = RFCIagg.all.RFindclassic==1 & RFCIagg.all.pRFclassic<0.05;
    case 'ctrCRF9'
        neuctr = RFCIagg.all.RFindclassic==1 & RFCIagg.all.Pkw_rfclassic<0.05;
    case 'ctrCRF9sigexcl'
        neuctr = RFCIagg.all.RFindclassic==1 & RFCIagg.all.RFsigexclclassic==1;
    case 'ctrCRF9exclsig'
        neuctr = RFCIagg.all.RFindclassic==1 & RFCIagg.all.RFexclsigclassic==1;
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
ttinds = [506 511; 107 110; 106 111];
ttcol = [0 0 1; 1 0.5 0; 0 0.7 0];
yl = [2 18];
figure
hold all
for ii = 1:numel(probes)
    probeind = find(visind==ii);
neuoi = neuvisaggall==probeind & neuctr;
subplot(2,3,ii)
hold all
for typi = 1:3
    temppsth_ICwcfg1 = [];
    temppsth_ICwcfg0 = [];
    for iprobe = 1:numel(probes)
        ttoi = ismember(ICtrialtypes, ttinds(typi,:));
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
xlabel('Time (ms)')
ylabel('Firing Rate (Hz)')
title(sprintf('%s center-CRF Neurons N=%d', visareas{probeind}, nnz(neuoi) ))
end

% isequal( mean(convn(temppsth,kergauss,'same'),2), conv(mean(temppsth,2),kergauss,'same') )
% max(abs( mean(convn(temppsth,kergauss,'same'),2)-conv(mean(temppsth,2),kergauss,'same') ))
