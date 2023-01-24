% uses variables aggregated in aggregate_psth -- but these are loaded here,
% so no need to run aggregate_psth before running this script

load(['G:\My Drive\DATA\ICexpts_submission22\openscope_popavg_agg.mat'])

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
%% V1: pref-ori distribution of IC/RC encoders
% RUN AFTER RUNNING FIRST TWO SECTIONS OF aggregate_psth
whichprobe = 'all';
neuinarea = true(size(neuprobeagg));

% whichprobe = 'C';
% iprobe = find(strcmp(probes, whichprobe));
% % neuinarea = contains(neulocagg{iprobe}, 'VISp');
% neuinarea = neuctxagg{iprobe}==1;
% neuinarea = strcmp(neulocagg{iprobe}, 'VISp2/3');
% if ~isequal(contains(neulocagg{iprobe}, 'VISp'), neuctxagg{iprobe}==1)
%     warning('PROBE AND VISUAL AREAS DO NOT ALWAYS MATCH UP!!')
% end

opt678 = true;
if opt678
prefidiragg = oriparamsagg.(whichprobe).prefiori678;
prefioriagg = ori4paramsagg.(whichprobe).prefiori4678;
sigori = ori4paramsagg.(whichprobe).Pkw_ori4678<0.05 & ori4paramsagg.(whichprobe).OP4678>0.5;
% Rori = ( oriparamsagg.(whichprobe).Rori678(:,1:4) + oriparamsagg.(whichprobe).Rori678(:,5:8) )/2;
% [~,prefioriagg] = max(Rori, [], 2);
Rori = ori4paramsagg.(whichprobe).Rori4678(:,1:4);
else
prefidiragg = oriparamsagg.(whichprobe).prefiori;
prefioriagg = ori4paramsagg.(whichprobe).prefiori4;
sigori = ori4paramsagg.(whichprobe).Pkw_ori4<0.05; % & ori4paramsagg.(whichprobe).OSI4>0.5;
Rori = ori4paramsagg.(whichprobe).Rori4(:,1:4);
end

normopt = false;
if normopt
    normRori = Rori./max(Rori,[],2);
else
    normRori = Rori;
end
fs = 12;
xtorder = [4 1 2 3]; xtl = [-45 0 45 90];
xtorder = 1:4; xtl = orivec;
figure
for ienc = 1:4
    subplot(2,2,ienc)
    hold all
    switch ienc
        case 1
            neuoi = ICsigagg.ICwcfg0_presentations.(whichprobe).ICencoder2==1;
            neuoi = ismember(ICsigagg.ICwcfg0_presentations.(whichprobe).sigmcBK, [0 0 0 1], 'rows');
        case 2
            neuoi = ICsigagg.ICwcfg1_presentations.(whichprobe).ICencoder2==1;
            neuoi = ismember(ICsigagg.ICwcfg1_presentations.(whichprobe).sigmcBK, [0 0 0 1], 'rows');
        case 3
            neuoi = ICsigagg.ICwcfg0_presentations.(whichprobe).ICencoder1==1;
            neuoi = ismember(ICsigagg.ICwcfg0_presentations.(whichprobe).sigmcBK, [1 0 0 0], 'rows');
        case 4
            neuoi = ICsigagg.ICwcfg1_presentations.(whichprobe).ICencoder1==1;
            neuoi = ismember(ICsigagg.ICwcfg1_presentations.(whichprobe).sigmcBK, [1 0 0 0], 'rows');
    end
    %     neuoi = neuoi & neuinarea & sigori;
    errorbar(xtl, nanmean(normRori(neuoi,xtorder),1), nanstd(normRori(neuoi,xtorder),0,1)/sqrt(nnz(neuoi)), 'ko-', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'LineWidth', 2)
    xlim([xtl(1) xtl(end)])
    set(gca, 'XTick', xtl, 'FontSize', fs)
    if normopt
        ylabel('Normalized Firing Rate', 'FontSize', fs)
    else
        ylabel('Firing Rate (Hz)', 'FontSize', fs)
    end
    xlabel('Orientation (deg)', 'FontSize', fs)
    title(sprintf('%d-deg IC-encoder (N=%d)', orivec(ienc), nnz(neuoi)), 'FontSize', fs)
end

fs = 14;
xtorder = [4 1 2 3];
xtl = [-45 0 45 90];
figure('Position', [100 100 300 300])
hold all
for ienc = 1%:4
    switch ienc
        case 1
            neuoi = ICsigagg.ICwcfg0_presentations.(whichprobe).ICencoder2==1;
            neuoi = ismember(ICsigagg.ICwcfg0_presentations.(whichprobe).sigmcBK, [0 0 0 1], 'rows');
        case 2
            neuoi = ICsigagg.ICwcfg1_presentations.(whichprobe).ICencoder2==1;
            neuoi = ismember(ICsigagg.ICwcfg1_presentations.(whichprobe).sigmcBK, [0 0 0 1], 'rows');
        case 3
            neuoi = ICsigagg.ICwcfg0_presentations.(whichprobe).ICencoder1==1;
            neuoi = ismember(ICsigagg.ICwcfg0_presentations.(whichprobe).sigmcBK, [1 0 0 0], 'rows');
        case 4
            neuoi = ICsigagg.ICwcfg1_presentations.(whichprobe).ICencoder1==1;
            neuoi = ismember(ICsigagg.ICwcfg1_presentations.(whichprobe).sigmcBK, [1 0 0 0], 'rows');
    end
%     neuoi = neuoi & neuinarea & sigori;
    errorbar(xtl, mean(Rori(neuoi,xtorder),1), std(Rori(neuoi,xtorder),0,1)/sqrt(nnz(neuoi)), 'ko-', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'LineWidth', 2)
end
xlim([xtl(1) xtl(end)])
set(gca, 'XTick', xtl, 'FontSize', fs)
ylabel('Firing Rate (Hz)', 'FontSize', fs)
xlabel('Orientation (deg)', 'FontSize', fs)
title(sprintf('0-deg IC-encoder (N=%d)', nnz(neuoi)), 'FontSize', fs)

% % sanity check
% figure; histogram2(ori4paramsagg.(whichprobe).prefiori4(sigori), ori4paramsagg.(whichprobe).prefiori4678(sigori), 'displaystyle', 'tile')

hcprefidir_wICenc = zeros(4, Ndirs);
hcprefiori_wICenc = zeros(4, Noris);
hcsigprefidir_wICenc = zeros(4, Ndirs);
hcsigprefiori_wICenc = zeros(4, Noris);
for ienc = 1:4
    switch ienc
        case 1
            neuoi = ICsigagg.ICwcfg0_presentations.(whichprobe).ICencoder2==1;
            neuoi = ismember(ICsigagg.ICwcfg0_presentations.(whichprobe).sigmcBK, [0 0 0 1], 'rows');
        case 2
            neuoi = ICsigagg.ICwcfg1_presentations.(whichprobe).ICencoder2==1;
            neuoi = ismember(ICsigagg.ICwcfg1_presentations.(whichprobe).sigmcBK, [0 0 0 1], 'rows');
        case 3
            neuoi = ICsigagg.ICwcfg0_presentations.(whichprobe).ICencoder1==1;
            neuoi = ismember(ICsigagg.ICwcfg0_presentations.(whichprobe).sigmcBK, [1 0 0 0], 'rows');
        case 4
            neuoi = ICsigagg.ICwcfg1_presentations.(whichprobe).ICencoder1==1;
            neuoi = ismember(ICsigagg.ICwcfg1_presentations.(whichprobe).sigmcBK, [1 0 0 0], 'rows');
    end
    neuoi = neuoi & neuinarea;
    hcprefidir_wICenc(ienc,:) = histcounts(prefidiragg(neuoi), 0.5:Ndirs+0.5);
    hcprefiori_wICenc(ienc,:) = histcounts(prefioriagg(neuoi), 0.5:Noris+0.5);
    hcsigprefidir_wICenc(ienc,:) = histcounts(prefidiragg(sigori & neuoi), 0.5:Ndirs+0.5);
    hcsigprefiori_wICenc(ienc,:) = histcounts(prefioriagg(sigori & neuoi), 0.5:Noris+0.5);
end


hcprefidir_kICenc = zeros(4, Ndirs);
hcprefiori_kICenc = zeros(4, Noris);
hcsigprefidir_kICenc = zeros(4, Ndirs);
hcsigprefiori_kICenc = zeros(4, Noris);
for ienc = 1:4
    switch ienc
        case 1
            neuoi = ICsigagg.ICkcfg0_presentations.(whichprobe).ICencoder2==1;
            neuoi = ismember(ICsigagg.ICkcfg0_presentations.(whichprobe).sigmcBK, [0 0 0 1], 'rows');
        case 2
            neuoi = ICsigagg.ICkcfg1_presentations.(whichprobe).ICencoder2==1;
            neuoi = ismember(ICsigagg.ICkcfg1_presentations.(whichprobe).sigmcBK, [0 0 0 1], 'rows');
        case 3
            neuoi = ICsigagg.ICkcfg0_presentations.(whichprobe).ICencoder1==1;
            neuoi = ismember(ICsigagg.ICkcfg0_presentations.(whichprobe).sigmcBK, [1 0 0 0], 'rows');
        case 4
            neuoi = ICsigagg.ICkcfg1_presentations.(whichprobe).ICencoder1==1;
            neuoi = ismember(ICsigagg.ICkcfg1_presentations.(whichprobe).sigmcBK, [1 0 0 0], 'rows');
    end
    neuoi = neuoi & neuinarea;
    hcprefidir_kICenc(ienc,:) = histcounts(prefidiragg(neuoi), 0.5:Ndirs+0.5);
    hcprefiori_kICenc(ienc,:) = histcounts(prefioriagg(neuoi), 0.5:Noris+0.5);
    hcsigprefidir_kICenc(ienc,:) = histcounts(prefidiragg(sigori & neuoi), 0.5:Ndirs+0.5);
    hcsigprefiori_kICenc(ienc,:) = histcounts(prefioriagg(sigori & neuoi), 0.5:Noris+0.5);
end

disp('hcsigprefiori_wICenc')
disp([hcsigprefiori_wICenc sum(hcsigprefiori_wICenc,2)])
disp('hcsigprefiori_kICenc')
disp([hcsigprefiori_kICenc sum(hcsigprefiori_kICenc,2)])

hcprefidir = histcounts(prefidiragg, 0.5:Ndirs+0.5);
hcprefiori = histcounts(prefioriagg, 0.5:Noris+0.5);
hcsigprefidir = histcounts(prefidiragg(sigori), 0.5:Ndirs+0.5);
hcsigprefiori = histcounts(prefioriagg(sigori), 0.5:Noris+0.5);

fs=10;
ytl = 0:45:180-1; ylab = 'IC Orientation';

figure('Position', [0 0 1000 500])
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'All Trials', 'edgecolor', 'none')
subplot(2,3,1); imagesc(hcprefiori_wICenc)
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
title('White Inducers Black Background')
subplot(2,3,2); imagesc(hcprefiori_wICenc./sum(hcprefiori_wICenc,2))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
title('rows sum to 1')
subplot(2,3,3); imagesc(hcprefiori_wICenc./sum(hcprefiori_wICenc,1))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
title('columns sum to 1')
subplot(2,3,4); imagesc(hcprefiori_kICenc)
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
title('Black Inducers White Background')
subplot(2,3,5); imagesc(hcprefiori_kICenc./sum(hcprefiori_kICenc,2))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
title('rows sum to 1')
subplot(2,3,6); imagesc(hcprefiori_kICenc./sum(hcprefiori_kICenc,1))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
title('columns sum to 1')

figure('Position', [0 500 1000 500])
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'All Trials: Significantly Tuned Neurons', 'edgecolor', 'none')
subplot(2,3,1); imagesc(hcsigprefiori_wICenc)
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
title('White Inducers Black Background')
subplot(2,3,2); imagesc(hcsigprefiori_wICenc./sum(hcsigprefiori_wICenc,2))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
title('rows sum to 1')
subplot(2,3,3); imagesc(hcsigprefiori_wICenc./sum(hcsigprefiori_wICenc,1))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
title('columns sum to 1')
subplot(2,3,4); imagesc(hcsigprefiori_kICenc)
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
title('Black Inducers White Background')
subplot(2,3,5); imagesc(hcsigprefiori_kICenc./sum(hcsigprefiori_kICenc,2))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
title('rows sum to 1')
subplot(2,3,6); imagesc(hcsigprefiori_kICenc./sum(hcsigprefiori_kICenc,1))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
title('columns sum to 1')

%% bootstrapped confidence interval: 
whichprobe = 'all';
neuinarea = true(size(neuprobeagg));

% whichprobe = 'C';
% iprobe = find(strcmp(probes, whichprobe));
% % neuinarea = strcmp(neulocagg{iprobe}, 'VISp2/3');
% % neuinarea = contains(neulocagg{iprobe}, 'VISp');
% % neuinarea = neuctxagg{iprobe}==1;
% neuinarea = true(size(neuctxagg{iprobe}));

opt678 = true;
if opt678
prefidiragg = oriparamsagg.(whichprobe).prefiori678;
prefioriagg = ori4paramsagg.(whichprobe).prefiori4678;
sigori = ori4paramsagg.(whichprobe).Pkw_ori4678<0.05;
% Rori = ( oriparamsagg.(whichprobe).Rori(:,1:4) + oriparamsagg.(whichprobe).Rori(:,5:8) )/2;
% [~,prefioriagg] = max(Rori, [], 2);
Rori = ori4paramsagg.(whichprobe).Rori4678;
else
prefidiragg = oriparamsagg.(whichprobe).prefiori;
prefioriagg = ori4paramsagg.(whichprobe).prefiori4;
sigori = ori4paramsagg.(whichprobe).Pkw_ori4<0.05;
Rori = ori4paramsagg.(whichprobe).Rori4;
end

NICenc = zeros(8,1);
NICenc_sigori = zeros(8,1);

hpICenc = cell(8,1);
hpICenc_sigori = cell(8,1);
propICenc = cell(8,1);
propICenc_sigori = cell(8,1);

% rows 1-3 are mean, negative errorbar, positive errorbar
Nboot = 1000;
hpbootCI = cell(8,1);
hpbootCI_sigori = cell(8,1);
propbootCI = cell(8,1);
propbootCI_sigori = cell(8,1);
for ienc = 1:8
switch ienc
    case 1
        neuoi = ICsigagg.ICwcfg0_presentations.(whichprobe).ICencoder2==1;
        neuoi = ismember(ICsigagg.ICwcfg0_presentations.(whichprobe).sigmcBK, [0 0 0 1], 'rows');
    case 2
        neuoi = ICsigagg.ICwcfg1_presentations.(whichprobe).ICencoder2==1;
        neuoi = ismember(ICsigagg.ICwcfg1_presentations.(whichprobe).sigmcBK, [0 0 0 1], 'rows');
    case 3
        neuoi = ICsigagg.ICwcfg0_presentations.(whichprobe).ICencoder1==1;
        neuoi = ismember(ICsigagg.ICwcfg0_presentations.(whichprobe).sigmcBK, [1 0 0 0], 'rows');
    case 4
        neuoi = ICsigagg.ICwcfg1_presentations.(whichprobe).ICencoder1==1;
        neuoi = ismember(ICsigagg.ICwcfg1_presentations.(whichprobe).sigmcBK, [1 0 0 0], 'rows');
    case 5
        neuoi = ICsigagg.ICkcfg0_presentations.(whichprobe).ICencoder2==1;
        neuoi = ismember(ICsigagg.ICkcfg0_presentations.(whichprobe).sigmcBK, [0 0 0 1], 'rows');
    case 6
        neuoi = ICsigagg.ICkcfg1_presentations.(whichprobe).ICencoder2==1;
        neuoi = ismember(ICsigagg.ICkcfg1_presentations.(whichprobe).sigmcBK, [0 0 0 1], 'rows');
    case 7
        neuoi = ICsigagg.ICkcfg0_presentations.(whichprobe).ICencoder1==1;
        neuoi = ismember(ICsigagg.ICkcfg0_presentations.(whichprobe).sigmcBK, [1 0 0 0], 'rows');
    case 8
        neuoi = ICsigagg.ICkcfg1_presentations.(whichprobe).ICencoder1==1;
        neuoi = ismember(ICsigagg.ICkcfg1_presentations.(whichprobe).sigmcBK, [1 0 0 0], 'rows');
end
    neuoi = neuoi & neuinarea;
neuoind = find(neuoi);
sigoriind = find(sigori);
sigprefioriagg = prefioriagg(sigori);

NICenc(ienc) = nnz(neuoi);
NICenc_sigori(ienc) = nnz(sigori & neuoi);

hpICenc{ienc} = histcounts(prefioriagg(neuoi), 0.5:Noris+0.5, 'normalization', 'probability');
hpICenc_sigori{ienc} = histcounts(prefioriagg(sigori & neuoi), 0.5:Noris+0.5, 'normalization', 'probability');

hcprefiori = histcounts(prefioriagg, 0.5:Noris+0.5);
hcprefiori_ICenc = histcounts(prefioriagg(neuoi), 0.5:Noris+0.5);
propICenc{ienc} = hcprefiori_ICenc./hcprefiori;

hcprefiori = histcounts(prefioriagg(sigori), 0.5:Noris+0.5);
hcprefiori_ICenc = histcounts(prefioriagg(sigori & neuoi), 0.5:Noris+0.5);
propICenc_sigori{ienc} = hcprefiori_ICenc./hcprefiori;

hpboot = NaN(Nboot,Noris);
hpboot_sigori = NaN(Nboot,Noris);
propboot = NaN(Nboot,Noris);
propboot_sigori = NaN(Nboot,Noris);
hcprefiori = histcounts(prefioriagg, 0.5:Noris+0.5);
hcprefiori_sigori = histcounts(prefioriagg(sigori), 0.5:Noris+0.5);
% tic
for iboot = 1:Nboot
    tempneu = false(size(sigori));
    tempneu( randperm(length(prefioriagg), nnz(neuoi)) ) = true;
    hpboot(iboot,:) = histcounts(prefioriagg(tempneu), 0.5:Noris+0.5, 'normalization', 'probability');
    
    tempneusig = false(size(sigprefioriagg));
    tempneusig( randperm(length(sigprefioriagg), nnz(neuoi(sigoriind))) ) = true;
    hpboot_sigori(iboot,:) = histcounts(sigprefioriagg(tempneusig), 0.5:Noris+0.5, 'normalization', 'probability');
    
    hcprefiori_ICenc = histcounts(prefioriagg(tempneu), 0.5:Noris+0.5);
    propboot(iboot,:) = hcprefiori_ICenc./hcprefiori;
    
    hcprefiori_ICenc = histcounts(sigprefioriagg(tempneusig), 0.5:Noris+0.5);
    propboot_sigori(iboot,:) = hcprefiori_ICenc./hcprefiori_sigori;
end
% toc
hpbootCI{ienc} = [nanmean(hpboot,1); nanmean(hpboot,1)-prctile(hpboot,2.5,1); prctile(hpboot,97.5,1)-nanmean(hpboot,1)];
hpbootCI_sigori{ienc} = [nanmean(hpboot_sigori,1); nanmean(hpboot_sigori,1)-prctile(hpboot_sigori,2.5,1); prctile(hpboot_sigori,97.5,1)-nanmean(hpboot_sigori,1)];

propbootCI{ienc} = [nanmean(propboot,1); nanmean(propboot,1)-prctile(propboot,2.5,1); prctile(propboot,97.5,1)-nanmean(propboot,1)];
propbootCI_sigori{ienc} = [nanmean(propboot_sigori,1); nanmean(propboot_sigori,1)-prctile(propboot_sigori,2.5,1); prctile(propboot_sigori,97.5,1)-nanmean(propboot_sigori,1)];

end

fs = 12;
enctitle = {'wICenc0', 'wICenc45', 'wICenc90', 'wICenc135', 'kICenc0', 'kICenc45', 'kICenc90', 'kICenc135'};

figure('Position', [0 0 1000 1000])
annotation('textbox', [0.1 0.9 0.9 0.1], 'string', 'All Trials', 'edgecolor', 'none', 'fontsize', fs)
for ienc = 1:8
subplot(4,4,ienc)
hold all
errorbar(1:Noris, hpbootCI{ienc}(1,:), hpbootCI{ienc}(2,:), hpbootCI{ienc}(3,:) , 'k')
plot(1:Noris, hpICenc{ienc}, 'r', 'LineWidth', 2)
indsigori = find(hpICenc{ienc}>hpbootCI{ienc}(1,:)+hpbootCI{ienc}(3,:));
plot(indsigori, hpICenc{ienc}(indsigori), 'ro', 'MarkerFaceColor', 'r', 'LineWidth', 2)
set(gca, 'FontSize', fs, 'XTick', 1:Noris, 'XTickLabel', orivec)
title(sprintf('%s N=%d\npref ori dist', enctitle{ienc}, NICenc(ienc) ))
ylabel('Probability')
xlim([0.5 Noris+0.5])
xlabel('Pref Ori')

subplot(4,4,ienc+8)
hold all
errorbar(1:Noris, propbootCI{ienc}(1,:), propbootCI{ienc}(2,:), propbootCI{ienc}(3,:) , 'k')
plot(1:Noris, propICenc{ienc}, 'r', 'LineWidth', 2)
indsigori = find(propICenc{ienc}>propbootCI{ienc}(1,:)+propbootCI{ienc}(3,:));
plot(indsigori, propICenc{ienc}(indsigori), 'ro', 'MarkerFaceColor', 'r', 'LineWidth', 2)
set(gca, 'FontSize', fs, 'XTick', 1:Noris, 'XTickLabel', orivec)
title(sprintf('%s N=%d\nprop for each pref ori', enctitle{ienc}, NICenc_sigori(ienc) ))
ylabel('Proportion')
xlim([0.5 Noris+0.5])
xlabel('Pref Ori')
end

figure('Position', [0 0 1000 1000])
annotation('textbox', [0.1 0.9 0.9 0.1], 'string', 'All Trials: Significantly Tuned Neurons', 'edgecolor', 'none', 'fontsize', fs)
for ienc = 1:8
subplot(4,4,ienc)
hold all
errorbar(1:Noris, hpbootCI_sigori{ienc}(1,:), hpbootCI_sigori{ienc}(2,:), hpbootCI_sigori{ienc}(3,:) , 'k')
plot(1:Noris, hpICenc_sigori{ienc}, 'r', 'LineWidth', 2)
indsigori = find(hpICenc_sigori{ienc}>hpbootCI_sigori{ienc}(1,:)+hpbootCI_sigori{ienc}(3,:));
plot(indsigori, hpICenc_sigori{ienc}(indsigori), 'ro', 'MarkerFaceColor', 'r', 'LineWidth', 2)
set(gca, 'FontSize', fs, 'XTick', 1:Noris, 'XTickLabel', orivec)
title(sprintf('%s N=%d\npref ori dist', enctitle{ienc}, NICenc(ienc) ))
ylabel('Probability')
xlim([0.5 Noris+0.5])
xlabel('Pref Ori')

subplot(4,4,ienc+8)
hold all
errorbar(1:Noris, propbootCI_sigori{ienc}(1,:), propbootCI_sigori{ienc}(2,:), propbootCI_sigori{ienc}(3,:) , 'k')
plot(1:Noris, propICenc_sigori{ienc}, 'r', 'LineWidth', 2)
indsigori = find(propICenc_sigori{ienc}>propbootCI_sigori{ienc}(1,:)+propbootCI_sigori{ienc}(3,:));
plot(indsigori, propICenc_sigori{ienc}(indsigori), 'ro', 'MarkerFaceColor', 'r', 'LineWidth', 2)
set(gca, 'FontSize', fs, 'XTick', 1:Noris, 'XTickLabel', orivec)
title(sprintf('%s N=%d\nprop for each pref ori', enctitle{ienc}, NICenc_sigori(ienc) ))
ylabel('Proportion')
xlim([0.5 Noris+0.5])
xlabel('Pref Ori')
end



%%
load('D:\OpenScopeData\000248v221123\postprocessed\psthavgagg.mat', 'ICblocks', 'ICtrialtypes', 'Ronavgagg')
for b = 1:numel(visblocks)
    Ronavgagg.(visblocks{b}).all = [];
    for iprobe = 1:numel(probes)
        Ronavgagg.(visblocks{b}).all = cat(2, Ronavgagg.(visblocks{b}).all, Ronavgagg.(visblocks{b}).(probes{iprobe}));
    end
end

%% orientation tuning of IC-encoders
whichprobe = 'all';
neuinarea = true(size(neuprobeagg));

% whichprobe = 'C';
% iprobe = find(strcmp(probes, whichprobe));
% % neuinarea = strcmp(neulocagg{iprobe}, 'VISp2/3');
% % neuinarea = contains(neulocagg{iprobe}, 'VISp');
% neuinarea = neuctxagg{iprobe}==1;

sigori = ori4paramsagg.(whichprobe).Pkw_ori4<0.05;

wICenc = neuinarea & (ICsigagg.ICwcfg0_presentations.(whichprobe).ICencoder | ICsigagg.ICwcfg1_presentations.(whichprobe).ICencoder);
kICenc = neuinarea & (ICsigagg.ICkcfg0_presentations.(whichprobe).ICencoder | ICsigagg.ICkcfg1_presentations.(whichprobe).ICencoder);

wRCenc = neuinarea & (ICsigagg.ICwcfg0_presentations.(whichprobe).RCencoder | ICsigagg.ICwcfg1_presentations.(whichprobe).RCencoder);
kRCenc = neuinarea & (ICsigagg.ICkcfg0_presentations.(whichprobe).RCencoder | ICsigagg.ICkcfg1_presentations.(whichprobe).RCencoder);

wemerge = neuinarea & ( (ICsigagg.ICwcfg0_presentations.(whichprobe).PkwBK<0.05 & ICsigagg.ICwcfg0_presentations.(whichprobe).PkwBI>0.05) ...
    | (ICsigagg.ICwcfg1_presentations.(whichprobe).PkwBK<0.05 & ICsigagg.ICwcfg1_presentations.(whichprobe).PkwBI>0.05) );
kemerge = neuinarea & ( (ICsigagg.ICkcfg0_presentations.(whichprobe).PkwBK<0.05 & ICsigagg.ICkcfg0_presentations.(whichprobe).PkwBI>0.05) ...
    | (ICsigagg.ICkcfg1_presentations.(whichprobe).PkwBK<0.05 & ICsigagg.ICkcfg1_presentations.(whichprobe).PkwBI>0.05) );

wsigkw = neuinarea & ( (ICsigagg.ICwcfg0_presentations.(whichprobe).PkwBK<0.05 & ICsigagg.ICwcfg0_presentations.(whichprobe).PkwBI<0.05) ...
    | (ICsigagg.ICwcfg1_presentations.(whichprobe).PkwBK<0.05 & ICsigagg.ICwcfg1_presentations.(whichprobe).PkwBI>0.05) );
ksigkw = neuinarea & ( (ICsigagg.ICkcfg0_presentations.(whichprobe).PkwBK<0.05 & ICsigagg.ICkcfg0_presentations.(whichprobe).PkwBI<0.05) ...
    | (ICsigagg.ICkcfg1_presentations.(whichprobe).PkwBK<0.05 & ICsigagg.ICkcfg1_presentations.(whichprobe).PkwBI>0.05) );

% 0 cfg0tt111 45 cfg1tt111 90 cfg0tt106 135 cfg1tt106
RonwICavgagg = [Ronavgagg.ICwcfg0_presentations.(whichprobe)(ICtrialtypes==111,:); ...
    Ronavgagg.ICwcfg1_presentations.(whichprobe)(ICtrialtypes==111,:); ...
    Ronavgagg.ICwcfg0_presentations.(whichprobe)(ICtrialtypes==106,:); ...
    Ronavgagg.ICwcfg1_presentations.(whichprobe)(ICtrialtypes==106,:)]';
RonkICavgagg = [Ronavgagg.ICkcfg0_presentations.(whichprobe)(ICtrialtypes==111,:); ...
    Ronavgagg.ICkcfg1_presentations.(whichprobe)(ICtrialtypes==111,:); ...
    Ronavgagg.ICkcfg0_presentations.(whichprobe)(ICtrialtypes==106,:); ...
    Ronavgagg.ICkcfg1_presentations.(whichprobe)(ICtrialtypes==106,:)]';
[~, prefwIC] = max(RonwICavgagg, [], 2);
[~, prefkIC] = max(RonkICavgagg, [], 2);

RonwICRCavgagg = [Ronavgagg.ICwcfg0_presentations.(whichprobe)(ICtrialtypes==111,:); ...
    Ronavgagg.ICwcfg1_presentations.(whichprobe)(ICtrialtypes==111,:); ...
    Ronavgagg.ICwcfg0_presentations.(whichprobe)(ICtrialtypes==106,:); ...
    Ronavgagg.ICwcfg1_presentations.(whichprobe)(ICtrialtypes==106,:); ...
    Ronavgagg.ICwcfg0_presentations.(whichprobe)(ICtrialtypes==110,:); ...
    Ronavgagg.ICwcfg1_presentations.(whichprobe)(ICtrialtypes==110,:); ...
    Ronavgagg.ICwcfg0_presentations.(whichprobe)(ICtrialtypes==107,:); ...
    Ronavgagg.ICwcfg1_presentations.(whichprobe)(ICtrialtypes==107,:)]';
RonkICRCavgagg = [Ronavgagg.ICkcfg0_presentations.(whichprobe)(ICtrialtypes==111,:); ...
    Ronavgagg.ICkcfg1_presentations.(whichprobe)(ICtrialtypes==111,:); ...
    Ronavgagg.ICkcfg0_presentations.(whichprobe)(ICtrialtypes==106,:); ...
    Ronavgagg.ICkcfg1_presentations.(whichprobe)(ICtrialtypes==106,:); ...
    Ronavgagg.ICkcfg0_presentations.(whichprobe)(ICtrialtypes==110,:); ...
    Ronavgagg.ICkcfg1_presentations.(whichprobe)(ICtrialtypes==110,:); ...
    Ronavgagg.ICkcfg0_presentations.(whichprobe)(ICtrialtypes==107,:); ...
    Ronavgagg.ICkcfg1_presentations.(whichprobe)(ICtrialtypes==107,:)]';
[~, prefwICRC] = max(RonwICRCavgagg, [], 2);
[~, prefkICRC] = max(RonkICRCavgagg, [], 2);

Ronwindinagg = cat(3,Ronavgagg.ICwcfg0_presentations.(whichprobe)(ismember(ICtrialtypes,[1302 1304]),:), ...
    Ronavgagg.ICwcfg1_presentations.(whichprobe)(ismember(ICtrialtypes,[1302 1304]),:), ...
    Ronavgagg.ICwcfg0_presentations.(whichprobe)(ismember(ICtrialtypes,[1301 1303]),:), ...
    Ronavgagg.ICwcfg1_presentations.(whichprobe)(ismember(ICtrialtypes,[1301 1303]),:));
Ronwindinagg = squeeze(mean(Ronwindinagg,1));
Ronkindinagg = cat(3,Ronavgagg.ICkcfg0_presentations.(whichprobe)(ismember(ICtrialtypes,[1302 1304]),:), ...
    Ronavgagg.ICkcfg1_presentations.(whichprobe)(ismember(ICtrialtypes,[1302 1304]),:), ...
    Ronavgagg.ICkcfg0_presentations.(whichprobe)(ismember(ICtrialtypes,[1301 1303]),:), ...
    Ronavgagg.ICkcfg1_presentations.(whichprobe)(ismember(ICtrialtypes,[1301 1303]),:));
Ronkindinagg = squeeze(mean(Ronkindinagg,1));
[~, prefwindin] = max(Ronwindinagg, [], 2);
[~, prefkindin] = max(Ronkindinagg, [], 2);

neuoi = wICenc & sigori;
figure; histogram2(prefkICRC(neuoi), ori4paramsagg.(whichprobe).prefiori4(neuoi), 'displaystyle', 'tile')
xlabel('Pref wIC Ori')
ylabel('Pref Grating Ori')

% hcall = histcounts2(prefwIC, ori4paramsagg.(whichprobe).prefiori4, 0.5:1:4.5, 0.5:1:4.5);
% hcwIC = histcounts2(prefwIC(wICenc), ori4paramsagg.(whichprobe).prefiori4(wICenc), 0.5:1:4.5, 0.5:1:4.5);
% figure; imagesc(hcwIC./hcall)

figure('Position', [0 0 2000 600])
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'All Trials: Significantly Tuned Neurons (rows sum to 1)', 'edgecolor', 'none')
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
hc = histcounts2(prefwIC(neuoi), ori4paramsagg.(whichprobe).prefiori4(neuoi), 0.5:4.5, 0.5:4.5);
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
hc = histcounts2(prefkIC(neuoi), ori4paramsagg.(whichprobe).prefiori4(neuoi), 0.5:4.5, 0.5:4.5);
imagesc( (hc./sum(hc,2))' ); colorbar
axis([0.5 4.5 0.5 4.5])
ylabel('Pref Grating Ori')
xlabel('Pref kIC')
title(sprintf('%s N=%d', neudesc, nnz(neuoi)))
end

figure('Position', [0 0 2000 600])
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'All Trials: Significantly Tuned Neurons (rows sum to 1)', 'edgecolor', 'none')
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
hc = histcounts2(prefwindin(neuoi), ori4paramsagg.(whichprobe).prefiori4(neuoi), 0.5:4.5, 0.5:4.5);
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
hc = histcounts2(prefkindin(neuoi), ori4paramsagg.(whichprobe).prefiori4(neuoi), 0.5:4.5, 0.5:4.5);
imagesc( (hc./sum(hc,2))' ); colorbar
axis([0.5 4.5 0.5 4.5])
ylabel('Pref Grating Ori')
xlabel('Pref kPac-In')
title(sprintf('%s N=%d', neudesc, nnz(neuoi)))
end

figure('Position', [0 0 2000 600])
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'All Trials: Significantly Tuned Neurons (rows sum to 1)', 'edgecolor', 'none')
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
hc = histcounts2(prefwICRC(neuoi), ori4paramsagg.(whichprobe).prefiori4(neuoi), 0.5:8.5, 0.5:4.5);
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
hc = histcounts2(prefkICRC(neuoi), ori4paramsagg.(whichprobe).prefiori4(neuoi), 0.5:8.5, 0.5:4.5);
imagesc( (hc./sum(hc,2))' ); colorbar
plot([4.5 4.5], [0.5 4.5], 'w-', 'Linewidth', 2)
axis([0.5 8.5 0.5 4.5])
ylabel('Pref Grating Ori')
xlabel('Pref kIC/RC')
title(sprintf('%s N=%d', neudesc, nnz(neuoi)))
end

% does OP correlate with SP_ICvsRC for "emergent" cells? -- no, not at all obvious
neuoi = wemerge & (ori4paramsagg.(whichprobe).prefiori4==1 | ori4paramsagg.(whichprobe).prefiori4==3);
corr(ori4paramsagg.(whichprobe).OP4678(neuoi), ICsigagg.ICwcfg0_presentations.(whichprobe).SP_ICvsRC(neuoi), 'type', 'spearman')
figure; plot(ori4paramsagg.(whichprobe).OP4678(neuoi), ICsigagg.ICwcfg0_presentations.(whichprobe).SP_ICvsRC(neuoi), 'o')
figure; histogram2(ori4paramsagg.(whichprobe).OP4678(neuoi), ICsigagg.ICwcfg0_presentations.(whichprobe).SP_ICvsRC(neuoi), 'displaystyle', 'tile')

neuoi = wemerge & (ori4paramsagg.(whichprobe).prefiori4==2 | ori4paramsagg.(whichprobe).prefiori4==4);
corr(ori4paramsagg.(whichprobe).OP4678(neuoi), ICsigagg.ICwcfg0_presentations.(whichprobe).SP_ICvsRC(neuoi), 'type', 'spearman')
figure; plot(ori4paramsagg.(whichprobe).OP4678(neuoi), ICsigagg.ICwcfg1_presentations.(whichprobe).SP_ICvsRC(neuoi), 'o')
figure; histogram2(ori4paramsagg.(whichprobe).OP4678(neuoi), ICsigagg.ICwcfg1_presentations.(whichprobe).SP_ICvsRC(neuoi), 'displaystyle', 'tile')
