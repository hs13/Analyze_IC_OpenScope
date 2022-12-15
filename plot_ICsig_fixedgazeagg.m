% uses variables aggregated in aggregate_fixedgaze -- but these are loaded here,
% so no need to run aggregate_fixedgaze before running this script

load(['G:\My Drive\DATA\ICexpts_submission22\openscope_popavg_fixedgaze_agg.mat'])
Nsessions = numel(nwbsessions);
Nrfs = RFCIfieldsize2(strcmp(RFCIfields, 'Rrfclassic'));
Nszs = sizeCIfieldsize2(strcmp(sizeCIfields, 'Rsizeclassic'));
Ndirs = oriparamsfieldsize2(strcmp(oriparamsfields, 'Rori'));
Noris = Ndirs/2;

%%
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

%% V1: pref-ori distribution of IC/RC encoders
% RUN AFTER RUNNING FIRST TWO SECTIONS OF aggregate_psth
whichprobe = 'all';
neuinarea = true(size(neuprobe_fixedgazeagg));

% whichprobe = 'C';
% iprobe = find(strcmp(probes, whichprobe));
% % neuinarea = strcmp(neuloc_fixedgazeagg{iprobe}, 'VISp2/3');
% % neuinarea = contains(neuloc_fixedgazeagg{iprobe}, 'VISp');
% neuinarea = neuctx_fixedgazeagg{iprobe}==1;
% if ~isequal(contains(neuloc_fixedgazeagg{iprobe}, 'VISp'), neuctx_fixedgazeagg{iprobe}==1)
%     warning('PROBE AND VISUAL AREAS DO NOT ALWAYS MATCH UP!!')
% end

opt678 = true;
if opt678
prefidiragg = oriparams_fixedgazeagg.(whichprobe).prefiori678;
prefioriagg = ori4params_fixedgazeagg.(whichprobe).prefiori4678;
sigori = ori4params_fixedgazeagg.(whichprobe).Pkw_ori4678<0.05 & ori4params_fixedgazeagg.(whichprobe).OP4678>0.5;
% Rori = ( oriparams_fixedgazeagg.(whichprobe).Rori678(:,1:4) + oriparams_fixedgazeagg.(whichprobe).Rori678(:,5:8) )/2;
% [~,prefiori_fixedgazeagg] = max(Rori, [], 2);
Rori = ori4params_fixedgazeagg.(whichprobe).Rori4678(:,1:4);
else
prefidiragg = oriparams_fixedgazeagg.(whichprobe).prefiori;
prefioriagg = ori4params_fixedgazeagg.(whichprobe).prefiori4;
sigori = ori4params_fixedgazeagg.(whichprobe).Pkw_ori4<0.05; % & ori4params_fixedgazeagg.(whichprobe).OSI4>0.5;
Rori = ori4params_fixedgazeagg.(whichprobe).Rori4(:,1:4);
end

% % sanity check
% figure; histogram2(ori4params_fixedgazeagg.(whichprobe).prefiori4(sigori), ori4params_fixedgazeagg.(whichprobe).prefiori4678(sigori), 'displaystyle', 'tile')

hcprefidir_wICenc = zeros(4, Ndirs);
hcprefiori_wICenc = zeros(4, Noris);
hcsigprefidir_wICenc = zeros(4, Ndirs);
hcsigprefiori_wICenc = zeros(4, Noris);
for ienc = 1:4
    switch ienc
        case 1
            neuoi = ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).ICencoder2==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).sigmcBK, [0 0 0 1], 'rows');
        case 2
            neuoi = ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).ICencoder2==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).sigmcBK, [0 0 0 1], 'rows');
        case 3
            neuoi = ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).ICencoder1==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).sigmcBK, [1 0 0 0], 'rows');
        case 4
            neuoi = ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).ICencoder1==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).sigmcBK, [1 0 0 0], 'rows');
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
            neuoi = ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).ICencoder2==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).sigmcBK, [0 0 0 1], 'rows');
        case 2
            neuoi = ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).ICencoder2==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).sigmcBK, [0 0 0 1], 'rows');
        case 3
            neuoi = ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).ICencoder1==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).sigmcBK, [1 0 0 0], 'rows');
        case 4
            neuoi = ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).ICencoder1==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).sigmcBK, [1 0 0 0], 'rows');
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

figure; 
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'All Trials', 'edgecolor', 'none')
subplot(2,2,1); imagesc(hcprefiori_wICenc)
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
title('White Inducers Black Background')
subplot(2,2,2); imagesc(hcprefiori_wICenc./sum(hcprefiori_wICenc,2))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
subplot(2,2,3); imagesc(hcprefiori_kICenc)
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
title('Black Inducers White Background')
subplot(2,2,4); imagesc(hcprefiori_kICenc./sum(hcprefiori_kICenc,2))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')

figure; 
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'All Trials: Significantly Tuned Neurons', 'edgecolor', 'none')
subplot(2,2,1); imagesc(hcsigprefiori_wICenc)
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
title('White Inducers Black Background')
subplot(2,2,2); imagesc(hcsigprefiori_wICenc./sum(hcsigprefiori_wICenc,2))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
subplot(2,2,3); imagesc(hcsigprefiori_kICenc)
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
title('Black Inducers White Background')
subplot(2,2,4); imagesc(hcsigprefiori_kICenc./sum(hcsigprefiori_kICenc,2))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')

%% bootstrapped confidence interval: 
whichprobe = 'all';
neuinarea = true(size(neuprobe_fixedgazeagg));

% whichprobe = 'C';
% iprobe = find(strcmp(probes, whichprobe));
% % neuinarea = strcmp(neuloc_fixedgazeagg{iprobe}, 'VISp2/3');
% % neuinarea = contains(neuloc_fixedgazeagg{iprobe}, 'VISp');
% neuinarea = neuctx_fixedgazeagg{iprobe}==1;
% % neuinarea = true(size(neuctx_fixedgazeagg{iprobe}));

opt678 = true;
if opt678
prefidiragg = oriparams_fixedgazeagg.(whichprobe).prefiori678;
prefioriagg = ori4params_fixedgazeagg.(whichprobe).prefiori4678;
sigori = ori4params_fixedgazeagg.(whichprobe).Pkw_ori4678<0.05;
% Rori = ( oriparams_fixedgazeagg.(whichprobe).Rori(:,1:4) + oriparams_fixedgazeagg.(whichprobe).Rori(:,5:8) )/2;
% [~,prefioriagg] = max(Rori, [], 2);
Rori = ori4params_fixedgazeagg.(whichprobe).Rori4678;
else
prefidiragg = oriparams_fixedgazeagg.(whichprobe).prefiori;
prefioriagg = ori4params_fixedgazeagg.(whichprobe).prefiori4;
sigori = ori4params_fixedgazeagg.(whichprobe).Pkw_ori4<0.05;
Rori = ori4params_fixedgazeagg.(whichprobe).Rori4;
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
        neuoi = ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).ICencoder2==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).sigmcBK, [0 0 0 1], 'rows');
    case 2
        neuoi = ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).ICencoder2==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).sigmcBK, [0 0 0 1], 'rows');
    case 3
        neuoi = ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).ICencoder1==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).sigmcBK, [1 0 0 0], 'rows');
    case 4
        neuoi = ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).ICencoder1==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).sigmcBK, [1 0 0 0], 'rows');
    case 5
        neuoi = ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).ICencoder2==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).sigmcBK, [0 0 0 1], 'rows');
    case 6
        neuoi = ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).ICencoder2==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).sigmcBK, [0 0 0 1], 'rows');
    case 7
        neuoi = ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).ICencoder1==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).sigmcBK, [1 0 0 0], 'rows');
    case 8
        neuoi = ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).ICencoder1==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).sigmcBK, [1 0 0 0], 'rows');
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

%% receptive field distribution
whichprobe = 'all';
neuinarea = true(size(neuprobe_fixedgazeagg));

whichprobe = 'C';
iprobe = find(strcmp(probes, whichprobe));
% neuinarea = strcmp(neuloc_fixedgazeagg{iprobe}, 'VISp2/3');
% neuinarea = contains(neuloc_fixedgazeagg{iprobe}, 'VISp');
neuinarea = neuctx_fixedgazeagg{iprobe}==1;
% neuinarea = true(size(neuctx_fixedgazeagg{iprobe}));

if strcmp(whichprobe, 'all')
    iprobe = [];
else
    iprobe = find( strcmp(whichprobe, probes) );
end

indcoldesc = 'White Pacmans';
sigBK = ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).PkwBK<0.05 | ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).PkwBK<0.05;
ICenc = ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).ICencoder==1 | ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).ICencoder==1;
RCenc = ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).RCencoder==1 | ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).RCencoder==1;

% indcoldesc = 'Black Pacmans';
% sigBK = ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).PkwBK<0.05 | ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).PkwBK<0.05;
% ICenc = ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).ICencoder==1 | ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).ICencoder==1;
% RCenc = ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).RCencoder==1 | ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).RCencoder==1;

RFdist_IC_acc = histcounts(RFCI_fixedgazeagg.(whichprobe).RFindclassic(ICenc), 0.5:1:9.5);
RFdist_sigBK_acc = histcounts(RFCI_fixedgazeagg.(whichprobe).RFindclassic(sigBK), 0.5:1:9.5);
RFdist_all_acc = histcounts(RFCI_fixedgazeagg.(whichprobe).RFindclassic, 0.5:1:9.5);

RFdist_IC_agg = zeros(Nsessions, 9);
RFdist_sigBK_agg = zeros(Nsessions, 9);
RFdist_all_agg = zeros(Nsessions, 9);
for ises = 1:Nsessions
    if strcmp(whichprobe, 'all')
    tempsesneu = cat(1,sesneu_fixedgazeagg{:})==ises;
    else
    tempsesneu = sesneu_fixedgazeagg{iprobe}==ises;
    end
RFdist_IC_agg(ises,:) = histcounts(RFCI_fixedgazeagg.(whichprobe).RFindclassic(ICenc(tempsesneu)), 0.5:1:9.5);
RFdist_sigBK_agg(ises,:) = histcounts(RFCI_fixedgazeagg.(whichprobe).RFindclassic(sigBK(tempsesneu)), 0.5:1:9.5);
RFdist_all_agg(ises,:) = histcounts(RFCI_fixedgazeagg.(whichprobe).RFindclassic(tempsesneu), 0.5:1:9.5);
end

figure; hold all
plot(1:9, RFdist_IC_agg./RFdist_sigBK_agg)
plot(1:9, RFdist_IC_acc./RFdist_sigBK_acc, 'k.-', 'LineWidth', 2)

figure; hold all
neudescs = cell(1,0);
for ii = 1:3
    switch ii
        case 1
            tempneu = true(size(ICenc));
            tempcol = [0 0 0];
            tempneudesc = 'All';
        case 2
            tempneu = RCenc;
            tempcol = [1 0.5 0];
            tempneudesc = 'RC-encoder';
        case 3
            tempneu = ICenc;
            tempcol = [0 0.7 0];
            tempneudesc = 'IC-encoder';
    end
    tempneu = tempneu & neuinarea;
    tempvec = RFCI_fixedgazeagg.(whichprobe).RFindclassic(tempneu);
    %neudescs{ii} = sprintf('%s N=%d/%d',tempneudesc, nnz(~isnan(tempvec)), nnz(tempneu));
    neudescs{ii} = sprintf('%s N=%d',tempneudesc, nnz(tempneu));
histogram(tempvec, 0.5:1:Nrfs+.5, 'FaceColor', tempcol, 'FaceAlpha', 0.5, 'normalization', 'probability')
end
legend(neudescs)
xlabel('CRF Position')
ylabel('Probability')
title(['Fixed Gaze Trials, ' indcoldesc])

figure; hold all
neudescs = cell(1,0);
for ii = 1:3
    switch ii
        case 1
            tempneu = true(size(ICenc));
            tempcol = [0 0 0];
            tempneudesc = 'All';
        case 2
            tempneu = RCenc;
            tempcol = [1 0.5 0];
            tempneudesc = 'RC-encoder';
        case 3
            tempneu = ICenc;
            tempcol = [0 0.7 0];
            tempneudesc = 'IC-encoder';
    end
    tempneu = tempneu & neuinarea;
    %neudescs{ii} = tempneudesc;
    neudescs{ii} = sprintf('%s N=%d',tempneudesc, nnz(tempneu));
histogram(RFCI_fixedgazeagg.(whichprobe).RFindinverse(tempneu), 0.5:1:Nrfs+.5, 'FaceColor', tempcol, 'FaceAlpha', 0.5, 'normalization', 'probability')
end
legend(neudescs)
xlabel('IRF Position')
ylabel('Probability')
title(['Fixed Gaze Trials, ' indcoldesc])

lw = 3;
figure; hold all
neudescs = cell(1,0);
for ii = 1:3
    switch ii
        case 1
            tempneu = true(size(ICenc));
            tempcol = [0 0 0];
            tempneudesc = 'All';
        case 2
            tempneu = RCenc;
            tempcol = [1 0.5 0];
            tempneudesc = 'RC-encoder';
        case 3
            tempneu = ICenc;
            tempcol = [0 0.7 0];
            tempneudesc = 'IC-encoder';
    end
    tempneu = tempneu & neuinarea;
    %neudescs{ii} = tempneudesc;
    neudescs{ii} = sprintf('%s N=%d',tempneudesc, nnz(tempneu));
errorbar(1:Nrfs, nanmean(RFCI_fixedgazeagg.(whichprobe).Rrfclassic(tempneu,:),1), nanstd(RFCI_fixedgazeagg.(whichprobe).Rrfclassic(tempneu,:),0,1)/sqrt(nnz(tempneu)), 'o-', 'Color', tempcol, 'MarkerFaceColor', tempcol, 'LineWidth', lw)
end
legend(neudescs)
xlabel('CRF Position')
ylabel('Rate (Hz)')

figure; hold all
neudescs = cell(1,0);
for ii = 1:3
    switch ii
        case 1
            tempneu = true(size(ICenc));
            tempcol = [0 0 0];
            tempneudesc = 'All';
        case 2
            tempneu = RCenc;
            tempcol = [1 0.5 0];
            tempneudesc = 'RC-encoder';
        case 3
            tempneu = ICenc;
            tempcol = [0 0.7 0];
            tempneudesc = 'IC-encoder';
    end
    tempneu = tempneu & neuinarea;
    %neudescs{ii} = tempneudesc;
    neudescs{ii} = sprintf('%s N=%d',tempneudesc, nnz(tempneu));
errorbar(1:Nrfs, nanmean(RFCI_fixedgazeagg.(whichprobe).Rrfinverse(tempneu,:),1), nanstd(RFCI_fixedgazeagg.(whichprobe).Rrfinverse(tempneu,:),0,1)/sqrt(nnz(tempneu)), 'o-', 'Color', tempcol, 'MarkerFaceColor', tempcol, 'LineWidth', lw)
end
legend(neudescs)
xlabel('IRF Position')
ylabel('Rate (Hz)')

tempneu = neuinarea;
figure; histogram2(RFCI_fixedgazeagg.(whichprobe).RFindclassic(tempneu), RFCI_fixedgazeagg.(whichprobe).RFindinverse(tempneu), 0.5:1:Nrfs+.5, 0.5:1:Nrfs+.5, 'displaystyle', 'tile')

% tempneu = neuinarea;
% tempneu = neuinarea & RFCI_fixedgazeagg.(whichprobe).Pkw_rfclassic<0.05;
tempneu = neuinarea & ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).ICencoder==1 | ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).ICencoder==1;
% tempneu = neuinarea & ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).ICencoder==1 | ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).ICencoder==1;
% tempneu = neuinarea & ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).ICencoder==1 | ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).ICencoder==1 ...
%     | ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).ICencoder==1 | ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).ICencoder==1;
hc = histcounts2(RFCI_fixedgazeagg.(whichprobe).RFindinverse(tempneu), RFCI_fixedgazeagg.(whichprobe).RFindclassic(tempneu), 0.5:1:Nrfs+.5, 0.5:1:Nrfs+.5);
figure; imagesc(hc./sum(hc,1)); colorbar
xlabel('CRF Position')
ylabel('IRF Position')
title(sprintf('Proportion out of each CRF\nDiagonal Sum %.2f', sum(diag(hc./sum(hc,2)))))

hc = histcounts2(RFCI_fixedgazeagg.(whichprobe).RFindclassic(tempneu), RFCI_fixedgazeagg.(whichprobe).RFindinverse(tempneu), 0.5:1:Nrfs+.5, 0.5:1:Nrfs+.5);
figure; imagesc(hc./sum(hc,2)); colorbar
xlabel('IRF Position')
ylabel('CRF Position')

%% size tuning for cells with receptive field in the center
whichprobe = 'all';
neuinarea = true(size(neuprobe_fixedgazeagg));

whichprobe = 'C';
iprobe = find(strcmp(probes, whichprobe));
% neuinarea = strcmp(neuloc_fixedgazeagg{iprobe}, 'VISp2/3');
% neuinarea = contains(neuloc_fixedgazeagg{iprobe}, 'VISp');
neuinarea = neuctx_fixedgazeagg{iprobe}==1;
% neuinarea = true(size(neuctx_fixedgazeagg{iprobe}));

szvec = [0, 4, 8, 16, 32, 64];

neuctrCRF = RFCI_fixedgazeagg.(whichprobe).RFindclassic==1 & RFCI_fixedgazeagg.(whichprobe).Pkw_rfclassic<0.05;
neuctrIRF = RFCI_fixedgazeagg.(whichprobe).RFindinverse==1 & RFCI_fixedgazeagg.(whichprobe).Pkw_rfinverse<0.05;

indcoldesc = 'White Pacmans';
sigBK = ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).PkwBK<0.05 | ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).PkwBK<0.05;
ICenc = ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).ICencoder==1 | ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).ICencoder==1;
RCenc = ICsig_fixedgazeagg.ICwcfg1_presentations.(whichprobe).RCencoder==1 | ICsig_fixedgazeagg.ICwcfg0_presentations.(whichprobe).RCencoder==1;

% indcoldesc = 'Black Pacmans';
% sigBK = ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).PkwBK<0.05 | ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).PkwBK<0.05;
% ICenc = ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).ICencoder==1 | ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).ICencoder==1;
% RCenc = ICsig_fixedgazeagg.ICkcfg1_presentations.(whichprobe).RCencoder==1 | ICsig_fixedgazeagg.ICkcfg0_presentations.(whichprobe).RCencoder==1;

lw = 3;
figure; hold all
neudescs = cell(1,0);
for ii = 1:3
    switch ii
        case 1
            tempneu = neuctrCRF;
            tempcol = [0 0 0];
            tempneudesc = 'ctrCRF';
        case 2
            tempneu = neuctrCRF & RCenc;
            tempcol = [1 0.5 0];
            tempneudesc = 'ctrCRF & RC-encoder';
        case 3
            tempneu = neuctrCRF & ICenc;
            tempcol = [0 0.7 0];
            tempneudesc = 'ctrCRF & IC-encoder';
    end
    tempneu = tempneu & neuinarea;
    neudescs{ii} = tempneudesc;
    neudescs{ii} = sprintf('%s N=%d',tempneudesc, nnz(tempneu));
errorbar(szvec, nanmean(sizeCI_fixedgazeagg.(whichprobe).Rsizeclassic(tempneu,:),1), nanstd(sizeCI_fixedgazeagg.(whichprobe).Rsizeclassic(tempneu,:),0,1)/sqrt(nnz(tempneu)), 'o-', 'Color', tempcol, 'MarkerFaceColor', tempcol, 'LineWidth', lw)
end
legend(neudescs, 'location', 'best')
set(gca, 'XTick', szvec)
xlabel('CG Size (Visual Degrees)')
ylabel('Rate (Hz)')
title(['Fixed Gaze Trials, ' indcoldesc])

figure; hold all
neudescs = cell(1,0);
for ii = 1:3
    switch ii
        case 1
            tempneu = neuctrCRF;
            tempcol = [0.5 0.5 0.5];
            tempneudesc = 'ctrCRF';
        case 2
            tempneu = neuctrCRF & RCenc;
            tempcol = [1 0.5 0];
            tempneudesc = 'ctrCRF & RC-encoder';
        case 3
            tempneu = neuctrCRF & ICenc;
            tempcol = [0 0.7 0];
            tempneudesc = 'ctrCRF & IC-encoder';
    end
    tempneu = tempneu & neuinarea;
    neudescs{ii} = tempneudesc;
    neudescs{ii} = sprintf('%s N=%d',tempneudesc, nnz(tempneu));
histogram(sizeCI_fixedgazeagg.(whichprobe).sizeindclassic(tempneu), 0.5:length(szvec)+0.5, 'FaceColor', tempcol, 'FaceAlpha', 0.5, 'normalization', 'probability')
end
legend(neudescs, 'location', 'best')
set(gca, 'XTick', 1:length(szvec), 'XTickLabel', szvec)
xlabel('Preferred CG Size (Visual Degrees)')
ylabel('Probability')
title(['Fixed Gaze Trials, ' indcoldesc])

lw = 3;
figure; hold all
neudescs = cell(1,0);
for ii = 1:3
    switch ii
        case 1
            tempneu = true(size(ICenc));
            tempcol = [0 0 0];
            tempneudesc = 'All';
        case 2
            tempneu = RCenc;
            tempcol = [1 0.5 0];
            tempneudesc = 'RC-encoder';
        case 3
            tempneu = ICenc;
            tempcol = [0 0.7 0];
            tempneudesc = 'IC-encoder';
    end
    tempneu = tempneu & neuinarea;
    neudescs{ii} = tempneudesc;
    neudescs{ii} = sprintf('%s N=%d',tempneudesc, nnz(tempneu));
errorbar(szvec, nanmean(sizeCI_fixedgazeagg.(whichprobe).Rsizeinverse(tempneu,:),1), nanstd(sizeCI_fixedgazeagg.(whichprobe).Rsizeinverse(tempneu,:),0,1)/sqrt(nnz(tempneu)), 'o-', 'Color', tempcol, 'MarkerFaceColor', tempcol, 'LineWidth', lw)
end
legend(neudescs)
set(gca, 'XTick', szvec)
xlabel('IG Size (Visual Degrees)')
ylabel('Rate (Hz)')
title(['Fixed Gaze Trials, ' indcoldesc])

lw = 3;
figure; hold all
neudescs = cell(1,0);
for ii = 1:3
    switch ii
        case 1
            tempneu = neuctrCRF;
            tempcol = [0 0 0];
            tempneudesc = 'ctrCRF';
        case 2
            tempneu = neuctrCRF & RCenc;
            tempcol = [1 0.5 0];
            tempneudesc = 'ctrCRF & RC-encoder';
        case 3
            tempneu = neuctrCRF & ICenc;
            tempcol = [0 0.7 0];
            tempneudesc = 'ctrCRF & IC-encoder';
    end
    tempneu = tempneu & neuinarea;
    neudescs{ii} = tempneudesc;
    neudescs{ii} = sprintf('%s N=%d',tempneudesc, nnz(tempneu));
errorbar(szvec, mean(sizeCI_fixedgazeagg.(whichprobe).Rsizeinverse(tempneu,:),1), std(sizeCI_fixedgazeagg.(whichprobe).Rsizeinverse(tempneu,:),0,1)/sqrt(nnz(tempneu)), 'o-', 'Color', tempcol, 'MarkerFaceColor', tempcol, 'LineWidth', lw)
end
legend(neudescs)
set(gca, 'XTick', szvec)
xlabel('IG Size (Visual Degrees)')
ylabel('Rate (Hz)')
title(['Fixed Gaze Trials, ' indcoldesc])

lw = 3;
figure; hold all
neudescs = cell(1,0);
for ii = 1:3
    switch ii
        case 1
            tempneu = neuctrIRF;
            tempcol = [0 0 0];
            tempneudesc = 'ctrIRF';
        case 2
            tempneu = neuctrIRF & RCenc;
            tempcol = [1 0.5 0];
            tempneudesc = 'ctrIRF & RC-encoder';
        case 3
            tempneu = neuctrIRF & ICenc;
            tempcol = [0 0.7 0];
            tempneudesc = 'ctrIRF & IC-encoder';
    end
    tempneu = tempneu & neuinarea;
    neudescs{ii} = tempneudesc;
    neudescs{ii} = sprintf('%s N=%d',tempneudesc, nnz(tempneu));
errorbar(szvec, mean(sizeCI_fixedgazeagg.(whichprobe).Rsizeinverse(tempneu,:),1), std(sizeCI_fixedgazeagg.(whichprobe).Rsizeinverse(tempneu,:),0,1)/sqrt(nnz(tempneu)), 'o-', 'Color', tempcol, 'MarkerFaceColor', tempcol, 'LineWidth', lw)
end
legend(neudescs)
set(gca, 'XTick', szvec)
xlabel('IG Size (Visual Degrees)')
ylabel('Rate (Hz)')
title(['Fixed Gaze Trials, ' indcoldesc])
