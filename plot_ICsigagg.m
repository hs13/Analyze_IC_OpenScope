%% IC-encoder depth (peak channel index)
iprobe = find( strcmp('C', probes) );

ICencagg = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder==1;% | ICsigagg.ICwcfg0_presentations(iprobe).ICencoder==1;
RCencagg = ICsigagg.ICwcfg1_presentations(iprobe).RCencoder==1;% | ICsigagg.ICwcfg0_presentations(iprobe).RCencoder==1;

figure; hold all
histogram(neupeakchagg{iprobe}, 'binwidth', 10, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.5, 'normalization', 'probability')
histogram(neupeakchagg{iprobe}(ICencagg), 'binwidth', 10, 'FaceColor', [0 0.7 0], 'FaceAlpha', 0.5, 'normalization', 'probability')
xlabel('higher numbes more superficial')

%% V1: pref-ori distribution of IC/RC encoders
% RUN AFTER RUNNING FIRST TWO SECTIONS OF aggregate_psth
iprobe = find( strcmp('C', probes) );
opt678 = true;
if opt678
prefidiragg = oriparamsagg(iprobe).prefiori678;
prefioriagg = ori4paramsagg(iprobe).prefiori4678;
sigori = ori4paramsagg(iprobe).Pkw_ori4678<0.05;
% Rori = ( oriparamsagg(iprobe).Rori678(:,1:4) + oriparamsagg(iprobe).Rori678(:,5:8) )/2;
% [~,prefioriagg] = max(Rori, [], 2);
Rori = ori4paramsagg(iprobe).Rori4678(:,1:4);
else
prefidiragg = oriparamsagg(iprobe).prefiori;
prefioriagg = ori4paramsagg(iprobe).prefiori4;
sigori = ori4paramsagg(iprobe).Pkw_ori4<0.05;
Rori = ori4paramsagg(iprobe).Rori4(:,1:4);
end

% % sanity check
% figure; histogram2(ori4paramsagg(iprobe).prefiori4(sigori), ori4paramsagg(iprobe).prefiori4678(sigori), 'displaystyle', 'tile')

% report sigmcBK of ICencoder1 and ICencoder2 
for ienc = 1:4
    switch ienc
        case 1
            neuoi = ICsigagg.ICwcfg0_presentations(iprobe).ICencoder2==1;
            tempsigmcBK = ICsigagg.ICwcfg0_presentations(iprobe).sigmcBK;
        case 2
            neuoi = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder2==1;
            tempsigmcBK = ICsigagg.ICwcfg1_presentations(iprobe).sigmcBK;
        case 3
            neuoi = ICsigagg.ICwcfg0_presentations(iprobe).ICencoder1==1;
            tempsigmcBK = ICsigagg.ICwcfg0_presentations(iprobe).sigmcBK;
        case 4
            neuoi = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder1==1;
            tempsigmcBK = ICsigagg.ICwcfg1_presentations(iprobe).sigmcBK;
    end
    [C,IA,IC] = unique(tempsigmcBK(neuoi,:), 'rows');
    rowcnt = zeros(size(C,1),1);
    for ii = 1:size(C,1)
        rowcnt(ii) = nnz(IC==ii);
    end
    disp([C rowcnt])
end

% % tuning curve across directions for each ICencoder group
% figure; hold all
% for ienc = 1:4
%     switch ienc
%         case 1
%             neuoi = ICsigagg.ICwcfg0_presentations(iprobe).ICencoder2==1;
%             neuoi = ismember(ICsigagg.ICwcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
%         case 2
%             neuoi = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder2==1;
%             neuoi = ismember(ICsigagg.ICwcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
%         case 3
%             neuoi = ICsigagg.ICwcfg0_presentations(iprobe).ICencoder1==1;
%             neuoi = ismember(ICsigagg.ICwcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
%         case 4
%             neuoi = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder1==1;
%             neuoi = ismember(ICsigagg.ICwcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
%     end
%     errorbar(1:Ndirs, nanmean(oriparamsagg(iprobe).Rori678(neuoi,:),1), ...
%         nanstd(oriparamsagg(iprobe).Rori678(neuoi,:),0,1)/sqrt(nnz(neuoi)), 'o-')
% end
% 
% % tuning curve across orientations for each ICencoder group
% figure; hold all
% for ienc = 1:4
%     switch ienc
%         case 1
%             neuoi = ICsigagg.ICwcfg0_presentations(iprobe).ICencoder2==1;
%             neuoi = ismember(ICsigagg.ICwcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
%         case 2
%             neuoi = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder2==1;
%             neuoi = ismember(ICsigagg.ICwcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
%         case 3
%             neuoi = ICsigagg.ICwcfg0_presentations(iprobe).ICencoder1==1;
%             neuoi = ismember(ICsigagg.ICwcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
%         case 4
%             neuoi = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder1==1;
%             neuoi = ismember(ICsigagg.ICwcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
%     end
% errorbar(1:Noris, nanmean(Rori(neuoi,:),1), nanstd(Rori(neuoi,:),0,1)/sqrt(nnz(neuoi)), 'o-')
% end

% % Rori friemdan test across orientations for each ICencoder group
% for ienc = 1:4
%     switch ienc
%         case 1
%             neuoi = ICsigagg.ICwcfg0_presentations(iprobe).ICencoder2==1;
%             neuoi = ismember(ICsigagg.ICwcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
%         case 2
%             neuoi = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder2==1;
%             neuoi = ismember(ICsigagg.ICwcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
%         case 3
%             neuoi = ICsigagg.ICwcfg0_presentations(iprobe).ICencoder1==1;
%             neuoi = ismember(ICsigagg.ICwcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
%         case 4
%             neuoi = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder1==1;
%             neuoi = ismember(ICsigagg.ICwcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
%     end
%     [p,tbl,stats] = friedman(Rori(neuoi,:));
%     [c,m,h] = multcompare(stats);
%     title(sprintf('%d p=%.4f', ienc, p))
% end

hcprefidir_wICenc = zeros(4, Ndirs);
hcprefiori_wICenc = zeros(4, Noris);
hcsigprefidir_wICenc = zeros(4, Ndirs);
hcsigprefiori_wICenc = zeros(4, Noris);
for ienc = 1:4
    switch ienc
        case 1
            neuoi = ICsigagg.ICwcfg0_presentations(iprobe).ICencoder2==1;
            neuoi = ismember(ICsigagg.ICwcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
        case 2
            neuoi = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder2==1;
            neuoi = ismember(ICsigagg.ICwcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
        case 3
            neuoi = ICsigagg.ICwcfg0_presentations(iprobe).ICencoder1==1;
            neuoi = ismember(ICsigagg.ICwcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
        case 4
            neuoi = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder1==1;
            neuoi = ismember(ICsigagg.ICwcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    end
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
            neuoi = ICsigagg.ICkcfg0_presentations(iprobe).ICencoder2==1;
            neuoi = ismember(ICsigagg.ICkcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
        case 2
            neuoi = ICsigagg.ICkcfg1_presentations(iprobe).ICencoder2==1;
            neuoi = ismember(ICsigagg.ICkcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
        case 3
            neuoi = ICsigagg.ICkcfg0_presentations(iprobe).ICencoder1==1;
            neuoi = ismember(ICsigagg.ICkcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
        case 4
            neuoi = ICsigagg.ICkcfg1_presentations(iprobe).ICencoder1==1;
            neuoi = ismember(ICsigagg.ICkcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    end
    hcprefidir_kICenc(ienc,:) = histcounts(prefidiragg(neuoi), 0.5:Ndirs+0.5);
    hcprefiori_kICenc(ienc,:) = histcounts(prefioriagg(neuoi), 0.5:Noris+0.5);
    hcsigprefidir_kICenc(ienc,:) = histcounts(prefidiragg(sigori & neuoi), 0.5:Ndirs+0.5);
    hcsigprefiori_kICenc(ienc,:) = histcounts(prefioriagg(sigori & neuoi), 0.5:Noris+0.5);
end

hcprefidir = histcounts(prefidiragg, 0.5:Ndirs+0.5);
hcprefiori = histcounts(prefioriagg, 0.5:Noris+0.5);
hcsigprefidir = histcounts(prefidiragg(sigori), 0.5:Ndirs+0.5);
hcsigprefiori = histcounts(prefioriagg(sigori), 0.5:Noris+0.5);

fs=10;
ytl = 0:45:180-1; ylab = 'IC Orientation';
figure; 
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'All Trials', 'edgecolor', 'none')
subplot(2,2,1); imagesc(hcprefidir_wICenc./sum(hcprefidir_wICenc,2))
colorbar; colormap jet
set(gca, 'FontSize', fs, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Direction')
title('White Inducers Black Background')
subplot(2,2,2); imagesc(hcprefiori_wICenc./sum(hcprefiori_wICenc,2))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
subplot(2,2,3); imagesc(hcprefidir_kICenc./sum(hcprefidir_kICenc,2))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Direction')
title('Black Inducers White Background')
subplot(2,2,4); imagesc(hcprefiori_kICenc./sum(hcprefiori_kICenc,2))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')

figure; 
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'All Trials: Significantly Tuned Neurons', 'edgecolor', 'none')
subplot(2,2,1); imagesc(hcsigprefidir_wICenc./sum(hcsigprefidir_wICenc,2))
colorbar; colormap jet
xlabel('Preferred Direction')
title('White Inducers Black Background')
subplot(2,2,2); imagesc(hcsigprefiori_wICenc./sum(hcsigprefiori_wICenc,2))
colorbar; colormap jet
xlabel('Preferred Orientation')
subplot(2,2,3); imagesc(hcsigprefidir_kICenc./sum(hcsigprefidir_kICenc,2))
colorbar; colormap jet
xlabel('Preferred Direction')
title('Black Inducers White Background')
subplot(2,2,4); imagesc(hcsigprefiori_kICenc./sum(hcsigprefiori_kICenc,2))
colorbar; colormap jet
xlabel('Preferred Orientation')

%{
figure; 
subplot(2,2,1); imagesc(hcprefidir_wICenc./hcprefidir)
colorbar; colormap jet
xlabel('Preferred Direction')
title('White Inducers Black Background')
subplot(2,2,2); imagesc(hcprefiori_wICenc./hcprefiori)
colorbar; colormap jet
xlabel('Preferred Orientation')
subplot(2,2,3); imagesc(hcprefidir_kICenc./hcprefidir)
colorbar; colormap jet
xlabel('Preferred Direction')
title('Black Inducers White Background')
subplot(2,2,4); imagesc(hcprefiori_kICenc./hcprefiori)
colorbar; colormap jet
xlabel('Preferred Orientation')

figure; 
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'Significantly Tuned Neurons', 'edgecolor', 'none')
subplot(2,2,1); imagesc(hcsigprefidir_wICenc./hcsigprefidir)
colorbar; colormap jet
xlabel('Preferred Direction')
title('White Inducers Black Background')
subplot(2,2,2); imagesc(hcsigprefiori_wICenc./hcsigprefiori)
colorbar; colormap jet
xlabel('Preferred Orientation')
subplot(2,2,3); imagesc(hcsigprefidir_kICenc./hcsigprefidir)
colorbar; colormap jet
xlabel('Preferred Direction')
title('Black Inducers White Background')
subplot(2,2,4); imagesc(hcsigprefiori_kICenc./hcsigprefiori)
colorbar; colormap jet
xlabel('Preferred Orientation')


figure; 
subplot(2,2,1); imagesc(hcprefidir_wICenc)
colorbar; colormap jet
xlabel('Preferred Direction')
title('White Inducers Black Background')
subplot(2,2,2); imagesc(hcprefiori_wICenc)
colorbar; colormap jet
xlabel('Preferred Orientation')
subplot(2,2,3); imagesc(hcprefidir_kICenc)
colorbar; colormap jet
xlabel('Preferred Direction')
title('Black Inducers White Background')
subplot(2,2,4); imagesc(hcprefiori_kICenc)
colorbar; colormap jet
xlabel('Preferred Orientation')

figure; 
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'Significantly Tuned Neurons', 'edgecolor', 'none')
subplot(2,2,1); imagesc(hcsigprefidir_wICenc)
colorbar; colormap jet
xlabel('Preferred Direction')
title('White Inducers Black Background')
subplot(2,2,2); imagesc(hcsigprefiori_wICenc)
colorbar; colormap jet
xlabel('Preferred Orientation')
subplot(2,2,3); imagesc(hcsigprefidir_kICenc)
colorbar; colormap jet
xlabel('Preferred Direction')
title('Black Inducers White Background')
subplot(2,2,4); imagesc(hcsigprefiori_kICenc)
colorbar; colormap jet
xlabel('Preferred Orientation')
%}

%% all areas: pref-ori distribution of IC/RC encoders
prefidiragg = cat(1,oriparamsagg.prefiori);
prefioriagg = cat(1,ori4paramsagg.prefiori4);
sigori = cat(1,ori4paramsagg.Pkw_ori4)<0.05;
Rori = cat(1,ori4paramsagg.Rori4);

hcprefidir_wICenc = zeros(4, Ndirs);
hcprefiori_wICenc = zeros(4, Noris);
hcsigprefidir_wICenc = zeros(4, Ndirs);
hcsigprefiori_wICenc = zeros(4, Noris);
for ienc = 1:4
    switch ienc
        case 1
            neuoi = cat(1,ICsigagg.ICwcfg0_presentations.ICencoder2)==1;
            neuoi = ismember(cat(1,ICsigagg.ICwcfg0_presentations.sigmcBK), [0 0 0 1], 'rows');
        case 2
            neuoi = cat(1,ICsigagg.ICwcfg1_presentations.ICencoder2)==1;
            neuoi = ismember(cat(1,ICsigagg.ICwcfg1_presentations.sigmcBK), [0 0 0 1], 'rows');
        case 3
            neuoi = cat(1,ICsigagg.ICwcfg0_presentations.ICencoder1)==1;
            neuoi = ismember(cat(1,ICsigagg.ICwcfg0_presentations.sigmcBK), [0 0 0 1], 'rows');
        case 4
            neuoi = cat(1,ICsigagg.ICwcfg1_presentations.ICencoder1)==1;
            neuoi = ismember(cat(1,ICsigagg.ICwcfg1_presentations.sigmcBK), [0 0 0 1], 'rows');
    end
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
            neuoi = cat(1,ICsigagg.ICkcfg0_presentations.ICencoder2)==1;
            neuoi = ismember(cat(1,ICsigagg.ICkcfg0_presentations.sigmcBK), [0 0 0 1], 'rows');
        case 2
            neuoi = cat(1,ICsigagg.ICkcfg1_presentations.ICencoder2)==1;
            neuoi = ismember(cat(1,ICsigagg.ICkcfg1_presentations.sigmcBK), [0 0 0 1], 'rows');
        case 3
            neuoi = cat(1,ICsigagg.ICkcfg0_presentations.ICencoder1)==1;
            neuoi = ismember(cat(1,ICsigagg.ICkcfg0_presentations.sigmcBK), [0 0 0 1], 'rows');
        case 4
            neuoi = cat(1,ICsigagg.ICkcfg1_presentations.ICencoder1)==1;
            neuoi = ismember(cat(1,ICsigagg.ICkcfg1_presentations.sigmcBK), [0 0 0 1], 'rows');
    end
    hcprefidir_kICenc(ienc,:) = histcounts(prefidiragg(neuoi), 0.5:Ndirs+0.5);
    hcprefiori_kICenc(ienc,:) = histcounts(prefioriagg(neuoi), 0.5:Noris+0.5);
    hcsigprefidir_kICenc(ienc,:) = histcounts(prefidiragg(sigori & neuoi), 0.5:Ndirs+0.5);
    hcsigprefiori_kICenc(ienc,:) = histcounts(prefioriagg(sigori & neuoi), 0.5:Noris+0.5);
end

hcprefidir = histcounts(prefidiragg, 0.5:Ndirs+0.5);
hcprefiori = histcounts(prefioriagg, 0.5:Noris+0.5);
hcsigprefidir = histcounts(prefidiragg(sigori), 0.5:Ndirs+0.5);
hcsigprefiori = histcounts(prefioriagg(sigori), 0.5:Noris+0.5);

fs=10;
ytl = 0:45:180-1; ylab = 'IC Orientation';
figure; 
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'All Trials', 'edgecolor', 'none')
subplot(2,2,1); imagesc(hcprefidir_wICenc./sum(hcprefidir_wICenc,2))
colorbar; colormap jet
set(gca, 'FontSize', fs, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Direction')
title('White Inducers Black Background')
subplot(2,2,2); imagesc(hcprefiori_wICenc./sum(hcprefiori_wICenc,2))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
subplot(2,2,3); imagesc(hcprefidir_kICenc./sum(hcprefidir_kICenc,2))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Direction')
title('Black Inducers White Background')
subplot(2,2,4); imagesc(hcprefiori_kICenc./sum(hcprefiori_kICenc,2))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')

figure; 
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'All Trials: Significantly Tuned Neurons', 'edgecolor', 'none')
subplot(2,2,1); imagesc(hcsigprefidir_wICenc./sum(hcsigprefidir_wICenc,2))
colorbar; colormap jet
xlabel('Preferred Direction')
title('White Inducers Black Background')
subplot(2,2,2); imagesc(hcsigprefiori_wICenc./sum(hcsigprefiori_wICenc,2))
colorbar; colormap jet
xlabel('Preferred Orientation')
subplot(2,2,3); imagesc(hcsigprefidir_kICenc./sum(hcsigprefidir_kICenc,2))
colorbar; colormap jet
xlabel('Preferred Direction')
title('Black Inducers White Background')
subplot(2,2,4); imagesc(hcsigprefiori_kICenc./sum(hcsigprefiori_kICenc,2))
colorbar; colormap jet
xlabel('Preferred Orientation')

%% shuffle CI
iprobe = find( strcmp('C', probes) );
opt678 = true;
if opt678
prefidiragg = oriparamsagg(iprobe).prefiori678;
prefioriagg = ori4paramsagg(iprobe).prefiori4678;
sigori = ori4paramsagg(iprobe).Pkw_ori4678<0.05;
% Rori = ( oriparamsagg(iprobe).Rori678(:,1:4) + oriparamsagg(iprobe).Rori678(:,5:8) )/2;
% [~,prefioriagg] = max(Rori, [], 2);
Rori = ori4paramsagg(iprobe).Rori4678(:,1:4);
else
prefidiragg = oriparamsagg(iprobe).prefiori;
prefioriagg = ori4paramsagg(iprobe).prefiori4;
sigori = ori4paramsagg(iprobe).Pkw_ori4<0.05;
Rori = ori4paramsagg(iprobe).Rori4(:,1:4);
end

whichcol = 'w';
switch whichcol
    case 'w'
ICpreforivec = [];
ICgroupvec = [];
for ienc = 1:4
    switch ienc
        case 1
            neuoi = ICsigagg.ICwcfg0_presentations(iprobe).ICencoder2==1;
            neuoi = ismember(ICsigagg.ICwcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
        case 2
            neuoi = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder2==1;
            neuoi = ismember(ICsigagg.ICwcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
        case 3
            neuoi = ICsigagg.ICwcfg0_presentations(iprobe).ICencoder1==1;
            neuoi = ismember(ICsigagg.ICwcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
        case 4
            neuoi = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder1==1;
            neuoi = ismember(ICsigagg.ICwcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    end
ICpreforivec = cat(1, ICpreforivec, prefioriagg(neuoi));
ICgroupvec = cat(1, ICgroupvec, ienc*ones(nnz(neuoi),1));
end
    case 'k'
ICpreforivec = [];
ICgroupvec = [];
for ienc = 1:4
    switch ienc
        case 1
            neuoi = ICsigagg.ICkcfg0_presentations(iprobe).ICencoder2==1;
            neuoi = ismember(ICsigagg.ICkcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
        case 2
            neuoi = ICsigagg.ICkcfg1_presentations(iprobe).ICencoder2==1;
            neuoi = ismember(ICsigagg.ICkcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
        case 3
            neuoi = ICsigagg.ICkcfg0_presentations(iprobe).ICencoder1==1;
            neuoi = ismember(ICsigagg.ICkcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
        case 4
            neuoi = ICsigagg.ICkcfg1_presentations(iprobe).ICencoder1==1;
            neuoi = ismember(ICsigagg.ICkcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    end
ICpreforivec = cat(1, ICpreforivec, prefioriagg(neuoi));
ICgroupvec = cat(1, ICgroupvec, ienc*ones(nnz(neuoi),1));
end
end

% kruskalwallis of preferred orientation across neuron group
% p = kruskalwallis(ICpreforivec, ICgroupvec);

% figure; histogram2(ICpreforivec, ICgroupvec, 'displaystyle', 'tile')

hc = histcounts2(ICgroupvec, ICpreforivec, 0.5:1:4.5, 0.5:1:Noris+.5);
figure; imagesc(hc./sum(hc,2)); colorbar
xlabel('Grating Pref Ori')
ylabel('IC Pref Ori')

Nshuf = 1000;
hc_shuf = NaN(4, Noris, Nshuf);
for ishuf = 1:Nshuf
hc_shuf(:,:,ishuf) = histcounts2(ICgroupvec, ICpreforivec(randperm(length(ICpreforivec))), 0.5:1:4.5, 0.5:1:Noris+.5);
end
hc_shufmean = squeeze(mean(hc_shuf,3));
hc_shufnegerr = squeeze(mean(hc_shuf,3) - prctile(hc_shuf,92.5,3));
hc_shufposerr = squeeze(prctile(hc_shuf,97.5,3) - mean(hc_shuf,3));

figure
for ienc = 1:4
    subplot(2,2,ienc)
    hold all
    errorbar(1:Noris, hc_shufmean(ienc,:), hc_shufnegerr(ienc,:), hc_shufposerr(ienc,:), 'k.-')
    plot(1:Noris, hc(ienc,:), 'r-', 'LineWidth', 3)
    title(sprintf('ICenc%d', ienc))
end

%% friedman across sessions: proportion of IC-encoders in each orientation
iprobe = find( strcmp('C', probes) );
opt678 = true;
if opt678
prefidiragg = oriparamsagg(iprobe).prefiori678;
prefioriagg = ori4paramsagg(iprobe).prefiori4678;
sigori = ori4paramsagg(iprobe).Pkw_ori4678<0.05;
% Rori = ( oriparamsagg(iprobe).Rori(:,1:4) + oriparamsagg(iprobe).Rori(:,5:8) )/2;
% [~,prefioriagg] = max(Rori, [], 2);
Rori = ori4paramsagg(iprobe).Rori4678;
else
prefidiragg = oriparamsagg(iprobe).prefiori;
prefioriagg = ori4paramsagg(iprobe).prefiori4;
sigori = ori4paramsagg(iprobe).Pkw_ori4<0.05;
Rori = ori4paramsagg(iprobe).Rori4;
end

hpICenc = cell(8,1);
hpICenc_sigori = cell(8,1);
hpfried = zeros(8,1);
hpfried_sigori = zeros(8,1);
propICenc = cell(8,1);
propICenc_sigori = cell(8,1);
ppfried = zeros(8,1);
ppfried_sigori = zeros(8,1);
for ienc = 1:8
switch ienc
    case 1
        neuoi = ICsigagg.ICwcfg0_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsigagg.ICwcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 2
        neuoi = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsigagg.ICwcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 3
        neuoi = ICsigagg.ICwcfg0_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsigagg.ICwcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    case 4
        neuoi = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsigagg.ICwcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    case 5
        neuoi = ICsigagg.ICkcfg0_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsigagg.ICkcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 6
        neuoi = ICsigagg.ICkcfg1_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsigagg.ICkcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 7
        neuoi = ICsigagg.ICkcfg0_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsigagg.ICkcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    case 8
        neuoi = ICsigagg.ICkcfg1_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsigagg.ICkcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
end

hpICenc{ienc} = NaN(Nsessions,Noris);
hpICenc_sigori{ienc} = NaN(Nsessions,Noris);
propICenc{ienc} = NaN(Nsessions,Noris);
propICenc_sigori{ienc} = NaN(Nsessions,Noris);

for ises = 1:Nsessions
    tempsesneu = sesneuoiagg{iprobe}==ises;
    hpICenc{ienc}(ises,:) = histcounts(prefioriagg(tempsesneu & neuoi), 0.5:Noris+0.5, 'normalization', 'probability');
    hpICenc_sigori{ienc}(ises,:) = histcounts(prefioriagg(tempsesneu & sigori & neuoi), 0.5:Noris+0.5, 'normalization', 'probability');
    
    hcprefiori = histcounts(prefioriagg(tempsesneu), 0.5:Noris+0.5);
    hcprefiori_ICenc = histcounts(prefioriagg(tempsesneu & neuoi), 0.5:Noris+0.5);
    propICenc{ienc}(ises,:) = hcprefiori_ICenc./hcprefiori;
    
    hcprefiori = histcounts(prefioriagg(tempsesneu & sigori), 0.5:Noris+0.5);
    hcprefiori_ICenc = histcounts(prefioriagg(tempsesneu & sigori & neuoi), 0.5:Noris+0.5);
    propICenc_sigori{ienc}(ises,:) = hcprefiori_ICenc./hcprefiori;
end

tempmat = hpICenc{ienc};
tempmat = tempmat(all(~isnan(tempmat),2),:);
[p,tbl,stats] = friedman(tempmat, 1,'off');
hpfried(ienc) = p;

tempmat = hpICenc_sigori{ienc};
tempmat = tempmat(all(~isnan(tempmat),2),:);
[p,tbl,stats] = friedman(tempmat, 1,'off');
hpfried_sigori(ienc) = p;

tempmat = propICenc{ienc};
tempmat = tempmat(all(~isnan(tempmat),2),:);
[p,tbl,stats] = friedman(tempmat, 1,'off');
ppfried(ienc) = p;

tempmat = propICenc_sigori{ienc};
tempmat = tempmat(all(~isnan(tempmat),2),:);
[p,tbl,stats] = friedman(tempmat, 1,'off');
ppfried_sigori(ienc) = p;

end
disp([hpfried hpfried_sigori ppfried ppfried_sigori])

for ienc = 1:4
tempmat = hpICenc{ienc};
tempmat = tempmat(all(~isnan(tempmat),2),:);
[p,tbl,stats] = friedman(tempmat);
hpfried(ienc) = p;
[c,m,h] = multcompare(stats);
title(sprintf('%d p=%.4f', ienc, p))
end

%% jackknife proportion of IC-encoders in each orientation
hpICenc = cell(8,1);
hpICenc_sigori = cell(8,1);
hpfried = zeros(8,1);
hpfried_sigori = zeros(8,1);
propICenc = cell(8,1);
propICenc_sigori = cell(8,1);
Pfried = zeros(8,1);
Pfried_sigori = zeros(8,1);
for ienc = 1:4
switch ienc
    case 1
        neuoi = ICsigagg.ICwcfg0_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsigagg.ICwcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 2
        neuoi = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsigagg.ICwcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 3
        neuoi = ICsigagg.ICwcfg0_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsigagg.ICwcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    case 4
        neuoi = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsigagg.ICwcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    case 5
        neuoi = ICsigagg.ICkcfg0_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsigagg.ICkcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 6
        neuoi = ICsigagg.ICkcfg1_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsigagg.ICkcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 7
        neuoi = ICsigagg.ICkcfg0_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsigagg.ICkcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    case 8
        neuoi = ICsigagg.ICkcfg1_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsigagg.ICkcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
end
neuoind = find(neuoi);


hpICenc{ienc} = NaN(numel(neuoind),Noris);
hpICenc_sigori{ienc} = NaN(numel(neuoind),Noris);
propICenc{ienc} = NaN(numel(neuoind),Noris);
propICenc_sigori{ienc} = NaN(numel(neuoind),Noris);

for ci = 1:numel(neuoind)
    tempexclneu = true(size(neuoi));
    tempexclneu(neuoind(ci)) = false;
    hpICenc{ienc}(ci,:) = histcounts(prefioriagg(tempexclneu & neuoi), 0.5:Noris+0.5, 'normalization', 'probability');
    hpICenc_sigori{ienc}(ci,:) = histcounts(prefioriagg(tempexclneu & sigori & neuoi), 0.5:Noris+0.5, 'normalization', 'probability');
    
    hcprefiori = histcounts(prefioriagg(tempexclneu), 0.5:Noris+0.5);
    hcprefiori_ICenc = histcounts(prefioriagg(tempexclneu & neuoi), 0.5:Noris+0.5);
    propICenc{ienc}(ci,:) = hcprefiori_ICenc./hcprefiori;
    
    hcprefiori = histcounts(prefioriagg(tempexclneu & sigori), 0.5:Noris+0.5);
    hcprefiori_ICenc = histcounts(prefioriagg(tempexclneu & sigori & neuoi), 0.5:Noris+0.5);
    propICenc_sigori{ienc}(ci,:) = hcprefiori_ICenc./hcprefiori;
end

tempmat = hpICenc{ienc};
tempmat = tempmat(all(~isnan(tempmat),2),:);
[p,tbl,stats] = friedman(tempmat, 1,'off');
hpfried(ienc) = p;

tempmat = hpICenc_sigori{ienc};
tempmat = tempmat(all(~isnan(tempmat),2),:);
[p,tbl,stats] = friedman(tempmat, 1,'off');
hpfried_sigori(ienc) = p;

tempmat = propICenc{ienc};
tempmat = tempmat(all(~isnan(tempmat),2),:);
[p,tbl,stats] = friedman(tempmat, 1,'off');
ppfried(ienc) = p;

tempmat = propICenc_sigori{ienc};
tempmat = tempmat(all(~isnan(tempmat),2),:);
[p,tbl,stats] = friedman(tempmat, 1,'off');
ppfried_sigori(ienc) = p;
end
disp([hpfried hpfried_sigori ppfried ppfried_sigori])

for ienc = 1:4
tempmat = hpICenc{ienc};
tempmat = tempmat(all(~isnan(tempmat),2),:);
[p,tbl,stats] = friedman(tempmat);
hpfried(ienc) = p;
[c,m,h] = multcompare(stats);
title(sprintf('%d p=%.4f', ienc, p))
end

% figure; plot(hpICenc{4}')

%% bootstrapped confidence interval: 
iprobe = find( strcmp('C', probes) );
opt678 = true;
if opt678
prefidiragg = oriparamsagg(iprobe).prefiori678;
prefioriagg = ori4paramsagg(iprobe).prefiori4678;
sigori = ori4paramsagg(iprobe).Pkw_ori4678<0.05;
% Rori = ( oriparamsagg(iprobe).Rori(:,1:4) + oriparamsagg(iprobe).Rori(:,5:8) )/2;
% [~,prefioriagg] = max(Rori, [], 2);
Rori = ori4paramsagg(iprobe).Rori4678;
else
prefidiragg = oriparamsagg(iprobe).prefiori;
prefioriagg = ori4paramsagg(iprobe).prefiori4;
sigori = ori4paramsagg(iprobe).Pkw_ori4<0.05;
Rori = ori4paramsagg(iprobe).Rori4;
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
        neuoi = ICsigagg.ICwcfg0_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsigagg.ICwcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 2
        neuoi = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsigagg.ICwcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 3
        neuoi = ICsigagg.ICwcfg0_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsigagg.ICwcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    case 4
        neuoi = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsigagg.ICwcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    case 5
        neuoi = ICsigagg.ICkcfg0_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsigagg.ICkcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 6
        neuoi = ICsigagg.ICkcfg1_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsigagg.ICkcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 7
        neuoi = ICsigagg.ICkcfg0_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsigagg.ICkcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    case 8
        neuoi = ICsigagg.ICkcfg1_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsigagg.ICkcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
end
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

%% NICencoder, NRCencoder, Ninducerencoder, Ninducerresponsive, NsigBK
ICblockgroups = {'ICwcfg', 'ICwcfg0', 'ICwcfg1', 'ICkcfg', 'ICkcfg0', 'ICkcfg1'};
encodernames = {'ICencoder', 'RCencoder', 'inducerencoder', 'inducerresponsive'};
for e = 1:numel(encodernames)
    for b = 1:numel(ICblockgroups)
        Nencoder.(ICblockgroups{b}).(encodernames{e}) = zeros(Nsessions,numel(probes));
        for iprobe = 1:numel(probes)
            switch ICblockgroups{b}
                case 'ICwcfg'
                    tempneuvec = ICsigagg.ICwcfg1_presentations(iprobe).(encodernames{e}) | ICsigagg.ICwcfg0_presentations(iprobe).(encodernames{e});
                case 'ICkcfg'
                    tempneuvec = ICsigagg.ICkcfg1_presentations(iprobe).(encodernames{e}) | ICsigagg.ICkcfg0_presentations(iprobe).(encodernames{e});
                otherwise
                    tempneuvec = ICsigagg.([ICblockgroups{b} '_presentations'])(iprobe).(encodernames{e});
            end
            for ises = 1:Nsessions
                tempsesneu = sesneuoiagg{iprobe}==ises;
                Nencoder.(ICblockgroups{b}).(encodernames{e})(ises,iprobe) = nnz(tempneuvec(tempsesneu));
            end
        end
    end
end

% compare number of IC vs RC encoders
for b = 1:numel(ICblockgroups)
    disp(ICblockgroups{b})
    for a = 1:numel(probes)
        iprobe = probehierarchy(a);
        tempNIC = Nencoder.(ICblockgroups{b}).ICencoder(:,iprobe);
        tempNRC = Nencoder.(ICblockgroups{b}).RCencoder(:,iprobe);
        p = signrank(tempNIC, tempNRC);
        fprintf('%s IC(%d) vs RC(%d) p=%.4f (IC>RC %d/%d)\n', visareas{iprobe}, ...
            sum(tempNIC), sum(tempNRC), p, nnz(tempNIC>tempNRC), Nsessions)
    end
end

%% report percentage of IC/RC/ind-enc
% white inducer block
disp('ICwcfg0|ICwcfg0')
for a = 1:numel(probes)
    iprobe = probehierarchy(a);
    fprintf('%s (cumulative, out of all) IC-enc %.2f%% RC-enc %.2f%% ind-enc %.2f%%\n', visareas{iprobe}, ...
        100*mean(ICsigagg.ICwcfg0_presentations(iprobe).ICencoder==1 | ICsigagg.ICwcfg1_presentations(iprobe).ICencoder==1), ...
        100*mean(ICsigagg.ICwcfg0_presentations(iprobe).RCencoder==1 | ICsigagg.ICwcfg1_presentations(iprobe).RCencoder==1), ...
        100*mean(ICsigagg.ICwcfg0_presentations(iprobe).inducerencoder==1 | ICsigagg.ICwcfg1_presentations(iprobe).inducerencoder==1))
end

for a = 1:numel(probes)
    iprobe = probehierarchy(a);
    tempNsigBK = nnz(ICsigagg.ICwcfg0_presentations(iprobe).PkwBK<0.05 | ICsigagg.ICwcfg1_presentations(iprobe).PkwBK<0.05);
    fprintf('%s (cumulative, out of PkwBK<0.05) IC-enc %.2f%% RC-enc %.2f%% ind-enc %.2f%%\n', visareas{iprobe}, ...
        100*nnz(ICsigagg.ICwcfg0_presentations(iprobe).ICencoder==1 | ICsigagg.ICwcfg1_presentations(iprobe).ICencoder==1)/tempNsigBK, ...
        100*nnz(ICsigagg.ICwcfg0_presentations(iprobe).RCencoder==1 | ICsigagg.ICwcfg1_presentations(iprobe).RCencoder==1)/tempNsigBK, ...
        100*nnz(ICsigagg.ICwcfg0_presentations(iprobe).inducerencoder==1 | ICsigagg.ICwcfg1_presentations(iprobe).inducerencoder==1)/tempNsigBK)
end

% black inducer block
disp('ICkcfg0|ICkcfg0')
for a = 1:numel(probes)
    iprobe = probehierarchy(a);
    fprintf('%s (cumulative, out of all) IC-enc %.2f%% RC-enc %.2f%% ind-enc %.2f%%\n', visareas{iprobe}, ...
        100*mean(ICsigagg.ICkcfg0_presentations(iprobe).ICencoder==1 | ICsigagg.ICkcfg1_presentations(iprobe).ICencoder==1), ...
        100*mean(ICsigagg.ICkcfg0_presentations(iprobe).RCencoder==1 | ICsigagg.ICkcfg1_presentations(iprobe).RCencoder==1), ...
        100*mean(ICsigagg.ICkcfg0_presentations(iprobe).inducerencoder==1 | ICsigagg.ICkcfg1_presentations(iprobe).inducerencoder==1))
end

for a = 1:numel(probes)
    iprobe = probehierarchy(a);
    tempNsigBK = nnz(ICsigagg.ICkcfg0_presentations(iprobe).PkwBK<0.05 | ICsigagg.ICkcfg1_presentations(iprobe).PkwBK<0.05);
    fprintf('%s (cumulative, out of PkwBK<0.05) IC-enc %.2f%% RC-enc %.2f%% ind-enc %.2f%%\n', visareas{iprobe}, ...
        100*nnz(ICsigagg.ICkcfg0_presentations(iprobe).ICencoder==1 | ICsigagg.ICkcfg1_presentations(iprobe).ICencoder==1)/tempNsigBK, ...
        100*nnz(ICsigagg.ICkcfg0_presentations(iprobe).RCencoder==1 | ICsigagg.ICkcfg1_presentations(iprobe).RCencoder==1)/tempNsigBK, ...
        100*nnz(ICsigagg.ICkcfg0_presentations(iprobe).inducerencoder==1 | ICsigagg.ICkcfg1_presentations(iprobe).inducerencoder==1)/tempNsigBK)
end

% for each block
whichvisblock = 'ICkcfg1_presentations';
disp(whichvisblock)
for a = 1:numel(probes)
    iprobe = probehierarchy(a);
    fprintf('%s (cumulative, out of all) IC-enc %.2f%% RC-enc %.2f%% ind-enc %.2f%%\n', visareas{iprobe}, ...
        100*mean(ICsigagg.(whichvisblock)(iprobe).ICencoder==1), ...
        100*mean(ICsigagg.(whichvisblock)(iprobe).RCencoder==1), ...
        100*mean(ICsigagg.(whichvisblock)(iprobe).inducerencoder==1))
end

for a = 1:numel(probes)
    iprobe = probehierarchy(a);
    fprintf('%s (cumulative, out of PkwBK<0.05) IC-enc %.2f%% RC-enc %.2f%% ind-enc %.2f%%\n', visareas{iprobe}, ...
        100*nnz(ICsigagg.(whichvisblock)(iprobe).ICencoder==1)/nnz(ICsigagg.(whichvisblock)(iprobe).PkwBK<0.05), ...
        100*nnz(ICsigagg.(whichvisblock)(iprobe).RCencoder==1)/nnz(ICsigagg.(whichvisblock)(iprobe).PkwBK<0.05), ...
        100*nnz(ICsigagg.(whichvisblock)(iprobe).inducerencoder==1)/nnz(ICsigagg.(whichvisblock)(iprobe).PkwBK<0.05))
end

%% orientation tuning
% ICRCorivec reorders images to: IC02 IC12 IC01 IC11 RC02 RC12 RC01 RC11
% orientation denotation is same as psychtoolbox (0 is 12h, 45 is 1h30m, clockwise)
ICRCorilabel = {'IC02', 'IC12', 'IC01', 'IC11', 'RC02', 'RC12', 'RC01', 'RC11'};
wICRCenc = [ICsigagg.ICwcfg0_presentations(iprobe).ICencoder2 ...
    ICsigagg.ICwcfg1_presentations(iprobe).ICencoder2 ...
    ICsigagg.ICwcfg0_presentations(iprobe).ICencoder1 ...
    ICsigagg.ICwcfg1_presentations(iprobe).ICencoder1 ...
    ICsigagg.ICwcfg0_presentations(iprobe).RCencoder2 ...
    ICsigagg.ICwcfg1_presentations(iprobe).RCencoder2 ...
    ICsigagg.ICwcfg0_presentations(iprobe).RCencoder1 ...
    ICsigagg.ICwcfg1_presentations(iprobe).RCencoder1];

kICRCenc = [ICsigagg.ICkcfg0_presentations(iprobe).ICencoder2 ...
    ICsigagg.ICkcfg1_presentations(iprobe).ICencoder2 ...
    ICsigagg.ICkcfg0_presentations(iprobe).ICencoder1 ...
    ICsigagg.ICkcfg1_presentations(iprobe).ICencoder1 ...
    ICsigagg.ICkcfg0_presentations(iprobe).RCencoder2 ...
    ICsigagg.ICkcfg1_presentations(iprobe).RCencoder2 ...
    ICsigagg.ICkcfg0_presentations(iprobe).RCencoder1 ...
    ICsigagg.ICkcfg1_presentations(iprobe).RCencoder1];

figure
for e = 1:size(wICRCenc,2)
    subplot(2,4,e); hold all
    plot(1:8, oriparamsagg(iprobe).Rori(wICRCenc(:,e)==1,:))
    plot(1:8, mean(oriparamsagg(iprobe).Rori(wICRCenc(:,e)==1,:),1), 'ko-', 'LineWidth', 2)
    title(ICRCorilabel{e})
end

figure
for e = 1:size(kICRCenc,2)
    subplot(2,4,e); hold all
    plot(1:8, oriparamsagg(iprobe).Rori(kICRCenc(:,e)==1,:))
    plot(1:8, mean(oriparamsagg(iprobe).Rori(kICRCenc(:,e)==1,:),1), 'ko-', 'LineWidth', 2)
    title(ICRCorilabel{e})
end

%%
dirvec = unique(vis.sizeCI_presentations.orientation);
Ndirs = numel(dirvec);
wICRCenc_prefidirdist = zeros(Ndirs, size(wICRCenc,2));
wICRCenc_prefioridist = zeros(Ndirs/2, size(wICRCenc,2));
for e = 1:size(wICRCenc,2)
    wICRCenc_prefidirdist(:,e) = histcounts(oriparamsagg(iprobe).prefiori(wICRCenc(:,e)==1), 0.5:1:Ndirs+0.5);%, 'normalization', 'probability');
    wICRCenc_prefioridist(:,e) = histcounts(mod(oriparamsagg(iprobe).prefiori(wICRCenc(:,e)==1)-1, Ndirs/2)+1, ...
        0.5:1:Ndirs/2+0.5);%, 'normalization', 'probability');
end

kICRCenc_prefidirdist = zeros(Ndirs, size(kICRCenc,2));
kICRCenc_prefioridist = zeros(Ndirs/2, size(kICRCenc,2));
for e = 1:size(wICRCenc,2)
    kICRCenc_prefidirdist(:,e) = histcounts(oriparamsagg(iprobe).prefiori(kICRCenc(:,e)==1), 0.5:1:Ndirs+0.5);%, 'normalization', 'probability');
    kICRCenc_prefioridist(:,e) = histcounts(mod(oriparamsagg(iprobe).prefiori(kICRCenc(:,e)==1)-1, Ndirs/2)+1, ...
        0.5:1:Ndirs/2+0.5);%, 'normalization', 'probability');
end

figure; 
subplot(2,2,1)
imagesc(wICRCenc_prefidirdist)
colorbar
set(gca, 'YTick', 1:Ndirs, 'YTickLabel', dirvec, 'XTick', 1:size(wICRCenc,2), 'XTickLabel', ICRCorilabel)
xlabel('IC/RC-encoders')
ylabel('Pref Dir')
title('White Inducers, Black Background')

subplot(2,2,2)
imagesc(kICRCenc_prefidirdist)
colorbar
set(gca, 'YTick', 1:Ndirs, 'YTickLabel', dirvec, 'XTick', 1:size(kICRCenc,2), 'XTickLabel', ICRCorilabel)
xlabel('IC/RC-encoders')
ylabel('Pref Dir')
title('Black Inducers, White Background')

subplot(2,2,3)
imagesc(wICRCenc_prefioridist)
colorbar
set(gca, 'YTick', 1:Ndirs/2, 'YTickLabel', dirvec(1:Ndirs/2), 'XTick', 1:size(wICRCenc,2), 'XTickLabel', ICRCorilabel)
xlabel('IC/RC-encoders')
ylabel('Pref Ori')
title('White Inducers, Black Background')

subplot(2,2,4)
imagesc(kICRCenc_prefioridist)
colorbar
set(gca, 'YTick', 1:Ndirs/2, 'YTickLabel', dirvec(1:Ndirs/2), 'XTick', 1:size(kICRCenc,2), 'XTickLabel', ICRCorilabel)
xlabel('IC/RC-encoders')
ylabel('Pref Ori')
title('Black Inducers, White Background')

%% receptive field distribution
iprobe = find( strcmp('C', probes) );
tempneusigBK = ICsigagg.ICwcfg0_presentations(iprobe).PkwBK<0.05;

ICencagg = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder==1 | ICsigagg.ICwcfg0_presentations(iprobe).ICencoder==1;
RCencagg = ICsigagg.ICwcfg1_presentations(iprobe).RCencoder==1 | ICsigagg.ICwcfg0_presentations(iprobe).RCencoder==1;

% ICencagg = ICsigagg.ICkcfg1_presentations(iprobe).ICencoder==1 | ICsigagg.ICkcfg0_presentations(iprobe).ICencoder==1;
% RCencagg = ICsigagg.ICkcfg1_presentations(iprobe).RCencoder==1 | ICsigagg.ICkcfg0_presentations(iprobe).RCencoder==1;

RFdist_IC_acc = histcounts(RFCIagg(iprobe).RFindclassic(tempneuic), 0.5:1:9.5);
RFdist_sigBK_acc = histcounts(RFCIagg(iprobe).RFindclassic(tempneusigBK), 0.5:1:9.5);
RFdist_all_acc = histcounts(RFCIagg(iprobe).RFindclassic, 0.5:1:9.5);

RFdist_IC_agg = zeros(Nsessions, 9);
RFdist_sigBK_agg = zeros(Nsessions, 9);
RFdist_all_agg = zeros(Nsessions, 9);
for ises = 1:Nsessions
    tempsesneu = sesneuoiagg{iprobe}==ises;
RFdist_IC_agg(ises,:) = histcounts(RFCIagg(iprobe).RFindclassic(tempneuic(tempsesneu)), 0.5:1:9.5);
RFdist_sigBK_agg(ises,:) = histcounts(RFCIagg(iprobe).RFindclassic(tempneusigBK(tempsesneu)), 0.5:1:9.5);
RFdist_all_agg(ises,:) = histcounts(RFCIagg(iprobe).RFindclassic(tempsesneu), 0.5:1:9.5);
end

figure; hold all
plot(1:9, RFdist_IC_agg./RFdist_sigBK_agg)
plot(1:9, RFdist_IC_acc./RFdist_sigBK_acc, 'k.-', 'LineWidth', 2)

figure; hold all
neudescs = cell(1,0);
for ii = 1:3
    switch ii
        case 1
            tempneu = true(size(ICencagg));
            tempcol = [0 0 0];
            tempneudesc = 'All';
        case 2
            tempneu = RCencagg;
            tempcol = [1 0.5 0];
            tempneudesc = 'RC-encoder';
        case 3
            tempneu = ICencagg;
            tempcol = [0 0.7 0];
            tempneudesc = 'IC-encoder';
    end
    neudescs{ii} = tempneudesc;
histogram(RFCIagg(iprobe).RFindclassic(tempneu), 0.5:1:Nrfs+.5, 'FaceColor', tempcol, 'FaceAlpha', 0.5, 'normalization', 'probability')
end
legend(neudescs)
xlabel('CRF Position')
ylabel('Probability')

figure; hold all
neudescs = cell(1,0);
for ii = 1:3
    switch ii
        case 1
            tempneu = true(size(ICencagg));
            tempcol = [0 0 0];
            tempneudesc = 'All';
        case 2
            tempneu = RCencagg;
            tempcol = [1 0.5 0];
            tempneudesc = 'RC-encoder';
        case 3
            tempneu = ICencagg;
            tempcol = [0 0.7 0];
            tempneudesc = 'IC-encoder';
    end
    neudescs{ii} = tempneudesc;
histogram(RFCIagg(iprobe).RFindinverse(tempneu), 0.5:1:Nrfs+.5, 'FaceColor', tempcol, 'FaceAlpha', 0.5, 'normalization', 'probability')
end
legend(neudescs)
xlabel('IRF Position')
ylabel('Probability')

lw = 3;
figure; hold all
neudescs = cell(1,0);
for ii = 1:3
    switch ii
        case 1
            tempneu = true(size(ICencagg));
            tempcol = [0 0 0];
            tempneudesc = 'All';
        case 2
            tempneu = RCencagg;
            tempcol = [1 0.5 0];
            tempneudesc = 'RC-encoder';
        case 3
            tempneu = ICencagg;
            tempcol = [0 0.7 0];
            tempneudesc = 'IC-encoder';
    end
    neudescs{ii} = tempneudesc;
errorbar(1:Nrfs, mean(RFCIagg(iprobe).Rrfclassic(tempneu,:),1), std(RFCIagg(iprobe).Rrfclassic(tempneu,:),0,1)/sqrt(nnz(tempneu)), 'o-', 'Color', tempcol, 'MarkerFaceColor', tempcol, 'LineWidth', lw)
end
legend(neudescs)
set(gca, 'XTick', szvec)
xlabel('CRF Position')
ylabel('Rate (Hz)')

figure; hold all
neudescs = cell(1,0);
for ii = 1:3
    switch ii
        case 1
            tempneu = true(size(ICencagg));
            tempcol = [0 0 0];
            tempneudesc = 'All';
        case 2
            tempneu = RCencagg;
            tempcol = [1 0.5 0];
            tempneudesc = 'RC-encoder';
        case 3
            tempneu = ICencagg;
            tempcol = [0 0.7 0];
            tempneudesc = 'IC-encoder';
    end
    neudescs{ii} = tempneudesc;
errorbar(1:Nrfs, mean(RFCIagg(iprobe).Rrfinverse(tempneu,:),1), std(RFCIagg(iprobe).Rrfinverse(tempneu,:),0,1)/sqrt(nnz(tempneu)), 'o-', 'Color', tempcol, 'MarkerFaceColor', tempcol, 'LineWidth', lw)
end
legend(neudescs)
set(gca, 'XTick', szvec)
xlabel('IRF Position')
ylabel('Rate (Hz)')

tempneu = true(size(RFCIagg(iprobe).RFindclassic,1),1);
figure; histogram2(RFCIagg(iprobe).RFindclassic(tempneu), RFCIagg(iprobe).RFindinverse(tempneu), 0.5:1:Nrfs+.5, 0.5:1:Nrfs+.5, 'displaystyle', 'tile')

tempneu = true(size(RFCIagg(iprobe).RFindclassic,1),1);
% tempneu = RFCIagg(iprobe).Pkw_rfclassic<0.05;
% tempneu = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder==1 | ICsigagg.ICwcfg0_presentations(iprobe).ICencoder==1;
% tempneu = ICsigagg.ICkcfg1_presentations(iprobe).ICencoder==1 | ICsigagg.ICkcfg0_presentations(iprobe).ICencoder==1;
% tempneu = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder==1 | ICsigagg.ICwcfg0_presentations(iprobe).ICencoder==1 ...
%     | ICsigagg.ICkcfg1_presentations(iprobe).ICencoder==1 | ICsigagg.ICkcfg0_presentations(iprobe).ICencoder==1;
hc = histcounts2(RFCIagg(iprobe).RFindinverse(tempneu), RFCIagg(iprobe).RFindclassic(tempneu), 0.5:1:Nrfs+.5, 0.5:1:Nrfs+.5);
figure; imagesc(hc./sum(hc,1)); colorbar
xlabel('CRF Position')
ylabel('IRF Position')
title(sprintf('Proportion out of each CRF\nDiagonal Sum %.2f', sum(diag(hc./sum(hc,2)))))

hc = histcounts2(RFCIagg(iprobe).RFindclassic(tempneu), RFCIagg(iprobe).RFindinverse(tempneu), 0.5:1:Nrfs+.5, 0.5:1:Nrfs+.5);
figure; imagesc(hc./sum(hc,2)); colorbar
xlabel('IRF Position')
ylabel('CRF Position')

%% size tuning for cells with receptive field in the center
iprobe = find( strcmp('C', probes) );
szvec = [0, 4, 8, 16, 32, 64];

ctrCRFneurons = RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).Pkw_rfclassic<0.05;
ctrIRFneurons = RFCIagg(iprobe).RFindinverse==1 & RFCIagg(iprobe).Pkw_rfinverse<0.05;

ICencagg = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder==1 | ICsigagg.ICwcfg0_presentations(iprobe).ICencoder==1;
RCencagg = ICsigagg.ICwcfg1_presentations(iprobe).RCencoder==1 | ICsigagg.ICwcfg0_presentations(iprobe).RCencoder==1;

% ICencagg = ICsigagg.ICkcfg1_presentations(iprobe).ICencoder==1 | ICsigagg.ICkcfg0_presentations(iprobe).ICencoder==1;
% RCencagg = ICsigagg.ICkcfg1_presentations(iprobe).RCencoder==1 | ICsigagg.ICkcfg0_presentations(iprobe).RCencoder==1;

lw = 3;
figure; hold all
errorbar(szvec, mean(sizeCIagg(iprobe).Rsizeclassic(ctrCRFneurons,:),1), std(sizeCIagg(iprobe).Rsizeclassic(ctrCRFneurons,:),0,1)/sqrt(nnz(ctrCRFneurons)), 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', lw)
errorbar(szvec, mean(sizeCIagg(iprobe).Rsizeclassic(ctrCRFneurons & RCencagg,:),1), std(sizeCIagg(iprobe).Rsizeclassic(ctrCRFneurons & RCencagg,:),0,1)/sqrt(nnz(ctrCRFneurons & RCencagg)), 'o-', 'Color', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0], 'LineWidth', lw)
errorbar(szvec, mean(sizeCIagg(iprobe).Rsizeclassic(ctrCRFneurons & ICencagg,:),1), std(sizeCIagg(iprobe).Rsizeclassic(ctrCRFneurons & ICencagg,:),0,1)/sqrt(nnz(ctrCRFneurons & ICencagg)), 'o-', 'Color', [0 0.7 0], 'MarkerFaceColor', [0 0.7 0], 'LineWidth', lw)
legend({'ctrCRF', 'ctrCRF & RC-encoder', 'ctrCRF & IC-encoder'})
set(gca, 'XTick', szvec)
xlabel('Size (Visual Degrees)')
ylabel('Rate (Hz)')

figure; hold all
histogram(sizeCIagg(iprobe).sizeindclassic(ctrCRFneurons), 0.5:length(szvec)+0.5, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.5, 'normalization', 'probability')
histogram(sizeCIagg(iprobe).sizeindclassic(ctrCRFneurons & RCencagg), 0.5:length(szvec)+0.5, 'FaceColor', [1 0.5 0], 'FaceAlpha', 0.5, 'normalization', 'probability')
histogram(sizeCIagg(iprobe).sizeindclassic(ctrCRFneurons & ICencagg), 0.5:length(szvec)+0.5, 'FaceColor', [0 0.7 0], 'FaceAlpha', 0.5, 'normalization', 'probability')
legend({'ctrCRF', 'ctrCRF & RC-encoder', 'ctrCRF & IC-encoder'})
set(gca, 'XTick', 1:length(szvec), 'XTickLabel', szvec)
xlabel('Preferred Size')
ylabel('Probability')

lw = 3;
figure; hold all
neudescs = cell(1,0);
for ii = 1:3
    switch ii
        case 1
            tempneu = true(size(tempneuic));
            tempcol = [0 0 0];
            tempneudesc = 'All';
        case 2
            tempneu = RCencagg;
            tempcol = [1 0.5 0];
            tempneudesc = 'RC-encoder';
        case 3
            tempneu = ICencagg;
            tempcol = [0 0.7 0];
            tempneudesc = 'IC-encoder';
    end
    neudescs{ii} = tempneudesc;
errorbar(szvec, mean(sizeCIagg(iprobe).Rsizeinverse(tempneu,:),1), std(sizeCIagg(iprobe).Rsizeinverse(tempneu,:),0,1)/sqrt(nnz(tempneu)), 'o-', 'Color', tempcol, 'MarkerFaceColor', tempcol, 'LineWidth', lw)
end
legend(neudescs)
set(gca, 'XTick', szvec)
xlabel('IG Size (Visual Degrees)')
ylabel('Rate (Hz)')

lw = 3;
figure; hold all
neudescs = cell(1,0);
for ii = 1:3
    switch ii
        case 1
            tempneu = ctrCRFneurons;
            tempcol = [0 0 0];
            tempneudesc = 'ctrCRF';
        case 2
            tempneu = ctrCRFneurons & RCencagg;
            tempcol = [1 0.5 0];
            tempneudesc = 'ctrCRF & RC-encoder';
        case 3
            tempneu = ctrCRFneurons & ICencagg;
            tempcol = [0 0.7 0];
            tempneudesc = 'ctrCRF & IC-encoder';
    end
    neudescs{ii} = tempneudesc;
errorbar(szvec, mean(sizeCIagg(iprobe).Rsizeinverse(tempneu,:),1), std(sizeCIagg(iprobe).Rsizeinverse(tempneu,:),0,1)/sqrt(nnz(tempneu)), 'o-', 'Color', tempcol, 'MarkerFaceColor', tempcol, 'LineWidth', lw)
end
legend(neudescs)
set(gca, 'XTick', szvec)
xlabel('IG Size (Visual Degrees)')
ylabel('Rate (Hz)')

lw = 3;
figure; hold all
neudescs = cell(1,0);
for ii = 1:3
    switch ii
        case 1
            tempneu = ctrIRFneurons;
            tempcol = [0 0 0];
            tempneudesc = 'ctrIRF';
        case 2
            tempneu = ctrIRFneurons & RCencagg;
            tempcol = [1 0.5 0];
            tempneudesc = 'ctrIRF & RC-encoder';
        case 3
            tempneu = ctrIRFneurons & ICencagg;
            tempcol = [0 0.7 0];
            tempneudesc = 'ctrIRF & IC-encoder';
    end
    neudescs{ii} = tempneudesc;
errorbar(szvec, mean(sizeCIagg(iprobe).Rsizeinverse(tempneu,:),1), std(sizeCIagg(iprobe).Rsizeinverse(tempneu,:),0,1)/sqrt(nnz(tempneu)), 'o-', 'Color', tempcol, 'MarkerFaceColor', tempcol, 'LineWidth', lw)
end
legend(neudescs)
set(gca, 'XTick', szvec)
xlabel('IG Size (Visual Degrees)')
ylabel('Rate (Hz)')

%% plot each IC blocks
whichvisblock = 'ICwcfg1_presentations';

tt2p = [106 107 110 111 1105 1109]';
tt2p = [106 107 110 111]';
ttcol = [0 0.5 0;
    0.75 0.5 0;
    1 0.75 0.25;
    0.25 0.75 0.25;
    0 0 0.5;
    0.25 0.25 0.75];
fs=14;
figure%('Position', [100 0 1800 1200])
annotation('textbox', [0.1 0.92 0.9 0.1], 'string', ['Neuropixels: IC-encoder1 ' whichvisblock], 'edgecolor', 'none', 'fontsize', fs)
for iprobe= 1:numel(visareas)
    subplot(2,3, visind(iprobe))
    neuoi = sesneuoiagg{iprobe}<=Nsessions & ICsigagg.ICwcfg1_presentations(iprobe).RCencoder;
    hold all
    for ii = 1:numel(tt2p)
        typi = ICtrialtypes==tt2p(ii);
        temppsth = convn(squeeze(psthavg_ctxagg.(whichvisblock){iprobe}(:,typi,neuoi)), kergauss, 'same');
        plot(psthtli/1000, squeeze(mean(temppsth,2)), '-', 'Color', ttcol(ii,:), 'LineWidth', 0.5)
        %shadedErrorBar(psthtli/1000, squeeze(mean(temppsth,2)), squeeze(std(temppsth,0,2)/sqrt(nnz(neuoi))), {'-', 'Color', ttcol(ii,:), 'LineWidth', 1}, 1)
    end
    xlim([-.1 .500])
%     ylim([0 25])
    set(gca, 'FontSize', fs)
    xlabel('Time (s)')
    ylabel('Rate (Hz)')
    title(sprintf('%s %s N=%d', whichvisblock, visareas{iprobe}, nnz(neuoi)), 'interpreter', 'none')
end

