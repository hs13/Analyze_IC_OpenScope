% NICencoder, NRCencoder, Ninducerencoder, Ninducerresponsive, NsigBK
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
iprobe = 3;
tempneuic = ICsigagg.ICwcfg0_presentations(iprobe).ICencoder==1;
tempneusigBK = ICsigagg.ICwcfg0_presentations(iprobe).PkwBK<0.05;
% neuoi = ICsigagg.ICwcfg1_presentations(iprobe).ICencoder==1 | ICsigagg.ICwcfg0_presentations(iprobe).ICencoder==1;
% neuoi = ICsigagg.ICkcfg1_presentations(iprobe).ICencoder==1 | ICsigagg.ICkcfg0_presentations(iprobe).ICencoder==1;

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
histogram(RFCIagg(iprobe).RFindclassic(tempneusigBK), 0.5:1:9.5, 'normalization', 'probability')
histogram(RFCIagg(iprobe).RFindclassic(tempneuic), 0.5:1:9.5, 'normalization', 'probability')

%% size tuning


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

