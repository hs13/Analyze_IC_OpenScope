datadir = 'S:/OpenScopeData/000248/';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));

ICblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations'};
visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
    'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
ICtrialtypes = [0 101 105 106 107 109 110 111 506 511 1105 1109 1201 1299 ...
    1301 1302 1303 1304 1305 1306 1307 1308];

%%
segcombos = {'BRTL_IC1_wcfg1', 'BRTR_LC1_wcfg1', 'BLTL_LC2_wcfg1', 'BLTR_IC2_wcfg1'};
Tcombos = {'BRTL_TC1_wcfg1', 'BRTR_TC1_wcfg1', 'BLTL_TC2_wcfg1', 'BLTR_TC2_wcfg1'};
encrespcombos = {'encresp_IC1_wcfg1', 'encresp_LC1_wcfg1', 'encresp_LC2_wcfg1', 'encresp_IC2_wcfg1'};
encsegcombos = {'encseg_IC1_wcfg1', 'encseg_LC1_wcfg1', 'encseg_LC2_wcfg1', 'encseg_IC2_wcfg1'};
encsegTcombos = {'IC1BR_TC1_wcfg1', 'IC1TL_TC1_wcfg1', 'IC1TR_TC1_wcfg1', 'IC2BL_TC2_wcfg1', 'IC2TR_TC2_wcfg1', 'IC2TL_TC2_wcfg1'};

Rscagg = struct();
for b = 1:numel(ICblocks)
    whichICblock = ICblocks{b};
Rscsegavg.(whichICblock) = zeros(numel(nwbsessions), numel(segcombos));
RscTavg.(whichICblock) = zeros(numel(nwbsessions), numel(segcombos));
Rscencrespavg.(whichICblock) = zeros(numel(nwbsessions), numel(encrespcombos));
Rscencsegavg.(whichICblock) = zeros(numel(nwbsessions), numel(encsegcombos));
RscencsegTavg.(whichICblock) = zeros(numel(nwbsessions), numel(encsegTcombos));
end

visttagg = struct();
Ravgagg = struct();
zRavgagg = struct();
neuinareaagg = cell(size(nwbsessions));
ICsigallagg = struct();

for ises = 1:numel(nwbsessions)
    clearvars Rscall Rall vis ICsigall
    tic
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'postprocessed.mat'])
    load([pathpp 'visresponses.mat'])

    neuinarea = contains(neuallloc, 'VISp') & ~contains(neuallloc, 'VISpm');
    neuinareaagg{ises} = neuinarea;

for b = 1:numel(visblocks)
    whichblock = visblocks{b};
tempzR = (Rall.(whichblock)-mean(Rall.(whichblock),1))./std(Rall.(whichblock),0,1);

    Nneurons = size(Rall.(whichblock),2);
    trialtypes = unique(vis.(whichblock).trialorder);
    Ravgblock = NaN(Nneurons, numel(trialtypes));
    zRavgblock = NaN(Nneurons, numel(trialtypes));
    for typi = 1:numel(trialtypes)
        trialsoi = vis.(whichblock).trialorder==trialtypes(typi);
        Ravgblock(:,typi) = mean(Rall.(whichblock)(trialsoi,:),1)';
        zRavgblock(:,typi) = mean(tempzR(trialsoi,:),1)';
    end
    Ravgagg(ises).(whichblock) = Ravgblock;
    zRavgagg(ises).(whichblock) = zRavgblock;
    visttagg.(whichblock) = trialtypes;
end

for b = 1:numel(ICblocks)
    whichICblock = ICblocks{b};

    ICsigallagg(ises).(whichICblock) = ICsigall.(whichICblock);

    Nneurons = size(Rall.(whichICblock),2);
    trialtypes = unique(vis.(whichICblock).trialorder);
    Rscall = NaN(Nneurons, Nneurons, numel(trialtypes));
    for typi = 1:numel(trialtypes)
        trialsoi = vis.(whichICblock).trialorder==trialtypes(typi);
        tempRsc = corr(Rall.(whichICblock)(trialsoi,:));
        tempRsc(eye(Nneurons)==1) = NaN;
        Rscall(:,:,typi) = tempRsc;
    end

    for f = 1:numel(segcombos)
        switch segcombos{f}
            case 'BRTL_IC1_wcfg1'
                typi = ICtrialtypes==106;
                tempneu1 = neuinarea & ICsigall.(whichICblock).indin1==1;
                tempneu2 = neuinarea & ICsigall.(whichICblock).indin3==1;
            case 'BRTR_LC1_wcfg1'
                typi = ICtrialtypes==107;
                tempneu1 = neuinarea & ICsigall.(whichICblock).indin1==1;
                tempneu2 = neuinarea & ICsigall.(whichICblock).indin4==1;
            case 'BLTL_LC2_wcfg1'
                typi = ICtrialtypes==110;
                tempneu1 = neuinarea & ICsigall.(whichICblock).indin2==1;
                tempneu2 = neuinarea & ICsigall.(whichICblock).indin3==1;
            case 'BLTR_IC2_wcfg1'
                typi = ICtrialtypes==111;
                tempneu1 = neuinarea & ICsigall.(whichICblock).indin2==1;
                tempneu2 = neuinarea & ICsigall.(whichICblock).indin4==1;
        end
        Rscagg.(whichICblock)(ises).(segcombos{f}) = reshape(Rscall(tempneu1, tempneu2, typi), [],1);
        Rscsegavg.(whichICblock)(ises, f) = squeeze( nanmean(Rscall(tempneu1, tempneu2, typi), [1,2]) );
    end

    for f = 1:numel(Tcombos)
        switch Tcombos{f}
            case 'BRTL_TC1_wcfg1'
                typi = ICtrialtypes==105;
                tempneu1 = neuinarea & ICsigall.(whichICblock).indin1==1;
                tempneu2 = neuinarea & ICsigall.(whichICblock).indin3==1;
            case 'BRTR_TC1_wcfg1'
                typi = ICtrialtypes==105;
                tempneu1 = neuinarea & ICsigall.(whichICblock).indin1==1;
                tempneu2 = neuinarea & ICsigall.(whichICblock).indin4==1;
            case 'BLTL_TC2_wcfg1'
                typi = ICtrialtypes==109;
                tempneu1 = neuinarea & ICsigall.(whichICblock).indin2==1;
                tempneu2 = neuinarea & ICsigall.(whichICblock).indin3==1;
            case 'BLTR_TC2_wcfg1'
                typi = ICtrialtypes==109;
                tempneu1 = neuinarea & ICsigall.(whichICblock).indin2==1;
                tempneu2 = neuinarea & ICsigall.(whichICblock).indin4==1;
        end
        Rscagg.(whichICblock)(ises).(Tcombos{f}) = reshape(Rscall(tempneu1, tempneu2, typi), [],1);
        RscTavg.(whichICblock)(ises, f) = squeeze( nanmean(Rscall(tempneu1, tempneu2, typi), [1,2]) );
    end

    for f = 1:numel(encrespcombos)
        switch encrespcombos{f}
            case 'encresp_IC1_wcfg1'
                typi = ICtrialtypes==106;
                tempneu1 = neuinarea & ICsigall.(whichICblock).ICencoder1==1;
                tempneu2 = neuinarea & (ICsigall.(whichICblock).ICresp1==1 & ICsigall.(whichICblock).ICencoder1==0);
            case 'encresp_LC1_wcfg1'
                typi = ICtrialtypes==107;
                tempneu1 = neuinarea & ICsigall.(whichICblock).RCencoder1==1;
                tempneu2 = neuinarea & (ICsigall.(whichICblock).RCresp1==1 & ICsigall.(whichICblock).RCencoder1==0);
            case 'encresp_LC2_wcfg1'
                typi = ICtrialtypes==110;
                tempneu1 = neuinarea & ICsigall.(whichICblock).RCencoder2==1;
                tempneu2 = neuinarea & (ICsigall.(whichICblock).RCresp2==1 & ICsigall.(whichICblock).RCencoder2==0);
            case 'encresp_IC2_wcfg1'
                typi = ICtrialtypes==111;
                tempneu1 = neuinarea & ICsigall.(whichICblock).ICencoder2==1;
                tempneu2 = neuinarea & (ICsigall.(whichICblock).ICresp2==1 & ICsigall.(whichICblock).ICencoder2==0);
        end
        Rscagg.(whichICblock)(ises).(encrespcombos{f}) = reshape(Rscall(tempneu1, tempneu2, typi), [],1);
        Rscencrespavg.(whichICblock)(ises, f) = squeeze( nanmean(Rscall(tempneu1, tempneu2, typi), [1,2]) );
    end

    for f = 1:numel(encsegcombos)
        switch encsegcombos{f}
            case 'encseg_IC1_wcfg1'
                typi = ICtrialtypes==106;
                tempneu1 = neuinarea & ICsigall.(whichICblock).ICencoder1==1;
                tempneu2 = neuinarea & (ICsigall.(whichICblock).indin1==1 | ICsigall.(whichICblock).indin3==1);
            case 'encseg_LC1_wcfg1'
                typi = ICtrialtypes==107;
                tempneu1 = neuinarea & ICsigall.(whichICblock).RCencoder1==1;
                tempneu2 = neuinarea & (ICsigall.(whichICblock).indin1==1 | ICsigall.(whichICblock).indin4==1);
            case 'encseg_LC2_wcfg1'
                typi = ICtrialtypes==110;
                tempneu1 = neuinarea & ICsigall.(whichICblock).RCencoder2==1;
                tempneu2 = neuinarea & (ICsigall.(whichICblock).indin2==1 | ICsigall.(whichICblock).indin3==1);
            case 'encseg_IC2_wcfg1'
                typi = ICtrialtypes==111;
                tempneu1 = neuinarea & ICsigall.(whichICblock).ICencoder2==1;
                tempneu2 = neuinarea & (ICsigall.(whichICblock).indin2==1 | ICsigall.(whichICblock).indin4==1);
        end
        Rscagg.(whichICblock)(ises).(encsegcombos{f}) = reshape(Rscall(tempneu1, tempneu2, typi), [],1);
        Rscencsegavg.(whichICblock)(ises, f) = squeeze( nanmean(Rscall(tempneu1, tempneu2, typi), [1,2]) );
    end

    for f = 1:numel(encsegTcombos)
        switch encsegTcombos{f}
            case 'IC1BR_TC1_wcfg1'
                typi = ICtrialtypes==105;
                tempneu1 = neuinarea & ICsigall.(whichICblock).ICencoder1==1;
                tempneu2 = neuinarea & ICsigall.(whichICblock).indin1==1;
            case 'IC1TL_TC1_wcfg1'
                typi = ICtrialtypes==105;
                tempneu1 = neuinarea & ICsigall.(whichICblock).ICencoder1==1;
                tempneu2 = neuinarea & ICsigall.(whichICblock).indin3==1;
            case 'IC1TR_TC1_wcfg1'
                typi = ICtrialtypes==105;
                tempneu1 = neuinarea & ICsigall.(whichICblock).ICencoder1==1;
                tempneu2 = neuinarea & ICsigall.(whichICblock).indin4==1;
            case 'IC2BL_TC2_wcfg1'
                typi = ICtrialtypes==109;
                tempneu1 = neuinarea & ICsigall.(whichICblock).ICencoder2==1;
                tempneu2 = neuinarea & ICsigall.(whichICblock).indin2==1;
            case 'IC2TR_TC2_wcfg1'
                typi = ICtrialtypes==109;
                tempneu1 = neuinarea & ICsigall.(whichICblock).ICencoder2==1;
                tempneu2 = neuinarea & ICsigall.(whichICblock).indin4==1;
            case 'IC2TL_TC2_wcfg1'
                typi = ICtrialtypes==109;
                tempneu1 = neuinarea & ICsigall.(whichICblock).ICencoder2==1;
                tempneu2 = neuinarea & ICsigall.(whichICblock).indin3==1;
        end
        Rscagg.(whichICblock)(ises).(encsegTcombos{f}) = reshape(Rscall(tempneu1, tempneu2, typi), [],1);
        RscencsegTavg.(whichICblock)(ises, f) = squeeze( nanmean(Rscall(tempneu1, tempneu2, typi), [1,2]) );
    end
end

    toc
end

ICsigagg = struct();
for b = 1:numel(ICblocks)
    whichICblock = ICblocks{b};
    ICsigagg.(whichICblock) = cat(1,ICsigallagg.(whichICblock));
end

%%
RavgsizeCI = cat(1,Ravgagg.sizeCI_presentations);
Ravgori = NaN(size(RavgsizeCI,1), 4);
for iori = 1:4
    typoi = ismember(visttagg.sizeCI_presentations, iori+[6000 6004 11000 11004 12000 12004]);
    % typoi = ismember(visttagg.sizeCI_presentations, iori+[6000 11000 12000]);
    % typoi = ismember(visttagg.sizeCI_presentations, iori+[6004 11004 12004]);
    Ravgori(:,iori) = mean(RavgsizeCI(:,typoi),2);
end
Ravgblank = mean(RavgsizeCI(:,ismember(visttagg.sizeCI_presentations,1001:1008)),2);
Ravgdir = NaN(size(RavgsizeCI,1), 8);
for idir = 1:8
    typoi = ismember(visttagg.sizeCI_presentations, idir+[6000 11000 12000]);
    Ravgdir(:,idir) = mean(RavgsizeCI(:,typoi),2);
end

zRavgsizeCI = cat(1,zRavgagg.sizeCI_presentations);
zRavgori = NaN(size(zRavgsizeCI,1), 4);
for iori = 1:4
    typoi = ismember(visttagg.sizeCI_presentations, iori+[6000 6004 11000 11004 12000 12004]);
    % typoi = ismember(visttagg.sizeCI_presentations, iori+[6000 11000 12000]);
    % typoi = ismember(visttagg.sizeCI_presentations, iori+[6004 11004 12004]);
    zRavgori(:,iori) = mean(zRavgsizeCI(:,typoi),2);
end
zRavgblank = mean(zRavgsizeCI(:,ismember(visttagg.sizeCI_presentations,1001:1008)),2);
zRavgdir = NaN(size(zRavgsizeCI,1), 8);
for idir = 1:8
    typoi = ismember(visttagg.sizeCI_presentations, idir+[6000 11000 12000]);
    zRavgdir(:,idir) = mean(zRavgsizeCI(:,typoi),2);
end

whichICblock = 'ICwcfg0_presentations';
ICenc1agg = cat(1,neuinareaagg{:}) & cat(1,ICsigagg.(whichICblock).ICencoder1)==1;
ICenc2agg = cat(1,neuinareaagg{:}) & cat(1,ICsigagg.(whichICblock).ICencoder2)==1;
% ICenc1agg = cat(1,ICsigagg.(whichICblock).ICencoder1)==1;
% ICenc2agg = cat(1,ICsigagg.(whichICblock).ICencoder2)==1;
tempRori = Ravgori;
tempRdir = Ravgdir;

fs=14;
ylab = 'Rate (Hz)';
ttcol = [0 .4 0; 0 1 0];
legs = {'I_C_1-enc', 'I_C_2-enc'};
yl = [-0.05 0.1];
yl = [-0.04 0.06];
% figure('Position', [100 100 270 240])
figure('Position', [100 100 250 220])
hold all
errorbar(0:45:180-1, mean(tempRori(ICenc1agg,:),1), std(tempRori(ICenc1agg,:),0,1)/sqrt(nnz(ICenc1agg)), 'o-', 'color', ttcol(1,:), 'markerfacecolor', ttcol(1,:), 'linewidth', 2)
errorbar(0:45:180-1, mean(tempRori(ICenc2agg,:),1), std(tempRori(ICenc2agg,:),0,1)/sqrt(nnz(ICenc2agg)), 'o-', 'color', ttcol(end,:), 'markerfacecolor', ttcol(end,:), 'linewidth', 2)
% legend(legs, 'location', 'northwest', 'FontSize', fs)
% legend('boxoff')
for itt = 1:numel(legs)
    if isequal(yl, [-0.05 0.1])
text(4.5, yl(2)-(itt-1)*.15*range(yl), legs{itt}, 'Color', ttcol(itt,:), 'FontSize', fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
    else
text(135-4.5, yl(1)+(numel(legs)-itt)*.15*range(yl), legs{itt}, 'Color', ttcol(itt,:), 'FontSize', fs, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')
    end
end
xlim([0 135])
% ylim(yl)
set(gca, 'XTick', 0:45:180-1, 'FontSize', fs)
xlabel('Orientation (°)', 'FontSize', fs)
ylabel(ylab, 'FontSize', fs)

ylab = 'Rate (Hz)';
ttcol = [0 .4 0; 0 1 0];
legs = {'I_C_1-enc', 'I_C_2-enc'};
yl = [-0.05 0.1];
yl = [-0.04 0.06];
% figure('Position', [100 100 270 240])
figure('Position', [100 100 250 220])
hold all
errorbar(0:45:360-1, mean(tempRdir(ICenc1agg,:),1), std(tempRdir(ICenc1agg,:),0,1)/sqrt(nnz(ICenc1agg)), 'o-', 'color', ttcol(1,:), 'markerfacecolor', ttcol(1,:), 'linewidth', 2)
errorbar(0:45:360-1, mean(tempRdir(ICenc2agg,:),1), std(tempRdir(ICenc2agg,:),0,1)/sqrt(nnz(ICenc2agg)), 'o-', 'color', ttcol(end,:), 'markerfacecolor', ttcol(end,:), 'linewidth', 2)
% legend(legs, 'location', 'northwest', 'FontSize', fs)
% legend('boxoff')
for itt = 1:numel(legs)
    if isequal(yl, [-0.05 0.1])
text(4.5, yl(2)-(itt-1)*.15*range(yl), legs{itt}, 'Color', ttcol(itt,:), 'FontSize', fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
    else
text(315-4.5, yl(1)+(numel(legs)-itt)*.15*range(yl), legs{itt}, 'Color', ttcol(itt,:), 'FontSize', fs, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')
    end
end
xlim([0 315])
% ylim(yl)
set(gca, 'XTick', 0:45:360-1, 'FontSize', fs)
xlabel('Orientation (°)', 'FontSize', fs)
ylabel(ylab, 'FontSize', fs)

%% static grating vs drifting grating?
Nori = size(Ravgori,2);
[Rprefori,iprefori]=max(Ravgori,[],2);
iorthori = mod(iprefori-1+Nori/2,Nori)+1;
isequal(Rprefori, Ravgori(sub2ind(size(Ravgori), (1:size(Ravgori,1))', iprefori)))
Rorthori = Ravgori(sub2ind(size(Ravgori), (1:size(Ravgori,1))', iorthori));
osi = (Rprefori-Rorthori)./(Rprefori+Rorthori);
figure; hold all
hc= histogram(osi, -10:0.05:10, 'Normalization', 'probability');
histogram(osi(ICenc1agg),hc.BinEdges, 'Normalization', 'probability')
histogram(osi(ICenc2agg),hc.BinEdges, 'Normalization', 'probability')

ranksum(osi, osi(ICenc1agg | ICenc2agg))

prctile(osi,[25 50 75])
prctile(osi(ICenc1agg | ICenc2agg),[25 50 75])

figure; hold all
histogram(iprefori(ICenc1agg), 0.5:1:4.5)
histogram(iprefori(ICenc2agg), 0.5:1:4.5)

figure; plot(iprefori(ICenc2agg), osi(ICenc2agg), 'o')
[p,t,s]=kruskalwallis(osi(ICenc2agg),iprefori(ICenc2agg));
figure; multcompare(s)

figure; imagesc(tempRori(ICenc1agg,:))

RavgICwcfg1 = cat(1,Ravgagg.ICwcfg1_presentations);
figure; plot(1:4, RavgICwcfg1(ICenc1agg,ismember(ICtrialtypes,[106 107 110 111])), 'o-')

Ndir = size(Ravgdir,2);
[Rprefdir,iprefdir]=max(Ravgdir,[],2);
iorthdir = mod(iprefdir-1+Ndir/2,Ndir)+1;
isequal(Rprefdir, Ravgdir(sub2ind(size(Ravgdir), (1:size(Ravgdir,1))', iprefdir)))
Rorthdir = Ravgdir(sub2ind(size(Ravgdir), (1:size(Ravgdir,1))', iorthdir));
dsi = (Rprefdir-Rorthdir)./(Rprefdir+Rorthdir);

ranksum(dsi, dsi(ICenc1agg | ICenc2agg))
prctile(dsi,[25 50 75])
prctile(dsi(ICenc1agg | ICenc2agg),[25 50 75])

%%
whichICblock = 'ICwcfg1_presentations';
ICenc1agg = cat(1,neuinareaagg{:}) & cat(1,ICsigagg.(whichICblock).ICencoder1)==1;
ICenc2agg = cat(1,neuinareaagg{:}) & cat(1,ICsigagg.(whichICblock).ICencoder2)==1;

RavgICwcfg1agg = cat(1,Ravgagg.(whichICblock));
RbsICwcfg1agg = RavgICwcfg1agg-RavgICwcfg1agg(:,ICtrialtypes==0);

tempR = RavgICwcfg1agg(:, ismember(ICtrialtypes,[106 107 110 111]) );
disp([range(mean(tempR(ICenc1agg,:),1)) range(mean(tempR(ICenc2agg,:),1))])
vec1 = tempR(ICenc1agg,1)-mean(tempR(ICenc1agg,[2,3]),2);
vec2 = tempR(ICenc2agg,4)-mean(tempR(ICenc2agg,[2,3]),2);
disp([mean(vec1) mean(vec2) mean([vec1;vec2])])

ttcol = [0 .4 0; 0 1 0];
legs = {'I_C_1-enc', 'I_C_2-enc'};
figure('Position', [100 100 250 220])
hold all
errorbar(1:4, mean(tempR(ICenc1agg,:),1), std(tempR(ICenc1agg,:),0,1)/sqrt(nnz(ICenc1agg)), 'o-', 'color', ttcol(1,:), 'markerfacecolor', ttcol(1,:), 'linewidth', 2)
errorbar(1:4, mean(tempR(ICenc2agg,:),1), std(tempR(ICenc2agg,:),0,1)/sqrt(nnz(ICenc2agg)), 'o-', 'color', ttcol(end,:), 'markerfacecolor', ttcol(end,:), 'linewidth', 2)


%%
whichICblock = 'ICwcfg1_presentations';

figure
plot(1:4, Rscsegavg.(whichICblock), 'o-')
set(gca, 'XTick', 1:numel(segcombos), 'XTickLabel', segcombos)

temp = Rscsegavg.(whichICblock);
temp(any(isnan(temp),2),:)=[];
[P,TABLE,STATS] = friedman(temp);
figure; multcompare(STATS)

ises = 1;
X= [];
G= [];
for f = 1:numel(segcombos)
    X = cat(1, X, Rscagg.(whichICblock)(ises).(segcombos{f}));
    G = cat(1, G, f*ones(size(Rscagg.(whichICblock)(ises).(segcombos{f}))));
end
[P,ANOVATAB,STATS] = kruskalwallis(X,G);
figure; multcompare(STATS)

X= [];
G= [];
for f = 1:numel(segcombos)
    X = cat(1, X, cat(1,Rscagg.(whichICblock).(segcombos{f})));
    G = cat(1, G, f*ones(size(cat(1,Rscagg.(whichICblock).(segcombos{f})))));
end
[P,ANOVATAB,STATS] = kruskalwallis(X,G);
figure; multcompare(STATS)

ranksum( X(ismember(G,[1,4])), X(ismember(G,[2,3])) ) % not significant

%%
figure
plot(1:4, Rscencrespavg.(whichICblock), 'o-')
set(gca, 'XTick', 1:numel(encrespcombos), 'XTickLabel', encrespcombos)

temp = Rscencrespavg.(whichICblock);
temp(any(isnan(temp),2),:)=[];
[P,TABLE,STATS] = friedman(temp);
figure; multcompare(STATS)

X= [];
G= [];
for f = 1:numel(encrespcombos)
    X = cat(1, X, cat(1,Rscagg.(whichICblock).(encrespcombos{f})));
    G = cat(1, G, f*ones(size(cat(1,Rscagg.(whichICblock).(encrespcombos{f})))));
end
[P,ANOVATAB,STATS] = kruskalwallis(X,G);
figure; multcompare(STATS)

ranksum( X(ismember(G,[1,4])), X(ismember(G,[2,3])) ) % not significant

%%
figure
plot(1:4, Rscencsegavg.(whichICblock), 'o-')
set(gca, 'XTick', 1:numel(encsegcombos), 'XTickLabel', encsegcombos)

temp = Rscencsegavg.(whichICblock);
temp(any(isnan(temp),2),:)=[];
[P,TABLE,STATS] = friedman(temp);
figure; multcompare(STATS)

X= [];
G= [];
for f = 1:numel(encsegcombos)
    X = cat(1, X, cat(1,Rscagg.(whichICblock).(encsegcombos{f})));
    G = cat(1, G, f*ones(size(cat(1,Rscagg.(whichICblock).(encsegcombos{f})))));
end
[P,ANOVATAB,STATS] = kruskalwallis(X,G);
figure; multcompare(STATS)

ranksum( X(ismember(G,[1,4])), X(ismember(G,[2,3])) ) % not significant

%% encsegTcombos
figure
plot(1:6, RscencsegTavg.(whichICblock), 'o-')
set(gca, 'XTick', 1:numel(encsegTcombos), 'XTickLabel', encsegTcombos)

X= [];
G= [];
for f = 1:numel(encsegTcombos)
    X = cat(1, X, cat(1,Rscagg.(whichICblock).(encsegTcombos{f})));
    G = cat(1, G, f*ones(size(cat(1,Rscagg.(whichICblock).(encsegTcombos{f})))));
end
[P,ANOVATAB,STATS] = kruskalwallis(X,G);
figure; multcompare(STATS)

ranksum( X(ismember(G,[1,2,4,5])), X(ismember(G,[3,6])) )

GIL=G;
GIL(ismember(G,[1,2,4,5]))=1;
GIL(ismember(G,[3,6]))=2;
figure; boxplot(X, GIL)

[P,ANOVATAB,STATS] = kruskalwallis(X,GIL);
figure; multcompare(STATS)

%%
ranksum( cat(1,Rscagg.(whichICblock).BRTL_TC1_wcfg1), cat(1,Rscagg.(whichICblock).BRTR_TC1_wcfg1) )
ranksum( cat(1,Rscagg.(whichICblock).BLTR_TC2_wcfg1), cat(1,Rscagg.(whichICblock).BLTL_TC2_wcfg1) )

%% R1C7
whichICblock = 'ICwcfg1_presentations';

RavgICblockagg = cat(1, Ravgagg.(whichICblock));
RbsICblockagg = RavgICblockagg-RavgICblockagg(:,ICtrialtypes==0);

Ravg_Ic = mean(RavgICblockagg(cat(1,neuinareaagg{:}), ismember(ICtrialtypes,[106 111]) ),2);
Ravg_Lc = mean(RavgICblockagg(cat(1,neuinareaagg{:}), ismember(ICtrialtypes,[107 110]) ),2);
signrank(Ravg_Ic, Ravg_Lc)
disp([mean(Ravg_Ic-Ravg_Lc) median(Ravg_Ic-Ravg_Lc)])

Rbs_Ic = mean(RbsICblockagg(cat(1,neuinareaagg{:}), ismember(ICtrialtypes,[106 111]) ),2);
Rbs_Lc = mean(RbsICblockagg(cat(1,neuinareaagg{:}), ismember(ICtrialtypes,[107 110]) ),2);
signrank(Rbs_Ic, Rbs_Lc)

figure; histogram(Ravg_Ic-Ravg_Lc)

figure; hold all
plot(Ravg_Ic, Ravg_Lc, 'o')
xl=xlim; plot(xl,xl,'r--')

figure; hold all
plot(Rbs_Ic, Rbs_Lc, 'o')
xl=xlim; plot(xl,xl,'r--')

Ravgindin_Ic = [];
Ravgindin_Lc = [];
for f = 1:4
    switch f
        case 1
            whichseg = 'indin1';
            typIc = ICtrialtypes==106;
            typLc = ICtrialtypes==107;
        case 2
            whichseg = 'indin2';
            typIc = ICtrialtypes==111;
            typLc = ICtrialtypes==110;
        case 3
            whichseg = 'indin3';
            typIc = ICtrialtypes==106;
            typLc = ICtrialtypes==110;
        case 4
            whichseg = 'indin4';
            typIc = ICtrialtypes==111;
            typLc = ICtrialtypes==107;
    end
    tempneu = cat(1,neuinareaagg{:}) & cat(1,ICsigagg.(whichICblock).(whichseg))==1;
    Ravgindin_Ic = cat(1, Ravgindin_Ic, RavgICblockagg(tempneu, typIc));
    Ravgindin_Lc = cat(1, Ravgindin_Lc, RavgICblockagg(tempneu, typLc));
end
signrank(Ravgindin_Ic, Ravgindin_Lc)
disp([mean(Ravgindin_Ic-Ravgindin_Lc) median(Ravgindin_Ic-Ravgindin_Lc)])

Ravgindout_Ic = [];
Ravgindout_Lc = [];
for f = 1:4
    switch f
        case 1
            whichseg = 'indout2';
            typIc = ICtrialtypes==106;
            typLc = ICtrialtypes==107;
        case 2
            whichseg = 'indout1';
            typIc = ICtrialtypes==111;
            typLc = ICtrialtypes==110;
        case 3
            whichseg = 'indout4';
            typIc = ICtrialtypes==106;
            typLc = ICtrialtypes==110;
        case 4
            whichseg = 'indout3';
            typIc = ICtrialtypes==111;
            typLc = ICtrialtypes==107;
    end
    tempneu = cat(1,neuinareaagg{:}) & cat(1,ICsigagg.(whichICblock).(whichseg))==1;
    Ravgindout_Ic = cat(1, Ravgindout_Ic, RavgICblockagg(tempneu, typIc));
    Ravgindout_Lc = cat(1, Ravgindout_Lc, RavgICblockagg(tempneu, typLc));
end
p = signrank(Ravgindout_Ic, Ravgindout_Lc)
disp([mean(Ravgindout_Ic-Ravgindout_Lc) median(Ravgindout_Ic-Ravgindout_Lc)])

figure
hold all
plot(Ravgindout_Ic, Ravgindout_Lc, 'o')
xl = xlim;
plot(xl, xl, 'r--')
xlabel('I_C')
ylabel('L_C')
title(sprintf('indout p=%.4f', p))


for f = 1:4
    switch f
        case 1
            whichseg = 'indin1';
            typIc = ICtrialtypes==106;
            typLc = ICtrialtypes==107;
            typT = ICtrialtypes==105;
        case 2
            whichseg = 'indin2';
            typIc = ICtrialtypes==111;
            typLc = ICtrialtypes==110;
            typT = ICtrialtypes==109;
        case 3
            whichseg = 'indin3';
            typIc = ICtrialtypes==106;
            typLc = ICtrialtypes==110;
            typT = ICtrialtypes==105;
        case 4
            whichseg = 'indin4';
            typIc = ICtrialtypes==111;
            typLc = ICtrialtypes==107;
            typT = ICtrialtypes==109;
    end
    tempneu = cat(1,neuinareaagg{:}) & cat(1,ICsigagg.(whichICblock).(whichseg))==1;
    p1 = signrank(RavgICblockagg(tempneu, typIc), RavgICblockagg(tempneu, typLc)); % none of these are significant
    p2 = signrank(RavgICblockagg(tempneu, typIc), RavgICblockagg(tempneu, typT));
    p3 = signrank(RavgICblockagg(tempneu, typT), RavgICblockagg(tempneu, typLc));
    disp([p1 p2 p3])
end

figure
for f = 1:2
    switch f
        case 1
            tempneu = cat(1,neuinareaagg{:}) & cat(1,ICsigagg.(whichICblock).indin3)==1;
            typIc = ICtrialtypes==106;
            typLc = ICtrialtypes==110;
            typT = ICtrialtypes==109;
        case 2
            tempneu = cat(1,neuinareaagg{:}) & cat(1,ICsigagg.(whichICblock).indin4)==1;
            typIc = ICtrialtypes==111;
            typLc = ICtrialtypes==107;
            typT = ICtrialtypes==105;
    end
    p1 = signrank(RavgICblockagg(tempneu, typIc), RavgICblockagg(tempneu, typLc));
    p2 = signrank(RavgICblockagg(tempneu, typIc), RavgICblockagg(tempneu, typT));
    p3 = signrank(RavgICblockagg(tempneu, typT), RavgICblockagg(tempneu, typLc));

    disp([p1 p2 p3])
    subplot(2,2,2*(f-1)+1)
    hold all
    plot(RavgICblockagg(tempneu, typIc), RavgICblockagg(tempneu, typT), 'o')
    xl=xlim;
    plot(xl,xl,'r-')
    xlabel('I_C')
    ylabel('T')
    subplot(2,2,2*f)
    hold all
    plot(RavgICblockagg(tempneu, typLc), RavgICblockagg(tempneu, typT), 'o')
    plot(xl,xl,'r-')
    xlabel('L_C')
    ylabel('T')
end

[P,TABLE,STATS] = friedman(RavgICblockagg(tempneu, ismember(ICtrialtypes,[106 107 110 111])));
figure; multcompare(STATS)

figure
for f = 1:4
    switch f
        case 1
            whichseg = 'indout2';
            typIc = ICtrialtypes==106;
            typLc = ICtrialtypes==107;
        case 2
            whichseg = 'indout1';
            typIc = ICtrialtypes==111;
            typLc = ICtrialtypes==110;
        case 3
            whichseg = 'indout4';
            typIc = ICtrialtypes==106;
            typLc = ICtrialtypes==110;
        case 4
            whichseg = 'indout3';
            typIc = ICtrialtypes==111;
            typLc = ICtrialtypes==107;
    end
    tempneu = cat(1,neuinareaagg{:}) & cat(1,ICsigagg.(whichseg))==1;
    p=signrank(RbsICblockagg(tempneu, typIc), RbsICblockagg(tempneu, typLc)); % none of these are significant
    subplot(2,2,f)
    hold all
    plot(RavgICblockagg(tempneu, typIc), RavgICblockagg(tempneu, typLc), 'o')
    xl = xlim;
    plot(xl, xl, 'r--')
    xlabel('I_C')
    ylabel('L_C')
    title(sprintf('%s p=%.4f', whichseg, p))
end

%%
SP_BKagg = cat(1,ICsigagg.ICwcfg1_presentations.SP_BK);
figure; hold all
hc = histogram(SP_BKagg(:,1));
histogram(SP_BKagg(:,2), hc.BinEdges)

figure; hold all
hc = histogram(SP_BKagg(:,4));
histogram(SP_BKagg(:,3), hc.BinEdges)

for typi = 1:4
    var(SP_BKagg(:,typi))
end

[h,p]=vartest2(SP_BKagg(:,1), SP_BKagg(:,2))
[h,p]=vartest2(SP_BKagg(:,1), SP_BKagg(:,3))
[h,p]=vartest2(SP_BKagg(:,4), SP_BKagg(:,2))
[h,p]=vartest2(SP_BKagg(:,4), SP_BKagg(:,3))

tt2p = [106 107 110 111];
ttcol = [0 0.5 0; 0.5 0.25 0; 1 0.5 0; 0 1 0];

figure; hold all
for typi = 1:4
c = cdfplot(SP_BKagg(:,typi)); c.Color=ttcol(typi,:);
end

figure; hold all
for typi = 1:4
[f,x] = ecdf(SP_BKagg(:,typi));
plot(x,f, 'Color', ttcol(typi,:))
% c.Color=ttcol(typi,:);
end
