datadir = 'S:/OpenScopeData/000248/';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

svmdesc = 'trainICRCtestRE';
preproc = 'zscore'; % '' is z-score train trials, '_zscoreall', or '_meancenter'
whichSVMkernel = 'Linear';

traintrialtypes = [106 107 110 111];
predlabel = {'IC1', 'RC1', 'RC2', 'IC2'};

whichICblock = 'ICwcfg1';
whichprobe = 'C';
whichvisarea = [whichprobe 'ctx'];

Nsplits = 10;
%%
ICsig_Cagg = struct();
RFCI_Cagg = struct();
RFCIspin_Cagg = struct();
ori4params_Cagg = struct();
oriparams_Cagg = struct();
sizeCI_Cagg = struct();

neuctxinprobe_Cagg = cell(1,Nsessions);
sesneuctxinprobe_Cagg = cell(1,Nsessions);
SVMcodingname_Cctxagg = cell(Nsplits,Nsessions);
Nlearners_Cctxagg = zeros(Nsplits,Nsessions);
betalearners_Cctxagg = cell(Nsplits,Nsessions);

for ises = 1:Nsessions
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    pathsv = [datadir 'postprocessed' filesep 'SVM' filesep 'SVM_' svmdesc filesep];
    pathsvm = [pathsv nwbsessions{ises} filesep];

    load([pathpp 'visresponses_probe' whichprobe '.mat'])
    if ises==1
        ICsig_Cagg = ICsig;
        RFCI_Cagg = RFCI;
        RFCIspin_Cagg = RFCIspin;
        ori4params_Cagg = ori4params;
        oriparams_Cagg = oriparams;
        sizeCI_Cagg = sizeCI;
    else
        ICsig_Cagg = cat(1, ICsig_Cagg, ICsig);
        RFCI_Cagg = cat(1, RFCI_Cagg, RFCI);
        RFCIspin_Cagg = cat(1, RFCIspin_Cagg, RFCIspin);
        ori4params_Cagg = cat(1, ori4params_Cagg, ori4params);
        oriparams_Cagg = cat(1, oriparams_Cagg, oriparams);
        sizeCI_Cagg = cat(1, sizeCI_Cagg, sizeCI);
    end

    load([pathsvm, 'SVMbeta_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'])
    neuctxinprobe_Cagg{ises} = neuctxinprobe;
    sesneuctxinprobe_Cagg{ises} = ises*ones(size(neuctxinprobe)) ;
    SVMcodingname_Cctxagg(:,ises) = SVMcodingname;
    Nlearners_Cctxagg(:,ises) = Nlearners;
    betalearners_Cctxagg(:,ises) = betalearners;
end

%% test accuracy

%%
comboonevsone = nchoosek(traintrialtypes, 2);

meanonevsone = cell(1,Nsessions);
meanonevsall = cell(1,Nsessions);

for ises = 1:Nsessions
    temponevsone = strcmp(SVMcodingname_Cctxagg(:,ises), 'onevsone');
    if nnz(temponevsone)>0
        meanonevsone{ises} = mean(cat(3, betalearners_Cctxagg{temponevsone,ises}),3);
    else
        meanonevsone{ises} = NaN( nnz(neuctxinprobe_Cagg{ises}), nchoosek(length(traintrialtypes),2));
    end

    temponevsall = strcmp(SVMcodingname_Cctxagg(:,ises), 'onevsall');
    if nnz(temponevsall)>0
        meanonevsall{ises} = mean(cat(3, betalearners_Cctxagg{temponevsall,ises}),3);
    else
        meanonevsall{ises} = NaN( nnz(neuctxinprobe_Cagg{ises}), length(traintrialtypes) );
    end
end

neuctxinprobeagg = cat(1, neuctxinprobe_Cagg{:});
sesneuctxinprobeagg = cat(1, sesneuctxinprobe_Cagg{:});

ICwcfg1agg = cat(1,ICsig_Cagg.([whichICblock '_presentations'])  );
IC1encagg = cat(1, ICwcfg1agg.ICencoder1);
IC2encagg = cat(1, ICwcfg1agg.ICencoder2);
indin1agg = cat(1, ICwcfg1agg.indin1);
indin2agg = cat(1, ICwcfg1agg.indin2);
indin3agg = cat(1, ICwcfg1agg.indin3);
indin4agg = cat(1, ICwcfg1agg.indin4);
RC1encagg = cat(1, ICwcfg1agg.RCencoder1);
RC2encagg = cat(1, ICwcfg1agg.RCencoder2);

ICencagg = cat(1, ICwcfg1agg.ICencoder);
indinagg = indin1agg | indin2agg | indin3agg | indin4agg;

meanonevsoneagg = cat(1,meanonevsone{:});
meanonevsallagg = cat(1,meanonevsall{:});

absmeanonevsallagg = abs(meanonevsallagg);
absmeanonevsoneagg = abs(meanonevsoneagg);

%% one vs all and one vs one decoder weights
absopt = false;
if absopt 
veconevsall = absmeanonevsallagg;
veconevsone = absmeanonevsoneagg;
else
veconevsall = meanonevsallagg;
veconevsone = meanonevsoneagg;
end

ICenc_betaonevsall = cat(1, veconevsall(IC1encagg(neuctxinprobeagg==1),1), veconevsall(IC2encagg(neuctxinprobeagg==1),4) );
segresp_betaonevsall = cat(1, veconevsall(indin1agg(neuctxinprobeagg==1) | indin3agg(neuctxinprobeagg==1),1), ...
    veconevsall(indin2agg(neuctxinprobeagg==1) | indin4agg(neuctxinprobeagg==1),4) );
p = ranksum(ICenc_betaonevsall, segresp_betaonevsall);
figure
hold all
hc = histogram(ICenc_betaonevsall, 'Normalization', 'probability');%, 'BinEdges', hc.BinEdges)
histogram(segresp_betaonevsall, hc.BinEdges, 'Normalization', 'probability')
title(sprintf('IC-encoders vs seg. resp. p=%.4f', p))
if absopt
    xlabel('SVM abs beta one vs all trials')
else
    xlabel('SVM beta one vs all trials')
end
ylabel('probability')


ICenc_betaonevsone = cat(1, mean( veconevsone(IC1encagg(neuctxinprobeagg==1),1:2), 2), ...
    -1*mean( veconevsone(IC2encagg(neuctxinprobeagg==1),5:6), 2) );
segresp_betaonevsone = cat(1, mean( veconevsone(indin1agg(neuctxinprobeagg==1) | indin3agg(neuctxinprobeagg==1),1:2), 2), ...
    -1*mean( veconevsone(indin2agg(neuctxinprobeagg==1) | indin4agg(neuctxinprobeagg==1),5:6), 2) );
p = ranksum(ICenc_betaonevsone, segresp_betaonevsone);
figure
hold all
hc = histogram(ICenc_betaonevsone, 'Normalization', 'probability');%, 'BinEdges', hc.BinEdges)
histogram(segresp_betaonevsone, hc.BinEdges, 'Normalization', 'probability')
title(sprintf('IC-encoders vs seg. resp. p=%.4f', p))
if absopt
    xlabel('SVM abs beta one vs one trials')
else
    xlabel('SVM beta one vs one trials')
end
ylabel('probability')

%% 
histnorm = 'Count';
plothistallneurons = false;

figure('Position', [100 100 800 600])
for ii = 1:4
    switch predlabel{ii}
        case 'IC1'
            neuic = IC1encagg(neuctxinprobeagg==1);
            neusr = indin1agg(neuctxinprobeagg==1) | indin3agg(neuctxinprobeagg==1);
        case 'LC1'
            neuic = RC1encagg(neuctxinprobeagg==1);
            neusr = indin1agg(neuctxinprobeagg==1) | indin4agg(neuctxinprobeagg==1);
        case 'LC2'
            neuic = RC2encagg(neuctxinprobeagg==1);
            neusr = indin2agg(neuctxinprobeagg==1) | indin3agg(neuctxinprobeagg==1);
        case 'IC2'
            neuic = IC2encagg(neuctxinprobeagg==1);
            neusr = indin2agg(neuctxinprobeagg==1) | indin4agg(neuctxinprobeagg==1);
    end
    subplot(2,2,ii)
    hold all
    if plothistallneurons
    hc = histogram(veconevsall(:,ii), 'Normalization', histnorm);
    hc = histogram(veconevsall(neusr,ii), 'Normalization', histnorm, 'BinEdges', hc.BinEdges);
    else
    hc = histogram(veconevsall(neusr,ii), 'Normalization', histnorm);
    end
    histogram(veconevsall(neuic,ii), hc.BinEdges, 'Normalization', histnorm)
    p = ranksum(veconevsall(neuic,ii), veconevsall(neusr,ii));
    pic = ranksum(veconevsall(neuic,ii), veconevsall(:,ii));
    psr = ranksum(veconevsall(neusr,ii), veconevsall(:,ii));
    if contains(predlabel{ii}, 'IC')
    title(sprintf('IC-encoders vs seg. resp. p=%.4f\nIC-enc vs all p=%.4f\nseg resp vs all p=%.4f', p, pic, psr))
    else
    title(sprintf('LC-encoders vs seg. resp. p=%.4f\nLC-enc vs all p=%.4f\nseg resp vs all p=%.4f', p, pic, psr))
    end
    if absopt
    xlabel(sprintf('SVM abs beta %s vs all trials', predlabel{ii}))
    else
    xlabel(sprintf('SVM beta %s vs all trials', predlabel{ii}))
    end
    ylabel(histnorm)
end

C= nchoosek(predlabel,2);
figure('Position', [1300 100 1200 600])
for ii = 1:6
    if strcmp(C{ii,1}, 'IC1') && contains(C{ii,2}, 'RC')
        neuic = IC1encagg(neuctxinprobeagg==1);
        neusr = indin1agg(neuctxinprobeagg==1) | indin3agg(neuctxinprobeagg==1);
    elseif strcmp(C{ii,2}, 'IC2') && contains(C{ii,1}, 'RC')
        neuic = IC2encagg(neuctxinprobeagg==1);
        neusr = indin2agg(neuctxinprobeagg==1) | indin4agg(neuctxinprobeagg==1);
    else % IC1 vs IC2, RC1 vs RC2
        neuic = ICencagg(neuctxinprobeagg==1);
        neusr = indinagg(neuctxinprobeagg==1);
    end
    subplot(2,3,ii)
    hold all
    if plothistallneurons
    hc = histogram(veconevsone(:,ii), 'Normalization', histnorm);
    hc = histogram(veconevsone(neusr,ii), 'Normalization', histnorm, 'BinEdges', hc.BinEdges);
    else
    hc = histogram(veconevsone(neusr,ii), 'Normalization', histnorm);
    end
    histogram(veconevsone(neuic,ii), hc.BinEdges, 'Normalization', histnorm)
    p = ranksum(veconevsone(neuic,ii), veconevsone(neusr,ii));
    pic = ranksum(veconevsone(neuic,ii), veconevsone(:,ii));
    psr = ranksum(veconevsone(neusr,ii), veconevsone(:,ii));
    title(sprintf('IC-encoders vs seg. resp. p=%.4f\nIC-enc vs all p=%.4f\nseg resp vs all p=%.4f', p, pic, psr))
    if absopt
    xlabel(sprintf('SVM abs beta %s vs %s trials', C{ii,1}, C{ii,2}))
    else
    xlabel(sprintf('SVM beta %s vs %s trials', C{ii,1}, C{ii,2}))
    end
    ylabel(histnorm)
end

%%
for ii = 1:6
    p1 = ranksum(meanonevsoneagg(IC1encagg(neuctxinprobeagg==1),ii), meanonevsoneagg(:,ii));
    p2 = ranksum(meanonevsoneagg(IC2encagg(neuctxinprobeagg==1),ii), meanonevsoneagg(:,ii));
    fprintf('%d %d : %.4f %.4f\n', comboonevsone(ii,1), comboonevsone(ii,2), p1, p2)
end

for ii = 1:6
    p1 = ranksum(meanonevsoneagg(indin1agg(neuctxinprobeagg==1),ii), meanonevsoneagg(:,ii));
    p2 = ranksum(meanonevsoneagg(indin2agg(neuctxinprobeagg==1),ii), meanonevsoneagg(:,ii));
    p3 = ranksum(meanonevsoneagg(indin3agg(neuctxinprobeagg==1),ii), meanonevsoneagg(:,ii));
    p4 = ranksum(meanonevsoneagg(indin4agg(neuctxinprobeagg==1),ii), meanonevsoneagg(:,ii));
    fprintf('%d %d : %.4f %.4f %.4f %.4f\n', comboonevsone(ii,1), comboonevsone(ii,2), p1, p2, p3, p4)
end

for ii = 1:6
    tempIC1 = meanonevsoneagg(IC1encagg(neuctxinprobeagg==1),ii);
    tempindin1 = meanonevsoneagg(indin1agg(neuctxinprobeagg==1),ii);
    tempindin3 = meanonevsoneagg(indin3agg(neuctxinprobeagg==1),ii);
    p1 = ranksum(tempIC1, tempindin1);
    p3 = ranksum(tempIC1, tempindin3);

    tempIC2 = meanonevsoneagg(IC2encagg(neuctxinprobeagg==1),ii);
    tempindin2 = meanonevsoneagg(indin2agg(neuctxinprobeagg==1),ii);
    tempindin4 = meanonevsoneagg(indin4agg(neuctxinprobeagg==1),ii);
    p2 = ranksum(tempIC2, tempindin2);
    p4 = ranksum(tempIC2, tempindin4);
    fprintf('%d %d : %.4f %.4f %.4f %.4f\n', comboonevsone(ii,1), comboonevsone(ii,2), p1, p2, p3, p4)
end

figure
hold all
% h = histogram(meanonevsoneagg(:,ii));
histogram(tempIC1, 'BinEdges', h.BinEdges)
histogram(tempindin13, 'BinEdges', h.BinEdges)

for ises = 1:Nsessions
    disp(sum(meanonevsone{ises},1))
end

%% sum of decoder weights
neugroups = {'ICencoder1', 'ICencoder2', 'indin1', 'indin2', 'indin3', 'indin4', 'all'};
onevsonesumbeta = struct();
onevsallsumbeta = struct();
onevsonemeanbeta = struct();
onevsallmeanbeta = struct();
for g = 1:numel(neugroups)
    onevsonesumbeta.(neugroups{g}) = NaN(Nsessions,6);
    onevsallsumbeta.(neugroups{g}) = NaN(Nsessions,4);
    onevsonemeanbeta.(neugroups{g}) = NaN(Nsessions,6);
    onevsallmeanbeta.(neugroups{g}) = NaN(Nsessions,4);
end
for ises = 1:Nsessions
for g = 1:numel(neugroups)
    switch neugroups{g}
        case 'ICencoder1'
            tempneu = ICwcfg1agg(ises).ICencoder1(neuctxinprobe_Cagg{ises}==1);
        case 'ICencoder2'
            tempneu = ICwcfg1agg(ises).ICencoder2(neuctxinprobe_Cagg{ises}==1);
        case 'indin1'
            tempneu = ICwcfg1agg(ises).indin1(neuctxinprobe_Cagg{ises}==1);
        case 'indin2'
            tempneu = ICwcfg1agg(ises).indin2(neuctxinprobe_Cagg{ises}==1);
        case 'indin3'
            tempneu = ICwcfg1agg(ises).indin3(neuctxinprobe_Cagg{ises}==1);
        case 'indin4'
            tempneu = ICwcfg1agg(ises).indin4(neuctxinprobe_Cagg{ises}==1);
        case 'all'
            tempneu = true(nnz(neuctxinprobe_Cagg{ises}==1),1);
    end
    if length(tempneu) ~= size(meanonevsone{ises},1)
        error('check tempneu')
    end
    onevsonesumbeta.(neugroups{g})(ises,:) = sum(meanonevsone{ises}(tempneu,:),1);
    onevsallsumbeta.(neugroups{g})(ises,:) = sum(meanonevsall{ises}(tempneu,:),1);
    onevsonemeanbeta.(neugroups{g})(ises,:) = mean(meanonevsone{ises}(tempneu,:),1);
    onevsallmeanbeta.(neugroups{g})(ises,:) = mean(meanonevsall{ises}(tempneu,:),1);
end
end

pmat = NaN(6, numel(neugroups)-1);
pmat0 = NaN(6, numel(neugroups)-1);
for g = 1:numel(neugroups)-1
    for ii = 1:6
    pmat(ii,g) = signrank(onevsonemeanbeta.(neugroups{g})(:,ii), onevsonemeanbeta.all(:,ii));
    pmat0(ii,g) = signrank(onevsonemeanbeta.(neugroups{g})(:,ii));
    end
end
isequal(pmat, pmat0) % not true
isequal(pmat<0.05, pmat0<0.05) % true

signrank(onevsallmeanbeta.ICencoder1(:,1), onevsallmeanbeta.indin1(:,1))
signrank(onevsallmeanbeta.ICencoder1(:,1), onevsallmeanbeta.indin3(:,1))

signrank(onevsallmeanbeta.ICencoder1(:,2), onevsallmeanbeta.indin1(:,2))
signrank(onevsallmeanbeta.ICencoder1(:,2), onevsallmeanbeta.indin3(:,2))

signrank(onevsallmeanbeta.ICencoder1(:,1), onevsallmeanbeta.indin3(:,1))

signrank(onevsallmeanbeta.ICencoder1(:,1), onevsallmeanbeta.indin3(:,1))
