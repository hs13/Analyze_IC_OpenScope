datadir = 'S:/OpenScopeData/000248/';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

svmdesc = 'trainRExtestICRC';
preproc = 'zscore'; % '' is z-score train trials, '_zscoreall', or '_meancenter'
whichSVMkernel = 'Linear';

whichICblock = 'ICwcfg1';
probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};

Twin = 50;
ises=1; pathsvm = [datadir 'postprocessed' filesep 'SVM' filesep 'SVM_' svmdesc filesep nwbsessions{ises} filesep];
load([pathsvm, 'SVMpsth' num2str(Twin) 'ms_', svmdesc, '_Cctx_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'])
psthtli = SVMpsth.psthtli;
tempTinds = SVMpsth.psthbinTinds( round(1+size(SVMpsth.psthbinTinds,1)/2),:);
SVMpsthbintli = psthtli(tempTinds);
Nsplits = size(SVMpsth.testtrialinds,2);
traintrialtypes = SVMpsth.traintrialtypes;
% RExtrialtypes = [1201 1299];
I_Ctrialtypes = [106, 111];
L_Ctrialtypes = [107, 110];

SVMtrainpsthbin = struct();
SVMtestpsthbin = struct();
SVMICasRExpsthbin = struct();
for iprobe = 1:numel(probes)
        whichvisarea = [probes{iprobe} 'ctx'];
SVMtrainpsthbin.(whichvisarea) = NaN(length(SVMpsthbintli), numel(traintrialtypes), Nsessions);
SVMtestpsthbin.(whichvisarea) = NaN(length(SVMpsthbintli), numel(traintrialtypes), Nsessions);
SVMICasRExpsthbin.(whichvisarea) = NaN(length(SVMpsthbintli), numel(I_Ctrialtypes), Nsessions);
end
tic
for ises = 1:Nsessions
    clearvars SVMpsth
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    pathsv = [datadir 'postprocessed' filesep 'SVM' filesep 'SVM_' svmdesc filesep];
    pathsvm = [pathsv nwbsessions{ises} filesep];
    for iprobe = 1:numel(probes)
        whichvisarea = [probes{iprobe} 'ctx'];
        svmpsthfn = [pathsvm, 'SVMpsth' num2str(Twin) 'ms_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'];
        if ~exist(svmpsthfn, 'file')
            fprintf('%d %s %s SVMpsth%dms does not exist\n', ises, nwbsessions{ises}, whichvisarea, Twin )
            continue
        end
        load(svmpsthfn)

        temptrainacc = cell(size(traintrialtypes));
        temptestacc = cell(size(traintrialtypes));
        for isplit = 1:Nsplits
            temptesttrials = SVMpsth.testtrialinds(:,isplit);
            temptraintrials = SVMpsth.traintrialinds(:,isplit);
            for itt = 1:numel(traintrialtypes)
                trialsoi = temptesttrials(SVMpsth.trialorder(temptesttrials)==traintrialtypes(itt));
                temp = squeeze(SVMpsth.label(:,trialsoi,isplit)==traintrialtypes(itt));
                temptestacc{itt} = cat(2, temptestacc{itt}, temp);

                trialsoi = temptraintrials(SVMpsth.trialorder(temptraintrials)==traintrialtypes(itt));
                temp = squeeze(SVMpsth.label(:,trialsoi,isplit)==traintrialtypes(itt));
                temptrainacc{itt} = cat(2, temptrainacc{itt}, temp);
            end

        end
        for itt = 1:numel(traintrialtypes)
            SVMtestpsthbin.(whichvisarea)(:,itt,ises) = mean(temptestacc{itt},2);
            SVMtrainpsthbin.(whichvisarea)(:,itt,ises) = mean(temptrainacc{itt},2);
        end

        for itt = 1:numel(I_Ctrialtypes)
            trialsoi = SVMpsth.trialorder==I_Ctrialtypes(itt);
            temp = SVMpsth.label(:,trialsoi,:)==traintrialtypes(itt);
            SVMICasRExpsthbin.(whichvisarea)(:,itt,ises) = mean(temp, [2,3]);
        end
    end
end
toc

%% peak timing statistics between areas
visareasoi = {'V1', 'LM', 'AL'};
probesoind = zeros(size(visareasoi));
for a= 1:numel(visareasoi)
    probesoind(a) = find(strcmp(visareas, visareasoi{a} ));
end

testpeaktiming = NaN(Nsessions, numel(probes));
testpeakperf = NaN(Nsessions, numel(probes));
ICasRExpeaktiming = NaN(Nsessions, numel(probes));
ICasRExpeakperf = NaN(Nsessions, numel(probes));
for ises = 1:Nsessions
    for iprobe = 1:numel(probes)
        whichvisarea = [probes{iprobe} 'ctx'];
        tempbintloi = SVMpsthbintli>=0 & SVMpsthbintli<400;
        tempbintli = SVMpsthbintli(tempbintloi);

        [mv,mi] = max(squeeze(mean(SVMtestpsthbin.(whichvisarea)(tempbintloi,:,ises),2)));
        if ~isnan(mv)
        testpeaktiming(ises, iprobe) = tempbintli(mi);
        end
        testpeakperf(ises, iprobe) = mv;

        [mv,mi] = max(squeeze(mean(SVMICasRExpsthbin.(whichvisarea)(tempbintloi,:,ises),2)));
        if ~isnan(mv)
        ICasRExpeaktiming(ises, iprobe) = tempbintli(mi);
        end
        ICasRExpeakperf(ises, iprobe) = mv;
        
    end
end

tempmat = testpeaktiming(:,probesoind);
[p,t,s]=friedman(tempmat(all(~isnan(tempmat),2),:));
figure; multcompare(s); title(sprintf('Friedman p=%.4f', p))
signrank(testpeaktiming(:,strcmp(visareas, 'V1')), testpeaktiming(:,strcmp(visareas, 'LM')))
signrank(testpeaktiming(:,strcmp(visareas, 'V1')), testpeaktiming(:,strcmp(visareas, 'AL')))

tempmat = ICasRExpeaktiming(:,probesoind);
[p,t,s]=friedman(tempmat(all(~isnan(tempmat),2),:));
figure; multcompare(s); title(sprintf('Friedman p=%.4f', p))
signrank(ICasRExpeaktiming(:,strcmp(visareas, 'V1')), ICasRExpeaktiming(:,strcmp(visareas, 'LM')))
signrank(ICasRExpeaktiming(:,strcmp(visareas, 'V1')), ICasRExpeaktiming(:,strcmp(visareas, 'AL')))

%%
addpath(genpath('C:\Users\USER\GitHub\helperfunctions'))

C = lines(numel(probes));
xl = [SVMpsthbintli(1) SVMpsthbintli(end)];

figure
for ises = 1:Nsessions
subplot(3,5,ises)
hold all
for ii = 1:numel(probesoind)
    iprobe = probesoind(ii);
    whichvisarea = [probes{iprobe} 'ctx'];
    plot(SVMpsthbintli, squeeze(mean(SVMtestpsthbin.(whichvisarea)(:,:,ises),2)), 'Color', C(ii,:))
end
plot(xl, (1/numel(traintrialtypes))*[1 1], 'k--')
xlabel('Time (ms)')
ylabel('Test Accuracy')
end

figure
for ises = 1:Nsessions
subplot(3,5,ises)
hold all
for ii = 1:numel(probesoind)
    iprobe = probesoind(ii);
    whichvisarea = [probes{iprobe} 'ctx'];
    plot(SVMpsthbintli, squeeze(mean(SVMICasRExpsthbin.(whichvisarea)(:,:,ises),2)), 'Color', C(ii,:))
end
plot(xl, [0 0], 'k--')
xlabel('Time (ms)')
ylabel('Inference Performance')
end

figure
hold all
for ii = 1:numel(probesoind)
    iprobe = probesoind(ii);
    whichvisarea = [probes{iprobe} 'ctx'];
    temppsth = SVMtrainpsthbin.(whichvisarea);
    valsessions = ~squeeze(all(isnan(temppsth),[1,2]));
    shadedErrorBar(SVMpsthbintli, squeeze(nanmean(mean(temppsth,2),3)), squeeze(nanstd(mean(temppsth,2),0,3))/sqrt(nnz(valsessions)), {'Color', C(ii,:), 'LineWidth', 1}, 1)
end
plot(xl, (1/numel(traintrialtypes))*[1 1], 'k--')
xlabel('Time (ms)')
ylabel('Train Accuracy')


figure('Position', [100 100 800 300])
subplot(1,2,1)
hold all
for ii = 1:numel(probesoind)
    iprobe = probesoind(ii);
    whichvisarea = [probes{iprobe} 'ctx'];
    temppsth = SVMtestpsthbin.(whichvisarea);
    valsessions = ~squeeze(all(isnan(temppsth),[1,2]));
    shadedErrorBar(SVMpsthbintli, squeeze(nanmean(mean(temppsth,2),3)), squeeze(nanstd(mean(temppsth,2),0,3))/sqrt(nnz(valsessions)), {'Color', C(ii,:), 'LineWidth', 1}, 1)
end
plot(xl, (1/numel(traintrialtypes))*[1 1], 'k--')
    xlim([-50 150])
xlabel('Time (ms)')
ylabel('Test Accuracy')
title(sprintf('trainREx psth %dms %s SVM (%s)', Twin, whichSVMkernel, preproc))

subplot(1,2,2)
hold all
for ii = 1:numel(probesoind)
    iprobe = probesoind(ii);
    whichvisarea = [probes{iprobe} 'ctx'];
    temppsth = SVMICasRExpsthbin.(whichvisarea);
    valsessions = ~squeeze(all(isnan(temppsth),[1,2]));
    shadedErrorBar(SVMpsthbintli, squeeze(nanmean(mean(temppsth,2),3)), squeeze(nanstd(mean(temppsth,2),0,3))/sqrt(nnz(valsessions)), {'Color', C(ii,:), 'LineWidth', 1}, 1)
    if ii==1
        yl = ylim;
    end
text(xl(2), yl(2)-(ii-1)*0.1*range(yl), visareas{iprobe}, 'Color', C(ii,:), 'HorizontalAlignment','right', 'VerticalAlignment','top')
end
plot(xl, (1/numel(traintrialtypes))*[1 1], 'k--')
    xlim([-50 150])
xlabel('Time (ms)')
ylabel('Inference Performance')


figure('Position', [100 100 1200 300])
for ii = 1:numel(probesoind)
    subplot(1,numel(probesoind),ii)
    hold all
    iprobe = probesoind(ii);
    whichvisarea = [probes{iprobe} 'ctx'];
    temppsth = SVMtestpsthbin.(whichvisarea);
    valsessions = ~squeeze(all(isnan(temppsth),[1,2]));
    shadedErrorBar(SVMpsthbintli, squeeze(nanmean(mean(temppsth,2),3)), squeeze(nanstd(mean(temppsth,2),0,3))/sqrt(nnz(valsessions)), {'Color', 'b', 'LineWidth', 1}, 1)

    temppsth = SVMICasRExpsthbin.(whichvisarea);
    valsessions = ~squeeze(all(isnan(temppsth),[1,2]));
    shadedErrorBar(SVMpsthbintli, squeeze(nanmean(mean(temppsth,2),3)), squeeze(nanstd(mean(temppsth,2),0,3))/sqrt(nnz(valsessions)), {'Color', 'r', 'LineWidth', 1}, 1)
    plot(xl, (1/numel(traintrialtypes))*[1 1], 'k--')

    if ii ==3
        yl = ylim;
        text(xl(2), yl(2), 'test', 'Color', 'b', 'HorizontalAlignment','right', 'VerticalAlignment','top')
        text(xl(2), yl(2)-0.1*range(yl), 'inference', 'Color', 'r', 'HorizontalAlignment','right', 'VerticalAlignment','top')
    end
    xlim([-50 150])
    xlabel('Time (ms)')
    ylabel('Test Accuracy')
    title(sprintf('%s psth %dms %s SVMtrainREx (%s)', visareas{iprobe}, Twin, whichSVMkernel, preproc))
end
