datadir = 'S:/OpenScopeData/000248/';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

svmdesc = 'trainICRCtestRE';
preproc = 'zscore'; % '' is z-score train trials, '_zscoreall', or '_meancenter'
whichSVMkernel = 'Linear';

whichICblock = 'ICwcfg1';
probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};

Twin = 100;
ises=1; pathsvm = [datadir 'postprocessed' filesep 'SVM' filesep 'SVM_' svmdesc filesep nwbsessions{ises} filesep];
load([pathsvm, 'SVMpsth' num2str(Twin) 'ms_', svmdesc, '_Cctx_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'])
psthtli = SVMpsth.psthtli;
tempTinds = SVMpsth.psthbinTinds( round(1+size(SVMpsth.psthbinTinds,1)/2),:);
SVMpsthbintli = psthtli(tempTinds);
Nsplits = size(SVMpsth.testtrialinds,2);
traintrialtypes = SVMpsth.traintrialtypes;
REttrialtypes = [1105, 1109];
I_Ctrialtypes = [106, 111];
L_Ctrialtypes = [107, 110];

SVMtrainpsthbin = struct();
SVMtestpsthbin = struct();
SVMREtasICpsthbin = struct();
SVMREtasLCpsthbin = struct();
for iprobe = 1:numel(probes)
        whichvisarea = [probes{iprobe} 'ctx'];
SVMtrainpsthbin.(whichvisarea) = NaN(length(SVMpsthbintli), numel(traintrialtypes), Nsessions);
SVMtestpsthbin.(whichvisarea) = NaN(length(SVMpsthbintli), numel(traintrialtypes), Nsessions);
SVMREtasICpsthbin.(whichvisarea) = NaN(length(SVMpsthbintli), numel(REttrialtypes), Nsessions);
SVMREtasLCpsthbin.(whichvisarea) = NaN(length(SVMpsthbintli), numel(REttrialtypes), Nsessions);
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

        for itt = 1:numel(REttrialtypes)
            trialsoi = SVMpsth.trialorder==REttrialtypes(itt);
            temp = SVMpsth.label(:,trialsoi,:)==I_Ctrialtypes(itt);
            SVMREtasICpsthbin.(whichvisarea)(:,itt,ises) = mean(temp, [2,3]);

            temp = SVMpsth.label(:,trialsoi,:)==L_Ctrialtypes(itt);
            SVMREtasLCpsthbin.(whichvisarea)(:,itt,ises) = mean(temp, [2,3]);
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
ICasREtpeaktiming = NaN(Nsessions, numel(probes));
ICasREtpeakperf = NaN(Nsessions, numel(probes));
infpeaktiming = NaN(Nsessions, numel(probes));
infpeakperf = NaN(Nsessions, numel(probes));
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

        [mv,mi] = max(squeeze(mean(SVMREtasICpsthbin.(whichvisarea)(tempbintloi,:,ises),2)));
        if ~isnan(mv)
        ICasREtpeaktiming(ises, iprobe) = tempbintli(mi);
        end
        ICasREtpeakperf(ises, iprobe) = mv;

        tempinf = SVMREtasICpsthbin.(whichvisarea) - SVMREtasLCpsthbin.(whichvisarea);
        [mv,mi] = max(squeeze(mean(tempinf(tempbintloi,:,ises),2)));
        if ~isnan(mv)
        infpeaktiming(ises, iprobe) = tempbintli(mi);
        end
        infpeakperf(ises, iprobe) = mv;
        
    end
end

figure; hold all
plot(testpeaktiming(:,strcmp(visareas, 'V1')), testpeaktiming(:,strcmp(visareas, 'LM')), 'o')
xl = xlim; 
plot(xl,xl,'-')

figure; hold all
plot(testpeakperf(:,strcmp(visareas, 'V1')), testpeaktiming(:,strcmp(visareas, 'V1')), 'o')
plot(testpeakperf(:,strcmp(visareas, 'LM')), testpeaktiming(:,strcmp(visareas, 'LM')), 'x')


tempmat = testpeaktiming(:,probesoind);
[p,t,s]=friedman(tempmat(all(~isnan(tempmat),2),:));
figure; multcompare(s); title(sprintf('Friedman p=%.4f', p))
signrank(testpeaktiming(:,strcmp(visareas, 'V1')), testpeaktiming(:,strcmp(visareas, 'LM'))) % 50ms bins: 0.1455; 100ms bins: 0.0586
signrank(testpeaktiming(:,strcmp(visareas, 'V1')), testpeaktiming(:,strcmp(visareas, 'AL')))

tempmat = ICasREtpeaktiming(:,probesoind);
[p,t,s]=friedman(tempmat(all(~isnan(tempmat),2),:));
figure; multcompare(s); title(sprintf('Friedman p=%.4f', p))
signrank(ICasREtpeaktiming(:,strcmp(visareas, 'V1')), ICasREtpeaktiming(:,strcmp(visareas, 'LM'))) % 50ms bins: 0.9824; 100ms bins: 0.5332
signrank(ICasREtpeaktiming(:,strcmp(visareas, 'V1')), ICasREtpeaktiming(:,strcmp(visareas, 'AL')))

tempmat = infpeaktiming(:,probesoind);
[p,t,s]=friedman(tempmat(all(~isnan(tempmat),2),:));
figure; multcompare(s); title(sprintf('Friedman p=%.4f', p))
signrank(infpeaktiming(:,strcmp(visareas, 'V1')), infpeaktiming(:,strcmp(visareas, 'LM'))) % 50ms bins: 0.8477; 100ms bins: 0.2354
signrank(infpeaktiming(:,strcmp(visareas, 'V1')), infpeaktiming(:,strcmp(visareas, 'AL')))

%%
addpath(genpath('C:\Users\USER\GitHub\helperfunctions'))

C = lines(numel(probes));

probesoind = find(ismember(visareas, {'V1', 'LM', 'AL'}));
figure
for ises = 1:Nsessions
subplot(3,5,ises)
hold all
for ii = 1:numel(probesoind)
    iprobe = probesoind(ii);
    whichvisarea = [probes{iprobe} 'ctx'];
    plot(SVMpsthbintli, squeeze(mean(SVMtestpsthbin.(whichvisarea)(:,:,ises),2)), 'Color', C(ii,:))
end
plot([SVMpsthbintli(1) SVMpsthbintli(end)], (1/numel(traintrialtypes))*[1 1], 'k--')
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
    plot(SVMpsthbintli, squeeze(mean(SVMREtasICpsthbin.(whichvisarea)(:,:,ises)-SVMREtasLCpsthbin.(whichvisarea)(:,:,ises),2)), 'Color', C(ii,:))
end
plot([SVMpsthbintli(1) SVMpsthbintli(end)], [0 0], 'k--')
xlabel('Time (ms)')
ylabel('Inference Performance')
end



xl = [SVMpsthbintli(1) SVMpsthbintli(end)];
figure('Position', [100 100 1200 300])
subplot(1,3,1)
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
title(sprintf('trainiCRC psth %dms %s SVM (%s)', Twin, whichSVMkernel, preproc))

subplot(1,3,2)
hold all
for ii = 1:numel(probesoind)
    iprobe = probesoind(ii);
    whichvisarea = [probes{iprobe} 'ctx'];
    temppsth = SVMREtasICpsthbin.(whichvisarea)-SVMREtasLCpsthbin.(whichvisarea);
    valsessions = ~squeeze(all(isnan(temppsth),[1,2]));
    shadedErrorBar(SVMpsthbintli, squeeze(nanmean(mean(temppsth,2),3)), squeeze(nanstd(mean(temppsth,2),0,3))/sqrt(nnz(valsessions)), {'Color', C(ii,:), 'LineWidth', 1}, 1)
end
plot(xl, [0 0], 'k--')
    xlim([-50 150])
xlabel('Time (ms)')
ylabel('Inference Performance')

subplot(1,3,3)
hold all
for ii = 1:numel(probesoind)
    iprobe = probesoind(ii);
    whichvisarea = [probes{iprobe} 'ctx'];
    temppsth = SVMREtasICpsthbin.(whichvisarea);
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
ylabel('T_R_E as I_C')


figure('Position', [100 100 1200 300])
for ii = 1:numel(probesoind)
    subplot(1,numel(probesoind),ii)
    hold all
    iprobe = probesoind(ii);
    whichvisarea = [probes{iprobe} 'ctx'];
    yyaxis left
    temppsth = SVMtestpsthbin.(whichvisarea);
    valsessions = ~squeeze(all(isnan(temppsth),[1,2]));
    shadedErrorBar(SVMpsthbintli, squeeze(nanmean(mean(temppsth,2),3)), squeeze(nanstd(mean(temppsth,2),0,3))/sqrt(nnz(valsessions)), {'Color', 'b', 'LineWidth', 1}, 1)
    ylabel('Test Accuracy')

    yyaxis right
    temppsth = SVMREtasICpsthbin.(whichvisarea)-SVMREtasLCpsthbin.(whichvisarea);
    valsessions = ~squeeze(all(isnan(temppsth),[1,2]));
    shadedErrorBar(SVMpsthbintli, squeeze(nanmean(mean(temppsth,2),3)), squeeze(nanstd(mean(temppsth,2),0,3))/sqrt(nnz(valsessions)), {'Color', 'r', 'LineWidth', 1}, 1)
    % plot(xl, (1/numel(traintrialtypes))*[1 1], 'k--')
    ylabel('Inference Performance')

    if ii ==3
        yl = ylim;
        text(xl(2), yl(2), 'test', 'Color', 'b', 'HorizontalAlignment','right', 'VerticalAlignment','top')
        text(xl(2), yl(2)-0.1*range(yl), 'inference', 'Color', 'r', 'HorizontalAlignment','right', 'VerticalAlignment','top')
    end
    xlim([-50 150])
    xlabel('Time (ms)')
    title(sprintf('%s psth %dms %s SVMtrainICRC (%s)', visareas{iprobe}, Twin, whichSVMkernel, preproc))
end
