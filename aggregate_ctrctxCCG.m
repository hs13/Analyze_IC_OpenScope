% FOUND SESSION 4 (sub_1174569641) TO BE AN EXACT DUPLICATE OF SESSION 1 (sub_1171903433)!!!!!!
Nsessions = numel(nwbsessions);

% 'neuctrctx', 'neulocctrctx', 'visarealabels', 'CCGtli', 'ctrctxCCG', 'ctrctxCCGweight'
neuctrctxagg = [];
neulocctrctxagg = [];
visarealabelagg = [];
sesctrctxagg = [];
ctrctxCCGweightagg = cell(Nsessions,1);
for ises = 1:Nsessions
pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load([pathpp 'spiketimes.mat'], 'ICblockstartend'); stlen = ICblockstartend(2)-ICblockstartend(1)+1;
load([pathpp 'ctrctxCCG.mat'])
neuctrctxagg = cat(1, neuctrctxagg, neuctrctx);
neulocctrctxagg = cat(1, neulocctrctxagg, neulocctrctx);
visarealabelagg = cat(1, visarealabelagg, visarealabels);
sesctrctxagg = cat(1, sesctrctxagg, ises*ones(size(visarealabels)) );
ctrctxCCGweightagg{ises} = ctrctxCCGweight*1000/(stlen-26);
end

%% aggregate visresponses: ICsig, RFCI, RFCIspin, sizeCI, oriparams
ises=1;
pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load(sprintf('%svisresponses_probeC.mat', pathpp))
load(sprintf('%spostprocessed_probeC.mat', pathpp), 'neuoind')
Nneurons = numel(neuoind);

ICblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations'};
% 'Palpha','BKtt','BKttpair','BItt','BIttpair','BICREl1tt','BICREl2tt','BICREl1ttpair','BICREl2ttpair',
allfields = fieldnames(ICsig.ICwcfg1_presentations);
validfields = true(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(ICsig.ICwcfg1_presentations.(allfields{f}),1)==Nneurons;
end
ICsigfields = allfields(validfields);
% {'SP_Ind','Pmww_Ind','SP_BK','sigmcBK','Pmww_BK', ...
%     'SP_ICvsRC','Pmww_ICvsRC','SP_BICREl','sigmcBICREl1','sigmcBICREl2', ...
%     'PkwBK','PmcBK', 'PkwBI','PmcBI', 'PkwBICREl1','PmcBICREl1', 'PkwBICREl2','PmcBICREl2', ...
%     'ICencoder','RCencoder','inducerencoder','inducerresponsive', ...
%     'indenc1','indenc2','indenc3','indenc4','indenc13','indenc14','indenc23','indenc24', ...
%     'ICencoder1','RCencoder1','RCencoder2','ICencoder2','indin1','indin2','indin3','indin4', ...
%     'indout1','indout2','indout3','indout4','RElfaith1','RElfaith2', ...
%     'ICresp1','RCresp1','RCresp2','ICresp2','ICtuned1','RCtuned1','RCtuned2','ICtuned2'};
ICsigallagg = struct();
    for b = 1:numel(ICblocks)
        for f= 1:numel(ICsigfields)
            ICsigallagg.(ICblocks{b}).(ICsigfields{f}) = [];
        end
    end

allfields = fieldnames(RFCI);
validfields = true(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(RFCI.(allfields{f}),1)==Nneurons;
end
RFCIfields = allfields(validfields);
% {'Rrfclassic','Rrfinverse','RFindclassic','RFindinverse', ...
%     'Pkw_rfclassic','Pkw_rfinverse','pRrfclassic','pRrfinverse','pRFclassic','pRFinverse'};
RFCIallagg = struct();
    for f= 1:numel(RFCIfields)
        RFCIallagg.(RFCIfields{f}) = [];
    end

allfields = fieldnames(RFCIspin);
validfields = true(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(RFCIspin.(allfields{f}),1)==Nneurons;
end
RFCIspinfields = allfields(validfields);
RFCIspinallagg = struct();
    for f= 1:numel(RFCIspinfields)
        RFCIspinallagg.(RFCIspinfields{f}) = [];
    end

allfields = fieldnames(sizeCI);
validfields = true(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(sizeCI.(allfields{f}),1)==Nneurons;
end
sizeCIfields = allfields(validfields);
% sizeCIfields = {'Rsizeclassic','Rsizeinverse','sizeindclassic','sizeindinverse','Pkw_sizeclassic','Pkw_sizeinverse'};
sizeCIallagg = struct();
    for f= 1:numel(sizeCIfields)
        sizeCIallagg.(sizeCIfields{f}) = [];
    end

allfields = fieldnames(oriparams);
validfields = true(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(oriparams.(allfields{f}),1)==Nneurons;
end
oriparamsfields = allfields(validfields);
% oriparamsfields = {'Rori','prefiori','orthiori','OSI','OP','Pmww_OP'};
oriparamsallagg = struct();
    for f= 1:numel(oriparamsfields)
        oriparamsallagg.(oriparamsfields{f}) = [];
    end

allfields = fieldnames(ori4params);
validfields = true(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(ori4params.(allfields{f}),1)==Nneurons;
end
ori4paramsfields = allfields(validfields);
% oriparamsfields = {'Rori','prefiori','orthiori','OSI','OP','Pmww_OP'};
ori4paramsallagg = struct();
    for f= 1:numel(ori4paramsfields)
        ori4paramsallagg.(ori4paramsfields{f}) = [];
    end

meanFRallagg = [];
sponFRallagg = [];

for ises = 1:Nsessions
    fprintf('Session %d/%d %s\n', ises, Nsessions, nwbsessions{ises} )
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'visresponses.mat'])
    
    meanFRallagg = cat(1, meanFRallagg, meanFRall');
    sponFRallagg = cat(1, sponFRallagg, sponFRall');
    
    for b = 1:numel(ICblocks)
        for f= 1:numel(ICsigfields)
            ICsigallagg.(ICblocks{b}).(ICsigfields{f}) = cat(1, ICsigallagg.(ICblocks{b}).(ICsigfields{f}), ICsigall.(ICblocks{b}).(ICsigfields{f}) );
        end
    end
    
    % RFCIfields = fieldnames(RFCI);
    for f= 1:numel(RFCIfields)
        RFCIallagg.(RFCIfields{f}) = cat(1, RFCIallagg.(RFCIfields{f}), RFCIall.(RFCIfields{f}) );
    end
    
    for f= 1:numel(RFCIspinfields)
        RFCIspinallagg.(RFCIspinfields{f}) = cat(1, RFCIspinallagg.(RFCIspinfields{f}), RFCIspinall.(RFCIspinfields{f}) );
    end
    
    % size vector [0, 4, 8, 16, 32, 64 ]
    for f= 1:numel(sizeCIfields)
        sizeCIallagg.(sizeCIfields{f}) = cat(1, sizeCIallagg.(sizeCIfields{f}), sizeCIall.(sizeCIfields{f}) );
    end
    
    for f= 1:numel(oriparamsfields)
        oriparamsallagg.(oriparamsfields{f}) = cat(1, oriparamsallagg.(oriparamsfields{f}), oriparamsall.(oriparamsfields{f}) );
    end
    
    for f= 1:numel(ori4paramsfields)
        ori4paramsallagg.(ori4paramsfields{f}) = cat(1, ori4paramsallagg.(ori4paramsfields{f}), ori4paramsall.(ori4paramsfields{f}) );
    end
    
end

%%
% sanity check
if ~( all(RFCIallagg.RFindclassic(neuctrctxagg==1)==1) && all(RFCIallagg.Pkw_rfclassic(neuctrctxagg==1)<0.05) )
    error('check that RFCIallagg has the same ordering as ctrctxCCG')
end

Nctrctx = nnz(neuctrctxagg);
ctrctxICenc = ICsigallagg.ICwcfg1_presentations.ICencoder(neuctrctxagg==1);
ctrctxsz = cellfun(@size, ctrctxCCGweightagg, 'uniformoutput', false);
ctrctxsz = cat(1, ctrctxsz{:});

ctrV1ICenc = ctrctxICenc==1 & visarealabelagg==1;

if sum(ctrctxsz(:,1),1) ~= nnz(neuctrctxagg)
    error('check ctrctx in analyze_CCG.m')
end

CCGweightsum = NaN(Nctrctx,1);
degdiverge = NaN(Nctrctx,1); % fraction positive
degconverge = NaN(Nctrctx,1); % fraction negative
deg2stddiverge = NaN(Nctrctx,1); % fraction positive
deg2stdconverge = NaN(Nctrctx,1); % fraction negative
for ises = 1:Nsessions
    neuses = sesctrctxagg==ises;
    CCGweightsum(neuses) = nanmean(ctrctxCCGweightagg{ises},2);
    degdiverge(neuses) = sum(ctrctxCCGweightagg{ises}>0,2) ./ sum(~isnan(ctrctxCCGweightagg{ises}),2);
    degconverge(neuses) = sum(ctrctxCCGweightagg{ises}<0,2) ./ sum(~isnan(ctrctxCCGweightagg{ises}),2);
    
    stdx2 = 2*std(ctrctxCCGweightagg{ises}(:));
    deg2stddiverge(neuses) = sum(ctrctxCCGweightagg{ises}>stdx2,2) ./ sum(~isnan(ctrctxCCGweightagg{ises}),2);
    deg2stdconverge(neuses) = sum(ctrctxCCGweightagg{ises}<stdx2,2) ./ sum(~isnan(ctrctxCCGweightagg{ises}),2);
end

CCGweightsum_area = NaN(Nctrctx,7);
degdiverge_area = NaN(Nctrctx,7); % fraction positive
degconverge_area = NaN(Nctrctx,7); % fraction negative
deg2stddiverge_area = NaN(Nctrctx,7); % fraction positive
deg2stdconverge_area = NaN(Nctrctx,7); % fraction negative
for ises = 1:Nsessions
    neuses = sesctrctxagg==ises;
    for a = 1:7
    sesneuinarea = visarealabelagg(sesctrctxagg==ises)==mod(a,7);
    CCGweightsum_area(neuses,a) = nanmean(ctrctxCCGweightagg{ises}(:,sesneuinarea),2);
    degdiverge_area(neuses,a) = sum(ctrctxCCGweightagg{ises}(:,sesneuinarea)>0,2) ./ sum(~isnan(ctrctxCCGweightagg{ises}(:,sesneuinarea)),2);
    degconverge_area(neuses,a) = sum(ctrctxCCGweightagg{ises}(:,sesneuinarea)<0,2) ./ sum(~isnan(ctrctxCCGweightagg{ises}(:,sesneuinarea)),2);
    
    stdx2 = 2*std(ctrctxCCGweightagg{ises}(:));
    deg2stddiverge_area(neuses,a) = sum(ctrctxCCGweightagg{ises}(:,sesneuinarea)>stdx2,2) ./ sum(~isnan(ctrctxCCGweightagg{ises}(:,sesneuinarea)),2);
    deg2stdconverge_area(neuses,a) = sum(ctrctxCCGweightagg{ises}(:,sesneuinarea)<stdx2,2) ./ sum(~isnan(ctrctxCCGweightagg{ises}(:,sesneuinarea)),2);
    end
end

CCGweightsum_V1HVA = NaN(Nctrctx,2);
degdiverge_V1HVA = NaN(Nctrctx,2); % fraction positive
degconverge_V1HVA = NaN(Nctrctx,2); % fraction negative
deg2stddiverge_V1HVA = NaN(Nctrctx,2); % fraction positive
deg2stdconverge_V1HVA = NaN(Nctrctx,2); % fraction negative
for ises = 1:Nsessions
    neuses = sesctrctxagg==ises;
    for a = 1:2
        switch a
            case 1
    sesneuinarea = visarealabelagg(sesctrctxagg==ises)==1;
            case 2
    sesneuinarea = visarealabelagg(sesctrctxagg==ises)>1;
        end
        
    CCGweightsum_V1HVA(neuses,a) = nanmean(ctrctxCCGweightagg{ises}(:,sesneuinarea),2);
    degdiverge_V1HVA(neuses,a) = sum(ctrctxCCGweightagg{ises}(:,sesneuinarea)>0,2) ./ sum(~isnan(ctrctxCCGweightagg{ises}(:,sesneuinarea)),2);
    degconverge_V1HVA(neuses,a) = sum(ctrctxCCGweightagg{ises}(:,sesneuinarea)<0,2) ./ sum(~isnan(ctrctxCCGweightagg{ises}(:,sesneuinarea)),2);
    
    stdx2 = 2*std(ctrctxCCGweightagg{ises}(:));
    deg2stddiverge_V1HVA(neuses,a) = sum(ctrctxCCGweightagg{ises}(:,sesneuinarea)>stdx2,2) ./ sum(~isnan(ctrctxCCGweightagg{ises}(:,sesneuinarea)),2);
    deg2stdconverge_V1HVA(neuses,a) = sum(ctrctxCCGweightagg{ises}(:,sesneuinarea)<stdx2,2) ./ sum(~isnan(ctrctxCCGweightagg{ises}(:,sesneuinarea)),2);
    end
end

%%
% divergence: should decrease from V1 to HvA
[P,ANOVATAB,STATS] = kruskalwallis(degdiverge(visarealabelagg>0), visarealabelagg(visarealabelagg>0));
multcompare(STATS)

% convergence: should increase from V1 to HvA
[P,ANOVATAB,STATS] = kruskalwallis(degconverge(visarealabelagg>0), visarealabelagg(visarealabelagg>0));
multcompare(STATS)

%%
addpath(genpath('/Users/hyeyoung/Documents/CODE/helperfunctions'))
neuoi = ctrV1ICenc;

tempagg = CCGweightsum_V1HVA(:,1);

figure; hold all
histogram(tempagg(visarealabelagg==1))
histogram(tempagg(neuoi))

figure
boxchart([zeros(nnz(visarealabelagg==1),1); ones(nnz(neuoi),1)], [tempagg(visarealabelagg==1); tempagg(neuoi)])
ranksum(tempagg(neuoi), tempagg(visarealabelagg==1))
ranksum(tempagg(neuoi), tempagg(visarealabelagg==1), 'tail', 'right') % if ICenc has bigger population coupling
ranksum(tempagg(neuoi), tempagg(visarealabelagg==1), 'tail', 'left')


temp_areaagg = CCGweightsum_area;
whicharea = 2;

figure; hold all
histogram(temp_areaagg(visarealabelagg==1,whicharea))
histogram(temp_areaagg(neuoi,whicharea))

figure
boxchart([zeros(nnz(visarealabelagg==1),1); ones(nnz(neuoi),1)], [temp_areaagg(visarealabelagg==1,whicharea); temp_areaagg(neuoi,whicharea)])
ranksum(temp_areaagg(neuoi,whicharea), temp_areaagg(visarealabelagg==1,whicharea))
ranksum(temp_areaagg(neuoi,whicharea), temp_areaagg(visarealabelagg==1,whicharea), 'tail', 'right') % if ICenc has bigger population coupling
ranksum(temp_areaagg(neuoi,whicharea), temp_areaagg(visarealabelagg==1,whicharea), 'tail', 'left')

recarealabels = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'};

ises = 1;
neuses = sesctrctxagg==ises;
[sv,si] = sortrows([visarealabelagg(sesctrctxagg==ises) nanmean(ctrctxCCGweightagg{ises},2)]);
figure
imagesc(ctrctxCCGweightagg{ises}(si,si))
hold on
for a = 0:6
    visareaborder = nnz(visarealabelagg(sesctrctxagg==ises)<=a);
    plot(0.5+[0 size(ctrctxCCGweightagg{ises},2)], 0.5+visareaborder+[0 0], 'k-', 'linewidth', 1)
    plot(0.5+visareaborder+[0 0], 0.5+[0 size(ctrctxCCGweightagg{ises},2)], 'k-', 'linewidth', 1)
    disp(visareaborder)
end
colormap redblue
caxis(4*10^-5*[-1 1])
colorbar


ises = 8;
neuses = sesctrctxagg==ises;
neuoi = ctrV1ICenc(neuses);
[sv,si] = sortrows([visarealabelagg(neuses) nanmean(ctrctxCCGweightagg{ises},2)]);
figure
imagesc(ctrctxCCGweightagg{ises}(neuoi,si))
hold on
for a = 0:6
    visareaborder = nnz(visarealabelagg(sesctrctxagg==ises)<=a);
    %plot(0.5+[0 size(ctrctxCCGweightagg{ises},2)], 0.5+visareaborder+[0 0], 'k-', 'linewidth', 1)
    plot(0.5+visareaborder+[0 0], 0.5+[0 size(ctrctxCCGweightagg{ises},2)], 'k-', 'linewidth', 1)
    disp(visareaborder)
end
colormap redblue
caxis(4*10^-5*[-1 1])
colorbar

% tempmat = ctrctxCCGweightagg{ises};
% tempmat = tempmat(si,:);
% tempmat = tempmat(:,si);
% isequaln(tempmat, ctrctxCCGweightagg{ises}(si,si))
