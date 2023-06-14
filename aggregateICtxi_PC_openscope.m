datadir = '/Users/hyeyoung/Documents/DATA/OpenScopeData/000248/';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

neulocagg = [];
neulocctxagg = [];
neuctxagg = [];

PCagg = [];
PCragg = [];
PCnormagg = [];
PCzagg = [];

PC_areaagg = [];
PCr_areaagg = [];
PCnorm_areaagg = [];
PCz_areaagg = [];

for ises = 1:Nsessions
    disp(nwbsessions{ises})
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'PC_openscope.mat'])
    % 'unit_peakch', 'neuloc', 'neuctx', 'neulocctx', ...
    % 'PC', 'PCr', 'PCnorm', 'PCz', 'PC_area', 'PCr_area', 'PCnorm_area', 'PCz_area'
    
    neulocagg = cat(1, neulocagg, neuloc);
    neulocctxagg = cat(1, neulocctxagg, neulocctx);
    neuctxagg = cat(1, neuctxagg, neuctx);
    
    PCagg = cat(1, PCagg, PC);
    PCragg = cat(1, PCragg, PCr);
    PCnormagg = cat(1, PCnormagg, PCnorm);
    PCzagg = cat(1, PCzagg, PCz);
    
    PC_areaagg = cat(1, PC_areaagg, PC_area);
    PCr_areaagg = cat(1, PCr_areaagg, PCr_area);
    PCnorm_areaagg = cat(1, PCnorm_areaagg, PCnorm_area);
    PCz_areaagg = cat(1, PCz_areaagg, PCz_area);
end


recarealabels = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'};
visarealabelagg = zeros(size(neulocctxagg));
for a = 1:numel(recarealabels)
    if strcmp(recarealabels{a}, 'VISp')
        neuinarea = contains(neulocctxagg, 'VISp') & ~contains(neulocctxagg, 'VISpm');
    else
        neuinarea = contains(neulocctxagg, recarealabels{a});
    end
    visarealabelagg(neuinarea) = a;
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
ICencagg = ICsigallagg.ICwcfg1_presentations.ICencoder(neuctxagg==1);
neuoi = ICencagg==1 & visarealabelagg==1;

tempPC_areaagg = PCr_areaagg;
tempPCagg = PCragg;

whicharea = 2;
% figure; hold all
% histogram(tempPC_areaagg(visarealabelagg==1,whicharea))
% histogram(tempPC_areaagg(neuoi,whicharea))
figure
boxchart([zeros(nnz(visarealabelagg==1),1); ones(nnz(neuoi),1)], [tempPC_areaagg(visarealabelagg==1,whicharea); tempPC_areaagg(neuoi,whicharea)])
ranksum(tempPC_areaagg(neuoi,whicharea), tempPC_areaagg(visarealabelagg==1,whicharea))
ranksum(tempPC_areaagg(neuoi,whicharea), tempPC_areaagg(visarealabelagg==1,whicharea), 'tail', 'right') % if ICenc has bigger population coupling
ranksum(tempPC_areaagg(neuoi,whicharea), tempPC_areaagg(visarealabelagg==1,whicharea), 'tail', 'left')

figure
boxchart([zeros(nnz(visarealabelagg==1),1); ones(nnz(neuoi),1)], [tempPCagg(visarealabelagg==1); tempPCagg(neuoi)])
ranksum(tempPCagg(neuoi), tempPCagg(visarealabelagg==1))
ranksum(tempPCagg(neuoi), tempPCagg(visarealabelagg==1), 'tail', 'right') % if ICenc has bigger population coupling
ranksum(tempPCagg(neuoi), tempPCagg(visarealabelagg==1), 'tail', 'left')

