addpath(genpath('C:\Users\USER\GitHub\helperfunctions'))

% FOUND SESSION 4 (sub_1174569641) TO BE AN EXACT DUPLICATE OF SESSION 1 (sub_1171903433)!!!!!!
nwbsessions = {'sub_1171903433','sub_1172968426','sub_1172969394', 'sub_1174569641', ...
    'sub_1175512783','sub_1176214862','sub_1177693342','sub_1181314060', ...
    'sub_1181585608','sub_1182593231','sub_1183369803','sub_1186544726', ...
    'sub_1189891322','sub_1194090570'};
Nsessions = numel(nwbsessions);
ses2anal = 1:Nsessions;
ses2anal(4) = [];

% ccgpath = '/Users/hyeyoung/Library/CloudStorage/GoogleDrive-hyeyoung_shin@berkeley.edu/My Drive/DATA/ICexpts_postprocessed_OpenScope/CCG/';
ccgpath = 'S:\OpenScopeData\000248\CCG\'

%{
for ises = ses2anal
    fprintf('%d/%d %s\n', ises, Nsessions, nwbsessions{ises})
    clearvars -except ises ses2anal Nsessions nwbsessions ccgpath
% pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
pathccg = [ccgpath nwbsessions{ises} filesep];
load([pathccg 'ctxCCG.mat'])
load([pathccg 'ctxCCGsm.mat'])
ctxCCGsm = double(ctxCCGsm);
ctxCCGjc = ctxCCG - ctxCCGsm;

tic
[excpeak, inhtrough] = computeCCGxconn(CCGtli, ctxCCGjc, 10);
toc
save([pathccg 'ctxCCGjcxconn.mat'], 'excpeak', 'inhtrough')
end
%}
%%
% 'neuctrctx', 'neulocctrctx', 'visarealabels', 'CCGtli', 'ctrctxCCG', 'ctrctxCCGweight'
neuctxagg = [];
neulocctxagg = [];
Nneurons_visarea = NaN(Nsessions,6);
visarealabelagg = [];
sesctxagg = [];
ctxCCGweightagg = cell(Nsessions,1);
ctxCCGjcweightagg = cell(Nsessions,1);
ctxCCGjcinagg = cell(Nsessions,1);
ctxCCGjcoutagg = cell(Nsessions,1);

ctxCCGsmt0agg = cell(Nsessions,1);
ctxCCGjcpostpkagg = cell(Nsessions,1);
ctxCCGjcpostTpkagg = cell(Nsessions,1);
ctxCCGjcprepkagg = cell(Nsessions,1);
ctxCCGjcpreTpkagg = cell(Nsessions,1);

ctxCCGsm0agg = cell(Nsessions,1);
tic
for ises = ses2anal
    fprintf('%d/%d %s\n', ises, Nsessions, nwbsessions{ises})
    clearvars ctxCCGsm ctxCCG ctxCCGweight visarealabels neuctx neulocctx
% pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
pathccg = [ccgpath nwbsessions{ises} filesep];
load([pathccg 'ctxCCG.mat'])
[v,c]=uniquecnt(visarealabels);
Nneurons_visarea(ises,ismember(1:6, v)) = c(ismember(v,1:6));
neuctxagg = cat(1, neuctxagg, neuctx);
neulocctxagg = cat(1, neulocctxagg, neulocctx);
visarealabelagg = cat(1, visarealabelagg, visarealabels);
sesctxagg = cat(1, sesctxagg, ises*ones(size(visarealabels)) );
ctxCCGweightagg{ises} = ctxCCGweight;

load([pathccg 'ctxCCGsm.mat'])
ctxCCGsm = double(ctxCCGsm);
ctxCCGjc = ctxCCG - ctxCCGsm;
ctxCCGjcweight = squeeze( sum(ctxCCGjc(:,:,CCGtli>0),3) - sum(ctxCCGjc(:,:,CCGtli<0),3) );
ctxCCGjcin = squeeze( mean(ctxCCGjc(:,:,CCGtli<0),3) );
ctxCCGjcout = squeeze( mean(ctxCCGjc(:,:,CCGtli>0),3) );
ctxCCGjcweightagg{ises} = ctxCCGjcweight;
ctxCCGjcinagg{ises} = ctxCCGjcin;
ctxCCGjcoutagg{ises} = ctxCCGjcout;
ctxCCGsmt0agg{ises} = squeeze(ctxCCGsm(:,:,CCGtli==0));

% posttli = CCGtli(CCGtli>=0 & CCGtli<=10);
% [postmv,postmi]= max(ctxCCGjc(:,:,CCGtli>=0 & CCGtli<=10), [],3);
% postmv(posttli(postmi)==0)=NaN;
% ctxCCGjcpostpkagg{ises} = postmv;
% ctxCCGjcpostTpkagg{ises} = posttli(postmi);
% 
% pretli = CCGtli(CCGtli<=0 & CCGtli>=-10);
% [premv,premi]= max(ctxCCGjc(:,:,CCGtli<=0 & CCGtli>=-10), [],3);
% premv(pretli(premi)==0)=NaN;
% ctxCCGjcprepkagg{ises} = premv;
% ctxCCGjcpreTpkagg{ises} = pretli(premi);

load([pathccg 'ctxCCGjcxconn.mat'])
ctxCCGjcpostpkagg{ises} = squeeze(excpeak(:,:,1));
ctxCCGjcpostTpkagg{ises} = squeeze(excpeak(:,:,2));

ctxCCGjcprepkagg{ises} = squeeze(excpeak(:,:,1))';
ctxCCGjcpreTpkagg{ises} = squeeze(excpeak(:,:,2))';

load([pathccg 'ctxCCGsm0.mat']) % good match with squeeze(ctxCCGsm(:,:,CCGtli==0))
ctxCCGsm0agg{ises} = ctxCCGsm0;

toc
end

%% report number of neurons
fprintf('V1 neurons mean %.4f sem %.4f\n', nanmean(Nneurons_visarea(:,1)), ...
    nanstd(Nneurons_visarea(:,1))/sqrt(size(Nneurons_visarea,1)) )

fprintf('HVA neurons mean %.4f sem %.4f\n', nanmean(sum(Nneurons_visarea(:,2:6),2)), ...
nanstd(sum(Nneurons_visarea(:,2:6),2))/sqrt(size(Nneurons_visarea,1)) )

%% sanity check
%{
mean(postmv==excpk, 'all')
%nnz(isnan(postmv) & ~isnan(excpk))

posttli = CCGtli(CCGtli>0 & CCGtli<=10);
[postmv,postmi]= max(ctxCCGjc(:,:,CCGtli>0 & CCGtli<=10), [],3);

pretli = CCGtli(CCGtli<0 & CCGtli>=-10);
[premv,premi]= max(ctxCCGjc(:,:,CCGtli<0 & CCGtli>=-10), [],3);

temppost = postmv;
temppost(eye(size(temppost))==1)= NaN;
temppre = premv;
temppre(eye(size(temppre))==1)= NaN;

isequaln(temppost, temppre')

mv = max(abs((temppost-temppre')./temppost),[],'all');
[r,c]=find(abs((temppost-temppre')./temppost)==mv);

figure; hold all
plot(CCGtli, squeeze(ctxCCGjc(c,r,:)), 'k-')
plot(CCGtli, flip(squeeze(ctxCCGjc(r,c,:))), 'r--')
plot(posttli(postmi(c,r)), postmv(c,r), 'ko')
plot(pretli(premi(c,r)), premv(c,r), 'ro')
%}

%% aggregate visresponses: ICsig, RFCI, RFCIspin, sizeCI, oriparams
ises=14;
datadir = 'S:\OpenScopeData\000248\';
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

for ises = ses2anal
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
% unique(sesctxagg(degdiverge<0.4 & degconverge<0.4))

ctxsz = cellfun(@size, ctxCCGjcweightagg, 'uniformoutput', false);
ctxsz = cat(1, ctxsz{:});
if sum(ctxsz(:,1),1) ~= nnz(neuctxagg)
    error('check neuctx in analyze_CCG_savio.m')
end
if length(ICsigallagg.ICwcfg1_presentations.ICencoder) ~= length(neuctxagg)
    error('check ctx in analyze_CCG_savio.m')
end

Nctrctx = nnz(neuctxagg);

% neuctx2incl = true( nnz(neuctxagg==1), 1);
neuctx2incl = ICsigallagg.ICwcfg1_presentations.PkwBK(neuctxagg==1)<0.05;
% neuctx2incl = ICsigallagg.ICwcfg1_presentations.ICresp1(neuctxagg==1)==1;


CCGjcweightsum = NaN(Nctrctx,1);
CCGjcinsum = NaN(Nctrctx,1);
CCGjcoutsum = NaN(Nctrctx,1);
degdiverge = NaN(Nctrctx,1); % fraction positive
degconverge = NaN(Nctrctx,1); % fraction negative
deg2stddiverge = NaN(Nctrctx,1); % fraction positive
deg2stdconverge = NaN(Nctrctx,1); % fraction negative
degin = NaN(Nctrctx,1);
degout = NaN(Nctrctx,1);
for ises = 1:Nsessions
    neuses = sesctxagg==ises;
noteye = ~eye(nnz(neuses));
    sesneuctx2incl = neuctx2incl(sesctxagg==ises);
    CCGjcweightsum(neuses) = nanmean(ctxCCGjcweightagg{ises}(:,sesneuctx2incl),2);
    CCGjcinsum(neuses) = nanmean(ctxCCGjcinagg{ises}(:,sesneuctx2incl),2);
    CCGjcoutsum(neuses) = nanmean(ctxCCGjcoutagg{ises}(:,sesneuctx2incl),2);
    degdiverge(neuses) = sum(ctxCCGjcweightagg{ises}(:,sesneuctx2incl)>0,2) ./ sum(noteye(:,sesneuctx2incl),2);
    degconverge(neuses) = sum(ctxCCGjcweightagg{ises}(:,sesneuctx2incl)<0,2) ./ sum(noteye(:,sesneuctx2incl),2);
    
    degin(neuses) = sum(ctxCCGjcinagg{ises}(:,sesneuctx2incl)>0,2) ./ sum(noteye(:,sesneuctx2incl),2);
    degout(neuses) = sum(ctxCCGjcoutagg{ises}(:,sesneuctx2incl)>0,2) ./ sum(noteye(:,sesneuctx2incl),2);
    
    stdx2 = 2*std( reshape(ctxCCGjcweightagg{ises}(sesneuctx2incl, sesneuctx2incl),[],1) );
    deg2stddiverge(neuses) = sum(ctxCCGjcweightagg{ises}(:,sesneuctx2incl)>stdx2,2) ./ sum(noteye(:,sesneuctx2incl),2);
    deg2stdconverge(neuses) = sum(ctxCCGjcweightagg{ises}(:,sesneuctx2incl)<stdx2,2) ./ sum(noteye(:,sesneuctx2incl),2);
end

CCGjcweightsum_area = NaN(Nctrctx,7);
CCGjcinsum_area = NaN(Nctrctx,7);
CCGjcoutsum_area = NaN(Nctrctx,7);
degdiverge_area = NaN(Nctrctx,7); % fraction positive
degconverge_area = NaN(Nctrctx,7); % fraction negative
deg2stddiverge_area = NaN(Nctrctx,7); % fraction positive
deg2stdconverge_area = NaN(Nctrctx,7); % fraction negative
degin_area = NaN(Nctrctx,7);
degout_area = NaN(Nctrctx,7);
for ises = 1:Nsessions
    neuses = sesctxagg==ises;
noteye = ~eye(nnz(neuses));
    for a = 1:7
    %sesneuinarea = visarealabelagg(sesctxagg==ises)==mod(a,7) ;
    sesneuinarea = visarealabelagg(sesctxagg==ises)==mod(a,7) & neuctx2incl(sesctxagg==ises);
    CCGjcweightsum_area(neuses,a) = nanmean(ctxCCGjcweightagg{ises}(:,sesneuinarea),2);
    CCGjcinsum_area(neuses,a) = nanmean(ctxCCGjcinagg{ises}(:,sesneuinarea),2);
    CCGjcoutsum_area(neuses,a) = nanmean(ctxCCGjcoutagg{ises}(:,sesneuinarea),2);
    degdiverge_area(neuses,a) = sum(ctxCCGjcweightagg{ises}(:,sesneuinarea)>0,2) ./ sum(noteye(:,sesneuinarea),2);
    degconverge_area(neuses,a) = sum(ctxCCGjcweightagg{ises}(:,sesneuinarea)<0,2) ./ sum(noteye(:,sesneuinarea),2);
    
    degin_area(neuses,a) = sum(ctxCCGjcinagg{ises}(:,sesneuinarea)>0,2) ./ sum(noteye(:,sesneuinarea),2);
    degout_area(neuses,a) = sum(ctxCCGjcoutagg{ises}(:,sesneuinarea)>0,2) ./ sum(noteye(:,sesneuinarea),2);
    
    stdx2 = 2*std(ctxCCGjcweightagg{ises}(:));
    deg2stddiverge_area(neuses,a) = sum(ctxCCGjcweightagg{ises}(:,sesneuinarea)>stdx2,2) ./ sum(noteye(:,sesneuinarea),2);
    deg2stdconverge_area(neuses,a) = sum(ctxCCGjcweightagg{ises}(:,sesneuinarea)<stdx2,2) ./ sum(noteye(:,sesneuinarea),2);
    end
end

CCGjcweightsum_V1HVA = NaN(Nctrctx,2);
CCGjcinsum_V1HVA = NaN(Nctrctx,2);
CCGjcoutsum_V1HVA = NaN(Nctrctx,2);
degdiverge_V1HVA = NaN(Nctrctx,2); % fraction positive
degconverge_V1HVA = NaN(Nctrctx,2); % fraction negative
deg2stddiverge_V1HVA = NaN(Nctrctx,2); % fraction positive
deg2stdconverge_V1HVA = NaN(Nctrctx,2); % fraction negative
degin_V1HVA = NaN(Nctrctx,2);
degout_V1HVA = NaN(Nctrctx,2);
for ises = 1:Nsessions
    neuses = sesctxagg==ises;
noteye = ~eye(nnz(neuses));
    for a = 1:2
        switch a
            case 1
    sesneuinarea = visarealabelagg(sesctxagg==ises)==1 & neuctx2incl(sesctxagg==ises);
            case 2
    sesneuinarea = visarealabelagg(sesctxagg==ises)>1 & neuctx2incl(sesctxagg==ises);
        end
        
    CCGjcweightsum_V1HVA(neuses,a) = nanmean(ctxCCGjcweightagg{ises}(:,sesneuinarea),2);
    CCGjcinsum_V1HVA(neuses,a) = nanmean(ctxCCGjcinagg{ises}(:,sesneuinarea),2);
    CCGjcoutsum_V1HVA(neuses,a) = nanmean(ctxCCGjcoutagg{ises}(:,sesneuinarea),2);
    
    degdiverge_V1HVA(neuses,a) = sum(ctxCCGjcweightagg{ises}(:,sesneuinarea)>0,2) ./ sum(noteye(:,sesneuinarea),2);
    degconverge_V1HVA(neuses,a) = sum(ctxCCGjcweightagg{ises}(:,sesneuinarea)<0,2) ./ sum(noteye(:,sesneuinarea),2);
    
    degin_V1HVA(neuses,a) = sum(ctxCCGjcinagg{ises}(:,sesneuinarea)>0,2) ./ sum(noteye(:,sesneuinarea),2);
    degout_V1HVA(neuses,a) = sum(ctxCCGjcoutagg{ises}(:,sesneuinarea)>0,2) ./ sum(noteye(:,sesneuinarea),2);
    
    stdx2 = 2*std(ctxCCGjcweightagg{ises}(:));
    deg2stddiverge_V1HVA(neuses,a) = sum(ctxCCGjcweightagg{ises}(:,sesneuinarea)>stdx2,2) ./ sum(noteye(:,sesneuinarea),2);
    deg2stdconverge_V1HVA(neuses,a) = sum(ctxCCGjcweightagg{ises}(:,sesneuinarea)<stdx2,2) ./ sum(noteye(:,sesneuinarea),2);
    end
end

%%
% neuctx2incl = true( nnz(neuctxagg==1), 1);
neuctx2incl = ICsigallagg.ICwcfg1_presentations.PkwBK(neuctxagg==1)<0.05;

areas = {'all', 'V1', 'HVA', 'LM', 'sigkwBK', 'V1sigkwBK', 'HVAsigkwBK', 'LMsigkwBK'};
xsm0thr = 0:0.2:8;
Pout = struct();
Pin = struct();
for a = 1:numel(areas)
    Pout.(areas{a}) = NaN(Nctrctx,length(xsm0thr));
    Pin.(areas{a}) = NaN(Nctrctx,length(xsm0thr));
end
for ises = ses2anal
    neuses = sesctxagg==ises;
    noteye = ~eye(nnz(neuses));
    for a = 1:numel(areas)
        switch areas{a}
            case 'all'
                sesneuinarea = true( nnz(sesctxagg==ises), 1);
            case 'V1'
                sesneuinarea = visarealabelagg(sesctxagg==ises)==1;
            case 'HVA'
                sesneuinarea = visarealabelagg(sesctxagg==ises)>1;
            case 'LM'
                sesneuinarea = visarealabelagg(sesctxagg==ises)==2;
            case 'sigkwBK'
                sesneuinarea = neuctx2incl(sesctxagg==ises);
            case 'V1sigkwBK'
                sesneuinarea = visarealabelagg(sesctxagg==ises)==1 & neuctx2incl(sesctxagg==ises);
            case 'HVAsigkwBK'
                sesneuinarea = visarealabelagg(sesctxagg==ises)>1 & neuctx2incl(sesctxagg==ises);
            case 'LMsigkwBK'
                sesneuinarea = visarealabelagg(sesctxagg==ises)==2 & neuctx2incl(sesctxagg==ises);
        end
        for ix = 1:length(xsm0thr)
            tempCCG = ctxCCGjcpostpkagg{ises} > xsm0thr(ix)*ctxCCGsm0agg{ises};
            Pout.(areas{a})(neuses,ix) = sum(tempCCG(:,sesneuinarea),2)./sum(noteye(:,sesneuinarea),2) ;
            
            tempCCG = ctxCCGjcprepkagg{ises} > xsm0thr(ix)*ctxCCGsm0agg{ises};
            Pin.(areas{a})(neuses,ix) = sum(tempCCG(:,sesneuinarea),2)./sum(noteye(:,sesneuinarea),2) ;
        end
    end
end

%% firing rate confound...
figure
subplot(2,2,1)
plot(meanFRallagg(neuctxagg==1), Pin.HVA(:,xsm0thr==2), '.')
xlabel('Firing Rate')
ylabel('Input from HVA')

subplot(2,2,2)
plot(meanFRallagg(neuctxagg==1), Pin.V1(:,xsm0thr==2), '.')
xlabel('Firing Rate')
ylabel('Input from V1')

subplot(2,2,3)
histogram2(meanFRallagg(neuctxagg==1), Pin.HVA(:,xsm0thr==2), 'DisplayStyle', 'tile')
xlabel('Firing Rate')
ylabel('Input from HVA')

subplot(2,2,4)
histogram2(meanFRallagg(neuctxagg==1), Pin.V1(:,xsm0thr==2), 'DisplayStyle', 'tile')
xlabel('Firing Rate')
ylabel('Input from V1')

%%
fs = 18;

ICenc = ICsigallagg.ICwcfg1_presentations.ICencoder==1;
indin = ICsigallagg.ICwcfg1_presentations.indin1==1 | ICsigallagg.ICwcfg1_presentations.indin2==1 ...
    | ICsigallagg.ICwcfg1_presentations.indin3==1 | ICsigallagg.ICwcfg1_presentations.indin4==1;
% ICenc = ICsigallagg.ICwcfg1_presentations.ICencoder1==1;
% indin = ICsigallagg.ICwcfg1_presentations.indin1==1 | ICsigallagg.ICwcfg1_presentations.indin3==1;

neuoi = ICenc(neuctxagg==1) & visarealabelagg==1;
neuctrl = indin(neuctxagg==1) & visarealabelagg==1;
% neuctrl = neuctx2incl & visarealabelagg==1;

V1neugroups = {'ICenc', 'indin'};
Pout_med = struct();
Pin_med = struct();
for n = 1:2
    switch V1neugroups{n}
        case 'ICenc'
            tempneu = ICenc(neuctxagg==1) & visarealabelagg==1;
        case 'indin'
            tempneu = indin(neuctxagg==1) & visarealabelagg==1;
    end
for a = 1:numel(areas)
    Pout_med.(areas{a}).(V1neugroups{n}) = median(Pout.(areas{a})(tempneu,:),1);
    Pin_med.(areas{a}).(V1neugroups{n}) = median(Pin.(areas{a})(tempneu,:),1);
    Pout_uq.(areas{a}).(V1neugroups{n}) = prctile(Pout.(areas{a})(tempneu,:),75,1);
    Pin_uq.(areas{a}).(V1neugroups{n}) = prctile(Pin.(areas{a})(tempneu,:),75,1);
    Pout_lq.(areas{a}).(V1neugroups{n}) = prctile(Pout.(areas{a})(tempneu,:),25,1);
    Pin_lq.(areas{a}).(V1neugroups{n}) = prctile(Pin.(areas{a})(tempneu,:),25,1);
    Pout_avg.(areas{a}).(V1neugroups{n}) = mean(Pout.(areas{a})(tempneu,:),1);
    Pin_avg.(areas{a}).(V1neugroups{n}) = mean(Pin.(areas{a})(tempneu,:),1);
    Pout_sem.(areas{a}).(V1neugroups{n}) = std(Pout.(areas{a})(tempneu,:),0,1)/sqrt(nnz(tempneu));
    Pin_sem.(areas{a}).(V1neugroups{n}) = std(Pin.(areas{a})(tempneu,:),0,1)/sqrt(nnz(tempneu));
end
end

%%
ICenc = ICsigallagg.ICwcfg1_presentations.ICencoder==1;
indin = ICsigallagg.ICwcfg1_presentations.indin1==1 | ICsigallagg.ICwcfg1_presentations.indin2==1 ...
    | ICsigallagg.ICwcfg1_presentations.indin3==1 | ICsigallagg.ICwcfg1_presentations.indin4==1;
% ICenc = ICsigallagg.ICwcfg1_presentations.ICencoder1==1;
% indin = ICsigallagg.ICwcfg1_presentations.indin1==1 | ICsigallagg.ICwcfg1_presentations.indin3==1;

neuoi = ICenc(neuctxagg==1) & visarealabelagg==1;
neuctrl = indin(neuctxagg==1) & visarealabelagg==1;
% neuctrl = neuctx2incl & visarealabelagg==1;
V1neugroups = {'ICenc', 'indin'};

tempCData = [0 0.7 0; 0.42 0 0.42];
figure('Position', [100 100 300 240]);
hold all
for n = 2:-1:1
errorbar(xsm0thr, Pin_avg.HVA.(V1neugroups{n}), Pin_sem.HVA.(V1neugroups{n}), 'LineWidth', 2, 'Color', tempCData(n,:))
end
xlim([0 8])
legend({'IC-enc.', 'Seg.'})
legend('boxoff')
set(gca, 'FontSize', fs)
xlabel('Factor of Jitter Mean', 'FontSize', fs)
ylabel('Input from HVA', 'FontSize', fs)

tempP = Pin.HVA;
tempPcat = cat(1, tempP(neuoi,:), tempP(neuctrl,:));
tempgroup = [ones(nnz(neuoi),1); 0.5+ones(nnz(neuctrl),1)];

figure; 
a = boxchart(tempP(neuoi,:));
hold on
b = boxchart(tempP(neuctrl,:));

% figure; 
% boxchart( reshape(repmat(xsm0thr,length(tempgroup),1),[],1), reshape(tempPin,[],1), 'GroupByColor', reshape(repmat(tempgroup,1,length(xsm0thr)),[],1));

figure; 
b = boxchart( reshape(repmat(1:length(xsm0thr),length(tempgroup),1),[],1), reshape(tempPcat,[],1), 'GroupByColor', reshape(repmat(tempgroup,1,length(xsm0thr)),[],1));
b(1).BoxFaceColor = [0 0.7 0];
b(1).MarkerColor = [0 0.7 0];
b(2).BoxFaceColor = [0.42 0 0.42];
b(2).MarkerColor = [0.42 0 0.42];
set(gca, 'XTick', 1:length(xsm0thr), 'XTickLabel', xsm0thr)

%%
yl = [0 1];
xtl = {'IC-enc.', 'Seg. Resp.'};

tempCData = [0 0.7 0; 0.42 0 0.42];
figure('Position', [100 100 270 240]);
hold all
for n = 1:2
errorbar(xsm0thr, Pin_avg.HVA.(V1neugroups{n}), Pin_sem.HVA.(V1neugroups{n}), 'LineWidth', 2, 'Color', tempCData(n,:))
text(4, yl(2)-(n-1)*0.1*range(yl), xtl{n}, 'Color', tempCData(n,:), 'VerticalAlignment', 'top', 'FontSize', fs)
end
plot(2*[1 1], yl, 'k--')
ylim(yl)
xlim([0 8])
% legend({'IC-enc.', 'Seg.'})
% legend('boxoff')
set(gca, 'FontSize', fs)
xlabel('Threshold', 'FontSize', fs)
ylabel('Input from HVA', 'FontSize', fs)

tempCData = [0 0.7 0; 0.42 0 0.42];
figure('Position', [100 100 270 240]);
hold all
for n = 1:2
errorbar(xsm0thr, Pout_avg.V1.(V1neugroups{n}), Pout_sem.V1.(V1neugroups{n}), 'LineWidth', 2, 'Color', tempCData(n,:))
text(4, yl(2)-(n-1)*0.1*range(yl), xtl{n}, 'Color', tempCData(n,:), 'VerticalAlignment', 'top', 'FontSize', fs)
end
plot(2*[1 1], yl, 'k--')
ylim(yl)
xlim([0 8])
% legend({'IC-enc.', 'Seg.'})
% legend('boxoff')
set(gca, 'FontSize', fs)
xlabel('Threshold', 'FontSize', fs)
ylabel('Output in V1', 'FontSize', fs)

%{
tempCData = [0 0.7 0; 0.42 0 0.42];
figure('Position', [100 100 270 240]);
hold all
for n = 2:-1:1
errorbar(xsm0thr, Pin_med.HVA.(V1neugroups{n}), Pin_med.HVA.(V1neugroups{n})-Pin_lq.HVA.(V1neugroups{n}), Pin_uq.HVA.(V1neugroups{n})-Pin_med.HVA.(V1neugroups{n}), 'LineWidth', 2, 'Color', tempCData(n,:))
end
xlim([0 8])
legend({'IC-enc.', 'Seg.'})
legend('boxoff')
set(gca, 'FontSize', fs)
xlabel('Threshold', 'FontSize', fs)
ylabel('Input from HVA', 'FontSize', fs)
%}

%% bar plot 
fs = 18;

ICenc = ICsigallagg.ICwcfg1_presentations.ICencoder==1;
indin = ICsigallagg.ICwcfg1_presentations.indin1==1 | ICsigallagg.ICwcfg1_presentations.indin2==1 ...
    | ICsigallagg.ICwcfg1_presentations.indin3==1 | ICsigallagg.ICwcfg1_presentations.indin4==1;
% ICenc = ICsigallagg.ICwcfg1_presentations.ICencoder1==1;
% indin = ICsigallagg.ICwcfg1_presentations.indin1==1 | ICsigallagg.ICwcfg1_presentations.indin3==1;

neuoi = ICenc(neuctxagg==1) & visarealabelagg==1;
neuctrl = indin(neuctxagg==1) & visarealabelagg==1;

for ifg = [1 3]
switch ifg
    case 1
tempP = Pin.HVA(:,xsm0thr==2);
ylab = 'Input from HVA';
    case 2
tempP = Pin.LM(:,xsm0thr==2);
ylab = 'Input from LM';
    case 3
tempP = Pout.V1(:,xsm0thr==2);
ylab = 'Output in V1';
    case 4
tempP = Pin.V1(:,xsm0thr==2);
ylab = 'Input in V1';
    case 5
tempP = Pin.V1(:,xsm0thr==2);
ylab = 'Output to LM';
end
yl = [0 0.55];
tempPcat = cat(1, tempP(neuoi,:), tempP(neuctrl,:));
tempgroup = [ones(nnz(neuoi),1); 2*ones(nnz(neuctrl),1)];

xtl = {'IC-enc.', 'Seg.'};

tempCData = [0 0.7 0; 0.42 0 0.42];
figure('Position', [100 100 180 240])
hold all
for ii = 1:2
b = boxchart(tempgroup(tempgroup==ii), tempPcat(tempgroup==ii), 'notch' , 'on', 'linewidth', 2);%, 'BoxFaceColor', tempCData(ii,:)); %, 'FontName', 'Arial')
b.BoxFaceColor = tempCData(ii,:);
b.MarkerColor = 0.15+tempCData(ii,:);
end
% yl = ylim;
p=ranksum(tempP(neuoi,:), tempP(neuctrl,:));
if p<0.05
    scatter(1.5, yl(2), 50, 'k*', 'LineWidth', 1)
end
set(gca, 'FontSize', fs, 'XTickLabelRotation', 0, 'XTick', 1:2, ...
    'XTickLabel', xtl);%, 'YTick', yl(1):ytd:yl(2))%, 'YTickLabelRotation', 45)
xlim([0.5 2.5])
ylim(yl)
xlabel('V1 Neurons', 'FontSize', fs)
ylabel(ylab, 'FontSize', fs)
% annotation('textbox', [0 0.92 1 0.1], 'string', figtit, 'fontsize', fs, 'edgecolor', 'none', 'horizontalalignment', 'center')
end

%%
figure; boxplot(degconverge, visarealabelagg)

%%
ICenc = ICsigallagg.ICwcfg1_presentations.ICencoder==1;
indin = ICsigallagg.ICwcfg1_presentations.indin1==1 | ICsigallagg.ICwcfg1_presentations.indin2==1 ...
    | ICsigallagg.ICwcfg1_presentations.indin3==1 | ICsigallagg.ICwcfg1_presentations.indin4==1;
% ICenc = ICsigallagg.ICwcfg1_presentations.ICencoder1==1;
% indin = ICsigallagg.ICwcfg1_presentations.indin1==1 | ICsigallagg.ICwcfg1_presentations.indin3==1;

neuoi = ICenc(neuctxagg==1) & visarealabelagg==1;
% neuctrl = indin(neuctxagg==1) & visarealabelagg==1;
neuctrl = neuctx2incl & visarealabelagg==1;

tempagg = degin_V1HVA(:,2);
% tempagg = CCGweightsum;

figure; hold all
hc = histogram(tempagg(neuctrl));
histogram(tempagg(neuoi), 'BinEdges', hc.BinEdges);

figure
boxchart([zeros(nnz(neuctrl),1); ones(nnz(neuoi),1)], [tempagg(neuctrl); tempagg(neuoi)])
ranksum(tempagg(neuoi), tempagg(neuctrl))
ranksum(tempagg(neuoi), tempagg(neuctrl), 'tail', 'right') % if ICenc has bigger population coupling
ranksum(tempagg(neuoi), tempagg(neuctrl), 'tail', 'left')


temp_areaagg = degin_area;
whicharea = 2;

figure; hold all
histogram(temp_areaagg(neuctrl,whicharea))
histogram(temp_areaagg(neuoi,whicharea))

figure
boxchart([zeros(nnz(neuctrl),1); ones(nnz(neuoi),1)], [temp_areaagg(neuctrl,whicharea); temp_areaagg(neuoi,whicharea)])
ranksum(temp_areaagg(neuoi,whicharea), temp_areaagg(neuctrl,whicharea))
ranksum(temp_areaagg(neuoi,whicharea), temp_areaagg(neuctrl,whicharea), 'tail', 'right') % if ICenc has bigger population coupling
ranksum(temp_areaagg(neuoi,whicharea), temp_areaagg(neuctrl,whicharea), 'tail', 'left')

recarealabels = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'};

ises = 12;
neuses = sesctxagg==ises;
tempCCG = ctxCCGjcinagg{ises};
[sv,si] = sortrows([visarealabelagg(sesctxagg==ises) nanmean(tempCCG,2)]);
%[sv,si] = sortrows([nanmean(ctxCCGweightagg{ises},2)]);
figure
imagesc(tempCCG(si,si))
hold on
for a = 0:6
    visareaborder = nnz(visarealabelagg(sesctxagg==ises)<=a);
    plot(0.5+[0 size(tempCCG,2)], 0.5+visareaborder+[0 0], 'k-', 'linewidth', 1)
    plot(0.5+visareaborder+[0 0], 0.5+[0 size(tempCCG,2)], 'k-', 'linewidth', 1)
    disp(visareaborder)
end
colormap redblue
colorbar
caxis(nanstd(tempCCG(:))*[-2 2])


for ises = ses2anal
    tempCCG = ctxCCGjcweightagg{ises};
neuses = sesctxagg==ises;
tempneu = neuoi(neuses);
[sv,si] = sortrows([visarealabelagg(neuses) nanmean(tempCCG,2)]);
figure
imagesc(tempCCG(tempneu,si))
hold on
for a = 0:6
    visareaborder = nnz(visarealabelagg(sesctxagg==ises)<=a);
    %plot(0.5+[0 size(ctrctxCCGweightagg{ises},2)], 0.5+visareaborder+[0 0], 'k-', 'linewidth', 1)
    plot(0.5+visareaborder+[0 0], 0.5+[0 size(tempCCG,2)], 'k-', 'linewidth', 1)
    disp(visareaborder)
end
colormap redblue
colorbar
caxis(nanstd(tempCCG(:))*[-2 2])
end

% tempmat = ctrctxCCGweightagg{ises};
% tempmat = tempmat(si,:);
% tempmat = tempmat(:,si);
% isequaln(tempmat, ctrctxCCGweightagg{ises}(si,si))


[v,c] = uniquecnt(sesctxagg(neuoi))

