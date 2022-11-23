addpath(genpath('C:\Users\Hyeyoung\Documents\matnwb'))
addpath(genpath('H:\CODE\Analyze_OpenScope'))
addpath(genpath('H:\CODE\helperfunctions'))

datadir = 'D:\OpenScopeData\000248\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(contains(nwbsessions, 'sub-'));
Nsessions = numel(nwbsessions)-1;

probes2agg = {'A', 'B', 'C', 'D', 'E', 'F'};
% probes2agg = {'C'};

visesponsesagg = cell(numel(probes2agg), Nsessions);
for ises = 1:Nsessions
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    for iprobe = 1:numel(probes2agg)
        visesponsesagg{iprobe, ises} = load(sprintf('%svisresponses_probe%s.mat', pathpp, probes2agg{iprobe}));
    end
end

%% accumulate neuctx
neuctxagg = cell(numel(probes2agg), Nsessions);

%%
ICblocknames = {'ICkcfg0', 'ICkcfg1', 'ICwcfg0', 'ICwcfg1'};
encagg = struct();
vr = cat(1, visesponsesagg{3,:});
vr = cat(1, vr.ICsig);
for b = 1:numel(ICblocknames)
    tempICsig = cat(1, vr.([ICblocknames{b} '_presentations']) );
    ICencagg.(ICblocknames{b}) = cat(1, tempICsig.ICencoder);
end

%% pref-ori distribution of IC/RC encoders

%% RF distribution of IC-encoders

%% size tuning of IC-encoders

%% prefori probably needs more repeats

%% szvec = [0, 4, 8, 16, 32, 64];

ises=4; 
pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load(sprintf('%spostprocessed_probeC.mat', pathpp))
whichblock = 'sizeCI_presentations';
neuoi = true(size(neuoind));

vistrialtypes.sizeCI_presentations = unique(vis.sizeCI_presentations.trialorder);
vistrialorder.sizeCI_presentations = vis.sizeCI_presentations.trialorder;

Ntt = numel(vistrialtypes.sizeCI_presentations);
vistrialrep.sizeCI_presentations = zeros(Ntt,1);
temppsth = zeros(length(psthtli), Ntt, nnz(neuoi));
for ii = 1:Ntt
    trialsoi = vis.sizeCI_presentations.trialorder==vistrialtypes.sizeCI_presentations(ii);
    vistrialrep.sizeCI_presentations(ii) = nnz(trialsoi);
    temppsth(:,ii,:) = mean(1000*psth.sizeCI_presentations(:,trialsoi,neuoi), 2);
end

tlon = psthtli>0 & psthtli<=250;
tloff = psthtli>250 & psthtli<=750;

tempRavg = squeeze(mean(temppsth(tlon,:,:),1));

trialtypes = reshape(vistrialtypes.sizeCI_presentations, Ndirs,Nszs*2);
tempRavg = reshape(tempRavg, Ndirs,Nszs*2,nnz(neuoi));

[~,iprefdirCG64] = max(tempRavg(:,6,:),[],1);
[~,iprefdirIG0] = max(tempRavg(:,7,:),[],1);
[~,iprefdirIG4] = max(tempRavg(:,8,:),[],1);

figure; histogram2(iprefdirCG64, iprefdirIG0, 'displaystyle', 'tile')
