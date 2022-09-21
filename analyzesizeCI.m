% % 10000s: which type (classic 0 vs inverse 1), 1000s which ctrsizes, 
% % 10-100s: which RFcenter, 1s: which direction

% 'gkcellind', 'kSLtunedneurons', 'kSLprefori', 'kSGtunedneurons', 'kSGprefori'

function [sizeCI, oriparams] = analyzesizeCI(R, trialorder)

Nneurons = size(R,1); 

% determine tuning with fullscreen grating
% fulltrials = floor(trialorder/1000)==11;

% TODO: ADD Pkw_ori
Noris = 8;
Rori = zeros(Nneurons, Noris);
for iori = 1:Noris
    Rori(:,iori) = mean(R(:,trialorder==11000+iori),2);
end

[Rpref, prefiori] = max(Rori,[],2);
if ~isequal(Rpref, Rori(sub2ind(size(Rori), (1:Nneurons)', prefiori)))
    error('check Rpref/prefiori')
end

orthiori = zeros(Nneurons,1);
Rorth = NaN(Nneurons,1);
for iori = 1:Noris
    iorth0 = mod(iori+2 -1, Noris)+1;
    iorth1 = mod(iori-2 -1, Noris)+1;
    neuoi0 = prefiori==iori & Rori(:,iorth0)<=Rori(:,iorth1);
    orthiori(neuoi0) = iorth0;
    Rorth(neuoi0) = Rori(neuoi0, iorth0);
    neuoi1 = prefiori==iori & Rori(:,iorth0)>Rori(:,iorth1);
    orthiori(neuoi1) = iorth1;
    Rorth(neuoi1) = Rori(neuoi1, iorth1);
end
if nnz(isnan(Rorth)) || nnz(orthiori==0)
    error('check orthiori and Rorth')
end
if ~isequal(Rorth, Rori(sub2ind(size(Rori), (1:Nneurons)', orthiori)))
    error('check Rorth/orthiori')
end

OSI = (Rpref-Rorth)./(Rpref+Rorth);

%%
trialiori = mod(trialorder,100);

OP = NaN(Nneurons,1);
Pmww_OP = NaN(Nneurons, 1);
for ci = 1:Nneurons
    preforitrials = floor(trialorder/1000)==11 & trialiori==prefiori(ci);
    orthoritrials = floor(trialorder/1000)==11 & trialiori==orthiori(ci);
    
    scores = [R(ci, orthoritrials) R(ci, preforitrials)];
    labels = [zeros(1,nnz(orthoritrials)) ones(1,nnz(preforitrials))];

    [X,Y,T,AUC,OPTROCPT] = perfcurve(labels,scores,'1');% , 'NBoot',Nshuf);
    OP(ci) = AUC;

    Pmww_OP(ci) = signrank(R(ci, orthoritrials), R(ci, preforitrials));
end

orifields = {'Rori', 'prefiori', 'orthiori', 'OSI', 'OP', 'Pmww_OP'};
oriparams = struct();
for f = 1:numel(orifields)
    oriparams.(orifields{f}) = eval(orifields{f});
end

%%
Nsz = 6;
Rsizeclassic = zeros(Nneurons, Nsz);
Rsizeinverse = zeros(Nneurons, Nsz);
% 10000s: which type (classic 0 vs inverse 1), 1000s which ctrsizes, 10-100s: which RFcenter, 1s: which direction'
for szi = 1:Nsz
    csztrials = floor(trialorder/10000) == 0 & floor(mod(trialorder, 10000) / 1000) == szi ;
    Rsizeclassic(:,szi) = mean(R(:, csztrials),2);
    
    isztrials = floor(trialorder/10000) == 1 & floor(mod(trialorder, 10000) / 1000) == szi ;
    Rsizeinverse(:,szi) = mean(R(:, isztrials),2);
end

[~, sizeindclassic] = max(Rsizeclassic, [], 2);
[~, sizeindinverse] = max(Rsizeinverse, [], 2);

indRFCIsize = floor(mod(trialorder, 10000)/1000);
indRFCIloc = floor(mod(trialorder, 1000)/10);

% for each neuron, choose the preferred orientation 
Pkw_sizeclassic = zeros(Nneurons, 1);
Pkw_sizeinverse = zeros(Nneurons, 1);
for ci = 1:Nneurons
    csztrials =  trialorder~=0 & floor(trialorder/10000) == 0;
    [p,tbl,stats] = kruskalwallis(R(ci,csztrials), indRFCIsize(csztrials), 'off');
    Pkw_sizeclassic(ci) = p;
    
    isztrials =  trialorder~=0 & floor(trialorder/10000) == 1;
    [p,tbl,stats] = kruskalwallis(R(ci,isztrials), indRFCIsize(isztrials), 'off');
    Pkw_sizeinverse(ci) = p;
end

sizeCIfields = {'Rsizeclassic', 'Rsizeinverse', 'sizeindclassic', 'sizeindinverse', ...
    'Pkw_sizeclassic', 'Pkw_sizeinverse'};%, 'pRsizeclassic', 'pRsizeinverse', 'psizeclassic', 'psizeinverse'};
sizeCI = struct();
for f = 1:numel(sizeCIfields)
    sizeCI.(sizeCIfields{f}) = eval(sizeCIfields{f});
end


end