% % 10000s: which type (classic 0 vs inverse 1), 1000s which ctrsizes, 
% % 10-100s: which RFcenter, 1s: which direction
% RFCI = struct();
% tloi = psthtli>0 & psthtli<=1000;
% tempR = squeeze(1000*mean(psth.RFCI(tloi,:,:), 1))';
% temptrialorder = vis.RFCI_presentations.trialorder(1:4:end);
% RFCI.(visblocknames{b}) = analyzeRFCI(tempR, temptrialorder);

% 'RFcentersrel', 'cRFsigneurons', 'cRFind', 'iRFsigneurons', 'iRFind', ...
% 'RFcentersrel9', 'RFcenterinds9', 'cRFind9', 'pRFclassic9', 'cRFsigexcl9', ...
% 'iRFind9', 'pRFinverse9', 'iRFsigexcl9')

function RFCI = analyzeRFCI(R, trialorder, sponFR)

Nneurons = size(R,1);
if numel(spnFR) ~= Nneurons
    error('check sponFR')
end

Nrfs = 9;
Rrfclassic = zeros(Nneurons, Nrfs);
Rrfinverse = zeros(Nneurons, Nrfs);
% 10000s: which type (classic 0 vs inverse 1), 1000s which ctrsizes, 10-100s: which RFcenter, 1s: which direction'
for rfi = 1:Nrfs
    crftrials = floor(trialorder/10000) == 0 & floor(mod(trialorder, 1000) / 10) == rfi-1 ;
    Rrfclassic(:,rfi) = mean(R(:, crftrials),2);
    
    irftrials = floor(trialorder/10000) == 1 & floor(mod(trialorder, 1000) / 10) == rfi-1 ;
    Rrfinverse(:,rfi) = mean(R(:, irftrials),2);
end

[~, RFindclassic] = max(Rrfclassic, [], 2);
[~, RFindinverse] = max(Rrfinverse, [], 2);

Pkw_rfclassic = zeros(Nneurons, 1);
Pkw_rfinverse = zeros(Nneurons, 1);
for ci = 1:Nneurons
    %     indRFCIsize = floor(mod(rectrialorder, 10000)/1000);
    %     indRFCIloc = floor(mod(rectrialorder, 1000)/10);
    crftrials =  floor(trialorder/10000) == 0;
    [p,tbl,stats] = kruskalwallis(R(ci,crftrials), trialorder(crftrials), 'off');
    Pkw_rfclassic(ci) = p;
    
    irftrials =  floor(trialorder/10000) == 1;
    [p,tbl,stats] = kruskalwallis(R(ci,irftrials), trialorder(irftrials), 'off');
    Pkw_rfinverse(ci) = p;
end


pRrfclassic = zeros(Nneurons, Nrfs);
pRrfinverse = zeros(Nneurons, Nrfs);
for rfi = 1:Nrfs
    crftrials = floor(trialorder/10000) == 0 & floor(mod(trialorder, 1000) / 10) == rfi-1 ;    
    irftrials = floor(trialorder/10000) == 1 & floor(mod(trialorder, 1000) / 10) == rfi-1 ;
    for ci = 1:Nneurons
        pRrfclassic(ci,rfi) = signrank(R(ci,crftrials)-sponFR(ci));
        pRrfinverse(ci,rfi) = signrank(R(ci,irftrials)-sponFR(ci));
    end
end

pRFclassic = pRrfclassic(sub2ind(size(pRrfclassic), (1:Nneurons)', RFindclassic));
pRFinverse = pRrfinverse(sub2ind(size(pRrfinverse), (1:Nneurons)', RFindinverse));

RFCIfields = {'Rrfclassic', 'Rrfinverse', 'RFindclassic', 'RFindinverse', ...
    'Pkw_rfclassic', 'Pkw_rfinverse', 'pRrfclassic', 'pRrfinverse', 'pRFclassic', 'pRFinverse'};
RFCI = struct();
for f = 1:numel(RFCIfields)
    RFCI.(RFCIfields{f}) = eval(RFCIfields{f});
end

end