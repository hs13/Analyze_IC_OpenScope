function BKshuf = analyzeBKshuf(R, trialorder)

Nneurons = size(R, 1);
if size(R,2)~=length(trialorder)
    error('check R dimensions')
end

Palpha = 0.05;
correctmultcomp = true;
% verbosity = false;

%% Kruskalwallis statistics
% perhaps do ranksum tests instead of

% PkwBKI and PmcBKI can be skipped if time is scarce (will save ~15 seconds)
BKtt = [0 106 107 110 111]; % Blank and Kanizsa trial type
itt = 1;
switch itt
    case 1
        temptt = BKtt;
    case 2
        temptt = BItt;
    case 3
        temptt = BICREltt;
    case 4
        temptt = BICREl1tt;
    case 5
        temptt = BICREl2tt;
end
temptrialsoi = ismember(trialorder, temptt);
temptt = unique(trialorder(temptrialsoi));

tempPkw = zeros(Nneurons,1);
tempPmc = NaN(Nneurons, nchoosek(numel(temptt),2));
for ci = 1:Nneurons
    [p,tbl,stats] = kruskalwallis(R(ci,temptrialsoi), ...
        trialorder(temptrialsoi), 'off');
    tempPkw(ci) = p;
    % if p<0.05
    c = multcompare(stats, 'Display','off');
    tempPmc(ci,:) = c(:,end);
    % end
end
tempttpair = [temptt(c(:,1))' temptt(c(:,2))'];

switch itt
    case 1
        PkwBK = tempPkw;
        PmcBK = tempPmc;
        BKttpair = tempttpair;
    case 2
        PkwBI = tempPkw;
        PmcBI = tempPmc;
        BIttpair = tempttpair;
    case 3
        %             PkwBKI = tempPkw;
        %             PmcBKI = tempPmc;
        %             BKIttpair = tempttpair;
        PkwBICREl = tempPkw;
        PmcBICREl = tempPmc;
        BICRElttpair = tempttpair;
    case 4
        PkwBICREl1 = tempPkw;
        PmcBICREl1 = tempPmc;
        BICREl1ttpair = tempttpair;
    case 5
        PkwBICREl2 = tempPkw;
        PmcBICREl2 = tempPmc;
        BICREl2ttpair = tempttpair;
end

blanktrials = trialorder==0;

BKttpairrows = find(BKttpair(:,1)==0);
SP_BK = NaN(Nneurons, numel(BKttpairrows) );
Pmww_BK = NaN(Nneurons, numel(BKttpairrows) );
for typi = 1:numel(BKttpairrows)
    Ktrials = trialorder==BKttpair(BKttpairrows(typi),2);
    if nnz(blanktrials)==0 || nnz(Ktrials)==0
        continue
    end
    labels = [zeros(1,nnz(blanktrials)) ones(1,nnz(Ktrials))];
    
    for ci = 1:Nneurons
        scores = [R(ci, blanktrials) R(ci, Ktrials)];
        
        [X,Y,T,AUC,OPTROCPT] = perfcurve(labels,scores,'1');% , 'NBoot',Nshuf);
        SP_BK(ci,typi) = AUC;
        
        Pmww_BK(ci, typi) = ranksum(R(ci, blanktrials), R(ci, Ktrials));
    end
end

% sigmcBK = false(size(PmcBK));
% sigmcBK(PmcBK<0.05) = true;
sigmcBK = zeros( Nneurons, size(SP_BK,2) );
sigmcBK(SP_BK>0.5) = 1;
sigmcBK(SP_BK<0.5) = -1;
sigmcBK(PmcBK(:, BKttpair(:,1)==0)>=Palpha) = 0;
sigmcBK(PkwBK>=Palpha, :)=0;

%% IC vs RC encoder neuron
ICencoder = (PkwBK<Palpha) & (sigmcBK(:,1)==1 | sigmcBK(:,4)==1) ...
    & (sigmcBK(:,2)==0) & (sigmcBK(:,3)==0);
RCencoder = (PkwBK<Palpha) & (sigmcBK(:,2)==1 | sigmcBK(:,3)==1) ...
    & (sigmcBK(:,1)==0) & (sigmcBK(:,4)==0);
inducerencoder = (PkwBK<Palpha) & xor(sigmcBK(:,1)==1, sigmcBK(:,4)==1) ...
    & xor(sigmcBK(:,2)==1, sigmcBK(:,3)==1);
inducerresponsive = inducerencoder | ( PkwBK<Palpha & all(sigmcBK(:,1:4)==1, 2) );

% if verbosity
% fprintf('IC-encoder %d, RC-encoder %d, inducer-encoder %d inducer resposive %d\n', ...
%     nnz(ICencoder), nnz(RCencoder), nnz(inducerencoder), nnz(inducerresponsive) )
% end

bootfields = {'trialorder', 'PkwBK', 'PmcBK', 'SP_BK', 'Pmww_BK', 'sigmcBK', ...
    'ICencoder', 'RCencoder', 'inducerencoder', 'inducerresponsive'};
for f = 1:numel(bootfields)
    BKshuf.(bootfields{f}) = eval(bootfields{f});
end

end