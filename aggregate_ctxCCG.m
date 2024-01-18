%% standardize across sessions
% CCG was normalized by the geometric mean of spike count across the entire spike train length


%% compare ctxCCGsm0 with geomean of spike count
%% compare CCGfft flank mean and ctxCCGsm0
figure; plot( reshape(ctxCCGfft(:,:,ismember(CCGtli_fft,CCGtli)),[],1), reshape(ctxCCG,[],1), '.')

recarealabels = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'};
Nctxneu = nnz(neuctx);
if size(ctxCCG,1)~=Nctxneu
    error('check # visual cortex neurons')
end
flanktli = abs(CCGtli_fft)>=50;
ctxCCGflankstd = std(ctxCCGfft(:,:,flanktli),0,3);
ctxCCGflankmean = mean(ctxCCGfft(:,:,flanktli),3);

figure; plot( reshape(mean(ctxCCGsm,3),[],1), reshape(ctxCCGsm0,[],1), '.') % almost identical

figure; plot( reshape(ctxCCGflankmean,[],1), reshape(ctxCCGsm0,[],1), '.')
corr(ctxCCGflankmean(triu(true(size(ctxCCGflankmean)),1)), ctxCCGsm0(triu(true(size(ctxCCGsm0)),1))) % 0.9763 : high correlation

figure; plot( reshape(ctxCCGflankmean,[],1), reshape(ctxCCGflankstd,[],1), '.')
corr(ctxCCGflankmean(triu(true(size(ctxCCGflankmean)),1)), ctxCCGflankstd(triu(true(size(ctxCCGflankstd)),1))) % 0.1856 : low correlation

% standardize with ctxCCGsm0 and CCGfft flank std
ctxCCGz = (ctxCCG-ctxCCGsm0)./ctxCCGflankstd;

% presyn as rows, postsyn as columns
CCGtliright = CCGtli(CCGtli>=0);
ctxCCGzright = ctxCCGz(:,:,CCGtli>=0);

tempctxCCG = reshape(ctxCCGzright, size(ctxCCGzright,1)*size(ctxCCGzright,2), size(ctxCCGzright,3));
[mv,mi]=max(tempctxCCG,[],2);
[sv,si]=sort(mv);
figure; imagesc(tempctxCCG(si,:)); colorbar; clim([-5 5])
set(gca, 'YDir', 'normal')

figure; histogram(mv, 0:0.1:20)
figure; histogram2(CCGtliright(mi)', mv, CCGtliright(1)-0.5:CCGtliright(end)+0.5, 0:0.2:20, 'displaystyle', 'tile')

% excitatory connection is significant if 
% 1) there is a peak between 1~4ms (inclusive)
findpkinds = find(CCGtliright>0 & CCGtliright<5);
isccgpeak = (ctxCCGzright(:,:,findpkinds)>=ctxCCGzright(:,:,findpkinds-1)) & (ctxCCGzright(:,:,findpkinds)>=ctxCCGzright(:,:,findpkinds+1));
ctxCCGpeak = zeros(Nctxneu, Nctxneu);
ctxCCGpeakT = zeros(Nctxneu, Nctxneu);
% largest peak value for each CCG
for ii = 1:length(findpkinds)
    CCGtp = squeeze(ctxCCGzright(:,:,findpkinds(ii)));
    pairs2overwrite = squeeze(isccgpeak(:,:,ii)) & CCGtp>ctxCCGpeak;
    ctxCCGpeak(pairs2overwrite) = CCGtp(pairs2overwrite);
    ctxCCGpeakT(pairs2overwrite) = CCGtliright(findpkinds(ii));
end

V1neurons = contains(neulocctx, 'VISp') & ~contains(neulocctx, 'VISpm');
if ~isequal(unique(visarealabels(V1neurons)), 1)
    error('check visarealabels')
end
connthresh = 10;
% proportion of feedforward outputs from each V1 neuron (out of all HVA neurons)
V1FFout = mean(ctxCCGpeak(V1neurons, ~V1neurons)>=connthresh, 2);

% proportion of feedback inputs to each V1 neuron
V1FBin = mean(ctxCCGpeak(~V1neurons, V1neurons)>=connthresh, 1)';

figure; plot(V1FFout, V1FBin, 'o')
figure; histogram2(V1FFout, V1FBin, 'BinWidth',1/(nnz(~V1neurons)), 'displaystyle', 'tile')

% >0 connections
mean(V1FFout>0)
mean(V1FBin>0)
mean(V1FFout>0)*mean(V1FBin>0) % expected overlap betweeen V1 neurons with FFout vs FBin 
mean(V1FFout>0 & V1FBin>0) % actual overlap between V1 neurons with FFout vs FBin

% spearman correlation between number of FF vs FB connections
V1connwHVA = ~(V1FFout==0 & V1FBin==0);
corr(V1FFout(V1connwHVA), V1FBin(V1connwHVA), 'Type', 'Spearman')

% what if we limit to V1 neurons that send OR receive local recurrent (within-V1) connections?
V1recurneu = sum(ctxCCGpeak(V1neurons, V1neurons)>=connthresh,2)>0 | sum(ctxCCGpeak(V1neurons, V1neurons)>=connthresh,1)'>0;
V1neuind = find(V1neurons);
rV1neuind = V1neuind(V1recurneu);
%% 
% feedforward connectivity: V1 to HVA
% feedback connectivity: HVA to V1

% thresholds = [1:9 10:10:100];
thresholds = 1:0.1:20;
V1FFsenderprop = zeros(size(thresholds));
V1FBreceiverprop = zeros(size(thresholds));
V1FFoutFBinoverlap = zeros(size(thresholds));
V1FFoutFBinrho = zeros(size(thresholds)); % spearman correlation
V1Rsenderprop = zeros(size(thresholds));
V1Rreceiverprop = zeros(size(thresholds));
V1RoutRinoverlap = zeros(size(thresholds));
V1RinFBinoverlap = zeros(size(thresholds));
V1RoutFFoutoverlap = zeros(size(thresholds));
V1RinFFoutoverlap = zeros(size(thresholds));
V1RoutFBinoverlap = zeros(size(thresholds));
for t = 1:length(thresholds)
    connthresh = thresholds(t);
    V1FFout = mean(ctxCCGpeak(V1neurons, ~V1neurons)>=connthresh, 2);
    V1FBin = mean(ctxCCGpeak(~V1neurons, V1neurons)>=connthresh, 1)';
    V1connwHVA = ~(V1FFout==0 & V1FBin==0);
    V1FFsenderprop(t) = mean(V1FFout>0);
    V1FBreceiverprop(t) = mean(V1FBin>0);
    V1FFoutFBinoverlap(t) = mean(V1FFout>0 & V1FBin>0);
    V1FFoutFBinrho(t) = corr(V1FFout(V1connwHVA), V1FBin(V1connwHVA), 'Type', 'Spearman');

    V1Rout = mean(ctxCCGpeak(V1neurons, V1neurons)>=connthresh, 2);
    V1Rin = mean(ctxCCGpeak(V1neurons, V1neurons)>=connthresh, 1)';
    V1Rsenderprop(t) = mean(V1Rout>0);
    V1Rreceiverprop(t) = mean(V1Rin>0);
    V1RoutRinoverlap(t) = mean(V1Rout>0 & V1Rin>0);
    V1RinFBinoverlap(t) = mean(V1Rin>0 & V1FBin>0);
    V1RoutFFoutoverlap(t) = mean(V1Rout>0 & V1FFout>0);
    V1RinFFoutoverlap(t) = mean(V1Rin>0 & V1FFout>0);
    V1RoutFBinoverlap(t) = mean(V1Rout>0 & V1FBin>0);
end
% overlap turns out to be stronger than expected...

figure; 
for ii = 1:6
    switch ii
        case 1
            a = V1Rsenderprop;
            b = V1Rreceiverprop;
            o = V1RoutRinoverlap;
            alab = 'R-senders';
            blab = 'R-receivers';
        case 2
            a = V1FFsenderprop;
            b = V1Rsenderprop;
            o = V1RoutFFoutoverlap;
            alab = 'FF-senders';
            blab = 'R-senders';
        case 3
            a = V1FFsenderprop;
            b = V1Rreceiverprop;
            o = V1RinFFoutoverlap;
            alab = 'FF-senders';
            blab = 'R-receivers';
        case 4
            a = V1FFsenderprop;
            b = V1FBreceiverprop;
            o = V1FFoutFBinoverlap;
            alab = 'FF-senders';
            blab = 'FB-receivers';
        case 5
            a = V1Rsenderprop;
            b = V1FBreceiverprop;
            o = V1RoutFBinoverlap;
            alab = 'R-senders';
            blab = 'FB-receivers';
        case 6
            a = V1Rreceiverprop;
            b = V1FBreceiverprop;
            o = V1RinFBinoverlap;
            alab = 'R-receivers';
            blab = 'FB-receivers';
    end
subplot(2,3,ii)
hold all
scatter(a.*b, o, 10, thresholds)
xl = xlim; plot(xl, xl)
cb = colorbar; cb.Label.String = 'CCG flank std threshold';
xlabel(['expected overlap between V1 ' alab ' vs ' blab])
ylabel(['actual overlap between V1 ' alab ' vs ' blab])
end

figure; hold all
plot(thresholds, V1FFoutFBinrho)
plot([thresholds(1) thresholds(end)], [0 0])

% what if we limit to V1 neurons that send OR receive local recurrent
% (within-V1) connections? -- same result...
rV1FFsenderprop = zeros(size(thresholds));
rV1FBreceiverprop = zeros(size(thresholds));
rV1FFoutFBinoverlap = zeros(size(thresholds));
rV1FFoutFBinrho = zeros(size(thresholds)); % spearman correlation
for t = 1:length(thresholds)
    connthresh = thresholds(t);
    V1recurneu = sum(ctxCCGpeak(V1neurons, V1neurons)>=connthresh,2)>0 | sum(ctxCCGpeak(V1neurons, V1neurons)>=connthresh,1)'>0;
    V1neuind = find(V1neurons);
    rV1neuind = V1neuind(V1recurneu);
    rV1FFout = mean(ctxCCGpeak(rV1neuind, ~V1neurons)>=connthresh, 2);
    rV1FBin = mean(ctxCCGpeak(~V1neurons, rV1neuind)>=connthresh, 1)';
    rV1connwHVA = ~(rV1FFout==0 & rV1FBin==0);
    rV1FFsenderprop(t) = mean(rV1FFout>0);
    rV1FBreceiverprop(t) = mean(rV1FBin>0);
    rV1FFoutFBinoverlap(t) = mean(rV1FFout>0 & rV1FBin>0);
    rV1FFoutFBinrho(t) = corr(rV1FFout(rV1connwHVA), rV1FBin(rV1connwHVA), 'Type', 'Spearman');
end

figure; hold all
scatter(rV1FFsenderprop.*rV1FBreceiverprop, rV1FFoutFBinoverlap, 10, thresholds)
xl = xlim; plot(xl, xl)
cb = colorbar; cb.Label.String = 'CCG flank std threshold';
xlabel('expected overlap between recurrent-V1 FF-senders vs FB-receivers')
ylabel('actual overlap between recurrent-V1 FF-senders vs FB-receivers')
% again, overlap turns out to be stronger than expected...

figure; hold all
scatter(V1FFoutFBinoverlap, rV1FFoutFBinoverlap, 10, thresholds)
xl = xlim; plot(xl, xl)

figure; hold all
plot(thresholds, rV1FFoutFBinrho)
plot([thresholds(1) thresholds(end)], [0 0])
