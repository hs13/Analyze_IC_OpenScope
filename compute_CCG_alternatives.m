%% alternative methods for computing CCG: use xcorr or use fft
datadir = '/Users/hyeyoung/Documents/DATA/OpenScopeData/000248/';

% FOUND SESSION 4 (sub_1174569641) TO BE AN EXACT DUPLICATE OF SESSION 1 (sub_1171903433)!!!!!!
nwbsessions = {'sub_1171903433','sub_1172968426','sub_1172969394', 'sub_1174569641', ...
    'sub_1175512783','sub_1176214862','sub_1177693342','sub_1181314060', ...
    'sub_1181585608','sub_1182593231','sub_1183369803','sub_1186544726', ...
    'sub_1189891322','sub_1194090570'};
Nsessions = numel(nwbsessions);
ses2anal = 1:Nsessions;
ses2anal(4) = [];


ises = 14;
pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load([pathpp 'postprocessed_probeC.mat'], 'vis')

%%
load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location'
load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur') %'unit_times_data'
load([pathpp 'spiketimes.mat'])
load([pathpp 'visresponses.mat'])
% isequal(ststartend, [floor(min(unit_times_data)/Tres)+1 floor(max(unit_times_data)/Tres)+1])
Tres = 0.001; % 1ms

elecid = electrode_id+1;
revmapelecid = NaN(max(elecid),1);
revmapelecid(elecid) = 1:numel(elecid);
electrode_location = cellstr(electrode_location);
neuloc = electrode_location(revmapelecid(unit_peakch+1));

neuctx = contains(neuloc, 'VIS');
neuctxind = find(neuctx) ;

%%
Nneuctx = numel(neuctxind);
neulocctx = neuloc(neuctxind);
% isequal(neulocctx, find(contains(neuloc, 'VIS')))

recarealabels = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'};
visarealabels = zeros(Nneuctx,1);
for a = 1:numel(recarealabels)
    if strcmp(recarealabels{a}, 'VISp')
        neuinarea = contains(neulocctx, 'VISp') & ~contains(neulocctx, 'VISpm');
    else
        neuinarea = contains(neulocctx, recarealabels{a});
    end
    visarealabels(neuinarea) = a;
end

%%
ctxspiketrain = false(Nneuctx, ststartend(end) );
for ii = 1:Nneuctx
    ci = neuctxind(ii);
    ctxspiketrain(ii, floor(spiketimes{ci}/Tres)+1) = true;
end

ICblockstartend = [vis.ICwcfg1_presentations.start_time(1) vis.ICwcfg1_presentations.stop_time(end)];
ICblockstartend = floor(ICblockstartend/Tres)+1;
ctxspiketrain = ctxspiketrain(:,ICblockstartend(1):ICblockstartend(2));
% ctxspiketrain = ctxspiketrain(:,ststartend(1):ststartend(2));

% kersm = ones(1,25);
% kersm = (kersm/sum(kersm));
% 
% ctxstsm = convn(ctxspiketrain, kersm, 'valid'); % this takes ~2min
% ctxstsm = single(ctxstsm);

stlen = size(ctxspiketrain, 2);
spkcntvec = sum(ctxspiketrain,2);

% save([pathpp 'spiketimes.mat'], 'spiketimes', 'neuloc', 'neuctx', 'neuctxind', 'neulocctx', ...
%     'recarealabels', 'visarealabels', 'ststartend', 'ICblockstartend', '-v7.3')

%%
ccgpath = '/Volumes/GoogleDrive-116160770365018316196/My Drive/DATA/ICexpts_postprocessed_OpenScope/CCG/';
pathccg = [ccgpath nwbsessions{ises} filesep];
load([pathccg 'ctxCCG.mat'])

%%
NFFT = 2^ceil(log2(stlen));
CCGffttli

% cij = randperm(Nneuctx,2);
cij = [761 846];
temp1 = ctxspiketrain(cij(1),:);
temp2 = ctxspiketrain(cij(2),:);

% takes 0.32s
tic
tempccgfft = fftshift(ifft(fft(temp2, NFFT) .* conj(fft(temp1, NFFT))));
toc

% figure; hold all
% plot(-NFFT/2:NFFT/2-1, (ifft(fft(temp2, NFFT) .* conj(fft(temp1, NFFT)))))
% plot(-NFFT/2:NFFT/2-1, fftshift(ifft(fft(temp2, NFFT) .* conj(fft(temp1, NFFT)))))

figure
hold all
plot(CCGtli, squeeze(ctxCCG(cij(1),cij(2),:)), 'k-')
plot(-NFFT/2:NFFT/2-1, tempccgfft/sqrt(spkcntvec(cij(1))*spkcntvec(cij(2))), 'r--')
xlim([-13 13])

% takes 2.7 seconds for 10 cells
cvec = cij(2)-9:cij(2);
tic
tempCCGfft = fftshift(ifft(fft(ctxspiketrain(cvec,:)', NFFT) .* conj(fft(temp1', NFFT))),1);
toc

figure
hold all
plot(-NFFT/2:NFFT/2-1, tempCCGfft/sqrt(spkcntvec(cij(1))*spkcntvec(cij(2))), 'k-')
plot(CCGtli, squeeze(ctxCCG(cij(1),cij(2),:)), 'r--')
xlim([-13 13])

ii = find(cvec==cij(2));
figure
hold all
plot(-NFFT/2:NFFT/2-1, tempCCGfft(:,ii)/sqrt(spkcntvec(cij(1))*spkcntvec(cij(2))), 'k-')
plot(CCGtli, squeeze(ctxCCG(cij(1),cij(2),:)), 'r--')
xlim([-13 13])

cvec = cij(1):cij(1)+9;
tic
tempCCGfft = fftshift(ifft(fft(temp2', NFFT) .* conj(fft(ctxspiketrain(cvec,:)', NFFT))),1);
toc

figure
hold all
plot(-NFFT/2:NFFT/2-1, tempCCGfft(:,cvec==cij(1))/sqrt(spkcntvec(cij(1))*spkcntvec(cij(2))), 'k-')
plot(CCGtli, squeeze(ctxCCG(cij(1),cij(2),:)), 'r--')
xlim([-13 13])

% takes 26s for 100 cells
tic
tempCCGfft = fftshift(ifft(fft(ctxspiketrain(1:100,:)', NFFT) .* conj(fft(temp1', NFFT))),1);
toc
% for 10^6 pairs, estimated it would take ~72 hours, ~36 hours if only
% calculating the upper right triangle (ci+1:end)

%% compare CCGsm0 with standard deviation of CCG flanks
% ises = 5; % 5 or 6
% subid = nwbsessions{ises};

% subid = 'sub_1175512783';
subid = 'sub_1176214862';
pathccg = ['/Volumes/GoogleDrive-116160770365018316196/My Drive/DATA/ICexpts_postprocessed_OpenScope/CCG/' subid filesep];
load([pathccg 'ctxCCG.mat'])
load([pathccg 'ctxCCGsm.mat'])
load([pathccg 'ctxCCGlong.mat'])
load([pathccg 'ctxCCGsm0.mat'])

datadir = '/Users/hyeyoung/Documents/DATA/OpenScopeData/000248/';
pathpp = [datadir 'postprocessed' filesep subid filesep];
load([pathpp 'visresponses.mat'])

flanktli = abs(CCGtli_long)>=50 & abs(CCGtli_long)<=100;
ctxCCGflankstd = std(ctxCCGlong(:,:,flanktli),0,3);

ICenc = ICsigall.ICwcfg1_presentations.ICencoder==1;
segresp = ICsigall.ICwcfg1_presentations.indin1==1 | ICsigall.ICwcfg1_presentations.indin2==1 | ...
    ICsigall.ICwcfg1_presentations.indin3==1 | ICsigall.ICwcfg1_presentations.indin4==1;

ctxICenc = ICenc(neuctx);
ctxsegresp = segresp(neuctx);
meanFRctx = meanFRall(neuctx);

ctxCCGjc = ctxCCG - ctxCCGsm;
[mv,mi]= max(ctxCCGjc(:,:,CCGtli>0), [],3);
sigpostsm0 = mv>2*ctxCCGsm0;
sigpoststd = mv>7*ctxCCGflankstd;

[mv,mi]= max(ctxCCGjc(:,:,CCGtli<0), [],3);
sigpresm0 = mv>2*ctxCCGsm0;
sigprestd = mv>7*ctxCCGflankstd;


figure; plot(ctxCCGsm0, ctxCCGflankstd, '.')
xlabel('ctxCCGsm0')
ylabel('ctxCCGflankstd')

figure; hold all
plot(mean(sigpostsm0,2), mean(sigpoststd,2), '.')
plot(mean(sigpostsm0(ctxICenc,:),2), mean(sigpoststd(ctxICenc,:),2), 'o')
xlabel('Outputs>2*ctxCCGsm0')
ylabel('Outputs>7*ctxCCGflankstd')
ranksum(mean(sigpostsm0(ctxICenc,:),2),  mean(sigpostsm0(ctxsegresp,:),2))
ranksum(mean(sigpoststd(ctxICenc,:),2),  mean(sigpoststd(ctxsegresp,:),2))

ranksum(mean(sigpostsm0(ctxICenc,visarealabels==1),2),  mean(sigpostsm0(ctxsegresp,visarealabels==1),2))
ranksum(mean(sigpoststd(ctxICenc,visarealabels==1),2),  mean(sigpoststd(ctxsegresp,visarealabels==1),2))
ranksum(mean(sigpoststd(ctxICenc,visarealabels==1),2),  mean(sigpoststd(:,visarealabels==1),2))
ranksum(mean(sigpoststd(:,visarealabels==1),2),  mean(sigpoststd(ctxsegresp,visarealabels==1),2))

figure
hold all
hc = histogram(mean(sigpoststd(:,visarealabels==1),2));
histogram(mean(sigpoststd(ctxsegresp,visarealabels==1),2), 'BinEdges', hc.BinEdges)
histogram(mean(sigpoststd(ctxICenc,visarealabels==1),2), 'BinEdges', hc.BinEdges)

figure; hold all
plot(mean(sigpresm0,2), mean(sigprestd,2), '.')
plot(mean(sigpresm0(ctxICenc,:),2), mean(sigprestd(ctxICenc,:),2), 'o')
xlabel('Inputs>2*ctxCCGsm0')
ylabel('Inputs>7*ctxCCGflankstd')
ranksum(mean(sigpresm0(ctxICenc,:),2),  mean(sigpresm0(ctxsegresp,:),2))
ranksum(mean(sigprestd(ctxICenc,:),2),  mean(sigprestd(ctxsegresp,:),2))

ranksum(mean(sigpresm0(ctxICenc,visarealabels>=2),2),  mean(sigpresm0(ctxsegresp,visarealabels>=2),2))
ranksum(mean(sigprestd(ctxICenc,visarealabels>=2),2),  mean(sigprestd(ctxsegresp,visarealabels>=2),2))

figure
hold all
hc = histogram(mean(sigprestd(:,visarealabels>=2),2));
% histogram(mean(sigprestd(ctxsegresp,visarealabels>=2),2), 'BinEdges', hc.BinEdges)
histogram(mean(sigprestd(ctxICenc,visarealabels>=2),2), 'BinEdges', hc.BinEdges)

figure; hold all
plot(mean(sigpresm0,2), meanFRctx, '.')
plot(mean(sigpresm0(ctxICenc,:),2), meanFRctx(ctxICenc), 'o')
xlabel('Inputs>2*ctxCCGsm0')
ylabel('meanFRall')

figure; hold all
plot( mean(sigprestd,2), meanFRctx, '.')
plot(mean(sigprestd(ctxICenc,:),2), meanFRctx(ctxICenc), 'o')
xlabel('Inputs>7*ctxCCGflankstd')
ylabel('meanFRall')
