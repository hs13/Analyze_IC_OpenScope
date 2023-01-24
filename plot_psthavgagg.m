load('D:\OpenScopeData\000248\postprocessed\psthavgagg.mat')
load(['G:\My Drive\DATA\ICexpts_submission22\openscope_popavg_agg.mat'])

kerwinhalf = 5; kersigma = 2;
kerwinhalf = 25; kersigma = 5;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));

ICexptidn = 'ICwcfg1_presentations';
tt2p = [106 107 110 111];
ttcol = [0 0.5 0; 0.7 0.4 0; 1 0.7 0.3; 0.3 1 0.3];

neuoi = ICsigagg.(ICexptidn).C.ICencoder2==1;
neuoi = ICsigagg.(ICexptidn).C.indin1==1;
figure
hold all
for typi = 1:numel(tt2p)
ttoi = ICtrialtypes==tt2p(typi);
plot(psthtli, squeeze(mean(convn(psthavgagg.(ICexptidn).C(:,ttoi,neuoi), kergauss, 'same'),3)), 'linewidth', 1, 'color', ttcol(typi,:))
end

