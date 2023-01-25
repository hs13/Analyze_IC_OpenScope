% load('D:\OpenScopeData\000248\postprocessed\psthavgagg.mat')
load(['G:\My Drive\DATA\ICexpts_submission22\openscope_popavg_agg.mat'])
load(['G:\My Drive\DATA\ICexpts_submission22\openscope_psthavgagg.mat'])

kerwinhalf = 5; kersigma = 2;
kerwinhalf = 25; kersigma = 5;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));

ICexptidn = 'ICwcfg1_presentations';
tt2p = [106 107 110 111];
ttcol = [0 0.5 0; 0.7 0.4 0; 1 0.7 0.3; 0.3 1 0.3];

figure
for isp = 1:4
    switch isp
        case 1
            neudesc = 'ICencoder1';
        case 2
            neudesc = 'ICencoder2';
        case 3
            neudesc = 'RCencoder1';
        case 4
            neudesc = 'RCencoder2';
    end
neuoi = ICsigagg.(ICexptidn).C.(neudesc)==1;
subplot(2,2,isp)
hold all
for typi = 1:numel(tt2p)
ttoi = ICtrialtypes==tt2p(typi);
plot(psthtli, squeeze(mean(convn(psthavgagg.(ICexptidn).C(:,ttoi,neuoi), kergauss, 'same'),3)), 'linewidth', 1, 'color', ttcol(typi,:))
end
title(sprintf('%s N=%d', neudesc, nnz(neuoi)))
end

figure
for isp = 1:4
    switch isp
        case 1
            neudesc = 'indin1';
        case 2
            neudesc = 'indin2';
        case 3
            neudesc = 'indin3';
        case 4
            neudesc = 'indin4';
    end
neuoi = ICsigagg.(ICexptidn).C.(neudesc)==1;
subplot(2,2,isp)
hold all
for typi = 1:numel(tt2p)
ttoi = ICtrialtypes==tt2p(typi);
plot(psthtli, squeeze(mean(convn(psthavgagg.(ICexptidn).C(:,ttoi,neuoi), kergauss, 'same'),3)), 'linewidth', 1, 'color', ttcol(typi,:))
end
title(sprintf('%s N=%d', neudesc, nnz(neuoi)))
end

tt2p = [106 107 110 111 506 511];
ttcol = [0 0.5 0; 0.7 0.4 0; 1 0.7 0.3; 0.3 1 0.3];
