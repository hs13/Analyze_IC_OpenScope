if exist('/Volumes/GoogleDrive-116160770365018316196/My Drive/', 'dir')
    drivepath = '/Volumes/GoogleDrive-116160770365018316196/My Drive/';
elseif exist('G:/My Drive/', 'dir')
    drivepath = 'G:/My Drive/';
else
    error('check drivepath in this computer')
end
load([drivepath, 'DATA/ICexpts_submission22/openscope_popavg_agg.mat'])
load([drivepath, 'DATA/ICexpts_submission22/openscope_psthavgagg.mat'])

probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
visind = [6 5 1 2 4 3];

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
neuoi = neuctxagg{3}==1 & ICsigagg.(ICexptidn).C.(neudesc)==1;
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
neuoi = neuctxagg{3}==1 & ICsigagg.(ICexptidn).C.(neudesc)==1;
subplot(2,2,isp)
hold all
for typi = 1:numel(tt2p)
ttoi = ICtrialtypes==tt2p(typi);
plot(psthtli, squeeze(mean(convn(psthavgagg.(ICexptidn).C(:,ttoi,neuoi), kergauss, 'same'),3)), 'linewidth', 1, 'color', ttcol(typi,:))
end
title(sprintf('%s N=%d', neudesc, nnz(neuoi)))
end

%%
ICexptidn = 'ICwcfg1_presentations';
tt2p = [106 107 110 111 506 511];
ttcol = [0 0.5 0; 0.7 0.4 0; 1 0.7 0.3; 0.3 1 0.3; 0 0 0.5; 0.3 0.3 1];

neudesc = 'ctrCRF9sigexcl';
figure
for iprobe = 1:numel(probes)
    probeind = find( strcmp(probes{iprobe}, {'A', 'B', 'C', 'D', 'E', 'F'}) );
    switch neudesc
        case 'ctrCRF9'
            neuoi = RFCIagg.(probes{iprobe}).RFindclassic==1 & RFCIagg.(probes{iprobe}).Pkw_rfclassic<0.05;
        case 'ctrCRF9sigexcl'
            neuoi = RFCIagg.(probes{iprobe}).RFindclassic==1 & RFCIagg.(probes{iprobe}).RFsigexclclassic==1;
        case 'ctrCRF9exclsig'
            neuoi = RFCIagg.(probes{iprobe}).RFindclassic==1 & RFCIagg.(probes{iprobe}).RFexclsigclassic==1;
    end
neuoi = neuctxagg{probeind}==1 & neuoi;
subplot(2,3,visind(iprobe))
hold all
for typi = 1:numel(tt2p)
ttoi = ICtrialtypes==tt2p(typi);
plot(psthtli, squeeze(mean(convn(psthavgagg.(ICexptidn).(probes{iprobe})(:,ttoi,neuoi), kergauss, 'same'),3)), 'linewidth', 1, 'color', ttcol(typi,:))
end
xlim([-100 400])
title(sprintf('%s N=%d', visareas{iprobe}, nnz(neuoi)))
end
