probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
visind = [6 5 1 2 4 3];


Nrfs = size(RFCIagg(1).C.Rrfclassic, 2);
Nszs = size(sizeCIagg(1).C.Rsizeclassic, 2);
Ndirs = size(oriparamsagg(1).C.Rori, 2);
Noris = Ndirs/2;
dirvec = vis.sizeCI_presentations.directions;
if length(dirvec)~=Ndirs
    error('check sizeCI_presentations directions')
end
orivec = vis.sizeCI_presentations.directions(1:Noris);

kerwinhalf = 5; kersigma = 2;
kerwinhalf = 25; kersigma = 5;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));

%%
whichICblock = 'ICwcfg1_presentations';
tt2p = [106 107 110 111 0];
ttcol = [0 0.5 0; 0.7 0.4 0; 1 0.7 0.3; 0.3 1 0.3; 0 0 0];

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
neuoi = neuctxagg{3}==1 & ICsigagg.(whichICblock).C.(neudesc)==1;
subplot(2,2,isp)
hold all
for typi = 1:numel(tt2p)
ttoi = ICtrialtypes==tt2p(typi);
plot(psthtli, squeeze(mean(convn(psthavgagg.(whichICblock).C(:,ttoi,neuoi), kergauss, 'same'),3)), 'linewidth', 1, 'color', ttcol(typi,:))
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
neuoi = neuctxagg{3}==1 & ICsigagg.(whichICblock).C.(neudesc)==1;
subplot(2,2,isp)
hold all
for typi = 1:numel(tt2p)
ttoi = ICtrialtypes==tt2p(typi);
plot(psthtli, squeeze(mean(convn(psthavgagg.(whichICblock).C(:,ttoi,neuoi), kergauss, 'same'),3)), 'linewidth', 1, 'color', ttcol(typi,:))
end
title(sprintf('%s N=%d', neudesc, nnz(neuoi)))
end

%% center-RF neurons in different areas
whichICblock = 'ICwcfg1_presentations';
tt2p = [106 107 110 111 506 511 0];
ttcol = [0 0.5 0; 0.7 0.4 0; 1 0.7 0.3; 0.3 1 0.3; 0 0 0.5; 0.3 0.3 1; 0 0 0];

neudesc = 'ctrCG9';
figure
for iprobe = 1:numel(probes)
    probeind = find( strcmp(probes{iprobe}, {'A', 'B', 'C', 'D', 'E', 'F'}) );
    switch neudesc
        case 'ctrCG9'
            neuoi = RFCIagg.(probes{iprobe}).RFindclassic==1 & RFCIagg.(probes{iprobe}).pRFclassic<0.05;
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
plot(psthtli, squeeze(mean(convn(psthavgagg.(whichICblock).(probes{iprobe})(:,ttoi,neuoi), kergauss, 'same'),3)), 'linewidth', 1, 'color', ttcol(typi,:))
end
xlim([-100 400])
title(sprintf('%s N=%d', visareas{iprobe}, nnz(neuoi)))
end

%% center-RF neurons orientation tuning
whichprobe = 'C';
iprobe = find(strcmp({'A', 'B', 'C', 'D', 'E', 'F'}, whichprobe));
% neuinarea = strcmp(neulocagg{iprobe}, 'VISp2/3');
% neuinarea = contains(neulocagg{iprobe}, 'VISp');
neuinarea = neuctxagg{iprobe}==1;

neudesc = 'ctrCRF9exclsig';
switch neudesc
    case 'ctrCG9'
        neuoi = RFCIagg.(whichprobe).RFindclassic==1 & RFCIagg.(whichprobe).pRFclassic<0.05;
    case 'ctrCRF9'
        neuoi = RFCIagg.(whichprobe).RFindclassic==1 & RFCIagg.(whichprobe).Pkw_rfclassic<0.05;
    case 'ctrCRF9sigexcl' % Bonferroni-Holms correction
        neuoi = RFCIagg.(whichprobe).RFindclassic==1 & RFCIagg.(whichprobe).RFsigexclclassic==1;
    case 'ctrCRF9exclsig' % only one RF with p<0.05
        neuoi = RFCIagg.(whichprobe).RFindclassic==1 & RFCIagg.(whichprobe).RFexclsigclassic==1;
end

xl = [-100 500];
figure
for iori = 1:4
    %neu2plot = neuinarea & neuoi & ori4paramsagg.(whichprobe).prefiori4==iori & ori4paramsagg.(whichprobe).Pkw_ori4<0.05;
    neu2plot = neuinarea & neuoi & ori4paramsagg.(whichprobe).prefiori4==iori;
    
    subplot(2,4,iori)
    hold all
    for b = 0:1
        switch b
            case 0
                whichICblock = 'ICwcfg0_presentations';
                tt2p = [106 111 506 511];
                ttcol = [0 3/4 0; 0 1/4 0; 0 0 3/4; 0 0 1/4];
            case 1
                whichICblock = 'ICwcfg1_presentations';
                tt2p = [106 111 506 511];
                ttcol = [0 4/4 0; 0 2/4 0; 0 0 4/4; 0 0 2/4];
        end
        for typi = 1:numel(tt2p)
            ttoi = ICtrialtypes==tt2p(typi);
            plot(psthtli, squeeze(mean(convn(psthavgagg.(whichICblock).(whichprobe)(:,ttoi,neu2plot), kergauss, 'same'),3)), 'linewidth', 1, 'color', ttcol(typi,:))
        end
    end
    xlim(xl)
    title(sprintf('White Inducers: Pref %d-deg N=%d', orivec(iori), nnz(neu2plot)))
    
    subplot(2,4,4+iori)
    hold all
    for b = 0:1
        switch b
            case 0
                whichICblock = 'ICkcfg0_presentations';
                tt2p = [106 111 506 511];
                ttcol = [0 3/4 0; 0 1/4 0; 0 0 3/4; 0 0 1/4];
            case 1
                whichICblock = 'ICkcfg1_presentations';
                tt2p = [106 111 506 511];
                ttcol = [0 4/4 0; 0 2/4 0; 0 0 4/4; 0 0 2/4];
        end
        for typi = 1:numel(tt2p)
            ttoi = ICtrialtypes==tt2p(typi);
            plot(psthtli, squeeze(mean(convn(psthavgagg.(whichICblock).(whichprobe)(:,ttoi,neu2plot), kergauss, 'same'),3)), 'linewidth', 1, 'color', ttcol(typi,:))
        end
    end
    xlim(xl)
    title(sprintf('Black Inducers: Pref %d-deg N=%d', orivec(iori), nnz(neu2plot)))    
end


figure
for iori = 1:4
    %neu2plot = neuinarea & neuoi & ori4paramsagg.(whichprobe).prefiori4==iori & ori4paramsagg.(whichprobe).Pkw_ori4<0.05;
    neu2plot = neuinarea & neuoi & ori4paramsagg.(whichprobe).prefiori4==iori;
    
    for icol = 1:2
        switch icol
            case 1
                blockprefix = 'ICw';
                indcoldesc = 'White';
            case 2
                blockprefix = 'ICk';
                indcoldesc = 'Black';
        end
    subplot(4,4,8*(icol-1)+iori)
    hold all
    for b = 0:1
        switch b
            case 0
                whichICblock = [blockprefix 'cfg0_presentations'];
                tt2p = [106 111];
                ttcol = [0 1 1; 1 0 0];
            case 1
                whichICblock = [blockprefix 'cfg1_presentations'];
                tt2p = [106 111];
                ttcol = [1 0 1; 0 1 0];
        end
        for typi = 1:numel(tt2p)
            ttoi = ICtrialtypes==tt2p(typi);
            plot(psthtli, squeeze(mean(convn(psthavgagg.(whichICblock).(whichprobe)(:,ttoi,neu2plot), kergauss, 'same'),3)), 'linewidth', 1, 'color', ttcol(typi,:))
        end
    end
    xlim(xl)
    title(sprintf('%s Inducers IC: Pref %d-deg N=%d', indcoldesc, orivec(iori), nnz(neu2plot)))
    
    subplot(4,4,8*(icol-1)+4+iori)
    hold all
    for b = 0:1
        switch b
            case 0
                whichICblock = [blockprefix 'cfg0_presentations'];
                tt2p = [506 511];
                ttcol = [0 1 1; 1 0 0];
            case 1
                whichICblock = [blockprefix 'cfg1_presentations'];
                tt2p = [506 511];
                ttcol = [1 0 1; 0 1 0];
        end
        for typi = 1:numel(tt2p)
            ttoi = ICtrialtypes==tt2p(typi);
            plot(psthtli, squeeze(mean(convn(psthavgagg.(whichICblock).(whichprobe)(:,ttoi,neu2plot), kergauss, 'same'),3)), 'linewidth', 1, 'color', ttcol(typi,:))
        end
    end
    xlim(xl)
    title(sprintf('%s Inducers REl: Pref %d-deg N=%d', indcoldesc, orivec(iori), nnz(neu2plot)))
    end  
end

%% cosine similarity between all trial types
whichICblock = 'ICwcfg1_presentations';

cossimmat = NaN(numel(ICtrialtypes));
for ii = 1:numel(ICtrialtypes)
for jj = 1:numel(ICtrialtypes)
    ivec = Ronavgagg.(whichICblock).C(ii,:)';
    jvec = Ronavgagg.(whichICblock).C(jj,:)';
    cossimmat(ii,jj) = dot(ivec, jvec)/(norm(ivec)*norm(jvec));
end
end

figure; imagesc(cossimmat)
set(gca, 'XTick',1:numel(ICtrialtypes),'XTickLabel', ICtrialtypes, 'YTick',1:numel(ICtrialtypes),'YTickLabel', ICtrialtypes)
colormap redblue

figure; 
subplot(2,1,1)
plot(cossimmat(ICtrialtypes==506,:), 'o-')
set(gca, 'XTick',1:numel(ICtrialtypes),'XTickLabel', ICtrialtypes, 'XGrid', 'on');
subplot(2,1,2)
plot(cossimmat(ICtrialtypes==511,:), 'o-')
set(gca, 'XTick',1:numel(ICtrialtypes),'XTickLabel', ICtrialtypes, 'XGrid', 'on');
