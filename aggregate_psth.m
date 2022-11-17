addpath(genpath('H:\CODE\helperfunctions'))
addpath(genpath('H:\CODE\Analyze_OpenScope'))

% kerwinhalf = 2; kersigma = 1;
% kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
% kergauss = (kergauss/sum(kergauss));

% A-AM, B-PM, C-V1, D-LM, E-AL, F-RL
probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
visind = [6 5 1 2 4 3];

datadir = 'D:\OpenScopeData\000248\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name}; 
nwbsessions = nwbsessions(contains(nwbsessions, 'sub'));

Nsessions = numel(nwbsessions)-1;
probeneuronsagg = cell(size(probes));
neuctxagg = cell(size(probes));
neulocagg = cell(size(probes));
sesneuoiagg = cell(size(probes));

visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
    'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
vistrialtypes = struct();
vistrialrep = struct();
vistrialorder = struct();
psthavg_ctxagg = struct();
Ron_ctxagg = struct();
Roff_ctxagg = struct();
for b = 1:numel(visblocks)
    vistrialtypes.(visblocks{b})=[];
    vistrialrep.(visblocks{b})=[];
    vistrialorder.(visblocks{b})=[];
psthavg_ctxagg.(visblocks{b}) = cell(size(probes));
Ron_ctxagg.(visblocks{b}) = cell(size(probes));
Roff_ctxagg.(visblocks{b}) = cell(size(probes));
end

ICblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations'};
% 'Palpha','BKtt','BKttpair','BItt','BIttpair','BICREl1tt','BICREl2tt','BICREl1ttpair','BICREl2ttpair',
ICsigfields = {'SP_Ind','Pmww_Ind','SP_BK','sigmcBK','Pmww_BK', ...
    'SP_ICvsRC','Pmww_ICvsRC','SP_BICREl','sigmcBICREl1','sigmcBICREl2', ...
    'PkwBK','PmcBK', 'PkwBI','PmcBI', 'PkwBICREl1','PmcBICREl1', 'PkwBICREl2','PmcBICREl2', ...
    'ICencoder','RCencoder','inducerencoder','inducerresponsive', ...
    'indenc1','indenc2','indenc3','indenc4','indenc13','indenc14','indenc23','indenc24', ...
    'ICencoder1','RCencoder1','RCencoder2','ICencoder2','indin1','indin2','indin3','indin4', ...
    'indout1','indout2','indout3','indout4','RElfaith1','RElfaith2', ...
    'ICresp1','RCresp1','RCresp2','ICresp2','ICtuned1','RCtuned1','RCtuned2','ICtuned2'};
ICsigagg = struct();
for iprobe = 1:numel(probes)
    for b = 1:numel(ICblocks)
        for f= 1:numel(ICsigfields)
            ICsigagg.(ICblocks{b})(iprobe).(ICsigfields{f}) = [];
        end
    end
end

RFCIfields = {'Rrfclassic','Rrfinverse','RFindclassic','RFindinverse', ...
    'Pkw_rfclassic','Pkw_rfinverse','pRrfclassic','pRrfinverse','pRFclassic','pRFinverse'};
RFCIagg = struct();
for iprobe = 1:numel(probes)
for f= 1:numel(RFCIfields)
RFCIagg(iprobe).(RFCIfields{f}) = [];
end
end

sizeCIfields = {'Rsizeclassic','Rsizeinverse','sizeindclassic','sizeindinverse','Pkw_sizeclassic','Pkw_sizeinverse'};
sizeCIagg = struct();
for iprobe = 1:numel(probes)
for f= 1:numel(sizeCIfields)
sizeCIagg(iprobe).(sizeCIfields{f}) = [];
end
end

oriparamsfields = {'Rori','prefiori','orthiori','OSI','OP','Pmww_OP'};
oriparamsagg = struct();
for iprobe = 1:numel(probes)
for f= 1:numel(oriparamsfields)
oriparamsagg(iprobe).(oriparamsfields{f}) = [];
end
end

meanFRvecagg = cell(size(probes));
sponFRvecagg = cell(size(probes));

Nsessions=4;
for ises = 1:Nsessions
    fprintf('Session %d/%d %s\n', ises, Nsessions, nwbsessions{ises} )
pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur'

whichneuctx = 1; % 1: correct definition, 2: >=230, 3: topmost-250um
elecid = electrode_id+1;
revmapelecid = NaN(max(elecid),1);
revmapelecid(elecid) = 1:numel(elecid);

for iprobe = 1:numel(probes)
    tic
    load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}))
    % 'neuoind', 'vis', 'Tres', 'psthtli', 'psth'
    load(sprintf('%svisresponses_probe%s.mat', pathpp, probes{iprobe}))
    % 'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'ICsig', 'RFCI', 'sizeCI', 'oriparams'
    
    if ~isequal(unique(floor(unit_peakch(neuoind)/1000)), iprobe-1)
        error('check neuoind')
    end
    
    % check whether CCF registration is correct
    probelocs = electrode_location(ismember(electrode_id, unit_peakch(neuoind)));
    
    neuloc = electrode_location(revmapelecid(unit_peakch(neuoind)+1));
    if ~isequal(unique(probelocs), unique(neuloc))
        disp(unique(neuloc)')
        error('check neuloc')
    end
    
    switch whichneuctx
        case 1
            neuctx = contains(neuloc, 'VIS');
        case 2
            % assume that >=230 is neocortex
            neuctx = mod(unit_peakch(neuoind), 1000)>=230;
            warning('for now, assuming that electrode_localid>=230 is neocortex (until CCF registration is corrected)')
            % 2 electrodes per depth, vertical spacing is 20um
            % typical electrode span in cortex 120 : (120/2)*20 = ~1200um
        case 3
            % to approximately get layer 2/3 cells, take units within 250um of the topmost unit
            neuctx = mod(unit_peakch(neuoind), 1000)>=max(mod(unit_peakch(neuoind), 1000))-25;
            warning('for now, to approximately get layer 2/3 cells, take units within 250um of the topmost unit')
    end

    fprintf('Probe %s Area %s: %d/%d\n', probes{iprobe}, visareas{iprobe}, nnz(neuctx), numel(neuoind) )
    disp(unique(probelocs)')
    
    probeneuronsagg{iprobe} = cat(1, probeneuronsagg{iprobe}, neuoind);
    neuctxagg{iprobe} = cat(1, neuctxagg{iprobe}, neuctx);
    neulocagg{iprobe} = cat(1, neulocagg{iprobe}, neuloc);
    
    % aggregate cortex neurons
    neuoi = neuctx;
    sesneuoiagg{iprobe} = cat(1, sesneuoiagg{iprobe}, ises*ones(nnz(neuoi),1));

meanFRvecagg{iprobe} = cat(1, meanFRvecagg{iprobe}, meanFRvec(neuoi)');
sponFRvecagg{iprobe} = cat(1, sponFRvecagg{iprobe}, sponFRvec(neuoi)');

for b = 1:numel(ICblocks)
    for f= 1:numel(ICsigfields)
        ICsigagg.(ICblocks{b})(iprobe).(ICsigfields{f}) = cat(1, ICsigagg.(ICblocks{b})(iprobe).(ICsigfields{f}), ICsig.(ICblocks{b}).(ICsigfields{f})(neuoi,:) );
    end
end

% RFCIfields = fieldnames(RFCI);
for f= 1:numel(RFCIfields)
RFCIagg(iprobe).(RFCIfields{f}) = cat(1, RFCIagg(iprobe).(RFCIfields{f}), RFCI.(RFCIfields{f})(neuoi,:) );
end

for f= 1:numel(sizeCIfields)
sizeCIagg(iprobe).(sizeCIfields{f}) = cat(1, sizeCIagg(iprobe).(sizeCIfields{f}), sizeCI.(sizeCIfields{f})(neuoi,:) );
end

for f= 1:numel(oriparamsfields)
oriparamsagg(iprobe).(oriparamsfields{f}) = cat(1, oriparamsagg(iprobe).(oriparamsfields{f}), oriparams.(oriparamsfields{f})(neuoi,:) );
end


% visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
%     'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
%ICblocks: each stim is 0.4s, inter-trial interval is 0.4s, static images
%RFCI: each stim is 0.25s, inter-trial interval is 0s, spinning drifting grating
%sizeCI: each stim is 0.25s, inter-trial interval is 0.5s, drifting grating
% orientation denotation is same as psychtoolbox (0 is 12h, 45 is 1h30m, clockwise)
for b = 1:numel(visblocks)
    vistrialtypes.(visblocks{b}) = unique(vis.(visblocks{b}).trialorder);
    vistrialorder.(visblocks{b}) = vis.(visblocks{b}).trialorder;
    if strcmp(visblocks{b}, 'RFCI_presentations')
        vistrialtypes.(visblocks{b}) = unique(vis.(visblocks{b}).trialorder(1:4:end));
    vistrialorder.(visblocks{b}) = vis.(visblocks{b}).trialorder(1:4:end);
    end
    Ntt = numel(vistrialtypes.(visblocks{b}));
    vistrialrep.(visblocks{b}) = zeros(Ntt,1);
    temppsth = zeros(length(psthtli), Ntt, nnz(neuoi));
    for ii = 1:Ntt
        trialsoi = vis.(visblocks{b}).trialorder==vistrialtypes.(visblocks{b})(ii);
        vistrialrep.(visblocks{b})(ii) = nnz(trialsoi);
        temppsth(:,ii,:) = mean(1000*psth.(visblocks{b})(:,trialsoi,neuoi), 2);
    end
psthavg_ctxagg.(visblocks{b}){iprobe} = cat(3, psthavg_ctxagg.(visblocks{b}){iprobe}, temppsth );
    
    if ismember(visblocks{b}, ICblocks)
        tlon = psthtli>0 & psthtli<=400;
        tloff = psthtli>400 & psthtli<=800;
    elseif strcmp(visblocks{b}, 'RFCI_presentations')
        tlon = psthtli>0 & psthtli<=1000;
        tloff = psthtli>1000 & psthtli<=1000;
    elseif strcmp(visblocks{b}, 'sizeCI_presentations')
        tlon = psthtli>0 & psthtli<=250;
        tloff = psthtli>250 & psthtli<=750;
    end
    tempRon = squeeze(mean(1000*psth.(visblocks{b})(tlon,:,neuoi), 1));
    tempRoff = squeeze(mean(1000*psth.(visblocks{b})(tloff,:,neuoi), 1));
Ron_ctxagg.(visblocks{b}){iprobe} = cat(2, Ron_ctxagg.(visblocks{b}){iprobe}, tempRon );
Roff_ctxagg.(visblocks{b}){iprobe} = cat(2, Roff_ctxagg.(visblocks{b}){iprobe}, tempRoff );
end

% size vector [0, 4, 8, 16, 32, 64 ]
toc
end

end

% %%
save('G:\My Drive\DATA\OpenScope\psth_ctxagg.mat', 'probes', 'visareas', 'visind', 'nwbsessions', ...
    'probeneuronsagg', 'neuctxagg', 'neulocagg', 'sesneuoiagg', 'RFCIagg', 'sizeCIagg', ...
    'oriparamsagg', 'ICsigagg', 'meanFRvecagg', 'sponFRvecagg', 'vistrialtypes', 'vistrialrep', 'vistrialorder', ...
    'ICblocks', 'ICtrialtypes', 'psthtli', 'psthavg_ctxagg', 'Ron_ctxagg', 'Roff_ctxagg', '-v7.3')

%% cosine similarity between REt-IC vs REt-RC
pathpp = [datadir nwbsessions{1} filesep];
load(sprintf('%sspikelatency_probe%s.mat', pathpp, probes{3}))

cossimfields = {'IC1REt1', 'RC1REt1', 'RC2REt2', 'IC2REt2'};
cossimacc = cell(size(probes));
cossimagg = cell(size(probes));
zcossimacc = cell(size(probes));
zcossimagg = cell(size(probes));
bscossimacc = cell(size(probes));
bscossimagg = cell(size(probes));
bszcossimacc = cell(size(probes));
bszcossimagg = cell(size(probes));
tloi = psthtli>0 & psthtli<=400;
for iprobe = 1:numel(probes)
    cossimacc{iprobe} = NaN(numel(ICblocks), numel(cossimfields));
    cossimagg{iprobe} = NaN(numel(ICblocks), numel(cossimfields), Nsessions);
    zcossimacc{iprobe} = NaN(numel(ICblocks), numel(cossimfields));
    zcossimagg{iprobe} = NaN(numel(ICblocks), numel(cossimfields), Nsessions);
    bscossimacc{iprobe} = NaN(numel(ICblocks), numel(cossimfields));
    bscossimagg{iprobe} = NaN(numel(ICblocks), numel(cossimfields), Nsessions);
    bszcossimacc{iprobe} = NaN(numel(ICblocks), numel(cossimfields));
    bszcossimagg{iprobe} = NaN(numel(ICblocks), numel(cossimfields), Nsessions);
    for b = 1:numel(ICblocks)
        for c=1:numel(cossimfields)
            switch cossimfields{c}
                case 'IC1REt1'
                    typi1 = ICtrialtypes==106;
                    typi2 = ICtrialtypes==1105;
                case 'RC1REt1'
                    typi1 = ICtrialtypes==107;
                    typi2 = ICtrialtypes==1105;
                case 'RC2REt2'
                    typi1 = ICtrialtypes==110;
                    typi2 = ICtrialtypes==1109;
                case 'IC2REt2'
                    typi1 = ICtrialtypes==111;
                    typi2 = ICtrialtypes==1109;
            end
            vec1 = squeeze(mean(psthavg_ctxagg.(ICblocks{b}){iprobe}(tloi,typi1,:),1));
            vec2 = squeeze(mean(psthavg_ctxagg.(ICblocks{b}){iprobe}(tloi,typi2,:),1));
            cossimacc{iprobe}(b,c) = dot(vec1, vec2)/(norm(vec1)*norm(vec2));
            
            tempRavg = squeeze(mean(Ron_ctxagg.(ICblocks{b}){iprobe},1));
            tempRstd = squeeze(std(Ron_ctxagg.(ICblocks{b}){iprobe},0,1));
            
            zvec1 = (vec1-tempRavg)./tempRstd;
            zvec2 = (vec2-tempRavg)./tempRstd;
            
            zvec1 = zvec1(tempRstd>0);
            zvec2 = zvec2(tempRstd>0);
            
            zcossimacc{iprobe}(b,c) = dot(zvec1, zvec2)/(norm(zvec1)*norm(zvec2));
            
            blankvec = squeeze(mean(psthavg_ctxagg.(ICblocks{b}){iprobe}(tloi,ICtrialtypes==0,:),1));
            vec1 = squeeze(mean(psthavg_ctxagg.(ICblocks{b}){iprobe}(tloi,typi1,:),1))-blankvec;
            vec2 = squeeze(mean(psthavg_ctxagg.(ICblocks{b}){iprobe}(tloi,typi2,:),1))-blankvec;
            bscossimacc{iprobe}(b,c) = dot(vec1, vec2)/(norm(vec1)*norm(vec2));
            
            zvec1 = (vec1-blankvec)./tempRstd;
            zvec2 = (vec2-blankvec)./tempRstd;
            zvec1 = zvec1(tempRstd>0);
            zvec2 = zvec2(tempRstd>0);
            bszcossimacc{iprobe}(b,c) = dot(zvec1, zvec2)/(norm(zvec1)*norm(zvec2));
            
            for ises = 1:Nsessions
                neuoi = sesneuoiagg{iprobe}==ises;
                vec1 = squeeze(mean(psthavg_ctxagg.(ICblocks{b}){iprobe}(tloi,typi1,neuoi),1));
                vec2 = squeeze(mean(psthavg_ctxagg.(ICblocks{b}){iprobe}(tloi,typi2,neuoi),1));
                cossimagg{iprobe}(b,c,ises) = dot(vec1, vec2)/(norm(vec1)*norm(vec2));
                
                %Rmat = squeeze(mean(psthavg_ctxagg.(ICblocks{b}){iprobe}(tloi,:,neuoi),1));
                %Ravg = sum(Rmat.*vistrialrep.(ICblocks{b}),1)/sum(vistrialrep.(ICblocks{b}));
                %figure; plot(Ravg', squeeze(mean(Ron_ctxagg.(ICblocks{b}){iprobe}(:,:,neuoi),2)), 'o')
                tempRavg = squeeze(mean(Ron_ctxagg.(ICblocks{b}){iprobe}(:,neuoi),1));
                tempRstd = squeeze(std(Ron_ctxagg.(ICblocks{b}){iprobe}(:,neuoi),0,1));
                
                zvec1 = (vec1-tempRavg)./tempRstd;
                zvec2 = (vec2-tempRavg)./tempRstd;
                
                zvec1 = zvec1(tempRstd>0);
                zvec2 = zvec2(tempRstd>0);
                
                zcossimagg{iprobe}(b,c,ises) = dot(zvec1, zvec2)/(norm(zvec1)*norm(zvec2));
            
            blankvec = squeeze(mean(psthavg_ctxagg.(ICblocks{b}){iprobe}(tloi,ICtrialtypes==0,neuoi),1));
            vec1 = squeeze(mean(psthavg_ctxagg.(ICblocks{b}){iprobe}(tloi,typi1,neuoi),1))-blankvec;
            vec2 = squeeze(mean(psthavg_ctxagg.(ICblocks{b}){iprobe}(tloi,typi2,neuoi),1))-blankvec;
            bscossimagg{iprobe}(b,c, ises) = dot(vec1, vec2)/(norm(vec1)*norm(vec2));
            
            zvec1 = (vec1-blankvec)./tempRstd;
            zvec2 = (vec2-blankvec)./tempRstd;
            zvec1 = zvec1(tempRstd>0);
            zvec2 = zvec2(tempRstd>0);
            bszcossimagg{iprobe}(b,c, ises) = dot(zvec1, zvec2)/(norm(zvec1)*norm(zvec2));
            end
        end
    end
end

ICblocknames = {'ICkcfg0', 'ICkcfg1', 'ICwcfg0', 'ICwcfg1'};
figure
for c = 1:4
for iprobe= 1:numel(visareas)
    subplot(4,6, 6*(c-1)+visind(iprobe))
    switch c
        case 1
    imagesc(cossimacc{iprobe})
    ylabel('Cosine Similarity')
        case 2
    imagesc(cossimacc{iprobe})
    ylabel('z-Cosine Similarity')
        case 3
    imagesc(bscossimacc{iprobe})
    ylabel('bs-Cosine Similarity')
        case 4
    imagesc(bszcossimacc{iprobe})
    ylabel('bsz-Cosine Similarity')
    end
    set(gca, 'XTick', 1:numel(cossimfields), 'XTickLabel', cossimfields, ...
        'YTick', 1:numel(ICblocks), 'YTickLabel', ICblocknames)
    colorbar
    title(visareas{iprobe})
end
end

figure
for iprobe= 1:numel(visareas)
    subplot(2,3, visind(iprobe))
    hold all
    for b = 1:numel(ICblocks)
        plot(cossimacc{iprobe}(b,:))
    end
    set(gca, 'XTick', 1:numel(cossimfields), 'XTickLabel', cossimfields)
    title(visareas{iprobe})
end

figure
for iprobe= 1:numel(visareas)
    subplot(2,3, visind(iprobe))
    imagesc(squeeze(nanmean(cossimagg{iprobe},3)))
    set(gca, 'XTick', 1:numel(cossimfields), 'XTickLabel', cossimfields, ...
        'YTick', 1:numel(ICblocks), 'YTickLabel', ICblocknames)
    colorbar
    title(visareas{iprobe})
end
% %%
figure
for b = 1:numel(ICblocks)
for iprobe= 1:numel(visareas)
    subplot(4,6, 6*(b-1)+visind(iprobe))
    hold all
    for ises = 1:Nsessions
    plot(squeeze(cossimagg{iprobe}(b,:,ises)))
    end
    plot(squeeze(nanmean(cossimagg{iprobe}(b,:,:), 3)), 'ko-', 'LineWidth', 2)
    plot(cossimacc{iprobe}(b,:), 'bo-', 'LineWidth', 2)
    set(gca, 'XTick', 1:numel(cossimfields), 'XTickLabel', cossimfields)
    ylabel('Cosine Similarity')
    title([visareas{iprobe} ' ' ICblocknames{b}])
end
end

figure
for b = 1:numel(ICblocks)
for iprobe= 1:numel(visareas)
    subplot(4,6, 6*(b-1)+visind(iprobe))
    hold all
    for ises = 1:Nsessions
    plot(squeeze(zcossimagg{iprobe}(b,:,ises)))
    end
    plot(squeeze(nanmean(zcossimagg{iprobe}(b,:,:), 3)), 'ko-', 'LineWidth', 2)
    plot(zcossimacc{iprobe}(b,:), 'bo-', 'LineWidth', 2)
    set(gca, 'XTick', 1:numel(cossimfields), 'XTickLabel', cossimfields)
    ylabel('z-Cosine Similarity')
    title([visareas{iprobe} ' ' ICblocknames{b}])
end
end

figure
for b = 1:numel(ICblocks)
for iprobe= 1:numel(visareas)
    subplot(4,6, 6*(b-1)+visind(iprobe))
    hold all
    for ises = 1:Nsessions
    plot(squeeze(bscossimagg{iprobe}(b,:,ises)))
    end
    plot(squeeze(nanmean(bscossimagg{iprobe}(b,:,:), 3)), 'ko-', 'LineWidth', 2)
    plot(bscossimacc{iprobe}(b,:), 'bo-', 'LineWidth', 2)
    set(gca, 'XTick', 1:numel(cossimfields), 'XTickLabel', cossimfields)
    ylabel('bs-Cosine Similarity')
    title([visareas{iprobe} ' ' ICblocknames{b}])
end
end

figure
for b = 1:numel(ICblocks)
for iprobe= 1:numel(visareas)
    subplot(4,6, 6*(b-1)+visind(iprobe))
    hold all
    for ises = 1:Nsessions
    plot(squeeze(bszcossimagg{iprobe}(b,:,ises)))
    end
    plot(squeeze(nanmean(bszcossimagg{iprobe}(b,:,:), 3)), 'ko-', 'LineWidth', 2)
    plot(bszcossimacc{iprobe}(b,:), 'bo-', 'LineWidth', 2)
    set(gca, 'XTick', 1:numel(cossimfields), 'XTickLabel', cossimfields)
    ylabel('bsz-Cosine Similarity')
    title([visareas{iprobe} ' ' ICblocknames{b}])
end
end

%% for grant
%sizeCI: each stim is 0.25s, inter-trial interval is 0.5s, drifting grating

addpath(genpath('H:\CODE\helperfunctions'))

% kerwinhalf = 2; kersigma = 1;
kerwinhalf = 12; kersigma = 5;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));


% plot all IC blocks
tt2p = [106 107 110 111 1105 1109]';
% tt2p = [106 107 110 111]';
ttcol = [0 0.5 0;
    0.75 0.5 0;
    1 0.75 0.25;
    0.25 0.75 0.25;
    0 0 0.5;
    0.25 0.25 0.75];
ICblockorder = {'ICwcfg1', 'ICwcfg0', 'ICkcfg1', 'ICkcfg0'};
fs=14;
figure('Position', [100 0 1800 1200])
% annotation('textbox', [0.1 0.9 0.9 0.1], 'string', 'Neuropixels: center-CRF neurons', 'edgecolor', 'none', 'fontsize', fs)
for b = 1:numel(ICblockorder)
    whichvisblock = [ICblockorder{b} '_presentations'];
if ~isequal(length(vistrialtypes.(whichvisblock)), length(ICtrialtypes))
    error('check ICtrialtypes')
end
for iprobe= 1:numel(visareas)
    subplot(numel(ICblockorder),numel(visareas), numel(visareas)*(b-1)+visind(iprobe))
    neuoi = sesneuoiagg{iprobe}<=Nsessions & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
    hold all
    for ii = 1:numel(tt2p)
        typi = ICtrialtypes==tt2p(ii);
        temppsth = convn(squeeze(psthavg_ctxagg.(whichvisblock){iprobe}(:,typi,neuoi)), kergauss, 'same');
        plot(psthtli/1000, squeeze(mean(temppsth,2)), '-', 'Color', ttcol(ii,:), 'LineWidth', 0.5)
        %shadedErrorBar(psthtli/1000, squeeze(mean(temppsth,2)), squeeze(std(temppsth,0,2)/sqrt(nnz(neuoi))), {'-', 'Color', ttcol(ii,:), 'LineWidth', 1}, 1)
    end
    xlim([-.05 .100])
%     xlim([-.1 .400])
%     ylim([0 25])
    set(gca, 'FontSize', fs)
    xlabel('Time (s)')
    ylabel('Rate (Hz)')
    title(sprintf('%s %s N=%d', ICblockorder{b}, visareas{iprobe}, nnz(neuoi)))
end
end

%% plot each IC blocks
whichvisblock = 'ICwcfg1_presentations';

tt2p = [106 107 110 111 1105 1109]';
% tt2p = [106 107 110 111]';
ttcol = [0 0.5 0;
    0.75 0.5 0;
    1 0.75 0.25;
    0.25 0.75 0.25;
    0 0 0.5;
    0.25 0.25 0.75];
fs=14;
figure%('Position', [100 0 1800 1200])
annotation('textbox', [0.1 0.92 0.9 0.1], 'string', ['Neuropixels: center-CRF neurons ' whichvisblock], 'edgecolor', 'none', 'fontsize', fs)
for iprobe= 1:numel(visareas)
    subplot(2,3, visind(iprobe))
    neuoi = sesneuoiagg{iprobe}<=Nsessions & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
    hold all
    for ii = 1:numel(tt2p)
        typi = ICtrialtypes==tt2p(ii);
        temppsth = convn(squeeze(psthavg_ctxagg.(whichvisblock){iprobe}(:,typi,neuoi)), kergauss, 'same');
        plot(psthtli/1000, squeeze(mean(temppsth,2)), '-', 'Color', ttcol(ii,:), 'LineWidth', 0.5)
        %shadedErrorBar(psthtli/1000, squeeze(mean(temppsth,2)), squeeze(std(temppsth,0,2)/sqrt(nnz(neuoi))), {'-', 'Color', ttcol(ii,:), 'LineWidth', 1}, 1)
    end
    xlim([-.05 .100])
%     ylim([0 25])
    set(gca, 'FontSize', fs)
    xlabel('Time (s)')
    ylabel('Rate (Hz)')
    title(sprintf('%s %s N=%d', whichvisblock, visareas{iprobe}, nnz(neuoi)))
end

%%
szvec = [0, 4, 8, 16, 32, 64];

figure('Position', [100 100 1800 300])
% annotation('textbox', [0.1 0.9 0.9 0.1], 'string', 'Neuropixels: center-CRF neurons', 'edgecolor', 'none', 'fontsize', fs)
for iprobe= 1:6
subplot(1,6, visind(iprobe))
hold all
neuoi = sesneuoiagg{iprobe}<=Nsessions & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
errorbar(1:length(szvec), mean(sizeCIagg(iprobe).Rsizeclassic(neuoi,:), 1), std(sizeCIagg(iprobe).Rsizeclassic(neuoi,:), 0,1)/sqrt(nnz(neuoi)), 'ko-', 'MarkerFaceColor', 'k', 'LineWidth', 2);
errorbar(1:length(szvec), mean(sizeCIagg(iprobe).Rsizeinverse(neuoi,:), 1), std(sizeCIagg(iprobe).Rsizeinverse(neuoi,:), 0,1)/sqrt(nnz(neuoi)), 'co-', 'MarkerFaceColor', 'c', 'LineWidth', 2);
set(gca, 'FontSize', fs, 'Xtick', 1:length(szvec), 'xticklabel', szvec)
xlim([1 length(szvec)])
    xlabel('Size (deg)')
    ylabel('Rate (Hz)')
    title(sprintf('%s N=%d', visareas{iprobe}, nnz(neuoi)))
end


%% response to full/classic/inverse screen gratings across areas (compare latency)
%sizeCI: each stim is 0.25s, inter-trial interval is 0.5s, drifting grating
% baseline period: -100 to 0ms (alternatively -200 to 0ms would work)
%szvec = [0, 4, 8, 16, 32, 64];

kerwinhalf = 2; kersigma = 1;
kerwinhalf = 4; kersigma = 2;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));

[~,probehierarchy]=sort(visind);
figure %('Position', [100 100 1000 400])
% annotation('textbox', [0.05 0.9 0.8 0.1], 'string', [nwbsessions{ises} ' All Neurons'], 'FontSize', 14, 'EdgeColor', 'none')
for pltopt = 1:3 % 1 for full grating all neurons, 2 for center grating center neurons, 3 for inverse grating center neurons
switch pltopt
    case 1
        figtitle = 'All Neurons Full Screen Gratings';
    case 2
        figtitle = 'center-CRF Neurons Center Classical Gratings';
    case 3
        figtitle = 'center-CRF Neurons Center Inverse Gratings';
end

subplot(2,3,pltopt)
hold all
for a = 1:numel(probes)
    iprobe = probehierarchy(a);
    
    switch pltopt
        case 1
            neuoi = sesneuoiagg{iprobe}<=Nsessions;
            ttoi = floor(vistrialtypes.sizeCI_presentations/1000)==11;
        case 2
            neuoi = sesneuoiagg{iprobe}<=Nsessions & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
            ttoi = floor(vistrialtypes.sizeCI_presentations/1000)==4;
        case 3
            neuoi = sesneuoiagg{iprobe}<=Nsessions & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
            ttoi = floor(vistrialtypes.sizeCI_presentations/1000)==14;
    end
    
    temp = psthavg_ctxagg.sizeCI_presentations{iprobe}(:,ttoi,neuoi);
    % plot(psthtli, smooth(squeeze(mean(temp, [2,3])), 5), '-')
    temp = convn(temp, kergauss, 'same');
    plot(psthtli, squeeze(mean(temp, [2,3])))
end
legend(visareas{probehierarchy})
xlim([-50 200])
% xlim([-500 750])
set(gca, 'XGrid', 'on')
title(figtitle)
% title(sprintf('%s %s neurons N=%d/%d\nsizeCI_presentations', whicharea, neudesc, nnz(neuoi), nnz(neuctx)), 'interpreter', 'none')
xlabel('Time (ms)')
ylabel('Rate (Hz)')

subplot(2,3,3+pltopt)
hold all
for a = 1:numel(probes)
    iprobe = probehierarchy(a);
    
    switch pltopt
        case 1
    neuoi = sesneuoiagg{iprobe}<=Nsessions;
    ttoi = floor(vistrialtypes.sizeCI_presentations/1000)==11;
        case 2
    neuoi = sesneuoiagg{iprobe}<=Nsessions & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
    ttoi = floor(vistrialtypes.sizeCI_presentations/1000)==4;
        case 3
    neuoi = sesneuoiagg{iprobe}<=Nsessions & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
    ttoi = floor(vistrialtypes.sizeCI_presentations/1000)==14;
    end
    temp = psthavg_ctxagg.sizeCI_presentations{iprobe}(:,ttoi,neuoi);
    % plot(psthtli, smooth(squeeze(mean(temp, [2,3])), 5), '-')
    temp = convn(temp, kergauss, 'same');
    tempbase = squeeze(mean(temp(psthtli>=-100 & psthtli<0,:,:), 'all'));
    temppeak = max(squeeze(mean(temp, [2,3])));
    plot(psthtli, (squeeze(mean(temp, [2,3]))-tempbase)/(temppeak-tempbase) )
end
legend(visareas{probehierarchy})
xlim([-50 200])
ylim([-0.2 1])
% xlim([-500 750])
set(gca, 'XGrid', 'on')
xlabel('Time (ms)')
ylabel('Normalized Activity')
end

%% response to IC vs RE across areas (compare latency) 
%sizeCI: each stim is 0.25s, inter-trial interval is 0.5s, drifting grating
% baseline period: -100 to 0ms (alternatively -200 to 0ms would work)
%szvec = [0, 4, 8, 16, 32, 64];

kerwinhalf = 2; kersigma = 1;
kerwinhalf = 4; kersigma = 2;
% kerwinhalf = 7; kersigma = 3;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));

whichvisblock = 'ICwcfg1_presentations';
neuopt = 0;
switch neuopt
    case 0
        neutitle = 'non-center-CRF neurons';
    case 1
        neutitle = 'center-CRF neurons';
    case 2
        neutitle = 'IC-encoders';
end

lw = 1;
[~,probehierarchy]=sort(visind);
figure %('Position', [100 100 1000 400])
% annotation('textbox', [0.05 0.9 0.8 0.1], 'string', [nwbsessions{ises} ' All Neurons'], 'FontSize', 14, 'EdgeColor', 'none')
for pltopt = 1:2 % 1 real edges, 2 illusory contours
switch pltopt
    case 1
        figtitle = [neutitle ' Real Edge'];
    case 2
%         figtitle = [neutitle ' REt'];
%     case 3
        figtitle = [neutitle ' Illusory Contour'];
end

subplot(2,2,pltopt)
%subplot(2,3,pltopt)
hold all
for a = 1:2%numel(probes)
    iprobe = probehierarchy(a);
    
    switch neuopt
        case 0
            neuoi = sesneuoiagg{iprobe}<=Nsessions & RFCIagg(iprobe).RFindclassic~=1 & RFCIagg(iprobe).pRFclassic<0.05;
        case 1
            neuoi = sesneuoiagg{iprobe}<=Nsessions & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
        case 2
            neuoi = sesneuoiagg{iprobe}<=Nsessions & ICsigagg.(whichvisblock)(iprobe).ICencoder1==1;
    end
    switch pltopt
        case 1
            ttoi = ismember(ICtrialtypes, [506 511]);
        case 2
%             ttoi = ismember(ICtrialtypes, [1105 1109]);
%         case 3
            ttoi = ismember(ICtrialtypes, [106 111]);
    end
    
    temp = psthavg_ctxagg.(whichvisblock){iprobe}(:,ttoi,neuoi);
    % plot(psthtli, smooth(squeeze(mean(temp, [2,3])), 5), '-')
    temp = convn(temp, kergauss, 'same');
    plot(psthtli, squeeze(mean(temp, [2,3])), 'LineWidth', lw)
%     temppsth = squeeze(mean(temp, 2));
%     shadedErrorBar(psthtli, squeeze(mean(temppsth,2)), squeeze(std(temppsth,0,2)/sqrt(nnz(neuoi))), {'LineWidth', lw}, 1)
end
legend(visareas{probehierarchy})
% xlim([0 400])
% ylim ([0 12])
xlim([-50 200])
set(gca, 'XGrid', 'on')
title(figtitle)
% title(sprintf('%s %s neurons N=%d/%d\nsizeCI_presentations', whicharea, neudesc, nnz(neuoi), nnz(neuctx)), 'interpreter', 'none')
xlabel('Time (ms)')
ylabel('Rate (Hz)')

subplot(2,2,2+pltopt)
%subplot(2,3,3+pltopt)
hold all
for a = 1:2%numel(probes)
    iprobe = probehierarchy(a);
    
    switch neuopt
        case 0
            neuoi = sesneuoiagg{iprobe}<=Nsessions & RFCIagg(iprobe).RFindclassic~=1 & RFCIagg(iprobe).pRFclassic<0.05;
        case 1
            neuoi = sesneuoiagg{iprobe}<=Nsessions & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
        case 2
            neuoi = sesneuoiagg{iprobe}<=Nsessions & ICsigagg.(whichvisblock)(iprobe).ICencoder1==1;
    end
    switch pltopt
        case 1
            ttoi = ismember(ICtrialtypes, [506 511]);
        case 2
%             ttoi = ismember(ICtrialtypes, [1105 1109]);
%         case 3
            ttoi = ismember(ICtrialtypes, [106 111]);
    end
    
    temp = psthavg_ctxagg.(whichvisblock){iprobe}(:,ttoi,neuoi);
    % plot(psthtli, smooth(squeeze(mean(temp, [2,3])), 5), '-')
    temp = convn(temp, kergauss, 'same');
    tempbase = squeeze(mean(temp(psthtli>=-200 & psthtli<0,:,:), 'all'));
    temppeak = max(squeeze(mean(temp(psthtli>=0 & psthtli<400,:,:), [2,3])));
    plot(psthtli, (squeeze(mean(temp, [2,3]))-tempbase)/(temppeak-tempbase), 'LineWidth', lw )
end
legend(visareas{probehierarchy})
xlim([-50 200])
ylim([-0.2 1])
% xlim([-500 750])
set(gca, 'XGrid', 'on')
xlabel('Time (ms)')
ylabel('Normalized Activity')
end

%% plot each vis area
whichvisblock = 'ICwcfg1_presentations';

kerwinhalf = 2; kersigma = 1;
% kerwinhalf = 4; kersigma = 2;
% kerwinhalf = 7; kersigma = 3;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));

neuopt = 0;
switch neuopt
    case 0
        neutitle = 'non-center-CRF neurons';
    case 1
        neutitle = 'center-CRF neurons';
    case 2
        neutitle = 'IC-encoders';
end

fs=14;
figure%('Position', [100 0 1800 1200])
annotation('textbox', [0.1 0.9 0.9 0.1], 'string', ['Neuropixels: ' neutitle ' ' whichvisblock], 'edgecolor', 'none', 'fontsize', fs, 'interpreter', 'none')
for iprobe= 1:numel(visareas)
    subplot(2,3, visind(iprobe))
    switch neuopt
        case 0
            neuoi = sesneuoiagg{iprobe}<=Nsessions & RFCIagg(iprobe).RFindclassic~=1 & RFCIagg(iprobe).pRFclassic<0.05;
        case 1
            neuoi = sesneuoiagg{iprobe}<=Nsessions & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
        case 2
            neuoi = sesneuoiagg{iprobe}<=Nsessions & ICsigagg.(whichvisblock)(iprobe).ICencoder1==1;
    end

    hold all
    for typi = 1:2
        switch typi
            case 1
        ttoi = ismember(ICtrialtypes, [506 511]);
        ttcol = [0 0 1];
            case 2
        ttoi = ismember(ICtrialtypes, [106 111]);
        ttcol = [0 0.7 0];
        end
        temppsth = convn(squeeze(psthavg_ctxagg.(whichvisblock){iprobe}(:,ttoi,neuoi)), kergauss, 'same');
        plot(psthtli/1000, squeeze(mean(temppsth,[2 3])), '-', 'Color', ttcol, 'LineWidth', 0.5)
        %shadedErrorBar(psthtli/1000, squeeze(mean(temppsth,2)), squeeze(std(temppsth,0,2)/sqrt(nnz(neuoi))), {'-', 'Color', ttcol, 'LineWidth', 1}, 1)
    end
    xlim([-.05 .100])
%     ylim([0 25])
    set(gca, 'FontSize', fs)
    xlabel('Time (s)')
    ylabel('Rate (Hz)')
    title(sprintf('%s N=%d', visareas{iprobe}, nnz(neuoi)), 'interpreter', 'none')
end

%% plot each vis area
whichvisblock = 'ICwcfg1_presentations';
neuopt = 1;

kerwinhalf = 2; kersigma = 1;
% kerwinhalf = 4; kersigma = 2;
% kerwinhalf = 7; kersigma = 3;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));

fs=14;
figure%('Position', [100 0 1800 1200])
annotation('textbox', [0.1 0.9 0.9 0.1], 'string', ['Neuropixels: ' neutitle ' ' whichvisblock], 'edgecolor', 'none', 'fontsize', fs, 'interpreter', 'none')
for iprobe= 1:numel(visareas)
    subplot(2,3, visind(iprobe))
    switch neuopt
        case 0
            neuoi = sesneuoiagg{iprobe}<=Nsessions & RFCIagg(iprobe).RFindclassic~=1 & RFCIagg(iprobe).pRFclassic<0.05;
        case 1
            neuoi = sesneuoiagg{iprobe}<=Nsessions & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
        case 2
            neuoi = sesneuoiagg{iprobe}<=Nsessions & ICsigagg.(whichvisblock)(iprobe).ICencoder1==1;
    end
    hold all
    for typi = 1:2
        switch typi
            case 1
        ttoi = ismember(ICtrialtypes, [506 511]);
        ttcol = [0 0 1];
            case 2
        ttoi = ismember(ICtrialtypes, [106 111]);
        ttcol = [0 0.7 0];
        end
        temp = convn(squeeze(psthavg_ctxagg.(whichvisblock){iprobe}(:,ttoi,neuoi)), kergauss, 'same');
    tempbase = squeeze(mean(temp(psthtli>=-200 & psthtli<0,:,:), 'all'));
    temppeak = max(squeeze(mean(temp(psthtli>=0 & psthtli<400,:,:), [2,3])));
    plot(psthtli/1000, (squeeze(mean(temp, [2,3]))-tempbase)/(temppeak-tempbase), '-', 'Color', ttcol, 'LineWidth', 1)
        %shadedErrorBar(psthtli/1000, squeeze(mean(temppsth,2)), squeeze(std(temppsth,0,2)/sqrt(nnz(neuoi))), {'-', 'Color', ttcol, 'LineWidth', 1}, 1)
    end
    xlim([-.05 .100])
%     ylim([0 25])
    set(gca, 'FontSize', fs)
    xlabel('Time (s)')
ylabel('Normalized Activity')
    title(sprintf('%s N=%d', visareas{iprobe}, nnz(neuoi)), 'interpreter', 'none')
end
