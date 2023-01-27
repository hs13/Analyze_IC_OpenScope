addpath(genpath('H:\CODE\helperfunctions'))
addpath(genpath('H:\CODE\Analyze_OpenScope'))

datadir = 'D:\OpenScopeData\000248\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name}; 
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

% A-AM, B-PM, C-V1, D-LM, E-AL, F-RL
probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
visind = [6 5 1 2 4 3];

% kerwinhalf = 2; kersigma = 1;
% kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
% kergauss = (kergauss/sum(kergauss));

%% 
probevisareas = cell(Nsessions,numel(probes));
for ises = 1:Nsessions
    fprintf('Session %d/%d %s\n', ises, Nsessions, nwbsessions{ises} )
    
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur'
    
    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    
    for iprobe = 1:numel(probes)
                
        if exist([pathpp, 'probes.mat'], 'file')
            probelist = load([pathpp, 'probes.mat']);
            warning('HS 230126: this was inserted to handle the exception case of sub_1183369803, can delete with the next nwb update')
        else
            probelist.probes = {'A', 'B', 'C', 'D', 'E', 'F'};
        end
        probeind = find( strcmp(probes{iprobe}, probelist.probes) );
        %probeind = find( strcmp(probes{iprobe}, {'A', 'B', 'C', 'D', 'E', 'F'}) );

%         if ~isequal(unique(floor(unit_peakch(neuoind)/1000)), probeind-1)
%             error('check neuoind')
%         end
        neuoind = find(floor(unit_peakch/1000)==probeind-1);
        
        % check whether CCF registration is correct
        probelocs = electrode_location(ismember(electrode_id, unit_peakch(neuoind)));
        
        neuloc = electrode_location(revmapelecid(unit_peakch(neuoind)+1));
        if ~isequal(unique(probelocs), unique(neuloc))
            disp(unique(neuloc)')
            error('check neuloc')
        end
%         disp(probes{iprobe})
%         disp(unique(neuloc(contains(neuloc, 'VIS')))')
        probevisareas{ises, iprobe} = sprintf('%s ',unique(neuloc(contains(neuloc, 'VIS'))));
    end
end

disp('Sessions with wrong CCF labels (all probes labeled as VISl)')
disp(nwbsessions(all(contains(probevisareas, 'VISl'),2))')

disp('Sessions with no layer information')
disp(nwbsessions( ~all(contains(probevisareas, 'VISl'),2) & ~any(contains(probevisareas, '2'),2) )')

disp('Sessions with seemingly correct CCF labels & layer info')
disp(nwbsessions(any(contains(probevisareas, '2'),2))')

% strcmp(probevisareas, 'VISl')

%% aggregate visresponses: ICsig, RFCI, RFCIspin, sizeCI, oriparams
ises=1;
pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load(sprintf('%svisresponses_probeC.mat', pathpp))
load(sprintf('%spostprocessed_probeC.mat', pathpp), 'neuoind')
Nneurons = numel(neuoind);

probeneuronsagg = cell(size(probes));
neulocagg = cell(size(probes));
neupeakchagg = cell(size(probes));
neuctxagg = cell(size(probes));
sesneuagg = cell(size(probes));
sesneuctxagg = cell(size(probes));

ICblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations'};
% 'Palpha','BKtt','BKttpair','BItt','BIttpair','BICREl1tt','BICREl2tt','BICREl1ttpair','BICREl2ttpair',
allfields = fieldnames(ICsig.ICwcfg1_presentations);
validfields = true(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(ICsig.ICwcfg1_presentations.(allfields{f}),1)==Nneurons;
end
ICsigfields = allfields(validfields);
% {'SP_Ind','Pmww_Ind','SP_BK','sigmcBK','Pmww_BK', ...
%     'SP_ICvsRC','Pmww_ICvsRC','SP_BICREl','sigmcBICREl1','sigmcBICREl2', ...
%     'PkwBK','PmcBK', 'PkwBI','PmcBI', 'PkwBICREl1','PmcBICREl1', 'PkwBICREl2','PmcBICREl2', ...
%     'ICencoder','RCencoder','inducerencoder','inducerresponsive', ...
%     'indenc1','indenc2','indenc3','indenc4','indenc13','indenc14','indenc23','indenc24', ...
%     'ICencoder1','RCencoder1','RCencoder2','ICencoder2','indin1','indin2','indin3','indin4', ...
%     'indout1','indout2','indout3','indout4','RElfaith1','RElfaith2', ...
%     'ICresp1','RCresp1','RCresp2','ICresp2','ICtuned1','RCtuned1','RCtuned2','ICtuned2'};
ICsigagg = struct();
for iprobe = 1:numel(probes)
    for b = 1:numel(ICblocks)
        for f= 1:numel(ICsigfields)
            ICsigagg.(ICblocks{b}).(probes{iprobe}).(ICsigfields{f}) = [];
        end
    end
end

allfields = fieldnames(RFCI);
validfields = true(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(RFCI.(allfields{f}),1)==Nneurons;
end
RFCIfields = allfields(validfields);
% {'Rrfclassic','Rrfinverse','RFindclassic','RFindinverse', ...
%     'Pkw_rfclassic','Pkw_rfinverse','pRrfclassic','pRrfinverse','pRFclassic','pRFinverse'};
RFCIagg = struct();
for iprobe = 1:numel(probes)
    for f= 1:numel(RFCIfields)
        RFCIagg.(probes{iprobe}).(RFCIfields{f}) = [];
    end
end

allfields = fieldnames(RFCIspin);
validfields = true(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(RFCIspin.(allfields{f}),1)==Nneurons;
end
RFCIspinfields = allfields(validfields);
RFCIspinagg = struct();
for iprobe = 1:numel(probes)
    for f= 1:numel(RFCIspinfields)
        RFCIspinagg.(probes{iprobe}).(RFCIspinfields{f}) = [];
    end
end

allfields = fieldnames(sizeCI);
validfields = true(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(sizeCI.(allfields{f}),1)==Nneurons;
end
sizeCIfields = allfields(validfields);
% sizeCIfields = {'Rsizeclassic','Rsizeinverse','sizeindclassic','sizeindinverse','Pkw_sizeclassic','Pkw_sizeinverse'};
sizeCIagg = struct();
for iprobe = 1:numel(probes)
    for f= 1:numel(sizeCIfields)
        sizeCIagg.(probes{iprobe}).(sizeCIfields{f}) = [];
    end
end

allfields = fieldnames(oriparams);
validfields = true(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(oriparams.(allfields{f}),1)==Nneurons;
end
oriparamsfields = allfields(validfields);
% oriparamsfields = {'Rori','prefiori','orthiori','OSI','OP','Pmww_OP'};
oriparamsagg = struct();
for iprobe = 1:numel(probes)
    for f= 1:numel(oriparamsfields)
        oriparamsagg.(probes{iprobe}).(oriparamsfields{f}) = [];
    end
end

allfields = fieldnames(ori4params);
validfields = true(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(ori4params.(allfields{f}),1)==Nneurons;
end
ori4paramsfields = allfields(validfields);
% oriparamsfields = {'Rori','prefiori','orthiori','OSI','OP','Pmww_OP'};
ori4paramsagg = struct();
for iprobe = 1:numel(probes)
    for f= 1:numel(ori4paramsfields)
        ori4paramsagg.(probes{iprobe}).(ori4paramsfields{f}) = [];
    end
end

meanFRvecagg = cell(size(probes));
sponFRvecagg = cell(size(probes));

for ises = 1:Nsessions
    fprintf('Session %d/%d %s\n', ises, Nsessions, nwbsessions{ises} )
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur'
    
    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    
    for iprobe = 1:numel(probes)
        
        if exist([pathpp, 'probes.mat'], 'file')
            probelist = load([pathpp, 'probes.mat']);
            warning('HS 230126: this was inserted to handle the exception case of sub_1183369803, can delete with the next nwb update')
        else
            probelist.probes = {'A', 'B', 'C', 'D', 'E', 'F'};
        end
        probeind = find( strcmp(probes{iprobe}, probelist.probes) );
        %probeind = find( strcmp(probes{iprobe}, {'A', 'B', 'C', 'D', 'E', 'F'}) );
        
        if nnz(floor(unit_peakch/1000)==probeind-1)==0
            fprintf('Probe %s Area %s: NO UNITS!!!\n', probes{iprobe}, visareas{iprobe} )
            continue
        end
%         tic
%         load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}))
%         % 'neuoind', 'vis', 'Tres', 'psthtli', 'psth'
        load(sprintf('%svisresponses_probe%s.mat', pathpp, probes{iprobe}))
        % 'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'ICsig', 'RFCI', 'sizeCI', 'oriparams'
        
%         if ~isequal(unique(floor(unit_peakch(neuoind)/1000)), probeind-1)
%             error('check neuoind')
%         end
        neuoind = find(floor(unit_peakch/1000)==probeind-1);
        
        % check whether CCF registration is correct
        probelocs = electrode_location(ismember(electrode_id, unit_peakch(neuoind)));
        
        neuloc = electrode_location(revmapelecid(unit_peakch(neuoind)+1));
        if ~isequal(unique(probelocs), unique(neuloc))
            disp(unique(neuloc)')
            error('check neuloc')
        end
        
        neuctx = contains(neuloc, 'VIS');
        fprintf('Probe %s Area %s: %d/%d\n', probes{iprobe}, visareas{iprobe}, nnz(neuctx), numel(neuoind) )
%         disp(unique(probelocs)')
        
        probeneuronsagg{iprobe} = cat(1, probeneuronsagg{iprobe}, neuoind);
        neulocagg{iprobe} = cat(1, neulocagg{iprobe}, neuloc);
        neupeakchagg{iprobe} = cat(1, neupeakchagg{iprobe}, unit_peakch(neuoind));
        neuctxagg{iprobe} = cat(1, neuctxagg{iprobe}, neuctx);
        
        sesneuagg{iprobe} = cat(1, sesneuagg{iprobe}, ises*ones(length(neuoind),1));
        sesneuctxagg{iprobe} = cat(1, sesneuctxagg{iprobe}, ises*ones(nnz(neuctx),1));
        
%         % aggregate cortex neurons
%         neuoi = neuctx;
        % aggregate all neurons
        neuoi = true(size(neuoind));
        
        meanFRvecagg{iprobe} = cat(1, meanFRvecagg{iprobe}, meanFRvec(neuoi)');
        sponFRvecagg{iprobe} = cat(1, sponFRvecagg{iprobe}, sponFRvec(neuoi)');
        
        for b = 1:numel(ICblocks)
            for f= 1:numel(ICsigfields)
                ICsigagg.(ICblocks{b}).(probes{iprobe}).(ICsigfields{f}) = cat(1, ICsigagg.(ICblocks{b}).(probes{iprobe}).(ICsigfields{f}), ICsig.(ICblocks{b}).(ICsigfields{f})(neuoi,:) );
            end
        end
        
        % RFCIfields = fieldnames(RFCI);
        for f= 1:numel(RFCIfields)
            RFCIagg.(probes{iprobe}).(RFCIfields{f}) = cat(1, RFCIagg.(probes{iprobe}).(RFCIfields{f}), RFCI.(RFCIfields{f})(neuoi,:) );
        end
        
        for f= 1:numel(RFCIspinfields)
            RFCIspinagg.(probes{iprobe}).(RFCIspinfields{f}) = cat(1, RFCIspinagg.(probes{iprobe}).(RFCIspinfields{f}), RFCIspin.(RFCIspinfields{f})(neuoi,:) );
        end
        
        % size vector [0, 4, 8, 16, 32, 64 ]
        for f= 1:numel(sizeCIfields)
            sizeCIagg.(probes{iprobe}).(sizeCIfields{f}) = cat(1, sizeCIagg.(probes{iprobe}).(sizeCIfields{f}), sizeCI.(sizeCIfields{f})(neuoi,:) );
        end
        
        for f= 1:numel(oriparamsfields)
            oriparamsagg.(probes{iprobe}).(oriparamsfields{f}) = cat(1, oriparamsagg.(probes{iprobe}).(oriparamsfields{f}), oriparams.(oriparamsfields{f})(neuoi,:) );
        end
        
        for f= 1:numel(ori4paramsfields)
            ori4paramsagg.(probes{iprobe}).(ori4paramsfields{f}) = cat(1, ori4paramsagg.(probes{iprobe}).(ori4paramsfields{f}), ori4params.(ori4paramsfields{f})(neuoi,:) );
        end
        
%         toc
    end
    
end

Nrfs = size(RFCI.Rrfclassic, 2);
Nszs = size(sizeCI.Rsizeclassic, 2);
Ndirs = size(oriparams.Rori, 2);
Noris = Ndirs/2;

ises=4;
pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load(sprintf('%spostprocessed_probeC.mat', pathpp), 'vis')
dirvec = vis.sizeCI_presentations.directions;
if length(dirvec)~=Ndirs
    error('check sizeCI_presentations directions')
end
orivec = vis.sizeCI_presentations.directions(1:Noris);

%%
save(['G:\My Drive\DATA\ICexpts_submission22\openscope_popavg_agg.mat'], ...
    'nwbsessions', 'probeneuronsagg', 'neulocagg', 'neupeakchagg', 'sesneuagg', 'neuctxagg', 'sesneuctxagg', ...
    'meanFRvecagg', 'sponFRvecagg', 'vis', ...
    'ICsigfields', 'ICsigagg', 'RFCIfields', 'RFCIagg', ...
    'RFCIspinfields', 'RFCIspinagg', 'sizeCIfields', 'sizeCIagg', ...
    'oriparamsfields', 'oriparamsagg', 'ori4paramsfields', 'ori4paramsagg', '-v7.3')

%% aggregate Ron Roff and psthagg
aggpsth = true;
if ~aggpsth
probesR = {'C'};
else
    probesR = probes;
end

visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
    'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
ICblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations'};
vistrialtypesagg = struct();
vistrialrepagg = struct();
vistrialorderagg = struct();
if aggpsth
psthavgagg = struct();
end
Ronavgagg = struct();
Roffavgagg = struct();
for b = 1:numel(visblocks)
    for iprobe = 1:numel(probesR)
        if aggpsth
            psthavgagg.(visblocks{b}).(probesR{iprobe})=[];
        end
        Ronavgagg.(visblocks{b}).(probesR{iprobe})=[];
        Roffavgagg.(visblocks{b}).(probesR{iprobe})=[];
    end
end

for ises = 1:Nsessions
    fprintf('Session %d/%d %s\n', ises, Nsessions, nwbsessions{ises} )
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur'
    
    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    
    for iprobe = 1:numel(probesR)
        tic
        
        if exist([pathpp, 'probes.mat'], 'file')
            probelist = load([pathpp, 'probes.mat']);
            warning('HS 230126: this was inserted to handle the exception case of sub_1183369803, can delete with the next nwb update')
        else
            probelist.probes = {'A', 'B', 'C', 'D', 'E', 'F'};
        end
        probeind = find( strcmp(probesR{iprobe}, probelist.probes) );
        %probeind = find( strcmp(probes{iprobe}, {'A', 'B', 'C', 'D', 'E', 'F'}) );

        if nnz(floor(unit_peakch/1000)==probeind-1)==0
        fprintf('Probe %s Area %s: NO UNITS!!!\n', probesR{iprobe}, visareas{iprobe} )
            continue
        end
        
        load(sprintf('%spostprocessed_probe%s.mat', pathpp, probesR{iprobe}))
        % 'neuoind', 'vis', 'Tres', 'psthtli', 'psth'
        %     load(sprintf('%svisresponses_probe%s.mat', pathpp, probesR{iprobe}))
        %     % 'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'ICsig', 'RFCI', 'sizeCI', 'oriparams'
        
        if ~isequal(unique(floor(unit_peakch(neuoind)/1000)), probeind-1)
            error('check neuoind')
        end
        
        % check whether CCF registration is correct
        probelocs = electrode_location(ismember(electrode_id, unit_peakch(neuoind)));
        
        neuloc = electrode_location(revmapelecid(unit_peakch(neuoind)+1));
        if ~isequal(unique(probelocs), unique(neuloc))
            disp(unique(neuloc)')
            error('check neuloc')
        end
        
        neuctx = contains(neuloc, 'VIS');
        
        fprintf('Probe %s Area %s: %d/%d\n', probesR{iprobe}, visareas{iprobe}, nnz(neuctx), numel(neuoind) )
        disp(unique(probelocs)')
        
%         probeneuronsagg{iprobe} = cat(1, probeneuronsagg{iprobe}, neuoind);
%         neulocagg{iprobe} = cat(1, neulocagg{iprobe}, neuloc);
%         neupeakchagg{iprobe} = cat(1, neupeakchagg{iprobe}, unit_peakch(neuoind));
%         neuctxagg{iprobe} = cat(1, neuctxagg{iprobe}, neuctx);
%         
%         sesneuagg{iprobe} = cat(1, sesneuagg{iprobe}, ises*ones(length(neuoind),1));
%         sesneuctxagg{iprobe} = cat(1, sesneuctxagg{iprobe}, ises*ones(nnz(neuctx),1));
        
%         % aggregate cortex neurons
%         neuoi = neuctx;
        % aggregate all neurons
        neuoi = true(size(neuoind));
        
        
        % visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
        %     'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
        %ICblocks: each stim is 0.4s, inter-trial interval is 0.4s, static images
        %RFCI: each stim is 0.25s, inter-trial interval is 0s, spinning drifting grating
        %sizeCI: each stim is 0.25s, inter-trial interval is 0.5s, drifting grating
        % orientation denotation is same as psychtoolbox (0 is 12h, 45 is 1h30m, clockwise)
        for b = 1:numel(visblocks)
            vistrialtypesagg(ises).(visblocks{b}) = unique(vis.(visblocks{b}).trialorder);
            vistrialorderagg(ises).(visblocks{b}) = vis.(visblocks{b}).trialorder;
            if strcmp(visblocks{b}, 'RFCI_presentations')
                vistrialtypesagg(ises).(visblocks{b}) = unique(vis.(visblocks{b}).trialorder(1:4:end));
                vistrialorderagg(ises).(visblocks{b}) = vis.(visblocks{b}).trialorder(1:4:end);
            end
            
            if ismember(visblocks{b}, ICblocks)
                tlon = psthtli>0 & psthtli<=400;
                tloff = psthtli>400 & psthtli<=800;
            elseif strcmp(visblocks{b}, 'RFCI_presentations')
                tlon = psthtli>0 & psthtli<=1000;
                tloff = psthtli>1000 & psthtli<=1000;
            elseif strcmp(visblocks{b}, 'sizeCI_presentations')
                tlon = psthtli>0 & psthtli<=250;
                tloff = psthtli>250 & psthtli<=750;
            else
                error('vis block not recognized')
            end
            
            Ntt = numel(vistrialtypesagg(ises).(visblocks{b}));
            vistrialrepagg(ises).(visblocks{b}) = zeros(Ntt,1);
            temppsth = zeros(length(psthtli), Ntt, nnz(neuoi));
            tempRonavg = zeros(Ntt, nnz(neuoi));
            tempRoffavg = zeros(Ntt, nnz(neuoi));
            for ii = 1:Ntt
                trialsoi = vis.(visblocks{b}).trialorder==vistrialtypesagg(ises).(visblocks{b})(ii);
                vistrialrepagg(ises).(visblocks{b})(ii) = nnz(trialsoi);
                temppsth(:,ii,:) = mean(1000*psth.(visblocks{b})(:,trialsoi,neuoi), 2);
            tempRonavg(ii,:) = squeeze(mean(1000*psth.(visblocks{b})(tlon,trialsoi,neuoi), [1 2]));
            tempRoffavg(ii,:) = squeeze(mean(1000*psth.(visblocks{b})(tloff,trialsoi,neuoi), [1 2]));
            end
            if aggpsth
                psthavgagg.(visblocks{b}).(probesR{iprobe}) = cat(3, psthavgagg.(visblocks{b}).(probesR{iprobe}), temppsth );
            end
            Ronavgagg.(visblocks{b}).(probesR{iprobe}) = cat(2, Ronavgagg.(visblocks{b}).(probesR{iprobe}), tempRonavg );
            Roffavgagg.(visblocks{b}).(probesR{iprobe}) = cat(2, Roffavgagg.(visblocks{b}).(probesR{iprobe}), tempRoffavg );
        end
        
        % size vector [0, 4, 8, 16, 32, 64 ]
        toc
    end
    
end

if aggpsth
    save('D:\OpenScopeData\000248\postprocessed\psthavgagg.mat', 'probes', 'visareas', 'visind', 'nwbsessions', ...
    	'probeneuronsagg', 'neulocagg', 'neupeakchagg', 'sesneuagg', 'neuctxagg', 'sesneuctxagg', ...
        'vistrialtypesagg', 'vistrialrepagg', 'vistrialorderagg', ...
        'ICblocks', 'ICtrialtypes', 'psthtli', 'psthavgagg', 'Ronavgagg', 'Roffavgagg', '-v7.3')
    save('G:\My Drive\DATA\ICexpts_submission22\openscope_psthavgagg.mat', 'probes', 'visareas', 'visind', 'nwbsessions', ...
    	'probeneuronsagg', 'neulocagg', 'neupeakchagg', 'sesneuagg', 'neuctxagg', 'sesneuctxagg', ...
        'vistrialtypesagg', 'vistrialrepagg', 'vistrialorderagg', ...
        'ICblocks', 'ICtrialtypes', 'psthtli', 'psthavgagg', 'Ronavgagg', 'Roffavgagg', '-v7.3')
end

%% report number of units in each area/session/probe
aggneuloc = cat(1,neulocagg{:});
aggsesneu = cat(1,sesneuagg{:});

[v,c]=uniquecnt(aggneuloc);
disp([v c])

areaunitsperses = zeros(length(v), Nsessions);
for ii = 1:length(v)
    neuoi = strcmp(aggneuloc, v(ii));
    hc=histcounts(aggsesneu(neuoi), 0.5:1:Nsessions+0.5);
    if ~( nnz(neuoi)==c(ii) && sum(hc)==c(ii) )
        error('mismatch between uniquecnt and strcmp -- check')
    end
    areaunitsperses(ii,:) = hc;
end

areaunitstab = table(v,c,areaunitsperses);
open areaunitstab

