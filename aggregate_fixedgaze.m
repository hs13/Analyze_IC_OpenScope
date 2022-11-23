addpath(genpath('H:\CODE\helperfunctions'))
addpath(genpath('H:\CODE\Analyze_OpenScope'))

datadir = 'D:\OpenScopeData\000248\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(contains(nwbsessions, 'sub-'));

gazedistthresh = 20;
Nsessions = numel(nwbsessions)-1;
whichneuctx = 0; % 0: all enurons, 1: correct definition, 2: >=230, 3: topmost-250um

% A-AM, B-PM, C-V1, D-LM, E-AL, F-RL
probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
visind = [6 5 1 2 4 3];

% kerwinhalf = 2; kersigma = 1;
% kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
% kergauss = (kergauss/sum(kergauss));

%% aggregate visresponses: ICsig, RFCI, RFCIspin, sizeCI, oriparams
ises=4;
pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load(sprintf('%svisresponses_fixedgaze%dpix_probeC.mat', pathpp, gazedistthresh))
load(sprintf('%spostprocessed_probeC.mat', pathpp), 'neuoind')
Nneurons = numel(neuoind);

probeneuronsagg = cell(size(probes));
neuctxagg = cell(size(probes));
neulocagg = cell(size(probes));
sesneuoiagg = cell(size(probes));

ICblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations'};
% 'Palpha','BKtt','BKttpair','BItt','BIttpair','BICREl1tt','BICREl2tt','BICREl1ttpair','BICREl2ttpair',
allfields = fieldnames(ICsig_fixedgaze.ICwcfg1_presentations);
validfields = true(size(allfields));
size2fields = zeros(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(ICsig_fixedgaze.ICwcfg1_presentations.(allfields{f}),1)==Nneurons;
    size2fields(f) = size(ICsig_fixedgaze.ICwcfg1_presentations.(allfields{f}),2);
end
ICsigfields = allfields(validfields);
ICsigfieldsize2 = size2fields(validfields);
% {'SP_Ind','Pmww_Ind','SP_BK','sigmcBK','Pmww_BK', ...
%     'SP_ICvsRC','Pmww_ICvsRC','SP_BICREl','sigmcBICREl1','sigmcBICREl2', ...
%     'PkwBK','PmcBK', 'PkwBI','PmcBI', 'PkwBICREl1','PmcBICREl1', 'PkwBICREl2','PmcBICREl2', ...
%     'ICencoder','RCencoder','inducerencoder','inducerresponsive', ...
%     'indenc1','indenc2','indenc3','indenc4','indenc13','indenc14','indenc23','indenc24', ...
%     'ICencoder1','RCencoder1','RCencoder2','ICencoder2','indin1','indin2','indin3','indin4', ...
%     'indout1','indout2','indout3','indout4','RElfaith1','RElfaith2', ...
%     'ICresp1','RCresp1','RCresp2','ICresp2','ICtuned1','RCtuned1','RCtuned2','ICtuned2'};
ICsig_fixedgazeagg = struct();
for iprobe = 1:numel(probes)
    for b = 1:numel(ICblocks)
        for f= 1:numel(ICsigfields)
            ICsig_fixedgazeagg.(ICblocks{b})(iprobe).(ICsigfields{f}) = [];
        end
    end
end

allfields = fieldnames(RFCI_fixedgaze);
validfields = true(size(allfields));
size2fields = zeros(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(RFCI_fixedgaze.(allfields{f}),1)==Nneurons;
    size2fields(f) = size(RFCI_fixedgaze.(allfields{f}),2);
end
RFCIfields = allfields(validfields);
RFCIfieldsize2 = size2fields(validfields);
% {'Rrfclassic','Rrfinverse','RFindclassic','RFindinverse', ...
%     'Pkw_rfclassic','Pkw_rfinverse','pRrfclassic','pRrfinverse','pRFclassic','pRFinverse'};
RFCI_fixedgazeagg = struct();
for iprobe = 1:numel(probes)
    for f= 1:numel(RFCIfields)
        RFCI_fixedgazeagg(iprobe).(RFCIfields{f}) = [];
    end
end

allfields = fieldnames(RFCIspin_fixedgaze);
validfields = true(size(allfields));
size2fields = zeros(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(RFCIspin_fixedgaze.(allfields{f}),1)==Nneurons;
    size2fields(f) = size(RFCIspin_fixedgaze.(allfields{f}),2);
end
RFCIspinfields = allfields(validfields);
RFCIspinfieldsize2 = size2fields(validfields);
RFCIspin_fixedgazeagg = struct();
for iprobe = 1:numel(probes)
    for f= 1:numel(RFCIspinfields)
        RFCIspin_fixedgazeagg(iprobe).(RFCIspinfields{f}) = [];
    end
end

allfields = fieldnames(sizeCI_fixedgaze);
validfields = true(size(allfields));
size2fields = zeros(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(sizeCI_fixedgaze.(allfields{f}),1)==Nneurons;
    size2fields(f) = size(sizeCI_fixedgaze.(allfields{f}),2);
end
sizeCIfields = allfields(validfields);
sizeCIfieldsize2 = size2fields(validfields);
% sizeCIfields = {'Rsizeclassic','Rsizeinverse','sizeindclassic','sizeindinverse','Pkw_sizeclassic','Pkw_sizeinverse'};
sizeCI_fixedgazeagg = struct();
for iprobe = 1:numel(probes)
    for f= 1:numel(sizeCIfields)
        sizeCI_fixedgazeagg(iprobe).(sizeCIfields{f}) = [];
    end
end

allfields = fieldnames(oriparams_fixedgaze);
validfields = true(size(allfields));
size2fields = zeros(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(oriparams_fixedgaze.(allfields{f}),1)==Nneurons;
    size2fields(f) = size(oriparams_fixedgaze.(allfields{f}),2);
end
oriparamsfields = allfields(validfields);
oriparamsfieldsize2 = size2fields(validfields);
% oriparamsfields = {'Rori','prefiori','orthiori','OSI','OP','Pmww_OP'};
oriparams_fixedgazeagg = struct();
for iprobe = 1:numel(probes)
    for f= 1:numel(oriparamsfields)
        oriparams_fixedgazeagg(iprobe).(oriparamsfields{f}) = [];
    end
end

allfields = fieldnames(ori4params_fixedgaze);
validfields = true(size(allfields));
size2fields = zeros(size(allfields));
for f = 1:numel(allfields)
    validfields(f) = size(ori4params_fixedgaze.(allfields{f}),1)==Nneurons;
    size2fields(f) = size(ori4params_fixedgaze.(allfields{f}),2);
end
ori4paramsfields = allfields(validfields);
ori4paramsfieldsize2 = size2fields(validfields);
% oriparamsfields = {'Rori','prefiori','orthiori','OSI','OP','Pmww_OP'};
ori4params_fixedgazeagg = struct();
for iprobe = 1:numel(probes)
    for f= 1:numel(ori4paramsfields)
        ori4params_fixedgazeagg(iprobe).(ori4paramsfields{f}) = [];
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
        tic
%         load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}))
%         % 'neuoind', 'vis', 'Tres', 'psthtli', 'psth'
        load(sprintf('%svisresponses_fixedgaze%dpix_probe%s.mat', pathpp, gazedistthresh, probes{iprobe}))
        % 'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'ICsig', 'RFCI', 'sizeCI', 'oriparams'
        
        probeind = find( strcmp(probes{iprobe}, {'A', 'B', 'C', 'D', 'E', 'F'}) );
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
        
        switch whichneuctx
            case 0
                neuctx = true(size(neuoind));
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
                if ~isfield(ICsig_fixedgaze.(ICblocks{b}), ICsigfields{f})
                    tempmat = NaN( nnz(neuoi), ICsigfieldsize2(f) );
                else
                    tempmat = ICsig_fixedgaze.(ICblocks{b}).(ICsigfields{f})(neuoi,:);
                end
                ICsig_fixedgazeagg.(ICblocks{b})(iprobe).(ICsigfields{f}) = ...
                    cat(1, ICsig_fixedgazeagg.(ICblocks{b})(iprobe).(ICsigfields{f}), tempmat );
            end
        end
        
        % RFCIfields = fieldnames(RFCI);
        for f= 1:numel(RFCIfields)
            if ~isfield(RFCI_fixedgaze, RFCIfields{f})
                tempmat = NaN( nnz(neuoi), RFCIfieldsize2(f) );
            else
                tempmat = RFCI_fixedgaze.(RFCIfields{f})(neuoi,:);
            end
            RFCI_fixedgazeagg(iprobe).(RFCIfields{f}) = cat(1, ...
                RFCI_fixedgazeagg(iprobe).(RFCIfields{f}), tempmat );
        end
        
        for f= 1:numel(RFCIspinfields)
            if ~isfield(RFCIspin_fixedgaze, RFCIspinfields{f})
                tempmat = NaN( nnz(neuoi), RFCIspinfieldsize2(f) );
            else
                tempmat = RFCIspin_fixedgaze.(RFCIspinfields{f})(neuoi,:);
            end
            RFCIspin_fixedgazeagg(iprobe).(RFCIspinfields{f}) = cat(1, ...
                RFCIspin_fixedgazeagg(iprobe).(RFCIspinfields{f}), tempmat );
        end
        
        % size vector [0, 4, 8, 16, 32, 64 ]
        for f= 1:numel(sizeCIfields)
            if isfield(sizeCI_fixedgaze, sizeCIfields{f})
                tempmat = NaN( nnz(neuoi), sizeCIfieldsize2(f) );
            else
                tempmat = sizeCI_fixedgaze.(sizeCIfields{f})(neuoi,:);
            end
            sizeCI_fixedgazeagg(iprobe).(sizeCIfields{f}) = cat(1, ...
                sizeCI_fixedgazeagg(iprobe).(sizeCIfields{f}), tempmat );
        end
        
        for f= 1:numel(oriparamsfields)
            if ~isfield(oriparams_fixedgaze, oriparamsfields{f})
                tempmat = NaN( nnz(neuoi), oriparamsfieldsize2(f) );
            else
                tempmat = oriparams_fixedgaze.(oriparamsfields{f})(neuoi,:);
            end
            oriparams_fixedgazeagg(iprobe).(oriparamsfields{f}) = cat(1, ...
                oriparams_fixedgazeagg(iprobe).(oriparamsfields{f}), tempmat );
        end
        
        for f= 1:numel(ori4paramsfields)
            if ~isfield(ori4params_fixedgaze, ori4paramsfields{f})
                tempmat = NaN( nnz(neuoi), ori4paramsfieldsize2(f) );
            else
                tempmat = ori4params_fixedgaze.(ori4paramsfields{f})(neuoi,:);
            end
            ori4params_fixedgazeagg(iprobe).(ori4paramsfields{f}) = cat(1, ...
                ori4params_fixedgazeagg(iprobe).(ori4paramsfields{f}), tempmat );
        end
        
        toc
    end
    
end

Nrfs = RFCIfieldsize2(strcmp(RFCIfields, 'Rrfclassic'));
Nszs = sizeCIfieldsize2(strcmp(sizeCIfields, 'Rsizeclassic'));
Ndirs = oriparamsfieldsize2(strcmp(oriparamsfields, 'Rori'));
Noris = Ndirs/2;

ises=4;
pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load(sprintf('%spostprocessed_probeC.mat', pathpp), 'vis')
dirvec = vis.sizeCI_presentations.directions;
if length(dirvec)~=Ndirs
    error('check sizeCI_presentations directions')
end
orivec = vis.sizeCI_presentations.directions(1:Noris);

%% V1: pref-ori distribution of IC/RC encoders
iprobe = find( strcmp('C', probes) );
opt678 = true;
if opt678
prefidiragg = oriparams_fixedgazeagg(iprobe).prefiori678;
prefioriagg = ori4params_fixedgazeagg(iprobe).prefiori4678;
sigori = ori4params_fixedgazeagg(iprobe).Pkw_ori4678<0.05;
% Rori = ( oriparams_fixedgazeagg(iprobe).Rori(:,1:4) + oriparams_fixedgazeagg(iprobe).Rori(:,5:8) )/2;
% [~,prefioriagg] = max(Rori, [], 2);
Rori = ori4params_fixedgazeagg(iprobe).Rori4678;
else
prefidiragg = oriparams_fixedgazeagg(iprobe).prefiori;
prefioriagg = ori4params_fixedgazeagg(iprobe).prefiori4;
sigori = ori4params_fixedgazeagg(iprobe).Pkw_ori4<0.05;
Rori = ori4params_fixedgazeagg(iprobe).Rori4;
end

% figure; hold all
% for ienc = 1:4
%     switch ienc
%         case 1
%             neuoi = ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).ICencoder2==1;
%             neuoi = ismember(ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
%         case 2
%             neuoi = ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).ICencoder2==1;
%             neuoi = ismember(ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
%         case 3
%             neuoi = ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).ICencoder1==1;
%             neuoi = ismember(ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
%         case 4
%             neuoi = ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).ICencoder1==1;
%             neuoi = ismember(ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
%     end
% errorbar(1:Noris, nanmean(Rori(neuoi,:),1), nanstd(Rori(neuoi,:),0,1)/sqrt(nnz(neuoi)), 'o-')
% end

hcprefidir_wICenc = zeros(4, Ndirs);
hcprefiori_wICenc = zeros(4, Noris);
hcsigprefidir_wICenc = zeros(4, Ndirs);
hcsigprefiori_wICenc = zeros(4, Noris);
for ienc = 1:4
    switch ienc
        case 1
            neuoi = ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).ICencoder2==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
        case 2
            neuoi = ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).ICencoder2==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
        case 3
            neuoi = ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).ICencoder1==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
        case 4
            neuoi = ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).ICencoder1==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    end
    hcprefidir_wICenc(ienc,:) = histcounts(prefidiragg(neuoi), 0.5:Ndirs+0.5);
    hcprefiori_wICenc(ienc,:) = histcounts(prefioriagg(neuoi), 0.5:Noris+0.5);
    hcsigprefidir_wICenc(ienc,:) = histcounts(prefidiragg(sigori & neuoi), 0.5:Ndirs+0.5);
    hcsigprefiori_wICenc(ienc,:) = histcounts(prefioriagg(sigori & neuoi), 0.5:Noris+0.5);
end

hcprefidir_kICenc = zeros(4, Ndirs);
hcprefiori_kICenc = zeros(4, Noris);
hcsigprefidir_kICenc = zeros(4, Ndirs);
hcsigprefiori_kICenc = zeros(4, Noris);
for ienc = 1:4
    switch ienc
        case 1
            neuoi = ICsig_fixedgazeagg.ICkcfg0_presentations(iprobe).ICencoder2==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICkcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
        case 2
            neuoi = ICsig_fixedgazeagg.ICkcfg1_presentations(iprobe).ICencoder2==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICkcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
        case 3
            neuoi = ICsig_fixedgazeagg.ICkcfg0_presentations(iprobe).ICencoder1==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICkcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
        case 4
            neuoi = ICsig_fixedgazeagg.ICkcfg1_presentations(iprobe).ICencoder1==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICkcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    end
    hcprefidir_kICenc(ienc,:) = histcounts(prefidiragg(neuoi), 0.5:Ndirs+0.5);
    hcprefiori_kICenc(ienc,:) = histcounts(prefioriagg(neuoi), 0.5:Noris+0.5);
    hcsigprefidir_kICenc(ienc,:) = histcounts(prefidiragg(sigori & neuoi), 0.5:Ndirs+0.5);
    hcsigprefiori_kICenc(ienc,:) = histcounts(prefioriagg(sigori & neuoi), 0.5:Noris+0.5);
end

hcprefidir = histcounts(prefidiragg, 0.5:Ndirs+0.5);
hcprefiori = histcounts(prefioriagg, 0.5:Noris+0.5);
hcsigprefidir = histcounts(prefidiragg(sigori), 0.5:Ndirs+0.5);
hcsigprefiori = histcounts(prefioriagg(sigori), 0.5:Noris+0.5);

fs=10;
ytl = 0:45:180-1; ylab = 'IC Orientation';
figure; 
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'Fixed Gaze', 'edgecolor', 'none')
subplot(2,2,1); imagesc(hcprefidir_wICenc./sum(hcprefidir_wICenc,2))
colorbar; colormap jet
set(gca, 'FontSize', fs, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Direction')
title('White Inducers Black Background')
subplot(2,2,2); imagesc(hcprefiori_wICenc./sum(hcprefiori_wICenc,2))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
subplot(2,2,3); imagesc(hcprefidir_kICenc./sum(hcprefidir_kICenc,2))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Direction')
title('Black Inducers White Background')
subplot(2,2,4); imagesc(hcprefiori_kICenc./sum(hcprefiori_kICenc,2))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')

figure; 
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'Fixed Gaze: Significantly Tuned Neurons', 'edgecolor', 'none')
subplot(2,2,1); imagesc(hcsigprefidir_wICenc./sum(hcsigprefidir_wICenc,2))
colorbar; colormap jet
xlabel('Preferred Direction')
title('White Inducers Black Background')
subplot(2,2,2); imagesc(hcsigprefiori_wICenc./sum(hcsigprefiori_wICenc,2))
colorbar; colormap jet
xlabel('Preferred Orientation')
subplot(2,2,3); imagesc(hcsigprefidir_kICenc./sum(hcsigprefidir_kICenc,2))
colorbar; colormap jet
xlabel('Preferred Direction')
title('Black Inducers White Background')
subplot(2,2,4); imagesc(hcsigprefiori_kICenc./sum(hcsigprefiori_kICenc,2))
colorbar; colormap jet
xlabel('Preferred Orientation')

%{
figure; 
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'Fixed Gaze', 'edgecolor', 'none')
subplot(2,2,1); imagesc(hcprefidir_wICenc./hcprefidir)
colorbar; colormap jet
xlabel('Preferred Direction')
title('White Inducers Black Background')
subplot(2,2,2); imagesc(hcprefiori_wICenc./hcprefiori)
colorbar; colormap jet
xlabel('Preferred Orientation')
subplot(2,2,3); imagesc(hcprefidir_kICenc./hcprefidir)
colorbar; colormap jet
xlabel('Preferred Direction')
title('Black Inducers White Background')
subplot(2,2,4); imagesc(hcprefiori_kICenc./hcprefiori)
colorbar; colormap jet
xlabel('Preferred Orientation')

figure; 
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'Fixed Gaze: Significantly Tuned Neurons', 'edgecolor', 'none')
subplot(2,2,1); imagesc(hcsigprefidir_wICenc./hcsigprefidir)
colorbar; colormap jet
xlabel('Preferred Direction')
title('White Inducers Black Background')
subplot(2,2,2); imagesc(hcsigprefiori_wICenc./hcsigprefiori)
colorbar; colormap jet
xlabel('Preferred Orientation')
subplot(2,2,3); imagesc(hcsigprefidir_kICenc./hcsigprefidir)
colorbar; colormap jet
xlabel('Preferred Direction')
title('Black Inducers White Background')
subplot(2,2,4); imagesc(hcsigprefiori_kICenc./hcsigprefiori)
colorbar; colormap jet
xlabel('Preferred Orientation')


figure; 
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'Fixed Gaze', 'edgecolor', 'none')
subplot(2,2,1); imagesc(hcprefidir_wICenc)
colorbar; colormap jet
xlabel('Preferred Direction')
title('White Inducers Black Background')
subplot(2,2,2); imagesc(hcprefiori_wICenc)
colorbar; colormap jet
xlabel('Preferred Orientation')
subplot(2,2,3); imagesc(hcprefidir_kICenc)
colorbar; colormap jet
xlabel('Preferred Direction')
title('Black Inducers White Background')
subplot(2,2,4); imagesc(hcprefiori_kICenc)
colorbar; colormap jet
xlabel('Preferred Orientation')

figure; 
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'Fixed Gaze: Significantly Tuned Neurons', 'edgecolor', 'none')
subplot(2,2,1); imagesc(hcsigprefidir_wICenc)
colorbar; colormap jet
xlabel('Preferred Direction')
title('White Inducers Black Background')
subplot(2,2,2); imagesc(hcsigprefiori_wICenc)
colorbar; colormap jet
xlabel('Preferred Orientation')
subplot(2,2,3); imagesc(hcsigprefidir_kICenc)
colorbar; colormap jet
xlabel('Preferred Direction')
title('Black Inducers White Background')
subplot(2,2,4); imagesc(hcsigprefiori_kICenc)
colorbar; colormap jet
xlabel('Preferred Orientation')
%}

%% all areas: pref-ori distribution of IC/RC encoders
prefidiragg = cat(1,oriparams_fixedgazeagg.prefiori);
prefioriagg = cat(1,ori4params_fixedgazeagg.prefiori4);
sigori = cat(1,ori4params_fixedgazeagg.Pkw_ori4)<0.05;
Rori = cat(1,ori4params_fixedgazeagg.Rori4);

hcprefidir_wICenc = zeros(4, Ndirs);
hcprefiori_wICenc = zeros(4, Noris);
hcsigprefidir_wICenc = zeros(4, Ndirs);
hcsigprefiori_wICenc = zeros(4, Noris);
for ienc = 1:4
    switch ienc
        case 1
            neuoi = cat(1,ICsig_fixedgazeagg.ICwcfg0_presentations.ICencoder2)==1;
            neuoi = ismember(cat(1,ICsig_fixedgazeagg.ICwcfg0_presentations.sigmcBK), [0 0 0 1], 'rows');
        case 2
            neuoi = cat(1,ICsig_fixedgazeagg.ICwcfg1_presentations.ICencoder2)==1;
            neuoi = ismember(cat(1,ICsig_fixedgazeagg.ICwcfg1_presentations.sigmcBK), [0 0 0 1], 'rows');
        case 3
            neuoi = cat(1,ICsig_fixedgazeagg.ICwcfg0_presentations.ICencoder1)==1;
            neuoi = ismember(cat(1,ICsig_fixedgazeagg.ICwcfg0_presentations.sigmcBK), [0 0 0 1], 'rows');
        case 4
            neuoi = cat(1,ICsig_fixedgazeagg.ICwcfg1_presentations.ICencoder1)==1;
            neuoi = ismember(cat(1,ICsig_fixedgazeagg.ICwcfg1_presentations.sigmcBK), [0 0 0 1], 'rows');
    end
    hcprefidir_wICenc(ienc,:) = histcounts(prefidiragg(neuoi), 0.5:Ndirs+0.5);
    hcprefiori_wICenc(ienc,:) = histcounts(prefioriagg(neuoi), 0.5:Noris+0.5);
    hcsigprefidir_wICenc(ienc,:) = histcounts(prefidiragg(sigori & neuoi), 0.5:Ndirs+0.5);
    hcsigprefiori_wICenc(ienc,:) = histcounts(prefioriagg(sigori & neuoi), 0.5:Noris+0.5);
end

hcprefidir_kICenc = zeros(4, Ndirs);
hcprefiori_kICenc = zeros(4, Noris);
hcsigprefidir_kICenc = zeros(4, Ndirs);
hcsigprefiori_kICenc = zeros(4, Noris);
for ienc = 1:4
    switch ienc
        case 1
            neuoi = cat(1,ICsig_fixedgazeagg.ICkcfg0_presentations.ICencoder2)==1;
            neuoi = ismember(cat(1,ICsig_fixedgazeagg.ICkcfg0_presentations.sigmcBK), [0 0 0 1], 'rows');
        case 2
            neuoi = cat(1,ICsig_fixedgazeagg.ICkcfg1_presentations.ICencoder2)==1;
            neuoi = ismember(cat(1,ICsig_fixedgazeagg.ICkcfg1_presentations.sigmcBK), [0 0 0 1], 'rows');
        case 3
            neuoi = cat(1,ICsig_fixedgazeagg.ICkcfg0_presentations.ICencoder1)==1;
            neuoi = ismember(cat(1,ICsig_fixedgazeagg.ICkcfg0_presentations.sigmcBK), [0 0 0 1], 'rows');
        case 4
            neuoi = cat(1,ICsig_fixedgazeagg.ICkcfg1_presentations.ICencoder1)==1;
            neuoi = ismember(cat(1,ICsig_fixedgazeagg.ICkcfg1_presentations.sigmcBK), [0 0 0 1], 'rows');
    end
    hcprefidir_kICenc(ienc,:) = histcounts(prefidiragg(neuoi), 0.5:Ndirs+0.5);
    hcprefiori_kICenc(ienc,:) = histcounts(prefioriagg(neuoi), 0.5:Noris+0.5);
    hcsigprefidir_kICenc(ienc,:) = histcounts(prefidiragg(sigori & neuoi), 0.5:Ndirs+0.5);
    hcsigprefiori_kICenc(ienc,:) = histcounts(prefioriagg(sigori & neuoi), 0.5:Noris+0.5);
end

hcprefidir = histcounts(prefidiragg, 0.5:Ndirs+0.5);
hcprefiori = histcounts(prefioriagg, 0.5:Noris+0.5);
hcsigprefidir = histcounts(prefidiragg(sigori), 0.5:Ndirs+0.5);
hcsigprefiori = histcounts(prefioriagg(sigori), 0.5:Noris+0.5);

fs=10;
ytl = 0:45:180-1; ylab = 'IC Orientation';
figure; 
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'Fixed Gaze', 'edgecolor', 'none')
subplot(2,2,1); imagesc(hcprefidir_wICenc./sum(hcprefidir_wICenc,2))
colorbar; colormap jet
set(gca, 'FontSize', fs, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Direction')
title('White Inducers Black Background')
subplot(2,2,2); imagesc(hcprefiori_wICenc./sum(hcprefiori_wICenc,2))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')
subplot(2,2,3); imagesc(hcprefidir_kICenc./sum(hcprefidir_kICenc,2))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Direction')
title('Black Inducers White Background')
subplot(2,2,4); imagesc(hcprefiori_kICenc./sum(hcprefiori_kICenc,2))
colorbar; colormap jet
set(gca, 'XTick', 1:Ndirs, 'XTickLabel', dirvec, 'YTick', 1:4, 'YTickLabel', ytl)
ylabel(ylab)
xlabel('Preferred Orientation')

figure; 
annotation('textbox', [0.1 0.91 0.9 0.1], 'string', 'Fixed Gaze: Significantly Tuned Neurons', 'edgecolor', 'none')
subplot(2,2,1); imagesc(hcsigprefidir_wICenc./sum(hcsigprefidir_wICenc,2))
colorbar; colormap jet
xlabel('Preferred Direction')
title('White Inducers Black Background')
subplot(2,2,2); imagesc(hcsigprefiori_wICenc./sum(hcsigprefiori_wICenc,2))
colorbar; colormap jet
xlabel('Preferred Orientation')
subplot(2,2,3); imagesc(hcsigprefidir_kICenc./sum(hcsigprefidir_kICenc,2))
colorbar; colormap jet
xlabel('Preferred Direction')
title('Black Inducers White Background')
subplot(2,2,4); imagesc(hcsigprefiori_kICenc./sum(hcsigprefiori_kICenc,2))
colorbar; colormap jet
xlabel('Preferred Orientation')

%% shuffle CI
iprobe = find( strcmp('C', probes) );
opt678 = true;
if opt678
prefidiragg = oriparams_fixedgazeagg(iprobe).prefiori678;
prefioriagg = ori4params_fixedgazeagg(iprobe).prefiori4678;
sigori = ori4params_fixedgazeagg(iprobe).Pkw_ori4678<0.05;
% Rori = ( oriparams_fixedgazeagg(iprobe).Rori678(:,1:4) + oriparams_fixedgazeagg(iprobe).Rori678(:,5:8) )/2;
% [~,prefioriagg] = max(Rori, [], 2);
Rori = ori4params_fixedgazeagg(iprobe).Rori4678(:,1:4);
else
prefidiragg = oriparams_fixedgazeagg(iprobe).prefiori;
prefioriagg = ori4params_fixedgazeagg(iprobe).prefiori4;
sigori = ori4params_fixedgazeagg(iprobe).Pkw_ori4<0.05;
Rori = ori4params_fixedgazeagg(iprobe).Rori4(:,1:4);
end

whichcol = 'w';
switch whichcol
    case 'w'
ICpreforivec = [];
ICgroupvec = [];
for ienc = 1:4
    switch ienc
        case 1
            neuoi = ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).ICencoder2==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
        case 2
            neuoi = ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).ICencoder2==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
        case 3
            neuoi = ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).ICencoder1==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
        case 4
            neuoi = ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).ICencoder1==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    end
ICpreforivec = cat(1, ICpreforivec, prefioriagg(neuoi));
ICgroupvec = cat(1, ICgroupvec, ienc*ones(nnz(neuoi),1));
end
    case 'k'
ICpreforivec = [];
ICgroupvec = [];
for ienc = 1:4
    switch ienc
        case 1
            neuoi = ICsig_fixedgazeagg.ICkcfg0_presentations(iprobe).ICencoder2==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICkcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
        case 2
            neuoi = ICsig_fixedgazeagg.ICkcfg1_presentations(iprobe).ICencoder2==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICkcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
        case 3
            neuoi = ICsig_fixedgazeagg.ICkcfg0_presentations(iprobe).ICencoder1==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICkcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
        case 4
            neuoi = ICsig_fixedgazeagg.ICkcfg1_presentations(iprobe).ICencoder1==1;
            neuoi = ismember(ICsig_fixedgazeagg.ICkcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    end
ICpreforivec = cat(1, ICpreforivec, prefioriagg(neuoi));
ICgroupvec = cat(1, ICgroupvec, ienc*ones(nnz(neuoi),1));
end
end

% kruskalwallis of preferred orientation across neuron group
% p = kruskalwallis(ICpreforivec, ICgroupvec);

% figure; histogram2(ICpreforivec, ICgroupvec, 'displaystyle', 'tile')

hc = histcounts2(ICgroupvec, ICpreforivec, 0.5:1:4.5, 0.5:1:Noris+.5);
figure; imagesc(hc./sum(hc,2)); colorbar
xlabel('Grating Pref Ori')
ylabel('IC Pref Ori')

Nshuf = 1000;
hc_shuf = NaN(4, Noris, Nshuf);
for ishuf = 1:Nshuf
hc_shuf(:,:,ishuf) = histcounts2(ICgroupvec, ICpreforivec(randperm(length(ICpreforivec))), 0.5:1:4.5, 0.5:1:Noris+.5);
end
hc_shufmean = squeeze(mean(hc_shuf,3));
hc_shufnegerr = squeeze(mean(hc_shuf,3) - prctile(hc_shuf,92.5,3));
hc_shufposerr = squeeze(prctile(hc_shuf,97.5,3) - mean(hc_shuf,3));

figure
for ienc = 1:4
    subplot(2,2,ienc)
    hold all
    errorbar(1:Noris, hc_shufmean(ienc,:), hc_shufnegerr(ienc,:), hc_shufposerr(ienc,:), 'k.-')
    plot(1:Noris, hc(ienc,:), 'r-', 'LineWidth', 3)
    title(sprintf('ICenc%d', ienc))
end

%% friedman across sessions: proportion of IC-encoders in each orientation
iprobe = find( strcmp('C', probes) );
opt678 = true;
if opt678
prefidiragg = oriparams_fixedgazeagg(iprobe).prefiori678;
prefioriagg = ori4params_fixedgazeagg(iprobe).prefiori4678;
sigori = ori4params_fixedgazeagg(iprobe).Pkw_ori4678<0.05;
% Rori = ( oriparams_fixedgazeagg(iprobe).Rori(:,1:4) + oriparams_fixedgazeagg(iprobe).Rori(:,5:8) )/2;
% [~,prefioriagg] = max(Rori, [], 2);
Rori = ori4params_fixedgazeagg(iprobe).Rori4678;
else
prefidiragg = oriparams_fixedgazeagg(iprobe).prefiori;
prefioriagg = ori4params_fixedgazeagg(iprobe).prefiori4;
sigori = ori4params_fixedgazeagg(iprobe).Pkw_ori4<0.05;
Rori = ori4params_fixedgazeagg(iprobe).Rori4;
end

hpICenc = cell(8,1);
hpICenc_sigori = cell(8,1);
hpfried = zeros(8,1);
hpfried_sigori = zeros(8,1);
propICenc = cell(8,1);
propICenc_sigori = cell(8,1);
ppfried = zeros(8,1);
ppfried_sigori = zeros(8,1);
for ienc = 1:8
switch ienc
    case 1
        neuoi = ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 2
        neuoi = ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 3
        neuoi = ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    case 4
        neuoi = ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    case 5
        neuoi = ICsig_fixedgazeagg.ICkcfg0_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICkcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 6
        neuoi = ICsig_fixedgazeagg.ICkcfg1_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICkcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 7
        neuoi = ICsig_fixedgazeagg.ICkcfg0_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICkcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    case 8
        neuoi = ICsig_fixedgazeagg.ICkcfg1_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICkcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
end

hpICenc{ienc} = NaN(Nsessions,Noris);
hpICenc_sigori{ienc} = NaN(Nsessions,Noris);
propICenc{ienc} = NaN(Nsessions,Noris);
propICenc_sigori{ienc} = NaN(Nsessions,Noris);

for ises = 1:Nsessions
    tempsesneu = sesneuoiagg{iprobe}==ises;
    hpICenc{ienc}(ises,:) = histcounts(prefioriagg(tempsesneu & neuoi), 0.5:Noris+0.5, 'normalization', 'probability');
    hpICenc_sigori{ienc}(ises,:) = histcounts(prefioriagg(tempsesneu & sigori & neuoi), 0.5:Noris+0.5, 'normalization', 'probability');
    
    hcprefiori = histcounts(prefioriagg(tempsesneu), 0.5:Noris+0.5);
    hcprefiori_ICenc = histcounts(prefioriagg(tempsesneu & neuoi), 0.5:Noris+0.5);
    propICenc{ienc}(ises,:) = hcprefiori_ICenc./hcprefiori;
    
    hcprefiori = histcounts(prefioriagg(tempsesneu & sigori), 0.5:Noris+0.5);
    hcprefiori_ICenc = histcounts(prefioriagg(tempsesneu & sigori & neuoi), 0.5:Noris+0.5);
    propICenc_sigori{ienc}(ises,:) = hcprefiori_ICenc./hcprefiori;
end

tempmat = hpICenc{ienc};
tempmat = tempmat(all(~isnan(tempmat),2),:);
[p,tbl,stats] = friedman(tempmat, 1,'off');
hpfried(ienc) = p;

tempmat = hpICenc_sigori{ienc};
tempmat = tempmat(all(~isnan(tempmat),2),:);
[p,tbl,stats] = friedman(tempmat, 1,'off');
hpfried_sigori(ienc) = p;

tempmat = propICenc{ienc};
tempmat = tempmat(all(~isnan(tempmat),2),:);
[p,tbl,stats] = friedman(tempmat, 1,'off');
ppfried(ienc) = p;

tempmat = propICenc_sigori{ienc};
tempmat = tempmat(all(~isnan(tempmat),2),:);
[p,tbl,stats] = friedman(tempmat, 1,'off');
ppfried_sigori(ienc) = p;

end
disp([hpfried hpfried_sigori ppfried ppfried_sigori])

for ienc = 4%1:4
tempmat = hpICenc{ienc};
tempmat = tempmat(all(~isnan(tempmat),2),:);
[p,tbl,stats] = friedman(tempmat);
hpfried(ienc) = p;
[c,m,h] = multcompare(stats);
title(sprintf('%d p=%.4f', ienc, p))
end

% figure; plot(hpICenc{4}')

%% jackknife proportion of IC-encoders in each orientation
hpICenc = cell(8,1);
hpICenc_sigori = cell(8,1);
hpfried = zeros(8,1);
hpfried_sigori = zeros(8,1);
propICenc = cell(8,1);
propICenc_sigori = cell(8,1);
Pfried = zeros(8,1);
Pfried_sigori = zeros(8,1);
for ienc = 1:4
switch ienc
    case 1
        neuoi = ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 2
        neuoi = ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 3
        neuoi = ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    case 4
        neuoi = ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    case 5
        neuoi = ICsig_fixedgazeagg.ICkcfg0_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICkcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 6
        neuoi = ICsig_fixedgazeagg.ICkcfg1_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICkcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 7
        neuoi = ICsig_fixedgazeagg.ICkcfg0_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICkcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    case 8
        neuoi = ICsig_fixedgazeagg.ICkcfg1_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICkcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
end
neuoind = find(neuoi);


hpICenc{ienc} = NaN(numel(neuoind),Noris);
hpICenc_sigori{ienc} = NaN(numel(neuoind),Noris);
propICenc{ienc} = NaN(numel(neuoind),Noris);
propICenc_sigori{ienc} = NaN(numel(neuoind),Noris);

for ci = 1:numel(neuoind)
    tempexclneu = true(size(neuoi));
    tempexclneu(neuoind(ci)) = false;
    hpICenc{ienc}(ci,:) = histcounts(prefioriagg(tempexclneu & neuoi), 0.5:Noris+0.5, 'normalization', 'probability');
    hpICenc_sigori{ienc}(ci,:) = histcounts(prefioriagg(tempexclneu & sigori & neuoi), 0.5:Noris+0.5, 'normalization', 'probability');
    
    hcprefiori = histcounts(prefioriagg(tempexclneu), 0.5:Noris+0.5);
    hcprefiori_ICenc = histcounts(prefioriagg(tempexclneu & neuoi), 0.5:Noris+0.5);
    propICenc{ienc}(ci,:) = hcprefiori_ICenc./hcprefiori;
    
    hcprefiori = histcounts(prefioriagg(tempexclneu & sigori), 0.5:Noris+0.5);
    hcprefiori_ICenc = histcounts(prefioriagg(tempexclneu & sigori & neuoi), 0.5:Noris+0.5);
    propICenc_sigori{ienc}(ci,:) = hcprefiori_ICenc./hcprefiori;
end

tempmat = hpICenc{ienc};
tempmat = tempmat(all(~isnan(tempmat),2),:);
[p,tbl,stats] = friedman(tempmat, 1,'off');
hpfried(ienc) = p;

tempmat = hpICenc_sigori{ienc};
tempmat = tempmat(all(~isnan(tempmat),2),:);
[p,tbl,stats] = friedman(tempmat, 1,'off');
hpfried_sigori(ienc) = p;

tempmat = propICenc{ienc};
tempmat = tempmat(all(~isnan(tempmat),2),:);
[p,tbl,stats] = friedman(tempmat, 1,'off');
ppfried(ienc) = p;

tempmat = propICenc_sigori{ienc};
tempmat = tempmat(all(~isnan(tempmat),2),:);
[p,tbl,stats] = friedman(tempmat, 1,'off');
ppfried_sigori(ienc) = p;
end
disp([hpfried hpfried_sigori ppfried ppfried_sigori])

for ienc = 1:4
tempmat = hpICenc{ienc};
tempmat = tempmat(all(~isnan(tempmat),2),:);
[p,tbl,stats] = friedman(tempmat);
hpfried(ienc) = p;
[c,m,h] = multcompare(stats);
title(sprintf('%d p=%.4f', ienc, p))
end

figure; plot(hpICenc{4}')


%% bootstrapped confidence interval: 
iprobe = find( strcmp('C', probes) );
opt678 = true;
if opt678
prefidiragg = oriparams_fixedgazeagg(iprobe).prefiori678;
prefioriagg = ori4params_fixedgazeagg(iprobe).prefiori4678;
sigori = ori4params_fixedgazeagg(iprobe).Pkw_ori4678<0.05;
% Rori = ( oriparams_fixedgazeagg(iprobe).Rori(:,1:4) + oriparams_fixedgazeagg(iprobe).Rori(:,5:8) )/2;
% [~,prefioriagg] = max(Rori, [], 2);
Rori = ori4params_fixedgazeagg(iprobe).Rori4678;
else
prefidiragg = oriparams_fixedgazeagg(iprobe).prefiori;
prefioriagg = ori4params_fixedgazeagg(iprobe).prefiori4;
sigori = ori4params_fixedgazeagg(iprobe).Pkw_ori4<0.05;
Rori = ori4params_fixedgazeagg(iprobe).Rori4;
end

NICenc = zeros(8,1);
NICenc_sigori = zeros(8,1);

hpICenc = cell(8,1);
hpICenc_sigori = cell(8,1);
propICenc = cell(8,1);
propICenc_sigori = cell(8,1);

% rows 1-3 are mean, negative errorbar, positive errorbar
Nboot = 1000;
hpbootCI = cell(8,1);
hpbootCI_sigori = cell(8,1);
propbootCI = cell(8,1);
propbootCI_sigori = cell(8,1);
for ienc = 1:8
switch ienc
    case 1
        neuoi = ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 2
        neuoi = ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 3
        neuoi = ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICwcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    case 4
        neuoi = ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICwcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    case 5
        neuoi = ICsig_fixedgazeagg.ICkcfg0_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICkcfg0_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 6
        neuoi = ICsig_fixedgazeagg.ICkcfg1_presentations(iprobe).ICencoder2==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICkcfg1_presentations(iprobe).sigmcBK, [0 0 0 1], 'rows');
    case 7
        neuoi = ICsig_fixedgazeagg.ICkcfg0_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICkcfg0_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
    case 8
        neuoi = ICsig_fixedgazeagg.ICkcfg1_presentations(iprobe).ICencoder1==1;
        neuoi = ismember(ICsig_fixedgazeagg.ICkcfg1_presentations(iprobe).sigmcBK, [1 0 0 0], 'rows');
end
neuoind = find(neuoi);
sigoriind = find(sigori);
sigprefioriagg = prefioriagg(sigori);

NICenc(ienc) = nnz(neuoi);
NICenc_sigori(ienc) = nnz(sigori & neuoi);

hpICenc{ienc} = histcounts(prefioriagg(neuoi), 0.5:Noris+0.5, 'normalization', 'probability');
hpICenc_sigori{ienc} = histcounts(prefioriagg(sigori & neuoi), 0.5:Noris+0.5, 'normalization', 'probability');

hcprefiori = histcounts(prefioriagg, 0.5:Noris+0.5);
hcprefiori_ICenc = histcounts(prefioriagg(neuoi), 0.5:Noris+0.5);
propICenc{ienc} = hcprefiori_ICenc./hcprefiori;

hcprefiori = histcounts(prefioriagg(sigori), 0.5:Noris+0.5);
hcprefiori_ICenc = histcounts(prefioriagg(sigori & neuoi), 0.5:Noris+0.5);
propICenc_sigori{ienc} = hcprefiori_ICenc./hcprefiori;

hpboot = NaN(Nboot,Noris);
hpboot_sigori = NaN(Nboot,Noris);
propboot = NaN(Nboot,Noris);
propboot_sigori = NaN(Nboot,Noris);
hcprefiori = histcounts(prefioriagg, 0.5:Noris+0.5);
hcprefiori_sigori = histcounts(prefioriagg(sigori), 0.5:Noris+0.5);
% tic
for iboot = 1:Nboot
    tempneu = false(size(sigori));
    tempneu( randperm(length(prefioriagg), nnz(neuoi)) ) = true;
    hpboot(iboot,:) = histcounts(prefioriagg(tempneu), 0.5:Noris+0.5, 'normalization', 'probability');
    
    tempneusig = false(size(sigprefioriagg));
    tempneusig( randperm(length(sigprefioriagg), nnz(neuoi(sigoriind))) ) = true;
    hpboot_sigori(iboot,:) = histcounts(sigprefioriagg(tempneusig), 0.5:Noris+0.5, 'normalization', 'probability');
    
    hcprefiori_ICenc = histcounts(prefioriagg(tempneu), 0.5:Noris+0.5);
    propboot(iboot,:) = hcprefiori_ICenc./hcprefiori;
    
    hcprefiori_ICenc = histcounts(sigprefioriagg(tempneusig), 0.5:Noris+0.5);
    propboot_sigori(iboot,:) = hcprefiori_ICenc./hcprefiori_sigori;
end
% toc
hpbootCI{ienc} = [nanmean(hpboot,1); nanmean(hpboot,1)-prctile(hpboot,2.5,1); prctile(hpboot,97.5,1)-nanmean(hpboot,1)];
hpbootCI_sigori{ienc} = [nanmean(hpboot_sigori,1); nanmean(hpboot_sigori,1)-prctile(hpboot_sigori,2.5,1); prctile(hpboot_sigori,97.5,1)-nanmean(hpboot_sigori,1)];

propbootCI{ienc} = [nanmean(propboot,1); nanmean(propboot,1)-prctile(propboot,2.5,1); prctile(propboot,97.5,1)-nanmean(propboot,1)];
propbootCI_sigori{ienc} = [nanmean(propboot_sigori,1); nanmean(propboot_sigori,1)-prctile(propboot_sigori,2.5,1); prctile(propboot_sigori,97.5,1)-nanmean(propboot_sigori,1)];

end

fs = 12;
enctitle = {'wICenc0', 'wICenc45', 'wICenc90', 'wICenc135', 'kICenc0', 'kICenc45', 'kICenc90', 'kICenc135'};

figure('Position', [0 0 1000 1000])
annotation('textbox', [0.1 0.9 0.9 0.1], 'string', 'Fixed Gaze', 'edgecolor', 'none', 'fontsize', fs)
for ienc = 1:8
subplot(4,4,ienc)
hold all
errorbar(1:Noris, hpbootCI{ienc}(1,:), hpbootCI{ienc}(2,:), hpbootCI{ienc}(3,:) , 'k')
plot(1:Noris, hpICenc{ienc}, 'r', 'LineWidth', 2)
indsigori = find(hpICenc{ienc}>hpbootCI{ienc}(1,:)+hpbootCI{ienc}(3,:));
plot(indsigori, hpICenc{ienc}(indsigori), 'ro', 'MarkerFaceColor', 'r', 'LineWidth', 2)
set(gca, 'FontSize', fs, 'XTick', 1:Noris, 'XTickLabel', orivec)
title(sprintf('%s N=%d\npref ori dist', enctitle{ienc}, NICenc(ienc) ))
ylabel('Probability')
xlim([0.5 Noris+0.5])
xlabel('Pref Ori')

subplot(4,4,ienc+8)
hold all
errorbar(1:Noris, propbootCI{ienc}(1,:), propbootCI{ienc}(2,:), propbootCI{ienc}(3,:) , 'k')
plot(1:Noris, propICenc{ienc}, 'r', 'LineWidth', 2)
indsigori = find(propICenc{ienc}>propbootCI{ienc}(1,:)+propbootCI{ienc}(3,:));
plot(indsigori, propICenc{ienc}(indsigori), 'ro', 'MarkerFaceColor', 'r', 'LineWidth', 2)
set(gca, 'FontSize', fs, 'XTick', 1:Noris, 'XTickLabel', orivec)
title(sprintf('%s N=%d\nprop for each pref ori', enctitle{ienc}, NICenc_sigori(ienc) ))
ylabel('Proportion')
xlim([0.5 Noris+0.5])
xlabel('Pref Ori')
end

figure('Position', [0 0 1000 1000])
annotation('textbox', [0.1 0.9 0.9 0.1], 'string', 'Fixed Gaze: Significantly Tuned Neurons', 'edgecolor', 'none', 'fontsize', fs)
for ienc = 1:8
subplot(4,4,ienc)
hold all
errorbar(1:Noris, hpbootCI_sigori{ienc}(1,:), hpbootCI_sigori{ienc}(2,:), hpbootCI_sigori{ienc}(3,:) , 'k')
plot(1:Noris, hpICenc_sigori{ienc}, 'r', 'LineWidth', 2)
indsigori = find(hpICenc_sigori{ienc}>hpbootCI_sigori{ienc}(1,:)+hpbootCI_sigori{ienc}(3,:));
plot(indsigori, hpICenc_sigori{ienc}(indsigori), 'ro', 'MarkerFaceColor', 'r', 'LineWidth', 2)
set(gca, 'FontSize', fs, 'XTick', 1:Noris, 'XTickLabel', orivec)
title(sprintf('%s N=%d\npref ori dist', enctitle{ienc}, NICenc(ienc) ))
ylabel('Probability')
xlim([0.5 Noris+0.5])
xlabel('Pref Ori')

subplot(4,4,ienc+8)
hold all
errorbar(1:Noris, propbootCI_sigori{ienc}(1,:), propbootCI_sigori{ienc}(2,:), propbootCI_sigori{ienc}(3,:) , 'k')
plot(1:Noris, propICenc_sigori{ienc}, 'r', 'LineWidth', 2)
indsigori = find(propICenc_sigori{ienc}>propbootCI_sigori{ienc}(1,:)+propbootCI_sigori{ienc}(3,:));
plot(indsigori, propICenc_sigori{ienc}(indsigori), 'ro', 'MarkerFaceColor', 'r', 'LineWidth', 2)
set(gca, 'FontSize', fs, 'XTick', 1:Noris, 'XTickLabel', orivec)
title(sprintf('%s N=%d\nprop for each pref ori', enctitle{ienc}, NICenc_sigori(ienc) ))
ylabel('Proportion')
xlim([0.5 Noris+0.5])
xlabel('Pref Ori')
end

%% RFspin distribution of IC-encoders
iprobe = find( strcmp('C', probes) );
CRFindagg = RFCI_fixedgazeagg(iprobe).RFindclassic;
sigCRFagg = RFCI_fixedgazeagg(iprobe).Pkw_rfclassic<0.05;

ICblocks2agg = {'ICwcfg1'};
ICblocks2agg = {'ICwcfg1', 'ICwcfg0'};
ICblocks2agg = {'ICkcfg1', 'ICkcfg0'};
ICblocks2agg = {'ICwcfg1', 'ICwcfg0', 'ICkcfg1', 'ICkcfg0'};
ICencoderagg = false(size(CRFindspinagg));
for b = 1:numel(ICblocks2agg)
    whichblock = [ICblocks2agg{b} '_presentations'];
    ICencoderagg = ICencoderagg | (ICsig_fixedgazeagg.(whichblock)(iprobe).ICencoder==1);
end

propICencoder = zeros(Nsessions,Nrfs);
propICencoder_sigCRF = zeros(Nsessions,Nrfs);
for ises = 1:Nsessions
    tempneuoi = sesneuoiagg{iprobe}==ises;
    hcCRF = histcounts(CRFindagg(tempneuoi), 0.5:Nrfs+0.5);
    hcCRF_ICenc = histcounts(CRFindagg(tempneuoi & ICencoderagg), 0.5:Nrfs+0.5);
    propICencoder(ises,:) = hcCRF_ICenc./hcCRF;
    
    hcsigCRF = histcounts(CRFindagg(tempneuoi & sigCRFagg), 0.5:Nrfs+0.5);
    hcsigCRF_ICenc = histcounts(CRFindagg(tempneuoi & sigCRFagg & ICencoderagg), 0.5:Nrfs+0.5);
    propICencoder_sigCRF(ises,:) = hcsigCRF_ICenc./hcsigCRF;
end

figure; hold all
histogram(CRFindagg, 'normalization', 'probability')
histogram(CRFindagg(ICencoderagg), 'normalization', 'probability')
legend({'All', 'IC-encoder'})
xlabel('RF position')
ylabel('Probability')

figure; hold all
histogram(CRFindagg(sigCRFagg), 'normalization', 'probability')
histogram(CRFindagg(sigCRFagg & ICencoderagg), 'normalization', 'probability')
legend({'All', 'IC-encoder'})
xlabel('RF position')
ylabel('Probability')
title('Neurons with Significant CRF')

tempmat =  propICencoder(all(~isnan(propICencoder),2), :);
[p,tbl,stats]=friedman(tempmat);
[h,m,c]=multcompare(stats);

figure; plot(propICencoder')
xlabel('RF position')
ylabel('Proportion IC-encoder')
title(sprintf('Friedman Test p=%.4f', p))

tempmat =  propICencoder_sigCRF(all(~isnan(propICencoder_sigCRF),2), :);
[p,tbl,stats]=friedman(tempmat);
[h,m,c]=multcompare(stats);

figure; plot(propICencoder_sigCRF')
xlabel('RF position')
ylabel('Proportion IC-encoder')
title(sprintf('Neurons with Significant CRF\nFriedman Test p=%.4f', p))

%% RFspin distribution of IC-encoders
iprobe = find( strcmp('C', probes) );
CRFindspinagg = RFCIspin_fixedgazeagg(iprobe).RFindclassic;
sigCRFspinagg = RFCIspin_fixedgazeagg(iprobe).Pfried_rfclassic<0.05;

ICblocks2agg = {'ICwcfg1'};
ICblocks2agg = {'ICwcfg1', 'ICwcfg0'};
ICblocks2agg = {'ICkcfg1', 'ICkcfg0'};
ICblocks2agg = {'ICwcfg1', 'ICwcfg0', 'ICkcfg1', 'ICkcfg0'};
ICencoderagg = false(size(CRFindspinagg));
for b = 1:numel(ICblocks2agg)
    whichblock = [ICblocks2agg{b} '_presentations'];
    ICencoderagg = ICencoderagg | (ICsig_fixedgazeagg.(whichblock)(iprobe).ICencoder==1);
end

propICencoder = zeros(Nsessions,Nrfs);
propICencoder_sigCRF = zeros(Nsessions,Nrfs);
for ises = 1:Nsessions
    tempneuoi = sesneuoiagg{iprobe}==ises;
    hcCRF = histcounts(CRFindspinagg(tempneuoi), 0.5:Nrfs+0.5);
    hcCRF_ICenc = histcounts(CRFindspinagg(tempneuoi & ICencoderagg), 0.5:Nrfs+0.5);
    propICencoder(ises,:) = hcCRF_ICenc./hcCRF;
    
    hcsigCRF = histcounts(CRFindspinagg(tempneuoi & sigCRFspinagg), 0.5:Nrfs+0.5);
    hcsigCRF_ICenc = histcounts(CRFindspinagg(tempneuoi & sigCRFspinagg & ICencoderagg), 0.5:Nrfs+0.5);
    propICencoder_sigCRF(ises,:) = hcsigCRF_ICenc./hcsigCRF;
end

figure; hold all
histogram(CRFindspinagg, 'normalization', 'probability')
histogram(CRFindspinagg(ICencoderagg), 'normalization', 'probability')
legend({'All', 'IC-encoder'})
xlabel('RF position')
ylabel('Probability')

figure; hold all
histogram(CRFindspinagg(sigCRFspinagg), 'normalization', 'probability')
histogram(CRFindspinagg(sigCRFspinagg & ICencoderagg), 'normalization', 'probability')
legend({'All', 'IC-encoder'})
xlabel('RF position')
ylabel('Probability')
title('Neurons with Significant CRF')

tempmat =  propICencoder(all(~isnan(propICencoder),2), :);
[p,tbl,stats]=friedman(tempmat);
[h,m,c]=multcompare(stats);

figure; plot(propICencoder')
xlabel('RF position')
ylabel('Proportion IC-encoder')
title(sprintf('Friedman Test p=%.4f', p))

tempmat =  propICencoder_sigCRF(all(~isnan(propICencoder_sigCRF),2), :);
[p,tbl,stats]=friedman(tempmat);
[h,m,c]=multcompare(stats);

figure; plot(propICencoder_sigCRF')
xlabel('RF position')
ylabel('Proportion IC-encoder')
title(sprintf('Neurons with Significant CRF\nFriedman Test p=%.4f', p))

%% size tuning of IC-encoders
