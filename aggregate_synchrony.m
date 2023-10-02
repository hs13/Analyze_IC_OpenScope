ICtrialtypes = [0	101	105	106	107	109	110	111	506	511	1105	1109	1201	1299	1301	1302	1303	1304	1305	1306	1307	1308];
nwbsessions = {'sub_1171903433','sub_1172968426','sub_1172969394','sub_1174569641', ...
    'sub_1175512783','sub_1176214862','sub_1177693342','sub_1181314060', ...
    'sub_1181585608','sub_1182593231','sub_1183369803','sub_1186544726', ...
    'sub_1189891322','sub_1194090570'};

ICsig_CctxRSagg = struct();
syncpsthavgagg = struct();
pairijagg = cell(numel(nwbsessions), 1);

whichblock = 'ICwcfg1_presentations';
syncpsthavgagg.(whichblock) = cell(numel(nwbsessions), 1);
tic
for ises = 1:numel(nwbsessions)
pathpp = ['S:\OpenScopeData\000248\postprocessed\' nwbsessions{ises} '\'];
load([pathpp 'visresponses_probeC.mat'])
load([pathpp 'syncpsthavg_probeC.mat'])
% 'CindV1RS', 'neuindV1RS', 'pairij', 'pairCind', 'pairneuind',
% 'blocktrialtypes', 'syncpsthavg'
pairijagg{ises} = pairij;

ICsig.(whichblock).indin13 = ICsig.(whichblock).indin1 | ICsig.(whichblock).indin3;
ICsig.(whichblock).indin14 = ICsig.(whichblock).indin1 | ICsig.(whichblock).indin4;
ICsig.(whichblock).indin23 = ICsig.(whichblock).indin2 | ICsig.(whichblock).indin3;
ICsig.(whichblock).indin24 = ICsig.(whichblock).indin2 | ICsig.(whichblock).indin4;

ICsig.(whichblock).indin34 = ICsig.(whichblock).indin3 | ICsig.(whichblock).indin4;

Nneurons = numel(meanFRvec);
ICsigfields = fieldnames(ICsig.(whichblock));
for f = 1:numel(ICsigfields)
    whichfield = ICsigfields{f};
    if size(ICsig.(whichblock).(whichfield),1)==Nneurons
        ICsig_CctxRSagg(ises).(whichblock).(whichfield) = ICsig.(whichblock).(whichfield)(CindV1RS,:);
    end
end

syncpsthavgagg.(whichblock){ises} = syncpsthavg.(whichblock);
toc
end


%% average all pairs IC vs RC trials
% kerwinhalf = 25; kersigma = 12/sqrt(2);
kerwinhalf = 15; kersigma = 3;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));

tempsync = cat(3, syncpsthavgagg.(whichblock){:});
synctloi = syncpsthavg.syncpsthtli>0 & syncpsthavg.syncpsthtli<=400;
Rsync = squeeze(mean(tempsync(synctloi,:,:),1));
[P,TABLE,STATS] = friedman( Rsync(ismember(ICtrialtypes, [106 107 110 111]),:)' );
figure
multcompare(STATS)
signrank( mean(Rsync(ismember(ICtrialtypes, [106 111]),:),1), mean(Rsync(ismember(ICtrialtypes, [107 110]),:),1) )

% all pairs
figure
hold all
plot(syncpsthavg.syncpsthtli, squeeze(mean( convn(tempsync(:,ismember(ICtrialtypes, [106 111]),:), kergauss, 'same') ,[2 3])), 'Color', [0 0.7 0])
plot(syncpsthavg.syncpsthtli, squeeze(mean( convn(tempsync(:,ismember(ICtrialtypes, [107 110]),:), kergauss, 'same'), [2 3])), 'Color', [1 0.5 0])


%% select pairs synchrony on IC vs RC trials
pltopt =0;
switch pltopt
    % between segment responders: fairly conssitent across sessions
    case 0 % slightly higher earlier (<200ms) -- fairly consistent across sessions
        group1 = {'indin13', 'indin14', 'indin23', 'indin24'};
        group2 = group1;

        % between segment responders
    case 1
        group1 = {'indin1', 'indin1', 'indin2', 'indin2'};
        group2 = {'indin3', 'indin4', 'indin3', 'indin4'};

        % between IC/RC encoders
    case 2
        group1 = {'ICencoder1', 'RCencoder1', 'RCencoder2', 'ICencoder2'};
        group2 = group1;

        % between IC/RC encoders and segment responders
    case 3 % slightly higher later (>200ms)
        group1 = {'ICencoder1', 'RCencoder1', 'RCencoder2', 'ICencoder2'};
        group2 = {'indin13', 'indin14', 'indin23', 'indin24'};

        % between inducerencoders
    case 4
        group1 = {'indenc13', 'indenc14', 'indenc23', 'indenc24'};
        group2 = group1;

        % between IC/RC responders
    case 5
        group1 = {'ICresp1', 'RCresp1', 'RCresp2', 'ICresp2'};
        group2 = group1;

        % between IC/RC encoders and IC/RC responders
    case 6 % late synchronization beginning 200ms? -- no, most sessions actually higher synchrony for RC-encoders
        group1 = {'ICencoder1', 'RCencoder1', 'RCencoder2', 'ICencoder2'};
        group2 = {'ICresp1', 'RCresp1', 'RCresp2', 'ICresp2'};
end

tempsynccell = cell(4,1);
for ises = 1:numel(nwbsessions)
    % if ises==4
    %     continue
    % end
    tempNneurons = length(ICsig_CctxRSagg(ises).ICwcfg1_presentations.ICencoder);
    if size(syncpsthavgagg.(whichblock){ises},3) ~= tempNneurons*(tempNneurons-1)/2
        error('check aggregation')
    end
for itt = 1:numel(tt2p)
    tempsr1 = find( ICsig_CctxRSagg(ises).(whichblock).(group1{itt}) );
    tempsr2 = find( ICsig_CctxRSagg(ises).(whichblock).(group2{itt}) );

    if max(pairijagg{ises}(:)) < max([tempsr1;tempsr2])
        error('check aggregation')
    end
    tempsrpair = sort( combvec(tempsr1', tempsr2')', 2);
    rows2rm = tempsrpair(:,1)==tempsrpair(:,2);
    tempsrpair(rows2rm,:)=[];
    if ~all(ismember(tempsrpair, pairijagg{ises}, 'rows'))
        error('check aggregation')
    end
    pairoi = ismember(pairijagg{ises}, tempsrpair, 'rows');

    typi = ICtrialtypes==tt2p(itt);
    tempsync = squeeze( syncpsthavgagg.(whichblock){ises}(:,typi,pairoi) );
    tempsynccell{itt} = cat(2, tempsynccell{itt}, tempsync);
end
end


tt2p = [106 107 110 111];
ttcol = [0 0.5 0; 0.5 0.25 0; 1 0.5 0; 0 1 0];
figure('Position', xywh)
annotation('textbox',[0 0.9 0.5 0.1], 'String', sprintf('%s %s %s %s\n%s %s %s %s', group1{:}, group2{:}), 'edgecolor', 'none', 'HorizontalAlignment', 'center')
subplot(1,2,1)
hold all
for itt = 1:numel(tt2p)
    tempsync = tempsynccell{itt} ;
    plot(syncpsthavg.syncpsthtli, squeeze(mean( convn(tempsync,kergauss,'same'), 2) ), 'Color', ttcol(itt,:), 'LineWidth', 1)
end
xlim([-400 800])
subplot(1,2,2)
hold all
tempsync = cat(2, tempsynccell{1}, tempsynccell{4});
plot(syncpsthavg.syncpsthtli, squeeze(mean( convn(tempsync,kergauss,'same'), 2) ), 'Color', [0 0.7 0], 'LineWidth', 1)
tempsync = cat(2, tempsynccell{2}, tempsynccell{3});
plot(syncpsthavg.syncpsthtli, squeeze(mean( convn(tempsync,kergauss,'same'), 2) ), 'Color', [1 0.5 0], 'LineWidth', 1)
xlim([-400 800])

synctloi = syncpsthavg.syncpsthtli>0 & syncpsthavg.syncpsthtli<=400;
ICsync = cat(2, tempsynccell{1}, tempsynccell{4});
RICsync = mean(ICsync(synctloi,:),1);
RCsync = cat(2, tempsynccell{2}, tempsynccell{3});
RRCsync = mean(RCsync(synctloi,:),1);
ranksum(RICsync, RRCsync)
ranksum(RICsync, RRCsync, 'tail', 'right')
ranksum(RICsync, RRCsync, 'tail', 'left')

%% T-shaped combination
pltopt = 4;
switch pltopt
    case 1
typi = ICtrialtypes==105;
group1 = {'indin1', 'indin1', 'indin3'}; 
group2 = {'indin3', 'indin4', 'indin4'}; 
    case 2
typi = ICtrialtypes==105;
group1 = {'indin13', 'indin14', 'indin34'}; 
group2 = group1;
    case 3
typi = ICtrialtypes==109;
group1 = {'indin2', 'indin2', 'indin3'}; 
group2 = {'indin4', 'indin3', 'indin4'}; 
    case 4
typi = ICtrialtypes==109;
group1 = {'indin24', 'indin23', 'indin34'}; 
group2 = group1;
end

tempsynccell = cell(3,1);
for ises = 1:numel(nwbsessions)
    for ipair = 1:3
        tempsr1 = find( ICsig_CctxRSagg(ises).(whichblock).(group1{ipair}) );
        tempsr2 = find( ICsig_CctxRSagg(ises).(whichblock).(group2{ipair}) );
        tempsrpair = sort( combvec(tempsr1', tempsr2')', 2);
        %rows2rm = tempsrpair(:,1)==tempsrpair(:,2);
        %tempsrpair(rows2rm,:)=[];
        pairoi = ismember(pairijagg{ises}, tempsrpair, 'rows');

        tempsync = squeeze( syncpsthavgagg.(whichblock){ises}(:,typi,pairoi) ) ;
        tempsynccell{ipair} = cat(2, tempsynccell{ipair}, tempsync);

    end
end


paircol = [1 0 0; 0 1 0; 0 0 1];
figure
hold all
for ipair = 1:3
    tempsync = squeeze(mean( convn(tempsynccell{ipair},kergauss,'same'), 2) ) ;
    plot(syncpsthavg.syncpsthtli, tempsync, 'Color', paircol(ipair,:), 'LineWidth', 1)
end
