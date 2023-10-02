
% calculate for ctx units in Probe C and D that pass quality control
% shuffleCI: for each unit pair, shuffle pairing within trial type.
% 1000 shuffles and get mean, std, 0.5%, 1%, 2.5%, 5%, 50%, 95%, 97.5%, 99%, 99.5% percentile CI
% smoothing bin size 1ms, 5ms, 25ms

nwbsessions = {'sub_1171903433','sub_1172968426','sub_1172969394','sub_1174569641', ...
    'sub_1175512783','sub_1176214862','sub_1177693342','sub_1181314060', ...
    'sub_1181585608','sub_1182593231','sub_1183369803','sub_1186544726', ...
    'sub_1189891322','sub_1194090570'};

for ises = 1:numel(nwbsessions)
clearvars -except ises nwbsessions 
pathpp = ['S:\OpenScopeData\000248\postprocessed\' nwbsessions{ises} '\'];

load([pathpp 'qc_units.mat'])
load([pathpp 'postprocessed.mat'], 'neuallloc')
load([pathpp 'visresponses_probeC.mat'])
load([pathpp 'postprocessed_probeC.mat'])
% neuoind = find(floor(unit_peakch/1000)==iprobe-1);

% neu2anal = contains(neuallloc, 'VISp') & ~contains(neuallloc, 'VISpm') ...
%     & unit_wfdur>0.4 & (unit_isi_violations<0.5 & unit_amplitude_cutoff<0.5 & unit_presence_ratio>0.9);

CindV1RS = find(contains(neuallloc(neuoind), 'VISp') ...
    & unit_wfdur(neuoind)>0.4 & (unit_isi_violations(neuoind)<0.5 & ...
    unit_amplitude_cutoff(neuoind)<0.5 & unit_presence_ratio(neuoind)>0.9));
neuindV1RS = neuoind(CindV1RS);

% nnz(ICsig.ICwcfg1_presentations.ICencoder)
% nnz(ICsig.ICwcfg1_presentations.ICencoder(contains(neuallloc(neuoind), 'VISp')))
% nnz(ICsig.ICwcfg1_presentations.ICencoder(ismember(neuoind, neuanalind)))

Npairs = numel(CindV1RS)*(numel(CindV1RS)-1)/2;
pairij = zeros(Npairs,2);
pairCind = zeros(Npairs,2);
pairneuind = zeros(Npairs,2);

disp([numel(CindV1RS) Npairs])
tloi = psthtli>=-400 & psthtli<800;
blocknames = fieldnames(psth);

syncpsthavg = struct();
syncpsthavg.syncpsthtli = psthtli(tloi);
blocktrialtypes = struct();
% whichblock = 'ICwcfg1_presentations';
for b = 1:numel(blocknames)
    whichblock = blocknames{b};

    % takes 1 minute for 67 V1RS units, takes up 0.55 GB space
    trialtypes = unique(vis.(whichblock).trialorder);
    blocktrialtypes.(whichblock) = trialtypes;
    syncpsthavg.(whichblock) = zeros(length(syncpsthavg.syncpsthtli), numel(trialtypes), Npairs);
    cnt = 0;
    paircnt = 0;
    tic
    for ii = 1:nnz(CindV1RS)-1
        for jj = ii+1:nnz(CindV1RS)
            cnt = cnt+1;
            pairij(cnt,:) = [ii jj];
            pairCind(cnt,:) = CindV1RS([ii jj]);
            pairneuind(cnt,:) = neuindV1RS([ii jj]);
        end
        tempsyncpsth = psth.ICwcfg1_presentations(tloi,:,CindV1RS(ii)) .* psth.ICwcfg1_presentations(tloi,:,CindV1RS(ii+1:nnz(CindV1RS)));
        for typi = 1:numel(trialtypes)
            temppairs = paircnt+1:paircnt+length(ii+1:nnz(CindV1RS));
            temptrials = vis.(whichblock).trialorder==trialtypes(typi);
            syncpsthavg.(whichblock)(:,typi,temppairs) = 1000*mean(tempsyncpsth(:,temptrials,:),2);
        end
        paircnt = paircnt+length(ii+1:nnz(CindV1RS));
    end
    toc

    % sanity check
    if ~( isequal(neuindV1RS(pairij), pairneuind) && isequal(CindV1RS(pairij), pairCind) )
        error('check pairij/pairCind/pairneuind')
    end
    if ~all(pairCind(:,1)<pairCind(:,2) & pairneuind(:,1)<pairneuind(:,2))
        error('check pairCind and pairneuind -- columns are not sorted')
    end
end
% syncshufpsthavg.(whichblock) = zeros(length(psthtli), numel(trialtypes), Npairs, Nshuf);


save([pathpp 'syncpsthavg_probeC.mat'], 'CindV1RS', 'neuindV1RS', ...
    'pairij', 'pairCind', 'pairneuind', 'blocktrialtypes', 'syncpsthavg', '-v7.3')

%% plotting
whichblock = 'ICwcfg1_presentations';
trialtypes = blocktrialtypes.(whichblock);

ICsig.(whichblock).indin13 = ICsig.(whichblock).indin1 | ICsig.(whichblock).indin3;
ICsig.(whichblock).indin14 = ICsig.(whichblock).indin1 | ICsig.(whichblock).indin4;
ICsig.(whichblock).indin23 = ICsig.(whichblock).indin2 | ICsig.(whichblock).indin3;
ICsig.(whichblock).indin24 = ICsig.(whichblock).indin2 | ICsig.(whichblock).indin4;

% kerwinhalf = 25; kersigma = 12/sqrt(2);
kerwinhalf = 15; kersigma = 3;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));

for pltopt = [0 3 6]
    switch pltopt
        % between segment responders: fairly conssitent across sessions
        case 0 % slightly higher earlier (<200ms) -- fairly consistent across sessions
            group1 = {'indin13', 'indin14', 'indin23', 'indin24'};
            group2 = group1;
            xywh = [50 500 1200 400];
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
            xywh = [50 50 1200 400];

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
            xywh = [700 50 1200 400];
    end
    tt2p = [106 107 110 111];
    ttcol = [0 0.5 0; 0.5 0.25 0; 1 0.5 0; 0 1 0];
    figure('Position', xywh)
    annotation('textbox',[0 0.9 0.5 0.1], 'String', sprintf('%s %s %s %s\n%s %s %s %s', group1{:}, group2{:}), 'edgecolor', 'none', 'HorizontalAlignment', 'center')
    tempsynccell = cell(4,1);
    subplot(1,2,1)
    hold all
    for itt = 1:numel(tt2p)
        tempsr1 = neuoind(ICsig.(whichblock).(group1{itt}));
        tempsr2 = neuoind(ICsig.(whichblock).(group2{itt}));
        tempsrpair = sort( combvec(tempsr1', tempsr2')', 2);
        %rows2rm = tempsrpair(:,1)==tempsrpair(:,2);
        %tempsrpair(rows2rm,:)=[];
        pairoi = ismember(pairneuind, tempsrpair, 'rows');

        typi = ICtrialtypes(trialtypes+1)==tt2p(itt);
        tempsync = convn(syncpsthavg.(whichblock)(:,typi,pairoi),kergauss,'same') ;
        tempsynccell{itt} = tempsync;
        plot(syncpsthavg.syncpsthtli, squeeze(mean( tempsync, 3) ), 'Color', ttcol(itt,:), 'LineWidth', 1)
    end
    xlim([-400 800])
    subplot(1,2,2)
    hold all
    tempsync = cat(3, tempsynccell{1}, tempsynccell{4});
    plot(syncpsthavg.syncpsthtli, squeeze(mean( tempsync, 3) ), 'Color', [0 0.7 0], 'LineWidth', 1)
    tempsync = cat(3, tempsynccell{2}, tempsynccell{3});
    plot(syncpsthavg.syncpsthtli, squeeze(mean( tempsync, 3) ), 'Color', [1 0.5 0], 'LineWidth', 1)
    xlim([-400 800])
end

%{
typi = ICtrialtypes(trialtypes+1)==105;
paircol = [1 0 0; 0 1 0; 0 0 1];
figure
hold all
for ipair = 1:3
    switch ipair
        case 1
            tempsr1 = neuoind(ICsig.(whichblock).indin1);
            tempsr2 = neuoind(ICsig.(whichblock).indin3);
        case 2
            tempsr1 = neuoind(ICsig.(whichblock).indin1);
            tempsr2 = neuoind(ICsig.(whichblock).indin4);
        case 3
            tempsr1 = neuoind(ICsig.(whichblock).indin3);
            tempsr2 = neuoind(ICsig.(whichblock).indin4);
    end
    tempsrpair = sort( combvec(tempsr1', tempsr2')', 2);
    %rows2rm = tempsrpair(:,1)==tempsrpair(:,2);
    %tempsrpair(rows2rm,:)=[];
    pairoi = ismember(pairneuind, tempsrpair, 'rows');

tempsync = squeeze(mean( convn(syncpsthavg.(whichblock)(:,typi,pairoi),kergauss,'same'), 3) ) ;
plot(syncpsthavg.syncpsthtli, tempsync, 'Color', paircol(ipair,:), 'LineWidth', 1)
end
%}

end