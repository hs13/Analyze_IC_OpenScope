datadir = 'D:\OpenScopeData\000248v221123\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name}; 
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

for ises = 1:Nsessions
clearvars -except ises nwbsessions Nsessions datadir
sesclk = tic;
fprintf('\nSession %d %s\n', ises, nwbsessions{ises})

pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load([pathpp 'info_units.mat'])
visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
    'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};

probes = {'A', 'B', 'C', 'D', 'E', 'F'};
psthtli = (-500:1000)';
neucnt = 0;
for iprobe = 1:numel(probes)
neuoind = find(floor(unit_peakch/1000)==iprobe-1);
neucnt = neucnt + numel(neuoind);
fprintf('Probe %s: %d\n', probes{iprobe}, numel(neuoind) )
tic
load(sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe}))
% load(sprintf('%svisresponses_probe%s.mat', pathpp, probes{iprobe})); %, 'meanFRvec', 'sponFRvec')

ICtrialtypes = [0 101 105 106 107 109 110 111 506 511 1105 1109 1201 1299 ...
    1301 1302 1303 1304 1305 1306 1307 1308];
ICblocks = 1:4;
BKshuf = struct();
for b = ICblocks
    disp(visblocks{b})
    tloi = psthtli>0 & psthtli<=400;
    tempR = squeeze(1000*mean(psth.(visblocks{b})(tloi,:,:), 1))';
    temptrialorder = ICtrialtypes( vis.(visblocks{b}).trialorder + 1);

    % shuffle trialorder only between Kanizsa trials [106 107 110 111]
    % note, blank trials are not randomized, therefor sigkwBK should not
    % change much between shuffles
    Nshuf = 1000;
    shuftt = [106 107 110 111];
    trialsoi = ismember(temptrialorder, shuftt);
    trials2shuf = temptrialorder(trialsoi);
    BKshuf.(visblocks{b}) = struct();
for ishuf = 1:Nshuf
    shuftrialorder = temptrialorder;
    shuftrialorder(trialsoi) = trials2shuf(randperm(numel(trials2shuf)));
    tempBKshuf = analyzeBKshuf(tempR, shuftrialorder);
    if ishuf==1
        BKshuf.(visblocks{b}) = tempBKshuf;
    else
        BKshuf.(visblocks{b}) = cat(1, BKshuf.(visblocks{b}), tempBKshuf);
    end
end

end


save(sprintf('%sBKshuf_probe%s.mat', pathpp, probes{iprobe}), 'BKshuf')
toc
end
toc(sesclk)
end

