datadir = 'S:/OpenScopeData/000248/';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));

% decoder weight
% weight_vector=c1.Beta;
% bais_vector=c1.Bias;
%%
svmdesc = 'trainRExtestICRC';
preproc = 'zscore'; % '' is z-score train trials, '_zscoreall', or '_meancenter'
whichSVMkernel = 'Linear';

whichICblock = 'ICwcfg1';
probes = {'A', 'B', 'C', 'D', 'E', 'F'};

for ises = 1:numel(nwbsessions)
    clearvars -except datadir nwbsessions ises svmdesc preproc whichSVMkernel whichICblock probes 
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    pathsv = [datadir 'postprocessed' filesep 'SVM' filesep 'SVM_' svmdesc filesep];
    pathsvm = [pathsv nwbsessions{ises} filesep];
    load([pathpp 'spiketimes.mat'], 'neuctx')
    for iprobe = 1:numel(probes)
        whichvisarea = [probes{iprobe} 'ctx'];
        fprintf('%d/%d %s %s\n', ises, numel(nwbsessions), svmdesc, whichvisarea)
        tic

        svmfn = strcat(pathsvm, 'SVM_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_silencesubsets_', whichICblock, '.mat');
        svmmdlfn = strcat(pathsvm, 'SVMmodels_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_silencesubsets_', whichICblock, '.mat');
        if ~exist(svmmdlfn, 'file')
            disp([svmmdlfn ' does not exist'])
            continue
        end
        load(svmfn)
        load(svmmdlfn)
        load([pathpp 'postprocessed_probe' probes{iprobe} '.mat'])


        switch svmdesc
            case 'trainICRCtestRE'
                SVMout = SVMtrainICRC;
                SVM_models = SVMtrainICRC_models;
            case 'trainRExtestICRC'
                SVMout = SVMtrainREx;
                SVM_models = SVMtrainREx_models;
        end

        ICtrialtypes = [0 101 105 106 107 109 110 111 506 511 1105 1109 1201 1299 ...
            1301 1302 1303 1304 1305 1306 1307 1308];
        traintrialtypes = SVMout.trialtypes;
        c = nchoosek(1:length(traintrialtypes),2);
        traintrialcombos = traintrialtypes(c);

        neuctxinprobe = neuctx(neuoind);
        Nclasses = length(traintrialtypes);
        Nbinarylearners = size(traintrialcombos,1);
        Nsplits = numel(SVM_models.spkcnt);

        SVMcodingname = cell(Nsplits,1);
        Nlearners = zeros(Nsplits,1);
        betalearners = cell(Nsplits,1);
        for isplit = 1:Nsplits
            Ylabs = str2num(cat(1,SVMout.spkcnt.Ylabs{:,isplit}))';
            splitcombo = sort(Ylabs(c),2);
            Nlearners(isplit) = numel(SVM_models.spkcnt{isplit}.BinaryLearners);
            SVMcodingname{isplit} = SVM_models.spkcnt{isplit}.CodingName;
            betalearners{isplit} = NaN(nnz(neuctxinprobe), Nlearners(isplit));
            for ibl = 1:Nlearners(isplit)
                switch SVM_models.spkcnt{isplit}.CodingName
                    case 'onevsone'
                        indbl = find(ismember(splitcombo, traintrialcombos(ibl,:), 'rows' ));
                        if ~isequal(sort(Ylabs(c(indbl,:))), traintrialcombos(ibl,:))
                            error('sanity check did not pass, check indbl')
                        end
                        if isequal(Ylabs(c(indbl,:)), traintrialcombos(ibl,:)) % same order
                            betagain = 1;
                        else % flipped order
                            betagain = -1;
                        end
                    case 'onevsall'
                        if length(traintrialtypes)==2
                            indbl = 1;
                            if isequal(Ylabs, traintrialtypes) % same order
                                betagain = 1;
                            else % flipped order
                                betagain = -1;
                            end
                        else
                            indbl = find(Ylabs==traintrialtypes(ibl));
                            betagain = 1;
                        end
                    otherwise
                        error('unexpected CodingName -- need to make another exception case')
                end
                betalearners{isplit}(:,ibl) = betagain * SVM_models.spkcnt{isplit}.BinaryLearners{indbl}.Beta;
            end
        end
        save([pathsvm, 'SVMbeta_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'], ...
            'ICtrialtypes', 'traintrialtypes', 'traintrialcombos', 'neuctxinprobe', 'SVMcodingname', 'Nlearners', 'betalearners', '-v7.3')

        % SVM PSTH
        trialorder = ICtrialtypes(vis.([whichICblock '_presentations']).trialorder + 1);
        temppsth = psth.([whichICblock '_presentations'])(:,:,neuctxinprobe==1);
        psth100msTinds = (0:99)'+(find(psthtli==-200):10:find(psthtli==600)-100);
        Nbins = size(psth100msTinds,2);
        psth100ms = NaN(Nbins, length(trialorder), nnz(neuctxinprobe==1) );
        for ibin = 1:Nbins
            temptli = psth100msTinds(:,ibin);
            psth100ms(ibin, :, :) = 1000*mean(temppsth(temptli, :, :), 1);
        end

        SVMpsth100ms = struct();
        SVMpsth100ms.trialorder = trialorder;
        SVMpsth100ms.traintrialinds = SVMout.spkcnt.traintrialinds;
        SVMpsth100ms.testtrialinds = SVMout.spkcnt.testtrialinds;
        SVMpsth100ms.traintrialtypes = SVMout.trialtypes;
        SVMpsth100ms.Ylabs = SVMout.spkcnt.Ylabs;
        SVMpsth100ms.psth100msTinds = psth100msTinds;
        SVMpsth100ms.psth100ms = psth100ms;
        SVMpsth100ms.label = NaN(Nbins, length(trialorder), Nsplits);
        SVMpsth100ms.score = NaN(Nbins, length(trialorder), Nclasses, Nsplits);
        for isplit = 1:Nsplits
            Ylabs = str2num(cat(1,SVMout.spkcnt.Ylabs{:,isplit}))';
            [~,reverseYlaborder]=sort(Ylabs);

            for ibin = 1:Nbins
                Xtemp = squeeze( psth100ms(ibin,:,:) );
                tempSVMmodel = SVM_models.spkcnt{isplit};
                [tempilabel,tempscore] = predict(tempSVMmodel,Xtemp); % Xtemp is Ntrials X Nneurons
                SVMpsth100ms.label(ibin,:,isplit) = Ylabs(tempilabel);
                SVMpsth100ms.score(ibin,:,:,isplit) = tempscore(:, reverseYlaborder);
            end
        end
        if exist([pathsvm, 'SVMpsth_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'], 'file')
            delete([pathsvm, 'SVMpsth_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'])
        end
        save([pathsvm, 'SVMpsth100ms_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'], ...
            'SVMpsth100ms', '-v7.3')
        toc
    end
end

%%

betabl = cat(3,betalearners{:});

for ibl = 1:Nbinarylearners
    corrbeta = corr(squeeze(betabl(:,ibl,Nlearners==Nbinarylearners)));
disp([ibl mean(corrbeta(triu(true(size(corrbeta)), 1))) mean(abs(corrbeta(triu(true(size(corrbeta)), 1))))])
end

figure; histogram( corrbeta(triu(true(size(corrbeta)))), -1:0.05:1 )
figure; imagesc(corrbeta); caxis([-1 1]); colorbar
figure; plot(betabl(:,1,1), betabl(:,1,6), '.')

Z = linkage(corrbeta,'complete','correlation');
figure
[H,T,outperm] = dendrogram(Z);

figure; imagesc(corrbeta(outperm,outperm)); caxis([-1 1]); colorbar

Z2 = linkage(squeeze(betabl(:,ibl,Nlearners=Nbinarylearners))','complete','correlation');
figure
[H2,T2,outperm2] = dendrogram(Z2);

figure; imagesc(corrbeta(outperm2,outperm2)); caxis([-1 1]); colorbar

% note, Z and Z2 are very similar

figure
hold all
for itt = 1:numel(traintrialtypes)
    typi = traintrialtypes(itt);
    trialsoi = trialorder==typi;
    temppsth = mean(SVMpsth100ms.label(:,trialsoi,:)==typi,[2,3]);
    plot(psthtli(psth100msTinds(51,:)), temppsth)
end
