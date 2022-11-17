datadir = 'D:\OpenScopeData\000248\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(contains(nwbsessions, 'sub-'));
Nsessions = 4; %numel(nwbsessions)-1;

fixgaze = true;

% visareas = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'};
% ICblocknames = {'ICkcfg0', 'ICkcfg1', 'ICwcfg0', 'ICwcfg1'};

whichvisarea = 'VISp';
whichICblock = 'ICwcfg1';

preproc = 'zscore'; % '' is z-score train trials, '_zscoreall', or '_meancenter'
svmdesc = 'trainICRCtestRE';
whichSVMkernel = 'Linear';

SVMsubICRCagg = struct();
tic
for ises = 1:Nsessions
    if fixgaze
    pathsvm = [datadir 'SVM_fixedgaze_' svmdesc '_subsets' filesep nwbsessions{ises} filesep];
    svmfn = strcat(pathsvm, 'SVM_fixedgaze_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_trainsubsets_', whichICblock, '.mat');
    else
    pathsvm = [datadir 'SVM_' svmdesc '_subsets' filesep nwbsessions{ises} filesep];
    svmfn = strcat(pathsvm, 'SVM_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_trainsubsets_', whichICblock, '.mat');
    end
    if exist(svmfn, 'file')
        load(svmfn, 'SVMtrainICRC')
        SVMsubICRCagg(ises).(whichvisarea) = SVMtrainICRC;
    else
        warning([svmfn ' does not exist'])
    end
end
toc

%%
subsetdescs = SVMtrainICRC.subsetdescs;
lmf = {'train', 'test', 'probe', 'REt', 'REx', 'X', 'blank'};

traintrialtypes = SVMtrainICRC.trialtypes';
probetrialtypes = unique(SVMtrainICRC.trialorder( SVMtrainICRC.probetrials ));

whichsubset = subsetdescs{1};
whichR = ['spkcnt_' whichsubset];
Nsplits=size(SVMtrainICRC.(whichR).traintrialinds,2);
Ntt = numel(traintrialtypes);
Nprobett = numel(probetrialtypes);

HR_SVMsubICRC = struct();
tic
for s = 1:numel(subsetdescs)
    whichsubset = subsetdescs{s};
    whichR = ['spkcnt_' whichsubset];

    HR_SVMsubICRC.(whichsubset) = struct();
    HR_SVMsubICRC.(whichsubset).train = NaN(Ntt,Ntt,Nsplits,Nsessions);
    HR_SVMsubICRC.(whichsubset).test = NaN(Ntt,Ntt,Nsplits,Nsessions);
    HR_SVMsubICRC.(whichsubset).probe = NaN(Nprobett,Ntt,Nsplits,Nsessions);
    HR_SVMsubICRC.(whichsubset).REt = NaN(2,Ntt,Nsplits,Nsessions);
    HR_SVMsubICRC.(whichsubset).T = NaN(2,Ntt,Nsplits,Nsessions);
    HR_SVMsubICRC.(whichsubset).REx = NaN(2,Ntt,Nsplits,Nsessions);
    HR_SVMsubICRC.(whichsubset).X = NaN(3,Ntt,Nsplits,Nsessions);
    HR_SVMsubICRC.(whichsubset).blank = NaN(1,Ntt,Nsplits,Nsessions);
    
    lmf = fieldnames(HR_SVMsubICRC.(whichsubset));
    for ises = 1:Nsessions
        if isempty( SVMsubICRCagg(ises).(whichvisarea) ) || isempty( fieldnames(SVMsubICRCagg(ises).(whichvisarea).(whichR)) )
            continue
        end
        for lm = 1:numel(lmf)
            for isplit = 1:Nsplits
                % tempyind: input trial indices, tempynn: output trial type, tempyu: input trial type labels
                switch lmf{lm}
                    case 'train'
                        tempyind = SVMsubICRCagg(ises).(whichvisarea).(whichR).traintrialinds(:,isplit);
                        tempysvm = SVMsubICRCagg(ises).(whichvisarea).(whichR).train.label(:,isplit);
                        tempyu = traintrialtypes;
                    case 'test'
                        tempyind = SVMsubICRCagg(ises).(whichvisarea).(whichR).testtrialinds(:,isplit);
                        tempysvm = SVMsubICRCagg(ises).(whichvisarea).(whichR).test.label(:,isplit);
                        tempyu = traintrialtypes;
                    case 'probe'
                        tempyind = SVMsubICRCagg(ises).(whichvisarea).probetrials;
                        tempysvm = SVMsubICRCagg(ises).(whichvisarea).(whichR).probe.label(:,isplit);
                        tempyu = probetrialtypes;
                    case 'REt'
                        tempyu = [1105, 1109];
                        probetrialinds = SVMsubICRCagg(ises).(whichvisarea).probetrials;
                        tempyind= find(ismember(SVMsubICRCagg(ises).(whichvisarea).trialorder , tempyu));
                        tempoi = ismember(probetrialinds, tempyind);
                        tempysvm = SVMsubICRCagg(ises).(whichvisarea).(whichR).probe.label(tempoi,isplit);
                    case 'T'
                        tempyu = [105, 109];
                        probetrialinds = SVMsubICRCagg(ises).(whichvisarea).probetrials;
                        tempyind= find(ismember(SVMsubICRCagg(ises).(whichvisarea).trialorder , tempyu));
                        tempoi = ismember(probetrialinds, tempyind);
                        tempysvm = SVMsubICRCagg(ises).(whichvisarea).(whichR).probe.label(tempoi,isplit);
                        if nnz(tempoi)==0
                            continue
                        end
                    case 'REx'
                        tempyu = [1201, 1299];
                        probetrialinds = SVMsubICRCagg(ises).(whichvisarea).probetrials;
                        tempyind= find(ismember(SVMsubICRCagg(ises).(whichvisarea).trialorder , tempyu));
                        tempoi = ismember(probetrialinds, tempyind);
                        tempysvm = SVMsubICRCagg(ises).(whichvisarea).(whichR).probe.label(tempoi,isplit);
                    case 'X'
                        tempyu = [101, 1201, 1299];
                        probetrialinds = SVMsubICRCagg(ises).(whichvisarea).probetrials;
                        tempyind= find(ismember(SVMsubICRCagg(ises).(whichvisarea).trialorder , tempyu));
                        tempoi = ismember(probetrialinds, tempyind);
                        tempysvm = SVMsubICRCagg(ises).(whichvisarea).(whichR).probe.label(tempoi,isplit);
                    case 'blank'
                        tempyind = SVMsubICRCagg(ises).(whichvisarea).alltrials;
                        tempysvm = SVMsubICRCagg(ises).(whichvisarea).(whichR).all.label(:,isplit);
                        tempyu = 0;
                end
                tempysvm = cellfun(@str2num, tempysvm);
                tempy = SVMsubICRCagg(ises).(whichvisarea).trialorder(tempyind);
                % rearrange ylabs order so that it is consistent across different sessions
                [sylabv,sylabi]=sort(SVMsubICRCagg(ises).(whichvisarea).(whichR).Ylabs(:,isplit));
                if ~isequal(sylabv, cellstr(num2str(traintrialtypes)))
                    error('check -- these should be the same')
                end
                tempylabs = sylabi;
                
                % calculate each element of matrix
                % rows: actual trial type, columns: categorized as
                for iy = 1:numel(tempyu)
                    trialsoy = tempy==tempyu(iy);
                    trialsnoty = tempy~=tempyu(iy);
                    %trialsnoty = true(size(tempy));
                    for jy = 1:numel(traintrialtypes)
                        HR_SVMsubICRC.(whichsubset).(lmf{lm})(iy,jy,isplit,ises) = nnz(tempysvm(trialsoy)==traintrialtypes(jy))/nnz(trialsoy);
                    end
                end
                                
            end
        end
    end
    
end
toc

%% adjusted performance metrics
if ~exist('dprime_SVMsubICRC', 'var')
    dprime_SVMsubICRC = struct();
    pcb_SVMsubICRC = struct();
    lognormHR_SVMsubICRC = struct();
end

for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    for b = 1:numel(ICblocknames)
        whichICblock = ICblocknames{b};
        for lm = 1:numel(lmf)
            if ~isequaln(HR_SVMsubICRC.(whichvisarea).probe(probetrialtypes==0,:,:,:), ...
                    HR_SVMsubICRC.(whichvisarea).blank)
                error('check base HR')
            end
            % add a tiny little offset
            tempHR = tempHR+2^-10;
            tempbase = HR_SVMsubICRC.(whichvisarea).blank+2^-10;
            dprime_SVMsubICRC.(whichvisarea).(lmf{lm}) = norminv(tempHR) - norminv(tempbase);
            % pcb: percent change from baseline, baseline being probability of any *test* trial being classified as a certain output class
            pcb_SVMsubICRC.(whichvisarea).(lmf{lm}) = 100*(tempHR./tempbase -1);
            % normHR: percent change from baseline, baseline being probability of any *test* trial being classified as a certain output class
            lognormHR_SVMsubICRC.(whichvisarea).(lmf{lm}) = log2(tempHR./tempbase);
        end
    end
end


%% subsets
fs = 10;
xtl = {'IC1', 'RC1', 'RC2', 'IC2'};
ytl = {'REt1', 'REt2'};
% whichsubset = 'allbutsigkwBI';
whichsubset = 'inducerencoder';

figure;
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', [preproc ' SVM V1 test accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for ises = 1:4
    subplot(2,2,ises)
    tempHR = squeeze(nanmean(HR_SVMsubICRC.(whichsubset).test(:,:,:,ises), 3 ));
    imagesc(tempHR)
    caxis([0 1]); colorbar
    set(gca, 'fontsize', fs, 'XTick', 1:4, 'XTickLabel', xtl, 'YTick', 1:4, 'YTickLabel', xtl)
    title(sprintf('Session%d %s %.4f', ises, whichICblock, mean(diag(tempHR))) )
end
colormap jet

figure
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', [preproc ' SVM V1 probe accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for ises = 1:4
subplot(2,2,ises)
tempHR = squeeze(nanmean(HR_SVMsubICRC.(whichsubset).REt(:,:,:,ises), 3 ));
imagesc(tempHR)
caxis([0 1]); colorbar
set(gca, 'fontsize', fs, 'XTick', 1:4, 'XTickLabel', xtl, 'YTick', 1:4, 'YTickLabel', ytl)
infscore = squeeze( (( tempHR(1,1,:)-tempHR(1,2,:) )+( tempHR(2,4,:)-tempHR(2,3,:) ))/2 );
title(sprintf('Session%d %s IC-RC %.4f', ises, whichICblock, mean(infscore)) )
colorbar
end
colormap jet
