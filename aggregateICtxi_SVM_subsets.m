
% datadir = 'D:\OpenScopeData\000248\';
% nwbdir = dir(datadir);
% nwbsessions = {nwbdir.name};
% nwbsessions = nwbsessions(contains(nwbsessions, 'sub-'));
% Nsessions = numel(nwbsessions)-1;

datadir = '/Users/hyeyoung/Documents/OpenScopeData/';
nwbsessions = {'sub-1171903426','sub-1172969386','sub-1174569632','sub-1175512776', ...
    'sub-1177693335','sub-1181585601','sub-1182593224','sub-1183369796','sub-1186544719'};
Nsessions = numel(nwbsessions);

fixgaze = false;
justctx = false;
subsetgroups = {'sigkw'}; % 'sigkw', 'encoder', 'CGIG'

whichvisarea = 'C';
whichICblock = 'ICkcfg1';

svmdesc = 'trainICRCtestRE';
preproc = 'zscore'; % '' is z-score train trials, '_zscoreall', or '_meancenter'
whichSVMkernel = 'Linear';

% visareas = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'};
% ICblocknames = {'ICkcfg0', 'ICkcfg1', 'ICwcfg0', 'ICwcfg1'};
%visareas = {'A', 'B', 'C', 'D', 'E', 'F'};
%probelabels = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};

SVMsubICRCagg = struct();
tic
for g = 1:numel(subsetgroups)
    whichsubsetgroup = subsetgroups{g};
    if justctx
        if fixgaze
            pathsv = [datadir 'SVM_fixedgaze_' svmdesc '_subsets_' whichsubsetgroup filesep];
        else
            pathsv = [datadir 'SVM_' svmdesc '_subsets_' whichsubsetgroup filesep];
        end
    else
        if fixgaze
            pathsv = [datadir 'SVM_fixedgaze_' svmdesc '_subsets_' whichsubsetgroup '_allunits' filesep];
        else
            pathsv = [datadir 'SVM_' svmdesc '_subsets_' whichsubsetgroup '_allunits' filesep];
        end
    end

for ises = 1:Nsessions
    if fixgaze
    pathsvm = [pathsv nwbsessions{ises} filesep];
    svmfn = strcat(pathsvm, 'SVM_fixedgaze_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_trainsubsets_', whichsubsetgroup, '_', whichICblock, '.mat');
    else
    pathsvm = [pathsv nwbsessions{ises} filesep];
    svmfn = strcat(pathsvm, 'SVM_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_trainsubsets_', whichsubsetgroup, '_', whichICblock, '.mat');
    end
    if exist(svmfn, 'file')
        load(svmfn, 'SVMtrainICRC')
        SVMsubICRCagg(ises).(whichsubsetgroup) = SVMtrainICRC;
    else
        warning([svmfn ' does not exist'])
    end
end
end
toc

%%
Nallneurons = zeros(Nsessions,1);
for ises = 1:Nsessions
    Nallneurons(ises) = SVMsubICRCagg(ises).(whichsubsetgroup).Nneurons;
end

Nsubsetneurons = struct();
for g = 1:numel(subsetgroups)
    whichsubsetgroup = subsetgroups{g};
subsetdescs = SVMsubICRCagg(1).(whichsubsetgroup).subsetdescs;
for s = 1:numel(subsetdescs)
    whichsubset = subsetdescs{s};
Nsubsetneurons.(whichsubsetgroup).(whichsubset) = zeros(Nsessions,1);
for ises = 1:Nsessions
    Nsubsetneurons.(whichsubsetgroup).(whichsubset)(ises) = nnz(SVMsubICRCagg(ises).(whichsubsetgroup).subsets2train{s});
end
end
end

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

    HR_SVMsubICRC.(whichsubsetgroup).(whichsubset) = struct();
    HR_SVMsubICRC.(whichsubsetgroup).(whichsubset).train = NaN(Ntt,Ntt,Nsplits,Nsessions);
    HR_SVMsubICRC.(whichsubsetgroup).(whichsubset).test = NaN(Ntt,Ntt,Nsplits,Nsessions);
    HR_SVMsubICRC.(whichsubsetgroup).(whichsubset).probe = NaN(Nprobett,Ntt,Nsplits,Nsessions);
    HR_SVMsubICRC.(whichsubsetgroup).(whichsubset).REt = NaN(2,Ntt,Nsplits,Nsessions);
    HR_SVMsubICRC.(whichsubsetgroup).(whichsubset).T = NaN(2,Ntt,Nsplits,Nsessions);
    HR_SVMsubICRC.(whichsubsetgroup).(whichsubset).REx = NaN(2,Ntt,Nsplits,Nsessions);
    HR_SVMsubICRC.(whichsubsetgroup).(whichsubset).X = NaN(3,Ntt,Nsplits,Nsessions);
    HR_SVMsubICRC.(whichsubsetgroup).(whichsubset).blank = NaN(1,Ntt,Nsplits,Nsessions);
    
    lmf = fieldnames(HR_SVMsubICRC.(whichsubsetgroup).(whichsubset));
    for ises = 1:Nsessions
        if isempty( SVMsubICRCagg(ises).(whichsubsetgroup) ) || isempty( fieldnames(SVMsubICRCagg(ises).(whichsubsetgroup).(whichR)) )
            continue
        end
        for lm = 1:numel(lmf)
            for isplit = 1:Nsplits
                % tempyind: input trial indices, tempynn: output trial type, tempyu: input trial type labels
                switch lmf{lm}
                    case 'train'
                        tempyind = SVMsubICRCagg(ises).(whichsubsetgroup).(whichR).traintrialinds(:,isplit);
                        tempysvm = SVMsubICRCagg(ises).(whichsubsetgroup).(whichR).train.label(:,isplit);
                        tempyu = traintrialtypes;
                    case 'test'
                        tempyind = SVMsubICRCagg(ises).(whichsubsetgroup).(whichR).testtrialinds(:,isplit);
                        tempysvm = SVMsubICRCagg(ises).(whichsubsetgroup).(whichR).test.label(:,isplit);
                        tempyu = traintrialtypes;
                    case 'probe'
                        tempyind = SVMsubICRCagg(ises).(whichsubsetgroup).probetrials;
                        tempysvm = SVMsubICRCagg(ises).(whichsubsetgroup).(whichR).probe.label(:,isplit);
                        tempyu = probetrialtypes;
                    case 'REt'
                        tempyu = [1105, 1109];
                        probetrialinds = SVMsubICRCagg(ises).(whichsubsetgroup).probetrials;
                        tempyind= find(ismember(SVMsubICRCagg(ises).(whichsubsetgroup).trialorder , tempyu));
                        tempoi = ismember(probetrialinds, tempyind);
                        tempysvm = SVMsubICRCagg(ises).(whichsubsetgroup).(whichR).probe.label(tempoi,isplit);
                    case 'T'
                        tempyu = [105, 109];
                        probetrialinds = SVMsubICRCagg(ises).(whichsubsetgroup).probetrials;
                        tempyind= find(ismember(SVMsubICRCagg(ises).(whichsubsetgroup).trialorder , tempyu));
                        tempoi = ismember(probetrialinds, tempyind);
                        tempysvm = SVMsubICRCagg(ises).(whichsubsetgroup).(whichR).probe.label(tempoi,isplit);
                        if nnz(tempoi)==0
                            continue
                        end
                    case 'REx'
                        tempyu = [1201, 1299];
                        probetrialinds = SVMsubICRCagg(ises).(whichsubsetgroup).probetrials;
                        tempyind= find(ismember(SVMsubICRCagg(ises).(whichsubsetgroup).trialorder , tempyu));
                        tempoi = ismember(probetrialinds, tempyind);
                        tempysvm = SVMsubICRCagg(ises).(whichsubsetgroup).(whichR).probe.label(tempoi,isplit);
                    case 'X'
                        tempyu = [101, 1201, 1299];
                        probetrialinds = SVMsubICRCagg(ises).(whichsubsetgroup).probetrials;
                        tempyind= find(ismember(SVMsubICRCagg(ises).(whichsubsetgroup).trialorder , tempyu));
                        tempoi = ismember(probetrialinds, tempyind);
                        tempysvm = SVMsubICRCagg(ises).(whichsubsetgroup).(whichR).probe.label(tempoi,isplit);
                    case 'blank'
                        tempyind = SVMsubICRCagg(ises).(whichsubsetgroup).alltrials;
                        tempysvm = SVMsubICRCagg(ises).(whichsubsetgroup).(whichR).all.label(:,isplit);
                        tempyu = 0;
                end
                tempysvm = cellfun(@str2num, tempysvm);
                tempy = SVMsubICRCagg(ises).(whichsubsetgroup).trialorder(tempyind);
                % rearrange ylabs order so that it is consistent across different sessions
                [sylabv,sylabi]=sort(SVMsubICRCagg(ises).(whichsubsetgroup).(whichR).Ylabs(:,isplit));
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
                        HR_SVMsubICRC.(whichsubsetgroup).(whichsubset).(lmf{lm})(iy,jy,isplit,ises) = ...
                            nnz(tempysvm(trialsoy)==traintrialtypes(jy))/nnz(trialsoy);
                    end
                end
                                
            end
        end
    end
    
end
toc

%% subsets
pltopt = 1;
switch pltopt
    case 1
        whichsubsetgroup = 'sigkw';
        subsets2plot = {'sigkwBI', 'sigkwBK', 'allbutsigkwBI', 'allbutsigkwBK'};
    case 2
        whichsubsetgroup = 'encoder';
        subsets2plot = {'ICRCencoder', 'inducerencoder', 'inducerresponsive', 'indin', 'allbutICRCencoder'};
    case 3
        whichsubsetgroup = 'CGIG';
        subsets2plot = {'sigkwCRF', 'ctrCRF9', 'allbutsigkwCRF', 'allbutctrCRF9'};
    case 4
        whichsubsetgroup = 'CGIG';
        subsets2plot = {'exclctrCRF9', 'exclCRF9', 'allbutctr9'};
    case 5
        whichsubsetgroup = 'CGIG';
        subsets2plot = {'ctrIRF9', 'ctrCRF9nonctrIRF9', 'faithCRF9', 'allbutctrIRF9'};
end
fs = 14;
xtl = {'IC1', 'RC1', 'RC2', 'IC2'};
ytl = {'REt1', 'REt2'};

figure('Position', [0 100 1500 400])
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', [preproc ' SVM Area ' whichvisarea ' ' whichICblock ' accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for s = 1:numel(subsets2plot)
    whichsubset = subsets2plot{s};
    %     subplot(2,ceil(numel(subsets2plot)/2),s)
    subplot(2,5,s)
    tempHR = squeeze(nanmean(HR_SVMsubICRC.(whichsubsetgroup).(whichsubset).test, 3 ));
    imagesc(squeeze(mean(tempHR,3)))
    caxis([0 0.5]); colorbar
    set(gca, 'fontsize', fs, 'XTick', 1:4, 'XTickLabel', xtl, 'YTick', 1:4, 'YTickLabel', xtl)
    testscore = squeeze( ( tempHR(1,1,:)+tempHR(2,2,:)+tempHR(3,3,:)+tempHR(4,4,:) )/4 );
    p= signrank(testscore-0.25);
    title(sprintf('%s N=%.0f test %.2f\np=%.4f', whichsubset, mean(Nsubsetneurons.(whichsubsetgroup).(whichsubset)), mean(testscore), p) )
    colorbar
    % end
    % colormap jet
    %
    % figure('Position', [800 100 800 400])
    % for s = 1:numel(subsets2plot)
    %     whichsubset = subsets2plot{s};
    %     subplot(2,ceil(numel(subsets2plot)/2),s)
    subplot(2,5,5+s)
    tempHR = squeeze(nanmean(HR_SVMsubICRC.(whichsubsetgroup).(whichsubset).REt, 3 ));
    imagesc(squeeze(mean(tempHR,3)))
    caxis([0 0.5]); colorbar
    set(gca, 'fontsize', fs, 'XTick', 1:4, 'XTickLabel', xtl, 'YTick', 1:2, 'YTickLabel', ytl)
    infscore = squeeze( (( tempHR(1,1,:)-tempHR(1,2,:) )+( tempHR(2,4,:)-tempHR(2,3,:) ))/2 );
    p= signrank(infscore);
    infscore1 = squeeze( tempHR(1,1,:)-tempHR(1,2,:) );
    infscore2 = squeeze( tempHR(1,1,:)-tempHR(1,2,:) );
    ppool= signrank([infscore1; infscore2]);
    title(sprintf('%s N=%.0f IC-RC %.2f\np=%.4f p_p_o_o_l=%.4f', whichsubset, mean(Nsubsetneurons.(whichsubsetgroup).(whichsubset)), mean(infscore), p, ppool) )
    colorbar
end
colormap jet

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
    tempHR = squeeze(nanmean(HR_SVMsubICRC.(whichsubsetgroup).(whichsubset).test(:,:,:,ises), 3 ));
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
tempHR = squeeze(nanmean(HR_SVMsubICRC.(whichsubsetgroup).(whichsubset).REt(:,:,:,ises), 3 ));
imagesc(tempHR)
caxis([0 1]); colorbar
set(gca, 'fontsize', fs, 'XTick', 1:4, 'XTickLabel', xtl, 'YTick', 1:4, 'YTickLabel', ytl)
infscore = squeeze( (( tempHR(1,1,:)-tempHR(1,2,:) )+( tempHR(2,4,:)-tempHR(2,3,:) ))/2 );
title(sprintf('Session%d %s IC-RC %.4f', ises, whichICblock, mean(infscore)) )
colorbar
end
colormap jet
