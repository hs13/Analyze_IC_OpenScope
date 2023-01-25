% aggregateICtxi_SVM_trainICRCtestRE_loadall.m requires too much RAM
% instead of aggregating SVMtrainICRCagg, compute HR_SVMtrainICRC while loading SVMtrainICRC

%%
datadir = 'D:\OpenScopeData\000248\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name}; 
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

svmdesc = 'trainICRCtestRE';
preproc = 'zscore'; % '' is z-score train trials, '_zscoreall', or '_meancenter'
whichSVMkernel = 'Linear';

justctx = true;
if justctx
    visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
    pathsv = [datadir 'postprocessed' filesep 'SVM' filesep 'SVM_' svmdesc filesep];
else
    visareas = {'A', 'B', 'C', 'D', 'E', 'F'};
    probelabels = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
    pathsv = [datadir 'postprocessed' filesep 'SVM' filesep 'SVM_' svmdesc '_allunits' filesep];
end
ICblocknames = {'ICkcfg0', 'ICkcfg1', 'ICwcfg0', 'ICwcfg1'};

% visareas = {'VISp'};
% ICblocknames = {'ICwcfg1'};

% % initialize
whichvisarea = visareas{1};
whichICblock = ICblocknames{1};
pathsvm = [pathsv nwbsessions{1} filesep];
svmfn = strcat(pathsvm, 'SVM_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_silencesubsets_', whichICblock, '.mat');
load(svmfn, 'SVMtrainICRC')

icf=1; lmf = {'train', 'test', 'probe', 'REt', 'T', 'REx', 'X', 'blank'};
traintrialtypes = SVMtrainICRC.trialtypes';
probetrialtypes = unique(SVMtrainICRC(icf).trialorder( SVMtrainICRC.probetrials ));
Ntt = numel(traintrialtypes);
Nprobett = numel(probetrialtypes);
Nsplits=size(SVMtrainICRC.spkcnt.traintrialinds,2);
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    for b = 1:numel(ICblocknames)
        whichICblock = ICblocknames{b};
        HR_SVMtrainICRC.(whichICblock).(whichvisarea) = struct();
        HR_SVMtrainICRC.(whichICblock).(whichvisarea).train = NaN(Ntt,Ntt,Nsplits,Nsessions);
        HR_SVMtrainICRC.(whichICblock).(whichvisarea).test = NaN(Ntt,Ntt,Nsplits,Nsessions);
        HR_SVMtrainICRC.(whichICblock).(whichvisarea).probe = NaN(Nprobett,Ntt,Nsplits,Nsessions);
        HR_SVMtrainICRC.(whichICblock).(whichvisarea).REt = NaN(2,Ntt,Nsplits,Nsessions);
        HR_SVMtrainICRC.(whichICblock).(whichvisarea).T = NaN(2,Ntt,Nsplits,Nsessions);
        HR_SVMtrainICRC.(whichICblock).(whichvisarea).REx = NaN(2,Ntt,Nsplits,Nsessions);
        HR_SVMtrainICRC.(whichICblock).(whichvisarea).X = NaN(3,Ntt,Nsplits,Nsessions);
        HR_SVMtrainICRC.(whichICblock).(whichvisarea).blank = NaN(1,Ntt,Nsplits,Nsessions);
    end
end

% SVMtrainICRCagg = struct();
tic
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    disp(whichvisarea)
    for b = 1:numel(ICblocknames)
        whichICblock = ICblocknames{b};
        disp(whichICblock)
        for ises = 1:Nsessions
            pathsvm = [pathsv nwbsessions{ises} filesep];
            svmfn = strcat(pathsvm, 'SVM_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_silencesubsets_', whichICblock, '.mat');
            if ~exist(svmfn, 'file')
                disp([svmfn ' does not exist'])
                continue
            end
            load(svmfn, 'SVMtrainICRC')
            %SVMtrainICRCagg(ises).(whichICblock).(whichvisarea) = SVMtrainICRC;
            for lm = 1:numel(lmf)
                for isplit = 1:Nsplits
                    % tempyind: input trial indices, tempynn: output trial type, tempyu: input trial type labels
                    switch lmf{lm}
                        case 'train'
                            tempyind = SVMtrainICRC.spkcnt.traintrialinds(:,isplit);
                            tempysvm = SVMtrainICRC.spkcnt.train.label(:,isplit);
                            tempyu = traintrialtypes;
                        case 'test'
                            tempyind = SVMtrainICRC.spkcnt.testtrialinds(:,isplit);
                            tempysvm = SVMtrainICRC.spkcnt.test.label(:,isplit);
                            tempyu = traintrialtypes;
                        case 'probe'
                            tempyind = SVMtrainICRC.probetrials;
                            tempysvm = SVMtrainICRC.spkcnt.probe.label(:,isplit);
                            tempyu = probetrialtypes;
                        case 'REt'
                            tempyu = [1105, 1109];
                            probetrialinds = SVMtrainICRC.probetrials;
                            tempyind= find(ismember(SVMtrainICRC.trialorder , tempyu));
                            tempoi = ismember(probetrialinds, tempyind);
                            tempysvm = SVMtrainICRC.spkcnt.probe.label(tempoi,isplit);
                        case 'T'
                            tempyu = [105, 109];
                            probetrialinds = SVMtrainICRC.probetrials;
                            tempyind= find(ismember(SVMtrainICRC.trialorder , tempyu));
                            tempoi = ismember(probetrialinds, tempyind);
                            tempysvm = SVMtrainICRC.spkcnt.probe.label(tempoi,isplit);
                            if nnz(tempoi)==0
                                continue
                            end
                        case 'REx'
                            tempyu = [1201, 1299];
                            probetrialinds = SVMtrainICRC.probetrials;
                            tempyind= find(ismember(SVMtrainICRC.trialorder , tempyu));
                            tempoi = ismember(probetrialinds, tempyind);
                            tempysvm = SVMtrainICRC.spkcnt.probe.label(tempoi,isplit);
                        case 'X'
                            tempyu = [101, 1201, 1299];
                            probetrialinds = SVMtrainICRC.probetrials;
                            tempyind= find(ismember(SVMtrainICRC.trialorder , tempyu));
                            tempoi = ismember(probetrialinds, tempyind);
                            tempysvm = SVMtrainICRC.spkcnt.probe.label(tempoi,isplit);
                        case 'blank'
                            tempyind = SVMtrainICRC.alltrials;
                            tempysvm = SVMtrainICRC.spkcnt.all.label(:,isplit);
                            tempyu = 0;
                    end
                    tempysvm = cellfun(@str2num, tempysvm);
                    tempy = SVMtrainICRC.trialorder(tempyind);
                    % rearrange ylabs order so that it is consistent across different sessions
                    [sylabv,sylabi]=sort(SVMtrainICRC.spkcnt.Ylabs(:,isplit));
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
                            HR_SVMtrainICRC.(whichICblock).(whichvisarea).(lmf{lm})(iy,jy,isplit,ises) = nnz(tempysvm(trialsoy)==traintrialtypes(jy))/nnz(trialsoy);
                        end
                    end
                    
                end % isplit
            end % lm
        end % ises
        toc
    end % b
end % a
toc

save([pathsv 'HR_SVMtrainICRC_' preproc '_agg.mat'], 'HR_SVMtrainICRC', '-v7.3')
save(['G:\My Drive\DATA\ICexpts_submission22\openscope_HR_SVMtrainICRC_' preproc '_agg.mat'], 'HR_SVMtrainICRC', '-v7.3')

%% adjusted performance metrics
if ~exist('dprime_SVMtrainICRC', 'var')
    dprime_SVMtrainICRC = struct();
    pcb_SVMtrainICRC = struct();
    lognormHR_SVMtrainICRC = struct();
end

for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    for b = 1:numel(ICblocknames)
        whichICblock = ICblocknames{b};
        for lm = 1:numel(lmf)
            if ~isequaln(HR_SVMtrainICRC.(whichICblock).(whichvisarea).probe(probetrialtypes==0,:,:,:), ...
                    HR_SVMtrainICRC.(whichICblock).(whichvisarea).blank)
                error('check base HR')
            end
            % add a tiny little offset
            tempHR = tempHR+2^-10;
            tempbase = HR_SVMtrainICRC.(whichICblock).(whichvisarea).blank+2^-10;
            dprime_SVMtrainICRC.(whichICblock).(whichvisarea).(lmf{lm}) = norminv(tempHR) - norminv(tempbase);
            % pcb: percent change from baseline, baseline being probability of any *test* trial being classified as a certain output class
            pcb_SVMtrainICRC.(whichICblock).(whichvisarea).(lmf{lm}) = 100*(tempHR./tempbase -1);
            % normHR: percent change from baseline, baseline being probability of any *test* trial being classified as a certain output class
            lognormHR_SVMtrainICRC.(whichICblock).(whichvisarea).(lmf{lm}) = log2(tempHR./tempbase);
        end
    end
end


%% compare blocks V1
fs = 14;
whichvisarea = 'V1';
xtl = {'IC1', 'RC1', 'RC2', 'IC2'};
ytl = {'REt1', 'REt2'};
Ntt = size(HR_SVMtrainICRC.(whichICblock).(whichvisarea).test,2);

figure;
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', [preproc ' SVM ' whichvisarea ' test accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for b = 1:numel(ICblocknames)
    whichICblock = ICblocknames{b};
    tempHR = squeeze(nanmean(HR_SVMtrainICRC.(whichICblock).(whichvisarea).test, 3 ));
    hrvec = reshape(tempHR, Ntt^2, size(tempHR,3));
    hrvec = mean(hrvec(find(eye(Ntt)),:), 1);
    p = signrank(hrvec - 1/Ntt);
    subplot(2,2,b)
    imagesc(mean(tempHR, 3))
    caxis([0 1]); colorbar
    set(gca, 'fontsize', fs, 'XTick', 1:numel(xtl), 'XTickLabel', xtl, 'YTick', 1:numel(xtl), 'YTickLabel', xtl)
    title(sprintf('%s %.2f\np=%.4f', whichICblock, mean(hrvec), p) )
end
colormap jet

figure;
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', [preproc ' SVM ' whichvisarea ' probe accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for b = 1:numel(ICblocknames)
    whichICblock = ICblocknames{b};
    subplot(2,2,b)
    tempHR = squeeze(nanmean(HR_SVMtrainICRC.(whichICblock).(whichvisarea).REt, 3 ));
    imagesc(mean(tempHR, 3))
    caxis([0 1]); colorbar
    set(gca, 'fontsize', fs, 'XTick', 1:numel(xtl), 'XTickLabel', xtl, 'YTick', 1:numel(ytl), 'YTickLabel', ytl)
    infscore = squeeze( (( tempHR(1,1,:)-tempHR(1,2,:) )+( tempHR(2,4,:)-tempHR(2,3,:) ))/2 );
    p = signrank(infscore);
    title(sprintf('%s IC-RC %.2f\np=%.4f', whichICblock, mean(infscore), p) )
end
colormap jet

%% compare areas
fs = 14;
whichICblock = 'ICwcfg1';
xtl = {'IC1', 'RC1', 'RC2', 'IC2'};
ytl = {'REt1', 'REt2'};

figure;
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', [preproc ' SVM ' whichICblock ' test accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    tempHR = squeeze(nanmean(HR_SVMtrainICRC.(whichICblock).(whichvisarea).test, 3 ));
    hrvec = reshape(tempHR, Ntt^2, size(tempHR,3));
    hrvec = mean(hrvec(find(eye(Ntt)),:), 1);
    p = signrank(hrvec - 1/Ntt);
    subplot(2,3,a)
    imagesc(nanmean(tempHR, 3))
    caxis([0 1]); colorbar
    set(gca, 'fontsize', fs, 'XTick', 1:numel(xtl), 'XTickLabel', xtl, 'YTick', 1:numel(xtl), 'YTickLabel', xtl)
    title(sprintf('%s %.2f\np=%.4f', whichvisarea, nanmean(hrvec), p) )
end
colormap jet

figure;
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', [preproc ' SVM ' whichICblock ' probe accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    subplot(2,3,a)
    tempHR = squeeze(nanmean(HR_SVMtrainICRC.(whichICblock).(whichvisarea).REt, 3 ));
    imagesc(nanmean(tempHR, 3))
    caxis([0 1]); colorbar
    set(gca, 'fontsize', fs, 'XTick', 1:numel(xtl), 'XTickLabel', xtl, 'YTick', 1:numel(ytl), 'YTickLabel', ytl)
    infscore = squeeze( (( tempHR(1,1,:)-tempHR(1,2,:) )+( tempHR(2,4,:)-tempHR(2,3,:) ))/2 );
    p = signrank(infscore);
    title(sprintf('%s IC-RC %.2f\np=%.4f', whichvisarea, nanmean(infscore), p) )
end
colormap jet

figure;
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', [preproc ' SVM ' whichICblock ' probe accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    subplot(2,3,a)
    tempHR = squeeze(nanmean(HR_SVMtrainICRC.(whichICblock).(whichvisarea).test, 3 ));
    hrvec = reshape(tempHR, Ntt^2, size(tempHR,3));
    hrvec = mean(hrvec(find(eye(Ntt)),:), 1)';
    tempHR = squeeze(nanmean(HR_SVMtrainICRC.(whichICblock).(whichvisarea).REt, 3 ));
    infscore = squeeze( (( tempHR(1,1,:)-tempHR(1,2,:) )+( tempHR(2,4,:)-tempHR(2,3,:) ))/2 );
    plot(hrvec, infscore, 'o')
    %axis([0 1 -0.5 0.5])
    set(gca, 'fontsize', fs, 'XGrid', 'on', 'YGrid', 'on')
    [rho, p] = corr(hrvec, infscore, 'type', 'spearman', 'rows', 'complete');
    title(sprintf('%s Spearman Corr %.2f\np=%.4f', whichvisarea, rho, p) )
end
