% aggregateICtxi_SVM_trainRExtestICRC_loadall.m requires too much RAM
% instead of aggregating SVMtrainRExagg, compute HR_SVMtrainREx while loading SVMtrainREx

%% takes 48 min
datadir = 'D:\OpenScopeData\000248\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name}; 
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

svmdesc = 'trainRExtestICRC';
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
load(svmfn, 'SVMtrainREx')

lmf = {'train', 'test', 'probe', 'blank'};
traintrialtypes = SVMtrainREx.trialtypes';
probetrialtypes = unique(SVMtrainREx.trialorder( SVMtrainREx.probetrials ));
Ntt = numel(traintrialtypes);
Nprobett = numel(probetrialtypes);
Nsplits=size(SVMtrainREx.spkcnt.traintrialinds,2);
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    for b = 1:numel(ICblocknames)
        whichICblock = ICblocknames{b};
        HR_SVMtrainREx.(whichICblock).(whichvisarea) = struct();
        HR_SVMtrainREx.(whichICblock).(whichvisarea).train = NaN(Ntt,Ntt,Nsplits,Nsessions);
        HR_SVMtrainREx.(whichICblock).(whichvisarea).test = NaN(Ntt,Ntt,Nsplits,Nsessions);
        HR_SVMtrainREx.(whichICblock).(whichvisarea).probe = NaN(Nprobett,Ntt,Nsplits,Nsessions);
        HR_SVMtrainREx.(whichICblock).(whichvisarea).blank = NaN(1,Ntt,Nsplits,Nsessions);
    end
end

% % load and compute
%SVMtrainRExagg = struct();
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
            load(svmfn, 'SVMtrainREx')
            %SVMtrainRExagg(ises).(whichICblock).(whichvisarea) = SVMtrainREx;
            
            for lm = 1:numel(lmf)
                for isplit = 1:Nsplits
                    % tempyind: input trial indices, tempynn: output trial type, tempyu: input trial type labels
                    switch lmf{lm}
                        case 'train'
                            tempyind = SVMtrainREx.spkcnt.traintrialinds(:,isplit);
                            tempysvm = SVMtrainREx.spkcnt.train.label(:,isplit);
                            tempyu = traintrialtypes;
                        case 'test'
                            tempyind = SVMtrainREx.spkcnt.testtrialinds(:,isplit);
                            tempysvm = SVMtrainREx.spkcnt.test.label(:,isplit);
                            tempyu = traintrialtypes;
                        case 'probe'
                            tempyind = SVMtrainREx.probetrials;
                            tempysvm = SVMtrainREx.spkcnt.probe.label(:,isplit);
                            tempyu = probetrialtypes;
                        case 'blank'
                            tempyind = SVMtrainREx.alltrials;
                            tempysvm = SVMtrainREx.spkcnt.all.label(:,isplit);
                            tempyu = 0;
                    end
                    tempysvm = cellfun(@str2num, tempysvm);
                    tempy = SVMtrainREx.trialorder(tempyind);
                    % rearrange ylabs order so that it is consistent across different sessions
                    [sylabv,sylabi]=sort(SVMtrainREx.spkcnt.Ylabs(:,isplit));
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
                            HR_SVMtrainREx.(whichICblock).(whichvisarea).(lmf{lm})(iy,jy,isplit,ises) = nnz(tempysvm(trialsoy)==traintrialtypes(jy))/nnz(trialsoy);
                        end
                    end
                end % isplit
            end % lm
        end % ises
        toc
    end % b
end % a
toc
save([pathsv 'HR_SVMtrainREx_' preproc '_agg.mat'], 'HR_SVMtrainREx', '-v7.3')
save(['G:\My Drive\DATA\ICexpts_submission22\openscope_HR_SVMtrainREx_' preproc '_agg.mat'], 'HR_SVMtrainREx', '-v7.3')

%% adjusted performance metrics
if ~exist('dprime_SVMtrainREx', 'var')
    dprime_SVMtrainREx = struct();
    pcb_SVMtrainREx = struct();
    lognormHR_SVMtrainREx = struct();
end

for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    for b = 1:numel(ICblocknames)
        whichICblock = ICblocknames{b};
        for lm = 1:numel(lmf)
            if ~isequaln(HR_SVMtrainREx.(whichICblock).(whichvisarea).probe(probetrialtypes==0,:,:,:), ...
                    HR_SVMtrainREx.(whichICblock).(whichvisarea).blank)
                error('check base HR')
            end
            % add a tiny little offset
            tempHR = tempHR+2^-10;
            tempbase = HR_SVMtrainREx.(whichICblock).(whichvisarea).blank+2^-10;
            dprime_SVMtrainREx.(whichICblock).(whichvisarea).(lmf{lm}) = norminv(tempHR) - norminv(tempbase);
            % pcb: percent change from baseline, baseline being probability of any *test* trial being classified as a certain output class
            pcb_SVMtrainREx.(whichICblock).(whichvisarea).(lmf{lm}) = 100*(tempHR./tempbase -1);
            % normHR: percent change from baseline, baseline being probability of any *test* trial being classified as a certain output class
            lognormHR_SVMtrainREx.(whichICblock).(whichvisarea).(lmf{lm}) = log2(tempHR./tempbase);
        end
    end
end


%% compare blocks V1
whichvisarea = 'V1';
fs = 12;
xtl = {'REx1', 'REx2'};
ytl = {'IC1', 'RC1', 'RC2', 'IC2'};

figure;
annotation('textbox', [0.1 0.92 0.8 0.1], 'string', [preproc ' SVM ' whichvisarea ' test accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for b = 1:numel(ICblocknames)
    whichICblock = ICblocknames{b};
    tempHR = squeeze(nanmean(HR_SVMtrainREx.(whichICblock).(whichvisarea).test, 3 ));
    hrvec = reshape(tempHR, length(traintrialtypes)^2, size(tempHR,3));
    hrvec = mean(hrvec(find(eye(length(traintrialtypes))),:), 1);
    p = signrank(hrvec - 1/length(traintrialtypes));
    subplot(2,2,b)
    imagesc(mean(tempHR, 3))
    caxis([0 1]); colorbar
    set(gca, 'fontsize', fs, 'XTick', 1:numel(xtl), 'XTickLabel', xtl, 'YTick', 1:numel(xtl), 'YTickLabel', xtl)
    title(sprintf('%s %.2f p=%.4f', whichICblock, mean(hrvec), p) )
end
colormap jet

figure;
annotation('textbox', [0.1 0.92 0.8 0.1], 'string', [preproc ' SVM ' whichvisarea ' probe accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for b = 1:numel(ICblocknames)
    whichICblock = ICblocknames{b};
    subplot(2,2,b)
    tempHR = squeeze(nanmean(HR_SVMtrainREx.(whichICblock).(whichvisarea).probe, 3 ));
    imagesc(mean(tempHR, 3))
    caxis([0 1]); colorbar
    set(gca, 'fontsize', fs, 'XTick', 1:numel(xtl), 'XTickLabel', xtl, 'YTick', 1:numel(ytl), 'YTickLabel', ytl)
    infscore = squeeze( (( tempHR(1,1,:)-0.5 )+( tempHR(4,2,:)-0.5 ))/2 );
    p = signrank(infscore);
    title(sprintf('%s IC-chance\nmean %.2f p=%.4f', whichICblock, mean(infscore), p) )
end
colormap jet

%% compare areas
whichICblock = 'ICwcfg1';
fs = 10;
xtl = {'REx1', 'REx2'};
ytl = {'IC1', 'RC1', 'RC2', 'IC2'};

figure;
annotation('textbox', [0.1 0.92 0.8 0.1], 'string', [preproc ' SVM ' whichICblock ' test accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    tempHR = squeeze(nanmean(HR_SVMtrainREx.(whichICblock).(whichvisarea).test, 3 ));
    hrvec = reshape(tempHR, length(traintrialtypes)^2, size(tempHR,3));
    hrvec = mean(hrvec(find(eye(length(traintrialtypes))),:), 1);
    p = signrank(hrvec - 1/length(traintrialtypes));
    subplot(2,3,a)
    imagesc(nanmean(tempHR, 3))
    caxis([0 1]); colorbar
    set(gca, 'fontsize', fs, 'XTick', 1:numel(xtl), 'XTickLabel', xtl, 'YTick', 1:numel(xtl), 'YTickLabel', xtl)
    title(sprintf('%s %.2f p=%.4f', whichvisarea, nanmean(hrvec), p) )
end
colormap jet

figure;
annotation('textbox', [0.1 0.92 0.8 0.1], 'string', [preproc ' SVM ' whichICblock ' probe accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    subplot(2,3,a)
    tempHR = squeeze(nanmean(HR_SVMtrainREx.(whichICblock).(whichvisarea).probe, 3 ));
    imagesc(nanmean(tempHR, 3))
    caxis([0 1]); colorbar
    set(gca, 'fontsize', fs, 'XTick', 1:numel(xtl), 'XTickLabel', xtl, 'YTick', 1:numel(ytl), 'YTickLabel', ytl)
    infscore = squeeze( (( tempHR(1,1,:)-0.5 )+( tempHR(4,2,:)-0.5 ))/2 );
    p = signrank(infscore);
    title(sprintf('%s IC-chance\nmean %.2f p=%.4f', whichvisarea, nanmean(infscore), p) )
end
colormap jet
