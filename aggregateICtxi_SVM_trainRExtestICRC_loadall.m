% datadir = 'D:\OpenScopeData\000248\';
% nwbdir = dir(datadir);
% nwbsessions = {nwbdir.name};
% nwbsessions = nwbsessions(contains(nwbsessions, 'sub-'));
% Nsessions = numel(nwbsessions)-1;

datadir = '/Users/hyeyoung/Documents/OpenScopeData/';
nwbsessions = {'sub-1171903426','sub-1172969386','sub-1174569632','sub-1175512776', ...
    'sub-1177693335','sub-1181585601','sub-1182593224','sub-1183369796','sub-1186544719'};
Nsessions = numel(nwbsessions);

preproc = 'zscore'; % '' is z-score train trials, '_zscoreall', or '_meancenter'
svmdesc = 'trainRExtestICRC';
whichSVMkernel = 'Linear';

justctx = false;
if justctx
    visareas = {'VISp', 'VISl', 'VISrl', 'VISal', 'VISpm', 'VISam'};
    pathsv = [datadir 'SVM_' svmdesc filesep];
else
    visareas = {'A', 'B', 'C', 'D', 'E', 'F'};
    pathsv = [datadir 'SVM_' svmdesc '_allunits' filesep];
end
ICblocknames = {'ICkcfg0', 'ICkcfg1', 'ICwcfg0', 'ICwcfg1'};

% visareas = {'VISp'};
% ICblocknames = {'ICwcfg1'};

SVMtrainRExagg = struct();
tic
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    for b = 1:numel(ICblocknames)
        whichICblock = ICblocknames{b};
        for ises = 1:Nsessions
            pathsvm = [pathsv nwbsessions{ises} filesep];
            svmfn = strcat(pathsvm, 'SVM_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_silencesubsets_', whichICblock, '.mat');
            if exist(svmfn, 'file')
                load(svmfn, 'SVMtrainREx')
                SVMtrainRExagg(ises).(whichICblock).(whichvisarea) = SVMtrainREx;
            end
        end
    end
end
toc

%%
SVMtrainREx = SVMtrainRExagg(1).(whichICblock).(whichvisarea);
Nsplits=size(SVMtrainREx.spkcnt.traintrialinds,2); icf=1;

traintrialtypes = SVMtrainREx.trialtypes';
probetrialtypes = unique(SVMtrainREx(icf).trialorder( SVMtrainREx.probetrials ));

Ntt = numel(traintrialtypes);
Nprobett = numel(probetrialtypes);

HR_SVMtrainREx = struct();
tic
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    for b = 1:numel(ICblocknames)
        whichICblock = ICblocknames{b};
        HR_SVMtrainREx.(whichICblock).(whichvisarea) = struct();
        HR_SVMtrainREx.(whichICblock).(whichvisarea).train = NaN(Ntt,Ntt,Nsplits,Nsessions);
        HR_SVMtrainREx.(whichICblock).(whichvisarea).test = NaN(Ntt,Ntt,Nsplits,Nsessions);
        HR_SVMtrainREx.(whichICblock).(whichvisarea).probe = NaN(Nprobett,Ntt,Nsplits,Nsessions);
        HR_SVMtrainREx.(whichICblock).(whichvisarea).blank = NaN(1,Ntt,Nsplits,Nsessions);
        
        lmf = fieldnames(HR_SVMtrainREx.(whichICblock).(whichvisarea));
        for ises = 1:Nsessions
            if ~isfield( SVMtrainRExagg(ises).(whichICblock), whichvisarea )
                continue
            end
            for lm = 1:numel(lmf)
                for isplit = 1:Nsplits
                    % tempyind: input trial indices, tempynn: output trial type, tempyu: input trial type labels
                    switch lmf{lm}
                        case 'train'
                            tempyind = SVMtrainRExagg(ises).(whichICblock).(whichvisarea).spkcnt.traintrialinds(:,isplit);
                            tempysvm = SVMtrainRExagg(ises).(whichICblock).(whichvisarea).spkcnt.train.label(:,isplit);
                            tempyu = traintrialtypes;
                        case 'test'
                            tempyind = SVMtrainRExagg(ises).(whichICblock).(whichvisarea).spkcnt.testtrialinds(:,isplit);
                            tempysvm = SVMtrainRExagg(ises).(whichICblock).(whichvisarea).spkcnt.test.label(:,isplit);
                            tempyu = traintrialtypes;
                        case 'probe'
                            tempyind = SVMtrainRExagg(ises).(whichICblock).(whichvisarea).probetrials;
                            tempysvm = SVMtrainRExagg(ises).(whichICblock).(whichvisarea).spkcnt.probe.label(:,isplit);
                            tempyu = probetrialtypes;
                        case 'blank'
                            tempyind = SVMtrainRExagg(ises).(whichICblock).(whichvisarea).alltrials;
                            tempysvm = SVMtrainRExagg(ises).(whichICblock).(whichvisarea).spkcnt.all.label(:,isplit);
                            tempyu = 0;
                    end
                    tempysvm = cellfun(@str2num, tempysvm);
                    tempy = SVMtrainRExagg(ises).(whichICblock).(whichvisarea).trialorder(tempyind);
                    % rearrange ylabs order so that it is consistent across different sessions
                    [sylabv,sylabi]=sort(SVMtrainRExagg(ises).(whichICblock).(whichvisarea).spkcnt.Ylabs(:,isplit));
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
                    
                end
            end
        end
        
    % %% print results
    if strcmp(whichICblock, 'ICwcfg1') && (strcmp(whichvisarea, 'VISp') || strcmp(whichvisarea, 'C'))
        disp(['SVMtrainREx ', whichICblock])
            disp('train HR')
            disp(squeeze(mean(mean(HR_SVMtrainREx.(whichICblock).(whichvisarea).train,3),4)))
            disp('test HR')
            disp(squeeze(mean(mean(HR_SVMtrainREx.(whichICblock).(whichvisarea).test,3),4)))
            disp('probe HR')
            disp(squeeze(mean(mean(HR_SVMtrainREx.(whichICblock).(whichvisarea).probe,3),4)))
        end
    end
end
toc
save([pathsv 'SVM_trainREx_' preproc '_agg.mat'], 'SVMtrainRExagg', 'HR_SVMtrainREx', '-v7.3')

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
fs = 14;
xtl = {'IC1', 'RC1', 'RC2', 'IC2'};
figure;
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', [preproc ' SVM V1 test accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for b = 1:numel(ICblocknames)
    whichICblock = ICblocknames{b};
    subplot(2,2,b)
    tempHR = squeeze(nanmean(HR_SVMtrainREx.(whichICblock).VISp.test, [3 4] ));
    imagesc(tempHR)
    caxis([0 1]); colorbar
    set(gca, 'fontsize', fs, 'XTick', 1:4, 'XTickLabel', xtl, 'YTick', 1:4, 'YTickLabel', xtl)
    title(sprintf('Session%d %s %.4f', ises, whichICblock, mean(diag(tempHR))) )
end
colormap jet

ytl = {'REt1', 'REt2'};
figure;
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', [preproc ' SVM V1 probe accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for b = 1:numel(ICblocknames)
    whichICblock = ICblocknames{b};
    subplot(2,2,b)
    tempHR = squeeze(nanmean(HR_SVMtrainREx.(whichICblock).VISp.REt, [3 4] ));
    imagesc(tempHR)
    caxis([0 1]); colorbar
    set(gca, 'fontsize', fs, 'XTick', 1:4, 'XTickLabel', xtl, 'YTick', 1:4, 'YTickLabel', ytl)
    infscore = squeeze( (( tempHR(1,1,:)-tempHR(1,2,:) )+( tempHR(2,4,:)-tempHR(2,3,:) ))/2 );
    title(sprintf('Session%d %s IC-RC %.4f', ises, whichICblock, mean(infscore)) )
end
colormap jet

%% compare areas
whichICblock = 'ICwcfg1';
fs = 12;
xtl = {'IC1', 'RC1', 'RC2', 'IC2'};
figure;
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', [preproc ' SVM ' whichICblock ' test accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for ises = 1:4
    for a = 1:numel(visareas)
        whichvisarea = visareas{a};
        subplot(4,numel(visareas),numel(visareas)*(ises-1)+a)
        tempHR = squeeze(nanmean(HR_SVMtrainREx.(whichICblock).(whichvisarea).test(:,:,:,ises), 3 ));
        imagesc(tempHR)
        caxis([0 1]); colorbar
        set(gca, 'fontsize', fs, 'XTick', 1:4, 'XTickLabel', xtl, 'YTick', 1:4, 'YTickLabel', xtl)
        title(sprintf('Session%d %s %.4f', ises, whichvisarea, mean(diag(tempHR))) )
    end
end
colormap jet

ytl = {'REt1', 'REt2'};
figure;
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', [preproc ' SVM ' whichICblock ' probe accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for ises = 1:4
    for a = 1:numel(visareas)
        whichvisarea = visareas{a};
        subplot(4,numel(visareas),numel(visareas)*(ises-1)+a)
        tempHR = squeeze(nanmean(HR_SVMtrainREx.(whichICblock).(whichvisarea).REt(:,:,:,ises), 3 ));
        imagesc(tempHR)
        caxis([0 1]); colorbar
        set(gca, 'fontsize', fs, 'XTick', 1:4, 'XTickLabel', xtl, 'YTick', 1:4, 'YTickLabel', ytl)
        infscore = squeeze( (( tempHR(1,1,:)-tempHR(1,2,:) )+( tempHR(2,4,:)-tempHR(2,3,:) ))/2 );
        title(sprintf('Session%d %s IC-RC %.4f', ises, whichvisarea, mean(infscore)) )
    end
end
colormap jet
