% datadir = 'D:\OpenScopeData\000248\';
% nwbdir = dir(datadir);
% nwbsessions = {nwbdir.name};
% nwbsessions = nwbsessions(contains(nwbsessions, 'sub-'));
% Nsessions = numel(nwbsessions)-1;

% datadir = '/Users/hyeyoung/Documents/OpenScopeData/';
% nwbsessions = {'sub-1171903426','sub-1172969386','sub-1174569632','sub-1175512776', ...
%     'sub-1177693335','sub-1181585601','sub-1182593224','sub-1183369796','sub-1186544719'};
% Nsessions = numel(nwbsessions);

preproc = 'zscore'; % '' is z-score train trials, '_zscoreall', or '_meancenter'
svmdesc = 'trainICRCtestRE';
whichSVMkernel = 'Linear';

justctx = true;
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

SVMtrainICRCagg = struct();
tic
for a = 1:numel(visareas)
    whichvisarea = visareas{a};
    for b = 1:numel(ICblocknames)
        whichICblock = ICblocknames{b};
        for ises = 1:Nsessions
            pathsvm = [pathsv nwbsessions{ises} filesep];
            svmfn = strcat(pathsvm, 'SVM_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_silencesubsets_', whichICblock, '.mat');
            if exist(svmfn, 'file')
                load(svmfn, 'SVMtrainICRC')
                SVMtrainICRCagg(ises).(whichICblock).(whichvisarea) = SVMtrainICRC;
            else
                warning([svmfn ' does not exist'])
            end
        end
    end
end
toc

%%
SVMtrainICRC = SVMtrainICRCagg(1).(whichICblock).(whichvisarea);
Nsplits=size(SVMtrainICRC.spkcnt.traintrialinds,2); icf=1;

traintrialtypes = SVMtrainICRC.trialtypes';
probetrialtypes = unique(SVMtrainICRC(icf).trialorder( SVMtrainICRC.probetrials ));

Ntt = numel(traintrialtypes);
Nprobett = numel(probetrialtypes);

HR_SVMtrainICRC = struct();
tic
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
    
    lmf = fieldnames(HR_SVMtrainICRC.(whichICblock).(whichvisarea));
    for ises = 1:Nsessions
        if ~isfield( SVMtrainICRCagg(ises).(whichICblock), whichvisarea )
            continue
        end
        for lm = 1:numel(lmf)
            for isplit = 1:Nsplits
                % tempyind: input trial indices, tempynn: output trial type, tempyu: input trial type labels
                switch lmf{lm}
                    case 'train'
                        tempyind = SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).spkcnt.traintrialinds(:,isplit);
                        tempysvm = SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).spkcnt.train.label(:,isplit);
                        tempyu = traintrialtypes;
                    case 'test'
                        tempyind = SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).spkcnt.testtrialinds(:,isplit);
                        tempysvm = SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).spkcnt.test.label(:,isplit);
                        tempyu = traintrialtypes;
                    case 'probe'
                        tempyind = SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).probetrials;
                        tempysvm = SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).spkcnt.probe.label(:,isplit);
                        tempyu = probetrialtypes;
                    case 'REt'
                        tempyu = [1105, 1109];
                        probetrialinds = SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).probetrials;
                        tempyind= find(ismember(SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).trialorder , tempyu));
                        tempoi = ismember(probetrialinds, tempyind);
                        tempysvm = SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).spkcnt.probe.label(tempoi,isplit);
                    case 'T'
                        tempyu = [105, 109];
                        probetrialinds = SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).probetrials;
                        tempyind= find(ismember(SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).trialorder , tempyu));
                        tempoi = ismember(probetrialinds, tempyind);
                        tempysvm = SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).spkcnt.probe.label(tempoi,isplit);
                        if nnz(tempoi)==0
                            continue
                        end
                    case 'REx'
                        tempyu = [1201, 1299];
                        probetrialinds = SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).probetrials;
                        tempyind= find(ismember(SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).trialorder , tempyu));
                        tempoi = ismember(probetrialinds, tempyind);
                        tempysvm = SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).spkcnt.probe.label(tempoi,isplit);
                    case 'X'
                        tempyu = [101, 1201, 1299];
                        probetrialinds = SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).probetrials;
                        tempyind= find(ismember(SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).trialorder , tempyu));
                        tempoi = ismember(probetrialinds, tempyind);
                        tempysvm = SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).spkcnt.probe.label(tempoi,isplit);
                    case 'blank'
                        tempyind = SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).alltrials;
                        tempysvm = SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).spkcnt.all.label(:,isplit);
                        tempyu = 0;
                end
                tempysvm = cellfun(@str2num, tempysvm);
                tempy = SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).trialorder(tempyind);
                % rearrange ylabs order so that it is consistent across different sessions
                [sylabv,sylabi]=sort(SVMtrainICRCagg(ises).(whichICblock).(whichvisarea).spkcnt.Ylabs(:,isplit));
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
                                
            end
        end
    end
    
    % %% print results
    if strcmp(whichICblock, 'ICwcfg1') && (strcmp(whichvisarea, 'VISp') || strcmp(whichvisarea, 'C'))
        disp(['SVMtrainICRC ', whichICblock])
        disp('train HR')
        disp(squeeze(mean(mean(HR_SVMtrainICRC.(whichICblock).(whichvisarea).train,3),4)))
        disp('test HR')
        disp(squeeze(mean(mean(HR_SVMtrainICRC.(whichICblock).(whichvisarea).test,3),4)))
        disp('REt HR')
        disp(squeeze(mean(mean(HR_SVMtrainICRC.(whichICblock).(whichvisarea).REt,3),4)))
    end
    end
end
toc
save([pathsv 'SVM_trainICRC_' preproc '_agg.mat'], 'SVMtrainICRCagg', 'HR_SVMtrainICRC', '-v7.3')

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
xtl = {'IC1', 'RC1', 'RC2', 'IC2'};
figure;
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', [preproc ' SVM V1 test accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for b = 1:numel(ICblocknames)
    whichICblock = ICblocknames{b};
    subplot(2,2,b)
    tempHR = squeeze(nanmean(HR_SVMtrainICRC.(whichICblock).VISp.test, [3 4] ));
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
    tempHR = squeeze(nanmean(HR_SVMtrainICRC.(whichICblock).VISp.REt, [3 4] ));
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
        tempHR = squeeze(nanmean(HR_SVMtrainICRC.(whichICblock).(whichvisarea).test(:,:,:,ises), 3 ));
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
        tempHR = squeeze(nanmean(HR_SVMtrainICRC.(whichICblock).(whichvisarea).REt(:,:,:,ises), 3 ));
        imagesc(tempHR)
        caxis([0 1]); colorbar
        set(gca, 'fontsize', fs, 'XTick', 1:4, 'XTickLabel', xtl, 'YTick', 1:4, 'YTickLabel', ytl)
        infscore = squeeze( (( tempHR(1,1,:)-tempHR(1,2,:) )+( tempHR(2,4,:)-tempHR(2,3,:) ))/2 );
        title(sprintf('Session%d %s IC-RC %.4f', ises, whichvisarea, mean(infscore)) )
    end
end
colormap jet
