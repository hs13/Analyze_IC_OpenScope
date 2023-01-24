datadir = 'D:\OpenScopeData\000248\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name}; 
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

% visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
% ICblocknames = {'ICkcfg0', 'ICkcfg1', 'ICwcfg0', 'ICwcfg1'};

whichvisarea = 'V1';
whichICblock = 'ICwcfg1';

svmdesc = 'trainICRCtestRE';
preproc = 'zscore'; % '' is z-score train trials, '_zscoreall', or '_meancenter'
whichSVMkernel = 'Linear';
pathsv = [datadir 'postprocessed' filesep 'SVM' filesep 'SVM_' svmdesc filesep];

SVMsilICRCagg = struct();
tic
for ises = 1:Nsessions
    pathsvm = [pathsv nwbsessions{ises} filesep];
    svmfn = strcat(pathsvm, 'SVM_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_silencesubsets_', whichICblock, '.mat');
    if exist(svmfn, 'file')
        load(svmfn, 'SVMtrainICRC')
        SVMsilICRCagg(ises).(whichvisarea) = SVMtrainICRC;
    else
        warning([svmfn ' does not exist'])
    end
end
toc

% %%
silencedescs = SVMtrainICRC.silencedescs;
lmf = {'train', 'test', 'probe', 'REt', 'REx', 'X', 'blank'};

Nsplits=size(SVMtrainICRC.spkcnt.traintrialinds,2); icf=1;

traintrialtypes = SVMtrainICRC.trialtypes';
probetrialtypes = unique(SVMtrainICRC(icf).trialorder( SVMtrainICRC.probetrials ));

Ntt = numel(traintrialtypes);
Nprobett = numel(probetrialtypes);

HR_SVMsilICRC.(whichvisarea) = struct();
tic
for isd = 0:numel(silencedescs)
    if isd == 0
        sil = '';
    else
        sil = ['_' silencedescs{isd}];
    end
    HR_SVMsilICRC.(whichvisarea).(['train' sil])= NaN(Ntt,Ntt,Nsplits,Nsessions);
    HR_SVMsilICRC.(whichvisarea).(['test' sil]) = NaN(Ntt,Ntt,Nsplits,Nsessions);
    HR_SVMsilICRC.(whichvisarea).(['probe' sil]) = NaN(Nprobett,Ntt,Nsplits,Nsessions);
    HR_SVMsilICRC.(whichvisarea).(['REt' sil]) = NaN(2,Ntt,Nsplits,Nsessions);
    HR_SVMsilICRC.(whichvisarea).(['T' sil]) = NaN(2,Ntt,Nsplits,Nsessions);
    HR_SVMsilICRC.(whichvisarea).(['REx' sil]) = NaN(2,Ntt,Nsplits,Nsessions);
    HR_SVMsilICRC.(whichvisarea).(['X' sil]) = NaN(3,Ntt,Nsplits,Nsessions);
    HR_SVMsilICRC.(whichvisarea).(['blank' sil]) = NaN(1,Ntt,Nsplits,Nsessions);
    
    for ises = 1:Nsessions
        if isempty(SVMsilICRCagg(ises).(whichvisarea))
            continue
        end
        for lm = 1:numel(lmf)
            for isplit = 1:Nsplits
                % tempyind: input trial indices, tempynn: output trial type, tempyu: input trial type labels
                switch lmf{lm}
                    case 'train'
                        tempyind = SVMsilICRCagg(ises).(whichvisarea).spkcnt.traintrialinds(:,isplit);
                        tempysvm = SVMsilICRCagg(ises).(whichvisarea).spkcnt.(['train' sil]).label(:,isplit);
                        tempyu = traintrialtypes;
                    case 'test'
                        tempyind = SVMsilICRCagg(ises).(whichvisarea).spkcnt.testtrialinds(:,isplit);
                        tempysvm = SVMsilICRCagg(ises).(whichvisarea).spkcnt.(['test' sil]).label(:,isplit);
                        tempyu = traintrialtypes;
                    case 'probe'
                        tempyind = SVMsilICRCagg(ises).(whichvisarea).probetrials;
                        tempysvm = SVMsilICRCagg(ises).(whichvisarea).spkcnt.(['probe' sil]).label(:,isplit);
                        tempyu = probetrialtypes;
                    case 'REt'
                        tempyu = [1105, 1109];
                        probetrialinds = SVMsilICRCagg(ises).(whichvisarea).probetrials;
                        tempyind= find(ismember(SVMsilICRCagg(ises).(whichvisarea).trialorder , tempyu));
                        tempoi = ismember(probetrialinds, tempyind);
                        tempysvm = SVMsilICRCagg(ises).(whichvisarea).spkcnt.(['probe' sil]).label(tempoi,isplit);
                    case 'T'
                        tempyu = [105, 109];
                        probetrialinds = SVMsilICRCagg(ises).(whichvisarea).probetrials;
                        tempyind= find(ismember(SVMsilICRCagg(ises).(whichvisarea).trialorder , tempyu));
                        tempoi = ismember(probetrialinds, tempyind);
                        tempysvm = SVMsilICRCagg(ises).(whichvisarea).spkcnt.(['probe' sil]).label(tempoi,isplit);
                        if nnz(tempoi)==0
                            continue
                        end
                    case 'REx'
                        tempyu = [1201, 1299];
                        probetrialinds = SVMsilICRCagg(ises).(whichvisarea).probetrials;
                        tempyind= find(ismember(SVMsilICRCagg(ises).(whichvisarea).trialorder , tempyu));
                        tempoi = ismember(probetrialinds, tempyind);
                        tempysvm = SVMsilICRCagg(ises).(whichvisarea).spkcnt.(['probe' sil]).label(tempoi,isplit);
                    case 'X'
                        tempyu = [101, 1201, 1299];
                        probetrialinds = SVMsilICRCagg(ises).(whichvisarea).probetrials;
                        tempyind= find(ismember(SVMsilICRCagg(ises).(whichvisarea).trialorder , tempyu));
                        tempoi = ismember(probetrialinds, tempyind);
                        tempysvm = SVMsilICRCagg(ises).(whichvisarea).spkcnt.(['probe' sil]).label(tempoi,isplit);
                    case 'blank'
                        tempyind = SVMsilICRCagg(ises).(whichvisarea).alltrials;
                        tempysvm = SVMsilICRCagg(ises).(whichvisarea).spkcnt.all.label(:,isplit);
                        tempyu = 0;
                end
                tempysvm = cellfun(@str2num, tempysvm);
                tempy = SVMsilICRCagg(ises).(whichvisarea).trialorder(tempyind);
                % rearrange ylabs order so that it is consistent across different sessions
                [sylabv,sylabi]=sort(SVMsilICRCagg(ises).(whichvisarea).spkcnt.Ylabs(:,isplit));
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
                        HR_SVMsilICRC.(whichvisarea).(lmf{lm})(iy,jy,isplit,ises) = nnz(tempysvm(trialsoy)==traintrialtypes(jy))/nnz(trialsoy);
                    end
                end
                
            end
        end
    end
end
toc

save([pathsv 'HR_SVMsilICRC_' preproc '_' whichICblock '_' whichvisarea '_agg.mat'], 'HR_SVMsilICRC', '-v7.3')
save(['G:\My Drive\DATA\ICexpts_submission22\openscope_HR_SVMsilICRC_' preproc '_' whichICblock '_' whichvisarea '_agg.mat'], 'HR_SVMsilICRC', '-v7.3')

% %% print results
disp(['SVMtrainICRC', whichSVMkernel])
disp('train HR')
disp(squeeze(mean(mean(HR_SVMsilICRC.(whichvisarea).train,3),4)))
disp('test HR')
disp(squeeze(mean(mean(HR_SVMsilICRC.(whichvisarea).test,3),4)))
disp('probe HR')
disp(squeeze(mean(mean(HR_SVMsilICRC.(whichvisarea).probe,3),4)))

disp('REt HR')
disp(squeeze(mean(mean(HR_SVMsilICRC.(whichvisarea).REt,3),4)))
toc

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
            if ~isequaln(HR_SVMsilICRC.(whichvisarea).probe(probetrialtypes==0,:,:,:), ...
                    HR_SVMsilICRC.(whichvisarea).blank)
                error('check base HR')
            end
            % add a tiny little offset
            tempHR = tempHR+2^-10;
            tempbase = HR_SVMsilICRC.(whichvisarea).blank+2^-10;
            dprime_SVMtrainICRC.(whichvisarea).(lmf{lm}) = norminv(tempHR) - norminv(tempbase);
            % pcb: percent change from baseline, baseline being probability of any *test* trial being classified as a certain output class
            pcb_SVMtrainICRC.(whichvisarea).(lmf{lm}) = 100*(tempHR./tempbase -1);
            % normHR: percent change from baseline, baseline being probability of any *test* trial being classified as a certain output class
            lognormHR_SVMtrainICRC.(whichvisarea).(lmf{lm}) = log2(tempHR./tempbase);
        end
    end
end


%% compare blocks V1
fs = 14;
xtl = {'IC1', 'RC1', 'RC2', 'IC2'};
figure;
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', [preproc ' SVM V1 test accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for ises = 1:4
    for isd = 0:numel(silencedescs)
        subplot(4,numel(silencedescs)+1,(numel(silencedescs)+1)*(ises-1)+isd)
        tempHR = squeeze(nanmean(HR_SVMsilICRC.VISp.test(:,:,:,ises), 3 ));
        imagesc(tempHR)
        caxis([0 1]); colorbar
        set(gca, 'fontsize', fs, 'XTick', 1:4, 'XTickLabel', xtl, 'YTick', 1:4, 'YTickLabel', xtl)
        title(sprintf('Session%d %s %.4f', ises, whichICblock, mean(diag(tempHR))) )
    end
end
colormap jet

ytl = {'REt1', 'REt2'};
figure;
annotation('textbox', [0.1 0.91 0.8 0.1], 'string', [preproc ' SVM V1 probe accuacy'], 'edgecolor', 'none', 'fontsize', fs)
for ises = 1:4
    for isd = 0:numel(silencedescs)
        subplot(4,numel(silencedescs)+1,(numel(silencedescs)+1)*(ises-1)+isd)
        tempHR = squeeze(nanmean(HR_SVMsilICRC.VISp.REt(:,:,:,ises), 3 ));
        imagesc(tempHR)
        caxis([0 1]); colorbar
        set(gca, 'fontsize', fs, 'XTick', 1:4, 'XTickLabel', xtl, 'YTick', 1:4, 'YTickLabel', ytl)
        infscore = squeeze( (( tempHR(1,1,:)-tempHR(1,2,:) )+( tempHR(2,4,:)-tempHR(2,3,:) ))/2 );
        title(sprintf('Session%d %s IC-RC %.4f', ises, whichICblock, mean(infscore)) )
    end
end
colormap jet
