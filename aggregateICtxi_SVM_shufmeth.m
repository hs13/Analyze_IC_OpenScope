whichSVMkernel = 'Linear';
whichvisarea = 'VISp'; preproc = 'zscore';
% whichvisarea = 'VISp'; preproc = 'meancenter';
% whichvisarea = 'VISpfiltRS'; preproc = 'meancenter';
whichICblock = 'ICwcfg1_presentations';

svmdesc = 'trainICRC_shuf';
traintrialtypes = [106, 107, 110, 111];
probetrialtypes = [1105, 1109, 1201, 1299, 0, 101, 105, 109];

Ntt = numel(traintrialtypes);
Nprobett = numel(probetrialtypes);
Nsplits = 10;

ses2agg = {'sub_1171903433','sub_1172968426','sub_1172969394','sub_1174569641', ...
    'sub_1175512783','sub_1176214862','sub_1177693342','sub_1181314060', ...
    'sub_1181585608','sub_1182593231','sub_1183369803','sub_1186544726', ...
    'sub_1189891322','sub_1194090570'};
Nsessions = numel(ses2agg);

% shufmeth = {'asis', 'shuftrials', 'shuf', 'shufall'};
shufmeth = {'_shuftrials', '_shuf', '_shufall'};%, '_shufallneurons', '_shufneurons'};
mdshufyn = {'', '_shuf', '_msdn', '_mnds'};
lmf = {'train', 'test', 'probe', 'REt', 'T', 'REx', 'X', 'blank'};

HR_SVMtrainICRCshuf = struct();
HR_SVMtrainICRCshuf.Nneurons = zeros(Nsessions,1);
for imeth = 1:numel(shufmeth)
    svmshuf = strcat('spkcnt', shufmeth{imeth});
    for imd = 1:numel(mdshufyn)
        whichmd = mdshufyn{imd};
        HR_SVMtrainICRCshuf.(svmshuf).(strcat('train', whichmd)) = NaN(Ntt,Ntt,Nsplits,Nsessions);
        HR_SVMtrainICRCshuf.(svmshuf).(strcat('test', whichmd)) = NaN(Ntt,Ntt,Nsplits,Nsessions);
        HR_SVMtrainICRCshuf.(svmshuf).(strcat('probe', whichmd)) = NaN(Nprobett,Ntt,Nsplits,Nsessions);
        HR_SVMtrainICRCshuf.(svmshuf).(strcat('REt', whichmd)) = NaN(2,Ntt,Nsplits,Nsessions);
        HR_SVMtrainICRCshuf.(svmshuf).(strcat('T', whichmd)) = NaN(2,Ntt,Nsplits,Nsessions);
        HR_SVMtrainICRCshuf.(svmshuf).(strcat('REx', whichmd)) = NaN(2,Ntt,Nsplits,Nsessions);
        HR_SVMtrainICRCshuf.(svmshuf).(strcat('X', whichmd)) = NaN(3,Ntt,Nsplits,Nsessions);
        HR_SVMtrainICRCshuf.(svmshuf).(strcat('blank', whichmd)) = NaN(1,Ntt,Nsplits,Nsessions);
    end
end

% ~30s per session
for ises = 1:Nsessions
    tic
    pathsvm = strcat('S:/OpenScopeData/000248/postprocessed/SVM/SVM_', svmdesc, '/', ses2agg{ises}, '/');
    svmfn = strcat(pathsvm, 'SVM_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat');
    load(svmfn)
    HR_SVMtrainICRCshuf.Nneurons(ises) = SVMtrainICRC.Nneurons;
    %SVMtrainICRCagg(ises).(whichICblock).(whichvisarea) = SVMtrainICRC;
    for imeth = 1:numel(shufmeth)
        svmshuf = strcat('spkcnt', shufmeth{imeth});
        for imd = 1:numel(mdshufyn)
            whichmd = mdshufyn{imd};
            for isplit = 1:Nsplits
                % tempyind: input trial indices, tempynn: output trial type, tempyu: input trial type labels
                for lm = 1:numel(lmf)
                    switch lmf{lm}
                        case 'train'
                            lmfmd = strcat('train', whichmd);
                            tempyind = SVMtrainICRC.traintrialinds(:,isplit);
                            tempysvm = SVMtrainICRC.(svmshuf).(lmfmd).label(:,isplit);
                            tempyu = traintrialtypes;
                        case 'test'
                            lmfmd = strcat('test', whichmd);
                            tempyind = SVMtrainICRC.testtrialinds(:,isplit);
                            tempysvm = SVMtrainICRC.(svmshuf).(lmfmd).label(:,isplit);
                            tempyu = traintrialtypes;
                        case 'probe'
                            lmfmd = strcat('probe', whichmd);
                            tempyind = SVMtrainICRC.probetrials;
                            tempysvm = SVMtrainICRC.(svmshuf).(lmfmd).label(:,isplit);
                            tempyu = probetrialtypes;
                        case 'REt'
                            lmfmd = strcat('probe', whichmd);
                            tempyu = [1105, 1109];
                            probetrialinds = SVMtrainICRC.probetrials;
                            tempyind= find(ismember(SVMtrainICRC.trialorder , tempyu));
                            tempoi = ismember(probetrialinds, tempyind);
                            tempysvm = SVMtrainICRC.(svmshuf).(lmfmd).label(tempoi,isplit);
                        case 'T'
                            lmfmd = strcat('probe', whichmd);
                            tempyu = [105, 109];
                            probetrialinds = SVMtrainICRC.probetrials;
                            tempyind= find(ismember(SVMtrainICRC.trialorder , tempyu));
                            tempoi = ismember(probetrialinds, tempyind);
                            tempysvm = SVMtrainICRC.(svmshuf).(lmfmd).label(tempoi,isplit);
                            if nnz(tempoi)==0
                                continue
                            end
                        case 'REx'
                            lmfmd = strcat('probe', whichmd);
                            tempyu = [1201, 1299];
                            probetrialinds = SVMtrainICRC.probetrials;
                            tempyind= find(ismember(SVMtrainICRC.trialorder , tempyu));
                            tempoi = ismember(probetrialinds, tempyind);
                            tempysvm = SVMtrainICRC.(svmshuf).(lmfmd).label(tempoi,isplit);
                        case 'X'
                            lmfmd = strcat('probe', whichmd);
                            tempyu = [101, 1201, 1299];
                            probetrialinds = SVMtrainICRC.probetrials;
                            tempyind= find(ismember(SVMtrainICRC.trialorder , tempyu));
                            tempoi = ismember(probetrialinds, tempyind);
                            tempysvm = SVMtrainICRC.(svmshuf).(lmfmd).label(tempoi,isplit);
                        case 'blank'
                            lmfmd = strcat('all', whichmd);
                            tempyind = SVMtrainICRC.alltrials;
                            tempysvm = SVMtrainICRC.(svmshuf).(lmfmd).label(:,isplit);
                            tempyu = 0;
                    end
                    tempysvm = cellfun(@str2num, tempysvm);
                    tempy = SVMtrainICRC.trialorder(tempyind);
                    % rearrange ylabs order so that it is consistent across different sessions
                    [sylabv,sylabi]=sort(SVMtrainICRC.Ylabs(:,isplit));
                    if ~isequal(sylabv, cellstr(num2str(traintrialtypes')))
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
                            HR_SVMtrainICRCshuf.(svmshuf).([lmf{lm} whichmd])(iy,jy,isplit,ises) = nnz(tempysvm(trialsoy)==traintrialtypes(jy))/nnz(trialsoy);
                        end
                    end

                end % lm
            end % isplit
        end
    end
    toc
end

%%
shufmeth = {'_shuftrials', '_shuf', '_shufall'};%, '_shufallneurons', '_shufneurons'};
mdshufyn = {'', '_shuf', '_msdn', '_mnds'};
mdshuf = {'asis', 'shuf', 'msdn', 'mnds'};
% HR_SVMtrainICRCshuf.(svmshuf).(strcat('test', whichmd)) % NaN(Ntt,Ntt,Nsplits,Nsessions);
% HR_SVMtrainICRCshuf.(svmshuf).(strcat('REt', whichmd)) % NaN(2,Ntt,Nsplits,Nsessions);

testHR = struct();
infHR = struct();
infICREt = struct();
infLCREt = struct();
for imeth = 1:numel(shufmeth)
    svmshuf = strcat('spkcnt', shufmeth{imeth});
    for imd = 1:numel(mdshufyn)
        tempmat = reshape( squeeze(nanmean(HR_SVMtrainICRCshuf.(svmshuf).(strcat('test', mdshufyn{imd})), 3)), Ntt*Ntt,Nsessions);
        testHR.(svmshuf).(mdshuf{imd}) = mean(tempmat(find(eye(Ntt)),:), 1)';

        tempmat = squeeze(nanmean(HR_SVMtrainICRCshuf.(svmshuf).(strcat('REt', mdshufyn{imd})), 3)) ;
        infHR.(svmshuf).(mdshuf{imd}) = squeeze( tempmat(1,1,:)-tempmat(1,2,:) + tempmat(2,4,:)-tempmat(2,3,:) )/2;
        infICREt.(svmshuf).(mdshuf{imd}) = squeeze( tempmat(1,1,:) + tempmat(2,4,:) )/2;
        infLCREt.(svmshuf).(mdshuf{imd}) = squeeze( tempmat(1,2,:) + tempmat(2,3,:) )/2;
    end
end

shufmethods = {'asis', 'shuftrials', 'shuf', 'shufall'};
testmat = [testHR.spkcnt_shuftrials.asis testHR.spkcnt_shuftrials.shuf testHR.spkcnt_shuf.shuf testHR.spkcnt_shufall.shuf];
infmat = [infHR.spkcnt_shuftrials.asis infHR.spkcnt_shuftrials.shuf infHR.spkcnt_shuf.shuf infHR.spkcnt_shufall.shuf];
figure; 
subplot(1,2,1)
plot(testmat', 'o-')
set(gca, 'XTick', 1:numel(shufmethods), 'XTickLabel', shufmethods)
title('Test')
subplot(1,2,2)
plot(infmat', 'o-')
set(gca, 'XTick', 1:numel(shufmethods), 'XTickLabel', shufmethods)
title('Inference')


[p,tbl,stats] = friedman(testmat);
figure
[c,m,h]=multcompare(stats);

[p,tbl,stats] = friedman(infmat);
figure
[c,m,h]=multcompare(stats);


signrank(testHR.spkcnt_shuftrials.asis, testHR.spkcnt_shuftrials.shuf)
signrank(infHR.spkcnt_shuftrials.asis, infHR.spkcnt_shuftrials.shuf) % not significant

signrank(infICREt.spkcnt_shuftrials.asis, infICREt.spkcnt_shuftrials.shuf) % trending for VISp-zscore, not significant for VISp-meancenter and VISpfiltRS-meancenter

signrank(testHR.spkcnt_shuf.asis, testHR.spkcnt_shuf.shuf)
signrank(infHR.spkcnt_shuf.asis, infHR.spkcnt_shuf.shuf) % significant

signrank(infHR.spkcnt_shuf.shuf, infHR.spkcnt_shufall.shuf) % significant for VISp-zscore and VISp-meancenter, not significant for VISpfiltRS-meancenter

% mnds (model trained on as-is data, input shuffle data)
shufmethods = {'asis', 'shuftrials', 'shuf', 'shufall'};
testmat = [testHR.spkcnt_shuftrials.asis testHR.spkcnt_shuftrials.mnds testHR.spkcnt_shuf.mnds testHR.spkcnt_shufall.mnds];
infmat = [infHR.spkcnt_shuftrials.asis infHR.spkcnt_shuftrials.mnds infHR.spkcnt_shuf.mnds infHR.spkcnt_shufall.mnds];
figure; 
subplot(1,2,1)
plot(testmat', 'o-')
set(gca, 'XTick', 1:numel(shufmethods), 'XTickLabel', shufmethods)
title('Test')
subplot(1,2,2)
plot(infmat', 'o-')
set(gca, 'XTick', 1:numel(shufmethods), 'XTickLabel', shufmethods)
title('Inference')

signrank(testHR.spkcnt_shuf.asis, testHR.spkcnt_shuf.mnds)
signrank(infHR.spkcnt_shuf.asis, infHR.spkcnt_shuf.mnds)


mcHR = [    0.8356
    0.8494
    0.5394
    0.8456
    0.9294
    0.6975
    0.7281
    0.7094
    0.7300
    0.7675
    0.8438
    0.8981
    0.9094
    0.6038];
mcfiltRSHR =    [0.8294
    0.8175
    0.5269
    0.8287
    0.9219
    0.6931
    0.6900
    0.6031
    0.6531
    0.7275
    0.8237
    0.8231
    0.8931
    0.5913];
signrank(mcHR, testHR.spkcnt_shuftrials.asis ) 
% p=0.1158 when comparing z-score vs mean-center, z-score has slightly
% higher mean and median
figure; hold all
plot(mcHR, testHR.spkcnt_shuftrials.asis, 'o')
xl= xlim; plot(xl, xl, 'r--')
