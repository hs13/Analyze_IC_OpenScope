%   CALCULATE _SHUF, _MSDN, _MNDS FOR ALL OF THESE METHODS

% svmdesc trainICRC_shuf / trainREx_shuf / trainIC1RC1_shuf / trainIC2RC2_shuf
svmdesc = 'trainICRC_shuf';
whichSVMkernel = 'Linear';

disp(['analyzeICtx_SVM_shufmeth ' svmdesc ' ' whichSVMkernel ' STARTED!'])
%%
% addpath('/Volumes/GoogleDrive/My Drive/CODE/Analyze_IC')

ses2agg = {'sub_1171903433','sub_1172968426','sub_1172969394','sub_1174569641', ...
    'sub_1175512783','sub_1176214862','sub_1177693342','sub_1181314060', ...
    'sub_1181585608','sub_1182593231','sub_1183369803','sub_1186544726', ...
    'sub_1189891322','sub_1194090570'};

%%
%shufmeth = {'asis', 'shuf', 'shufall', 'shufallneurons', 'shuftrials', 'shufneurons'};
shufmeth = {'_shuftrials', '_shuf', '_shufall'}; %, '_shufallneurons', '_shufneurons'};

% takes ~50hrs per session...
rng('shuffle')
for ises = 1:numel(ses2agg)
    clearvars -except ses2agg ises shufmeth SVMkernels ii whichSVMkernel recordingmethod svmdesc
    sesclk = tic;

    computeSVM = true;
    optimizeSVM = true;

    % NOTE -- WHEN Z-SCORING FOR NEUROPIXELS DATA need to change NaN values to 0
    % because some units have 0 firing rate for all trials (these units should normally be filtered out...)
    preproc = 'meancenter'; % zscore or meancenter
    %     whichSVMkernel = 'Linear';

    %Nsplits = 10;
    Nshuf = 1;

    % mousedate = 'HS_CamKIIGC6s_53/201204/'; % example session
    mousedate = ses2agg{ises};
    fprintf(strcat('%d  ', mousedate, '\n'), ises)

    pathpp = strcat('S:/OpenScopeData/000248/postprocessed/', mousedate, '/');
    pathsvm = strcat('S:/OpenScopeData/000248/postprocessed/SVM/SVM_', svmdesc, '/', mousedate, '/');

    if ~exist(pathsvm, 'dir')
        mkdir(pathsvm)
    end


    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur'
    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    neuallloc = electrode_location(revmapelecid(unit_peakch+1));
    load([pathpp 'qc_units.mat'])

    % Siegle et al's single unit filter criteria: isi_violations < 0.5 & amplitude_cutoff < 0.1 & presence_ratio > 0.9
    % RELAX amplitude_cutoff criterion
    whichvisarea = 'VISp';
    if strcmp(whichvisarea, 'VISp')
        neu2anal = contains(neuallloc, 'VISp') & ~contains(neuallloc, 'VISpm');
    elseif strcmp(whichvisarea, 'VISpfiltRS')
        neu2anal = contains(neuallloc, 'VISp') & ~contains(neuallloc, 'VISpm') ...
            & unit_wfdur>0.4 & (unit_isi_violations<0.5 & unit_amplitude_cutoff<0.5 & unit_presence_ratio>0.9);
    end
    

    whichICblock = 'ICwcfg1_presentations';
    ICtrialtypes = [0 101 105 106 107 109 110 111 506 511 1105 1109 1201 1299 ...
        1301 1302 1303 1304 1305 1306 1307 1308];

    svmfn = strcat(pathsvm, 'SVM_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat');
    svmmdlfn = strcat(pathsvm, 'SVMmodels_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat');

    load(sprintf('%spostprocessed.mat', pathpp ))
    load(sprintf('%svisresponses.mat', pathpp ))

    whichR = 'spkcnt';
    if size(Rall.(whichICblock),2)~=length(neu2anal)
        error('check neu2anal')
    end
    tempR = Rall.(whichICblock)(:,neu2anal)';
    trialorder = ICtrialtypes( vis.(whichICblock).trialorder + 1);

    Nneurons = nnz(neu2anal);
    numrectrials = size(tempR,2);

    %% discriminability index
    SVM_models = struct();
    SVM_shuf_models = struct();
    SVMout = struct();

    alltrialtypes = unique(trialorder);
    switch svmdesc
        case 'trainICRC_shuf'
            traintrialtypes = [106, 107, 110, 111];
            probetrialtypes = [1105, 1109, 1201, 1299, 0, 101, 105, 109];
            TRAINRATIO = 0.9;
        case 'trainREx_shuf'
            traintrialtypes = [1201, 1299];
            probetrialtypes = [106, 107, 110, 111];
            TRAINRATIO = 0.9;
        case 'trainIC1RC1_shuf'
            traintrialtypes = [106, 107];
            probetrialtypes = [1105];
            TRAINRATIO = 0.9;
        case 'trainIC2RC2_shuf'
            traintrialtypes = [111 110];
            probetrialtypes = [1109];
            TRAINRATIO = 0.9;
        case 'trainind1_shuf'
            traintrialtypes = [1301 1305];
            probetrialtypes = [1303 1307];
            TRAINRATIO = 0.8;
        case 'trainind2_shuf'
            traintrialtypes = [1302 1306];
            probetrialtypes = [1304 1308];
            TRAINRATIO = 0.8;
        case 'trainind3_shuf'
            traintrialtypes = [1303 1307];
            probetrialtypes = [1301 1305];
            TRAINRATIO = 0.8;
        case 'trainind4_shuf'
            traintrialtypes = [1304 1308];
            probetrialtypes = [1302 1306];
            TRAINRATIO = 0.8;
        otherwise
            error('svmdesc %s not recognized', svmdesc)
    end
    Nsplits = floor(1/(1-TRAINRATIO));

    Ntt = numel(traintrialtypes);
    Nprobett = numel(probetrialtypes);
    Nalltt = numel(alltrialtypes);

    SVMout.Nneurons = Nneurons;
    SVMout.exptid = whichICblock;
    SVMout.ICtrialtypes = ICtrialtypes;
    SVMout.trialtypes = traintrialtypes;

    SVMout.numtrials = zeros(Ntt,1);
    % DI.numtrialpairs = zeros(Ntt,Ntt);
    for typi1 = 1:Ntt
        SVMout.numtrials(typi1) = nnz(trialorder==SVMout.trialtypes(typi1));
    end


    % balance trials
    cftttrials = ismember(trialorder, SVMout.trialtypes);
    Ntrialspertype = min(SVMout.numtrials);
    if all(SVMout.numtrials==Ntrialspertype)
        trials2anal = cftttrials;
        SVMout.analtrials = find(trials2anal);
    else
        warning('balancing number of trials')
        trials2anal = false(numrectrials,1);
        for typi1 = 1:Ntt
            trialsintype = find(trialorder==SVMout.trialtypes(typi1));
            trialsintype = trialsintype(1:Ntrialspertype);
            trials2anal(trialsintype) = true;
        end
        if all(cftttrials(trials2anal)) && ~any(trials2anal(~cftttrials))
            SVMout.analtrials = find(trials2anal);
        else
            error('trials to analyze was not selected correctly')
        end
    end
    SVMout.analtriallabels = trialorder(trials2anal);

    Ntraintrialspertype = round(TRAINRATIO*Ntrialspertype);
    Ntesttrialspertype = Ntrialspertype - Ntraintrialspertype;
    SVMout.Ntt = Ntt;
    SVMout.Ntrialspertype = Ntrialspertype;
    SVMout.Ntraintrialspertype = Ntraintrialspertype;

    Ntraintrials = Ntt*Ntraintrialspertype;
    Ntesttrials = Ntt*(Ntrialspertype-Ntraintrialspertype);

    % probe trials
    probetrials = ismember(trialorder, probetrialtypes);
    SVMout.probetrials = find(probetrials);

    alltrials = true(size(trialorder));
    SVMout.alltrials = find(alltrials);

    SVMout.trialorder = trialorder;

    randtrialorder=randperm(numrectrials);
    SVMout.randtrialorder = randtrialorder;

    SVM_models.(whichR) = cell(1, Nsplits);

    SVMout.traintrialinds = zeros(Ntraintrials, Nsplits);
    SVMout.testtrialinds = zeros(Ntesttrials, Nsplits);

    SVMout.Ylabs = cell(Ntt, Nsplits);

    for imeth = 1:numel(shufmeth)
        SVM_shuf_models.(strcat(whichR, shufmeth{imeth})) = cell(1, Nsplits);
        for ts = 1:16
            switch ts
                case 1
                    svmmd = 'train';
                    tempNtrials = Ntraintrials;
                case 2
                    svmmd = 'train_shuf';
                    tempNtrials = Ntraintrials;
                case 3
                    svmmd = 'train_msdn';
                    tempNtrials = Ntraintrials;
                case 4
                    svmmd = 'train_mnds';
                    tempNtrials = Ntraintrials;
                case 5
                    svmmd = 'test';
                    tempNtrials = Ntesttrials;
                case 6
                    svmmd = 'test_shuf';
                    tempNtrials = Ntesttrials;
                case 7
                    svmmd = 'test_msdn';
                    tempNtrials = Ntesttrials;
                case 8
                    svmmd = 'test_mnds';
                    tempNtrials = Ntesttrials;
                case 9
                    svmmd = 'probe';
                    tempNtrials = nnz(probetrials);
                case 10
                    svmmd = 'probe_shuf';
                    tempNtrials = nnz(probetrials);
                case 11
                    svmmd = 'probe_msdn';
                    tempNtrials = nnz(probetrials);
                case 12
                    svmmd = 'probe_mnds';
                    tempNtrials = nnz(probetrials);
                case 13
                    svmmd = 'all';
                    tempNtrials = nnz(alltrials);
                case 14
                    svmmd = 'all_shuf';
                    tempNtrials = nnz(alltrials);
                case 15
                    svmmd = 'all_msdn';
                    tempNtrials = nnz(alltrials);
                case 16
                    svmmd = 'all_mnds';
                    tempNtrials = nnz(alltrials);
            end
            SVMout.(strcat(whichR, shufmeth{imeth})).(svmmd).label = cell(tempNtrials, Nsplits);
            SVMout.(strcat(whichR, shufmeth{imeth})).(svmmd).score = NaN(tempNtrials,Ntt, Nsplits);
        end
    end


    switch preproc
        case 'none'
            Tp = tempR';
        case 'zscore'
            % Z-score
            tempRmean = mean(tempR,2);
            tempRstd = std(tempR,0,2);

            Tp = ( (tempR-tempRmean)./tempRstd )';
            Tp(isnan(Tp))=0;
        case 'minmax'
            tempRmin = min(tempR,[],2);
            tempRrange = range(tempR,2);

            Tp = ( (tempR-tempRmin)./tempRrange )';
        case 'meancenter'
            tempRmean = mean(tempR,2);
            Tp = ( tempR-tempRmean )';
    end
    if size(Tp,2) ~= Nneurons
        error('Tp size not expected')
    end

    Tpttmean = zeros(size(Tp));
    for typi1 = 1:Nalltt
        trialsintype = find(trialorder==alltrialtypes(typi1));
        Tpttmean(trialsintype,:) = repmat(mean(Tp(trialsintype,:), 1), numel(trialsintype),1);
    end
    Tpres = Tp-Tpttmean;

    % shuf: for each trial type, shuffle responses across trials
    Tpshuf = Tp;
    for typi1 = 1:Nalltt
        trialsintype = find(trialorder==alltrialtypes(typi1));
        for ci = 1:Nneurons
            tempshuforder = trialsintype(randperm(numel(trialsintype)));
            Tpshuf(trialsintype,ci) = Tpshuf(tempshuforder,ci);
        end
    end

    % shuftrials: shuffle residuals across trials same way
    % across all neurons
    tempshuforder = randperm(size(Tp,1));
    Tpresshuftrials = Tpres(tempshuforder,:);
    Tpshuftrials = Tpresshuftrials + Tpttmean;

    % shufall: shuffle residuals across all trials
    Tpresshufall = Tpres;
    for ci = 1:Nneurons
        tempshuforder = randperm(size(Tp,1));
        Tpresshufall(:,ci) = Tpres(tempshuforder,ci);
    end
    Tpshufall = Tpresshufall + Tpttmean;

    % shufneurons: shuffle residuals across neurons same way
    % across all trials
    tempshuforder = randperm(Nneurons);
    Tpresshufneurons = Tpres(:,tempshuforder);
    Tpshufneurons = Tpresshufneurons + Tpttmean;

    % shufallneurons: shuffle residuals across all neurons
    Tpresshufallneurons = Tpres;
    for itrial = 1:size(Tp,1)
        tempshuforder = randperm(Nneurons);
        Tpresshufallneurons(itrial,:) = Tpres(itrial,tempshuforder);
    end
    Tpshufallneurons = Tpresshufallneurons + Tpttmean;



    % takes 20 min per trial type pair. 2000 min per session (33 hr)
    trials2anal = randtrialorder(ismember(randtrialorder, SVMout.analtrials));
    for isplit = 1:Nsplits
        close all
        ttclk = tic;

        % trials2anal = randtrialorder(ismember(randtrialorder, SVMtestRE.analtrials));
        testtrialinds = zeros(Ntesttrials,1);
        traintrialinds = zeros(Ntraintrials,1);
        for typi1 = 1:Ntt
            tempinds = trials2anal( trialorder(trials2anal)==SVMout.trialtypes(typi1) );
            tempinds = reshape(tempinds,[],1);
            if size(tempinds,1) ~= Ntrialspertype
                error('Ntrialspertype not consistent between trial types? check')
            end
            temptestintype = false(Ntrialspertype,1);
            temptestintype((isplit-1)*Ntesttrialspertype+1:isplit*Ntesttrialspertype) = true;
            temptrainintype = true(Ntrialspertype,1);
            temptrainintype((isplit-1)*Ntesttrialspertype+1:isplit*Ntesttrialspertype) = false;
            testtrialinds((typi1-1)*Ntesttrialspertype+1:typi1*Ntesttrialspertype) = tempinds(temptestintype);
            traintrialinds((typi1-1)*Ntraintrialspertype+1:typi1*Ntraintrialspertype) = tempinds(temptrainintype);
        end
        testtrialinds = trials2anal(ismember(trials2anal, testtrialinds));
        traintrialinds = trials2anal(ismember(trials2anal, traintrialinds));

        if any(ismember(traintrialinds, testtrialinds))
            error('train and test trials should not overlap')
        end
        if ~( all(ismember(trialorder(testtrialinds), SVMout.trialtypes)) && all(ismember(trialorder(traintrialinds), SVMout.trialtypes)) )
            error('train and test trials of incorrect type detected')
        end

        SVMout.traintrialinds(:,isplit) = traintrialinds;
        SVMout.testtrialinds(:, isplit) = testtrialinds;


        X = Tp(traintrialinds,:);
        Y = strsplit(sprintf('%d\n',trialorder(traintrialinds)), '\n')';
        Y = Y(1:end-1);

        Xtest = Tp(testtrialinds,:);
        Ytest = strsplit(sprintf('%d\n',trialorder(testtrialinds)), '\n')';
        Ytest = Ytest(1:end-1);

        Xprobe = Tp(probetrials,:);
        Xall = Tp(alltrials,:);

        % t is an SVM template. Most of its properties are empty.
        % When the software trains the ECOC classifier, it sets the applicable properties to their default values.
        % Train the ECOC classifier using the SVM template.
        % Transform classification scores to class posterior probabilities
        % (which are returned by predict or resubPredict) using the 'FitPosterior' name-value pair argument.
        % Specify the class order using the 'ClassNames' name-value pair argument.
        % Display diagnostic messages during training by using the 'Verbose' name-value pair argument.

        Ylabs = unique(Y);
        tempYlaborder = randperm(numel(Ylabs));
        [~,reverseYlaborder]=sort(tempYlaborder);
        Ylabs = Ylabs(tempYlaborder);
        SVMout.Ylabs(:,isplit) = Ylabs;

        switch whichSVMkernel
            case 'RBF'
                t = templateSVM('Standardize',true,'KernelFunction', 'rbf');
            case 'Linear'
                t = templateSVM('Standardize',true,'KernelFunction', 'linear');
            case 'Poly2'
                t = templateSVM('Standardize',true,'KernelFunction', 'polynomial' , 'PolynomialOrder', 2);
        end
        if optimizeSVM
            SVMModel = fitcecoc(X,Y,'Learners',t,'FitPosterior',false, ...
                'ClassNames', Ylabs, 'Verbose',0, 'OptimizeHyperparameters', 'auto', ...
                'HyperparameterOptimizationOptions', struct('UseParallel',true, 'ShowPlots', false));
        else
            SVMModel = fitcecoc(X,Y,'Learners',t,'FitPosterior',false, 'ClassNames', Ylabs, 'Verbose',0);
        end
        %                 CVMdl = crossval(SVMModel);


        SVM_models.(whichR){isplit} = SVMModel;

        %                     tempHR_shuf_train = NaN(1,Nshuf);
        %                     tempDIsq_shuf_train = NaN(1,Nshuf);
        %                     tempHR_shuf_test = NaN(1,Nshuf);
        %                     tempDIsq_shuf_test = NaN(1,Nshuf);
        %                     for ishuf = 1:Nshuf

        for imeth = 1:numel(shufmeth)
            switch shufmeth{imeth}
                case '_shuftrials'
                    tempTp = Tpshuftrials;
                case '_shuf'
                    tempTp = Tpshuf;
                case '_shufall'
                    tempTp = Tpshufall;
                case '_shufneurons'
                    tempTp = Tpshufneurons;
                case '_shufallneurons'
                    tempTp = Tpshufallneurons;
            end

            %                         trainRshufmean = mean(tempRshuf(:, traintrialsshuf ),2);
            %                         trainRshufstd = std(tempRshuf(:, traintrialsshuf ),0,2);
            %                         tempRshufsubz = (tempRshuf-trainRshufmean)./trainRshufstd;

            Xshuf = tempTp(traintrialinds, :);
            Xshuftest = tempTp(testtrialinds, :);
            Xshufprobe = tempTp(probetrials,:);
            Xshufall = tempTp(alltrials,:);

            switch whichSVMkernel
                case 'RBF'
                    t = templateSVM('Standardize',true,'KernelFunction', 'rbf');
                case 'Linear'
                    t = templateSVM('Standardize',true,'KernelFunction', 'linear');
                case 'Poly2'
                    t = templateSVM('Standardize',true,'KernelFunction', 'polynomial' , 'PolynomialOrder', 2);
            end
            if optimizeSVM
                SVMModelshuf = fitcecoc(Xshuf,Y, 'Learners',t,'FitPosterior',false, ...
                    'ClassNames', Ylabs, 'Verbose',0, 'OptimizeHyperparameters', 'auto', ...
                    'HyperparameterOptimizationOptions', struct('UseParallel',true, 'ShowPlots', false));
            else
                SVMModelshuf = fitcecoc(Xshuf,Y, 'Learners',t,'FitPosterior',false, 'ClassNames', Ylabs, 'Verbose',0);
            end


            SVM_shuf_models.(strcat(whichR, shufmeth{imeth})){isplit} = SVMModelshuf;

            for t = 1:16
                switch t
                    case 1
                        Xtemp = X;
                        tempSVMmodel = SVMModel;
                        svmmd = 'train';
                    case 2
                        Xtemp = Xshuf;
                        tempSVMmodel = SVMModelshuf;
                        svmmd = 'train_shuf';
                    case 3
                        Xtemp = X;
                        tempSVMmodel = SVMModelshuf;
                        svmmd = 'train_msdn';
                    case 4
                        Xtemp = Xshuf;
                        tempSVMmodel = SVMModel;
                        svmmd = 'train_mnds';
                    case 5
                        Xtemp = Xtest;
                        tempSVMmodel = SVMModel;
                        svmmd = 'test';
                    case 6
                        Xtemp = Xshuftest;
                        tempSVMmodel = SVMModelshuf;
                        svmmd = 'test_shuf';
                    case 7
                        Xtemp = Xtest;
                        tempSVMmodel = SVMModelshuf;
                        svmmd = 'test_msdn';
                    case 8
                        Xtemp = Xshuftest;
                        tempSVMmodel = SVMModel;
                        svmmd = 'test_mnds';
                    case 9
                        Xtemp = Xprobe;
                        tempSVMmodel = SVMModel;
                        svmmd = 'probe';
                    case 10
                        Xtemp = Xshufprobe;
                        tempSVMmodel = SVMModelshuf;
                        svmmd = 'probe_shuf';
                    case 11
                        Xtemp = Xprobe;
                        tempSVMmodel = SVMModelshuf;
                        svmmd = 'probe_msdn';
                    case 12
                        Xtemp = Xshufprobe;
                        tempSVMmodel = SVMModel;
                        svmmd = 'probe_mnds';
                    case 13
                        Xtemp = Xall;
                        tempSVMmodel = SVMModel;
                        svmmd = 'all';
                    case 14
                        Xtemp = Xshufall;
                        tempSVMmodel = SVMModelshuf;
                        svmmd = 'all_shuf';
                    case 15
                        Xtemp = Xall;
                        tempSVMmodel = SVMModelshuf;
                        svmmd = 'all_msdn';
                    case 16
                        Xtemp = Xshufall;
                        tempSVMmodel = SVMModel;
                        svmmd = 'all_mnds';
                end
                [templabel,tempscore] = predict(tempSVMmodel,Xtemp);
                SVMout.(strcat(whichR, shufmeth{imeth})).(svmmd).label(:,isplit) = templabel;
                SVMout.(strcat(whichR, shufmeth{imeth})).(svmmd).score(:,:,isplit) = tempscore(:, reverseYlaborder);
            end
        end
        fprintf('%s %s %d/%d\n', mousedate, whichSVMkernel, isplit, Nsplits)
        toc(ttclk)
    end

    switch svmdesc
        case 'trainICRC_shuf'
            SVMtrainICRC = SVMout;
            SVMtrainICRC_models = SVM_models;
            SVMtrainICRC_shuf_models = SVM_shuf_models;
            save(svmfn, 'preproc', 'whichSVMkernel', 'SVMtrainICRC', '-v7.3')
            save(svmmdlfn, 'preproc', 'whichSVMkernel', 'SVMtrainICRC_models', 'SVMtrainICRC_shuf_models', '-v7.3')
        case 'trainREx_shuf'
            SVMtrainREx = SVMout;
            SVMtrainREx_models = SVM_models;
            SVMtrainREx_shuf_models = SVM_shuf_models;
            save(svmfn, 'preproc', 'whichSVMkernel', 'SVMtrainREx', '-v7.3')
            save(svmmdlfn, 'preproc', 'whichSVMkernel', 'SVMtrainREx_models', 'SVMtrainREx_shuf_models', '-v7.3')
        case 'trainIC1RC1_shuf'
            SVMtrainIC1RC1 = SVMout;
            SVMtrainIC1RC1_models = SVM_models;
            SVMtrainIC1RC1_shuf_models = SVM_shuf_models;
            save(svmfn, 'preproc', 'whichSVMkernel', 'SVMtrainIC1RC1', '-v7.3')
            save(svmmdlfn, 'preproc', 'whichSVMkernel', 'SVMtrainIC1RC1_models', 'SVMtrainIC1RC1_shuf_models', '-v7.3')
        case 'trainIC2RC2_shuf'
            SVMtrainIC2RC2 = SVMout;
            SVMtrainIC2RC2_models = SVM_models;
            SVMtrainIC2RC2_shuf_models = SVM_shuf_models;
            save(svmfn, 'preproc', 'whichSVMkernel', 'SVMtrainIC2RC2', '-v7.3')
            save(svmmdlfn, 'preproc', 'whichSVMkernel', 'SVMtrainIC2RC2_models', 'SVMtrainIC2RC2_shuf_models', '-v7.3')
        case 'trainind1_shuf'
            SVMtrainind1 = SVMout;
            SVMtrainind1_models = SVM_models;
            SVMtrainind1_shuf_models = SVM_shuf_models;
            save(svmfn, 'preproc', 'whichSVMkernel', 'SVMtrainind1', '-v7.3')
            save(svmmdlfn, 'preproc', 'whichSVMkernel', 'SVMtrainind1_models', 'SVMtrainind1_shuf_models', '-v7.3')
        case 'trainind2_shuf'
            SVMtrainind2 = SVMout;
            SVMtrainind2_models = SVM_models;
            SVMtrainind2_shuf_models = SVM_shuf_models;
            save(svmfn, 'preproc', 'whichSVMkernel', 'SVMtrainind2', '-v7.3')
            save(svmmdlfn, 'preproc', 'whichSVMkernel', 'SVMtrainind2_models', 'SVMtrainind2_shuf_models', '-v7.3')
        case 'trainind3_shuf'
            SVMtrainind3 = SVMout;
            SVMtrainind3_models = SVM_models;
            SVMtrainind3_shuf_models = SVM_shuf_models;
            save(svmfn, 'preproc', 'whichSVMkernel', 'SVMtrainind3', '-v7.3')
            save(svmmdlfn, 'preproc', 'whichSVMkernel', 'SVMtrainind3_models', 'SVMtrainind3_shuf_models', '-v7.3')
        case 'trainind4_shuf'
            SVMtrainind4 = SVMout;
            SVMtrainind4_models = SVM_models;
            SVMtrainind4_shuf_models = SVM_shuf_models;
            save(svmfn, 'preproc', 'whichSVMkernel', 'SVMtrainind4', '-v7.3')
            save(svmmdlfn, 'preproc', 'whichSVMkernel', 'SVMtrainind4_models', 'SVMtrainind4_shuf_models', '-v7.3')
        otherwise
            error('svmdesc %s not recognized', svmdesc)
    end

    %end
    disp([mousedate ' ' svmdesc ' ' whichSVMkernel ' done!'])
    toc(sesclk)
end

disp(['analyzeICtx_SVM_shufmeth ' svmdesc ' ' whichSVMkernel ' FINISHED! READY TO MOVE DATA'])

