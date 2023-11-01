datadir = 'S:/OpenScopeData/000248/';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

svmdesc = 'trainRExtestICRC';
preproc = 'zscore'; % '' is z-score train trials, '_zscoreall', or '_meancenter'
whichSVMkernel = 'Linear';

whichICblock = 'ICwcfg1';
whichprobe = 'C';
whichvisarea = [whichprobe 'ctx'];

Nsplits = 10;
%%
ICsig_Cagg = struct();
RFCI_Cagg = struct();
RFCIspin_Cagg = struct();
ori4params_Cagg = struct();
oriparams_Cagg = struct();
sizeCI_Cagg = struct();

neuctxinprobe_Cagg = cell(1,Nsessions);
sesneuctxinprobe_Cagg = cell(1,Nsessions);
SVMcodingname_Cctxagg = cell(Nsplits,Nsessions);
Nlearners_Cctxagg = zeros(Nsplits,Nsessions);
betalearners_Cctxagg = cell(Nsplits,Nsessions);

for ises = 1:Nsessions
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    pathsv = [datadir 'postprocessed' filesep 'SVM' filesep 'SVM_' svmdesc filesep];
    pathsvm = [pathsv nwbsessions{ises} filesep];

    load([pathpp 'visresponses_probe' whichprobe '.mat'])
    if ises==1
        ICsig_Cagg = ICsig;
        RFCI_Cagg = RFCI;
        RFCIspin_Cagg = RFCIspin;
        ori4params_Cagg = ori4params;
        oriparams_Cagg = oriparams;
        sizeCI_Cagg = sizeCI;
    else
        ICsig_Cagg = cat(1, ICsig_Cagg, ICsig);
        RFCI_Cagg = cat(1, RFCI_Cagg, RFCI);
        RFCIspin_Cagg = cat(1, RFCIspin_Cagg, RFCIspin);
        ori4params_Cagg = cat(1, ori4params_Cagg, ori4params);
        oriparams_Cagg = cat(1, oriparams_Cagg, oriparams);
        sizeCI_Cagg = cat(1, sizeCI_Cagg, sizeCI);
    end

    load([pathsvm, 'SVMbeta_', svmdesc, '_', whichvisarea, '_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'])
    neuctxinprobe_Cagg{ises} = neuctxinprobe;
    sesneuctxinprobe_Cagg{ises} = ises*ones(size(neuctxinprobe)) ;
    SVMcodingname_Cctxagg(:,ises) = SVMcodingname;
    Nlearners_Cctxagg(:,ises) = Nlearners;
    betalearners_Cctxagg(:,ises) = betalearners;
end

%% test accuracy

%%
comboonevsone = nchoosek(traintrialtypes, 2);

meanonevsone = cell(1,Nsessions);
meanonevsall = cell(1,Nsessions);
bestonevsone = cell(1,Nsessions); % best test accuracy
bestonevsall = cell(1,Nsessions);

for ises = 1:Nsessions
    temponevsone = strcmp(SVMcodingname_Cctxagg(:,ises), 'onevsone');
    if nnz(temponevsone)>0
        meanonevsone{ises} = mean(cat(3, betalearners_Cctxagg{temponevsone,ises}),3);
    else
        meanonevsone{ises} = NaN( nnz(neuctxinprobe_Cagg{ises}), nchoosek(length(traintrialtypes),2));
    end

    temponevsall = strcmp(SVMcodingname_Cctxagg(:,ises), 'onevsall');
    if nnz(temponevsall)>0
        meanonevsall{ises} = mean(cat(3, betalearners_Cctxagg{temponevsall,ises}),3);
    else
        meanonevsall{ises} = NaN( nnz(neuctxinprobe_Cagg{ises}), length(traintrialtypes) );
    end
end

neuctxinprobeagg = cat(1, neuctxinprobe_Cagg{:});
sesneuctxinprobeagg = cat(1, sesneuctxinprobe_Cagg{:});

ICwcfg1agg = cat(1,ICsig_Cagg.([whichICblock '_presentations'])  );
IC1encagg = cat(1, ICwcfg1agg.ICencoder1);
IC2encagg = cat(1, ICwcfg1agg.ICencoder2);
indin1agg = cat(1, ICwcfg1agg.indin1);
indin2agg = cat(1, ICwcfg1agg.indin2);
indin3agg = cat(1, ICwcfg1agg.indin3);
indin4agg = cat(1, ICwcfg1agg.indin4);

meanonevsoneagg = cat(1,meanonevsone{:});
meanonevsallagg = cat(1,meanonevsall{:});

figure
hold all
histogram(meanonevsoneagg)
histogram(meanonevsoneagg(IC1encagg(neuctxinprobeagg==1)))
histogram(meanonevsoneagg(IC2encagg(neuctxinprobeagg==1)))

p1 = ranksum(meanonevsoneagg(IC1encagg(neuctxinprobeagg==1)), meanonevsoneagg);
p2 = ranksum(meanonevsoneagg(IC2encagg(neuctxinprobeagg==1)), meanonevsoneagg);
fprintf('%d %d : %.4f %.4f\n', comboonevsone(1), comboonevsone(2), p1, p2) % not significant

figure
hold all
histogram(meanonevsallagg)
histogram(meanonevsallagg(IC1encagg(neuctxinprobeagg==1)))
histogram(meanonevsallagg(IC2encagg(neuctxinprobeagg==1)))

p1 = ranksum(meanonevsallagg(IC1encagg(neuctxinprobeagg==1)), meanonevsallagg);
p2 = ranksum(meanonevsallagg(IC2encagg(neuctxinprobeagg==1)), meanonevsallagg);
fprintf('%d %d : %.4f %.4f\n', comboonevsone(1), comboonevsone(2), p1, p2) % not significant
