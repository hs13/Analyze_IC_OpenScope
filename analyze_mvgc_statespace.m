%% code taken from mvgc_demo_statespace.m from MVGCtoolbox
% also refer to https://users.sussex.ac.uk/~lionelb/MVGC/html/mvgchelp.html
% format X as Nneurons X Ntimesteps X Ntrials
addpath(genpath('C:\Users\USER\GitHub\MVGCtoolbox'))

datadir = 'S:/OpenScopeData/000248/';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions(~contains(nwbsessions, 'Placeholder') & ...
    ( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') ));
Nsessions = numel(nwbsessions);

svmdesc = 'trainICRCtestRE';
preproc = 'zscore'; % '' is z-score train trials, '_zscoreall', or '_meancenter'
whichSVMkernel = 'Linear';

whichICblock = 'ICwcfg1';
probes = {'A', 'B', 'C', 'D', 'E', 'F'};
visareas = {'AM', 'PM', 'V1', 'LM', 'AL', 'RL'};
ICtrialtypes = [0 101 105 106 107 109 110 111 506 511 1105 1109 1201 1299 ...
    1301 1302 1303 1304 1305 1306 1307 1308];

Twin = 100;
ises=2; 
pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
pathsvm = [datadir 'postprocessed' filesep 'SVM' filesep 'SVM_' svmdesc filesep nwbsessions{ises} filesep];

load([pathpp 'visresponses_probeC.mat'], 'ICsig')
load([pathsvm, 'SVMbeta_', svmdesc, '_Cctx_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'], 'neuctxinprobe')
load([pathsvm, 'SVMpsth' num2str(Twin) 'ms_', svmdesc, '_Cctx_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'])
if ~(length(neuctxinprobe)==length(ICsig.([whichICblock '_presentations']).ICresp1) && nnz(neuctxinprobe)==size(SVMpsth.psthbin,3))
    error('check neuctxinprobe/ICsig/SVMpsth')
end
ICsigC = ICsig;
neuinprobeCctx = neuctxinprobe;
psthbinCctx = SVMpsth.psthbin;
invalidCctx = reshape(all(psthbinCctx==0, [1,2]), [],1);

load([pathpp 'visresponses_probeD.mat'], 'ICsig')
load([pathsvm, 'SVMbeta_', svmdesc, '_Dctx_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'], 'neuctxinprobe')
load([pathsvm, 'SVMpsth' num2str(Twin) 'ms_', svmdesc, '_Dctx_', whichSVMkernel, '_', preproc, '_', whichICblock, '.mat'])
if ~(length(neuctxinprobe)==length(ICsig.([whichICblock '_presentations']).ICresp1) && nnz(neuctxinprobe)==size(SVMpsth.psthbin,3))
    error('check neuctxinprobe/ICsig/SVMpsth')
end
ICsigD = ICsig;
neuinprobeDctx = neuctxinprobe;
psthbinDctx = SVMpsth.psthbin;
invalidDctx = reshape(all(psthbinDctx==0, [1,2]), [],1);

psthbinCDctx = cat(3,psthbinCctx, psthbinDctx);
psthbinCDctx = permute(psthbinCDctx, [3 1 2]);

psthtli = SVMpsth.psthtli;
tempTinds = SVMpsth.psthbinTinds( round(1+size(SVMpsth.psthbinTinds,1)/2),:);
psthbintli = psthtli(tempTinds);

neuoiCctx = ~invalidCctx & ICsigC.([whichICblock '_presentations']).ICresp1(neuinprobeCctx==1);
neuoiDctx = ~invalidDctx & ICsigD.([whichICblock '_presentations']).ICresp1(neuinprobeDctx==1);
neuoi = cat(1,neuoiCctx, neuoiDctx);

trialsoi = SVMpsth.trialorder==106;
tloi = psthbintli>=0 & psthbintli<400;
X = psthbinCDctx(neuoi, tloi, trialsoi);

neugroups = {'ICencoder','RCencoder','inducerencoder','inducerresponsive', ...
    'indenc1','indenc2','indenc3','indenc4','indenc13','indenc14','indenc23','indenc24', ...
    'ICencoder1','RCencoder1','RCencoder2','ICencoder2','indin1','indin2','indin3','indin4', ...
    'indout1','indout2','indout3','indout4','RElfaith1','RElfaith2', ...
    'ICresp1','RCresp1','RCresp2','ICresp2','ICtuned1','RCtuned1','RCtuned2','ICtuned2'};
neuCDctxind = struct();
neuXind = struct();
for g = 1:numel(neugroups)
    neuCDctxind.(neugroups{g}) = cat(1, ICsigC.([whichICblock '_presentations']).(neugroups{g})(neuinprobeCctx==1), ICsigD.([whichICblock '_presentations']).(neugroups{g})(neuinprobeDctx==1));
    neuXind.(neugroups{g}) = neuCDctxind.(neugroups{g})(neuoi);
end
neuCDctxind.Cctx = false(length(neuoi),1); neuCDctxind.Cctx(1:length(neuoiCctx)) = true;
neuCDctxind.Dctx = false(length(neuoi),1); neuCDctxind.Dctx(length(neuoiCctx)+1:length(neuoiCctx)+length(neuoiDctx)) = true;
neuXind.Cctx = false(nnz(neuoi),1); neuXind.Cctx(1:nnz(neuoiCctx)) = true;
neuXind.Dctx = false(nnz(neuoi),1); neuXind.Dctx(nnz(neuoiCctx)+1:nnz(neuoiCctx)+nnz(neuoiDctx)) = true;

all(neuXind.ICresp1)
isequal(neuCDctxind.Cctx, ~neuCDctxind.Dctx)
isequal(neuXind.Cctx, ~neuXind.Dctx)

%% Parameters

regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 30;     % maximum model order for model order estimation

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = 'F';    % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDRD'; % multiple hypothesis test correction (see routine 'significance')

fs        = 1000/mean(diff(psthbintli));    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

seed      = 0;      % random seed (0 for unseeded)

%% Model order estimation (<mvgc_schema.html#3 |A2|>)
% Matrix must be positive definite error with 'chol' 

% Calculate information criteria up to specified maximum model order.

ptic('\n*** tsdata_to_infocrit\n');
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
ptoc('*** tsdata_to_infocrit took ');

% Plot information criteria.

figure(1); clf;
plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
title('Model order estimation');

amo = size(AT,3); % actual model order

fprintf('\nbest model order (AIC) = %d\n',moAIC);
fprintf('best model order (BIC) = %d\n',moBIC);
fprintf('actual model order     = %d\n',amo);

% Select model order.

if     strcmpi(morder,'actual')
    morder = amo;
    fprintf('\nusing actual model order = %d\n',morder);
elseif strcmpi(morder,'AIC')
    morder = moAIC;
    fprintf('\nusing AIC best model order = %d\n',morder);
elseif strcmpi(morder,'BIC')
    morder = moBIC;
    fprintf('\nusing BIC best model order = %d\n',morder);
else
    fprintf('\nusing specified model order = %d\n',morder);
end

%% VAR model estimation (<mvgc_schema.html#3 |A2|>)

% Estimate VAR model of selected order from data.

ptic('\n*** tsdata_to_var... ');
[A,SIG] = tsdata_to_var(X,morder,regmode);
assert(~isbad(A),'VAR estimation failed - bailing out');
ptoc;

% Report information on the estimated VAR, and check for errors.
%
% _IMPORTANT:_ We check the VAR model for stability and symmetric
% positive-definite residuals covariance matrix. _THIS CHECK SHOULD ALWAYS BE
% PERFORMED!_ - subsequent routines may fail if there are errors here. If there
% are problems with the data (e.g. non-stationarity, colinearity, etc.) there's
% also a good chance they'll show up at this point - and the diagnostics may
% supply useful information as to what went wrong.

info = var_info(A,SIG);
assert(~info.error,'VAR error(s) found - bailing out');

%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Calculate time-domain pairwise-conditional causalities from VAR model parameters
% by state-space method [4]. The VAR model is transformed into an equivalent state-
% space model for computation. Also return p-values for specified test (F-test or
% likelihood-ratio test; this is optional - if p-values are not required, then it
% is not necessary to supply time series |X|, regression mode |regmode|, or test
% specification |tstat|).

ptic('*** var_to_pwcgc... ');
[F,pval] = var_to_pwcgc(A,SIG,X,regmode,tstat);
ptoc;

% Check for failed GC calculation

assert(~isbad(F,false),'GC calculation failed - bailing out');

% Significance-test p-values, correcting for multiple hypotheses.

sig = significance(pval,alpha,mhtc);

% Plot time-domain causal graph, p-values and significance.

figure(2); clf;
sgtitlex('Pairwise-conditional Granger causality - time domain');
subplot(1,3,1);
plot_pw(F);
title('Pairwise-conditional GC');
subplot(1,3,2);
plot_pw(pval);
title(['p-values (' tstat '-test)']);
subplot(1,3,3);
plot_pw(sig);
title(['Significant at \alpha = ' num2str(alpha)]);


%% F/pval/sig rows "to", columns "from
[r,c]=find(sig==1);
figure
    hold all
plot(psthbintli(tloi), squeeze(mean(X(c, :, :),[1,3])), 'b-', 'linewidth',1)
plot(psthbintli(tloi), squeeze(mean(X(r, :, :),[1,3])), 'r-', 'linewidth',2)

randord = randperm(length(r));
figure
for ii = 1:5*8
    tocell = r(randord(ii));
    fromcell = c(randord(ii));
    subplot(5,8,ii)
    hold all
    yyaxis left
plot(psthbintli(tloi), squeeze(mean(X(fromcell, :, :),3)), 'b-', 'linewidth',1)
    yyaxis right
plot(psthbintli(tloi), squeeze(mean(X(tocell, :, :),3)), 'r-', 'linewidth',2)
end

%%
neuCctxinds = 1:nnz(neuoiCctx);
neuDctxinds = nnz(neuoiCctx)+1:nnz(neuoiCctx)+nnz(neuoiDctx);

figure; hold all
histogram(nanmean(sig(neuCctxinds,neuCctxinds),2))
histogram(nanmean(sig(neuDctxinds,neuCctxinds),2)) % C->D
histogram(nanmean(sig(neuCctxinds,neuDctxinds),2)) % D->C
histogram(nanmean(sig(neuDctxinds,neuDctxinds),2))

nanmean(sig(neuCctxinds,neuCctxinds),[1,2])
nanmean(sig(neuDctxinds,neuCctxinds),[1,2]) % C->D
nanmean(sig(neuCctxinds,neuDctxinds),[1,2]) % D->C
nanmean(sig(neuDctxinds,neuDctxinds),[1,2])


% proportion of inputs from LM
figure; hold all
hc = histogram(nanmean(sig(neuXind.Cctx, neuXind.Dctx),2));
histogram(nanmean(sig( (neuXind.indin1|neuXind.indin3) & neuXind.Cctx, neuXind.Dctx),2), 'BinEdges', hc.BinEdges)
histogram(nanmean(sig(neuXind.ICencoder1 & neuXind.Cctx, neuXind.Dctx),2), 'BinEdges', hc.BinEdges)

% proportion of outputs to LM
figure; hold all
hc = histogram(nanmean(sig(neuXind.Dctx, neuXind.Cctx),1));
histogram(nanmean(sig(neuXind.Dctx, (neuXind.indin1|neuXind.indin3) & neuXind.Cctx),1), 'BinEdges', hc.BinEdges)
histogram(nanmean(sig(neuXind.Dctx, neuXind.ICencoder1 & neuXind.Cctx),1), 'BinEdges', hc.BinEdges)

% proportion of inputs from within V1
figure; hold all
hc = histogram(nanmean(sig(neuXind.Cctx, neuXind.Cctx),2));
histogram(nanmean(sig( (neuXind.indin1|neuXind.indin3) & neuXind.Cctx, neuXind.Cctx),2), 'BinEdges', hc.BinEdges)
histogram(nanmean(sig(neuXind.ICencoder1 & neuXind.Cctx, neuXind.Cctx),2), 'BinEdges', hc.BinEdges)

% proportion of outputs to within V1
figure; hold all
hc = histogram(nanmean(sig(neuXind.Cctx, neuXind.Cctx),1));
histogram(nanmean(sig(neuXind.Cctx, (neuXind.indin1|neuXind.indin3) & neuXind.Cctx),1), 'BinEdges', hc.BinEdges)
histogram(nanmean(sig(neuXind.Cctx, neuXind.ICencoder1 & neuXind.Cctx),1), 'BinEdges', hc.BinEdges)

figure; hold all
hc = histogram(nanmean(sig(neuXind.Cctx, neuXind.Cctx),1));
histogram(nanmean(sig(neuXind.ICencoder1 & neuXind.Cctx, (neuXind.indin1|neuXind.indin3) & neuXind.Cctx),1), 'BinEdges', hc.BinEdges)
histogram(nanmean(sig((neuXind.indin1|neuXind.indin3) & neuXind.Cctx, neuXind.ICencoder1 & neuXind.Cctx),1), 'BinEdges', hc.BinEdges)

figure; hold all
histogram(nanmean(sig(neuXind.ICencoder1 & neuXind.Cctx, (neuXind.indin1|neuXind.indin3) & neuXind.Cctx),1))
histogram(nanmean(sig((neuXind.indin1|neuXind.indin3) & neuXind.Cctx, neuXind.ICencoder1 & neuXind.Cctx),1))

figure; 
for isp = 1:3
subplot(1,3,isp)
switch isp
    case 1
imshow( sig )
    case 2
imshow( sig(neuXind.ICencoder1 & neuXind.Cctx, :) )
    case 3
imshow( sig(:,neuXind.ICencoder1 & neuXind.Cctx) )
end
colormap(flipud(gray))
xlabel('from')
ylabel('to')
end

nanmean(sig(neuXind.ICencoder1 & neuXind.Cctx,neuXind.ICencoder1 & neuXind.Cctx), 'all')
nanmean(sig(neuXind.Cctx,neuXind.Cctx), 'all')
ranksum( reshape(pval(neuXind.ICencoder1 & neuXind.Cctx,neuXind.ICencoder1 & neuXind.Cctx),[],1), reshape(pval(neuXind.Cctx,neuXind.Cctx),[],1))

figure; imagesc( sig(neuXind.ICencoder1 & neuXind.Cctx,:) )
colormap gray

ranksum(nanmean(sig(neuXind.ICencoder1 & neuXind.Cctx, neuXind.Dctx),2), nanmean(sig(~neuXind.ICencoder1 & neuXind.Cctx, neuXind.Dctx),2))

ranksum(nanmean(sig( (neuXind.indin1|neuXind.indin3) & neuXind.Cctx, neuXind.Dctx),2), nanmean(sig(neuXind.ICencoder1 & neuXind.Cctx, neuXind.Dctx),2))

ranksum(nanmean(sig( (neuXind.indin1|neuXind.indin3) & neuXind.Cctx, neuXind.Dctx),2), nanmean(sig( ~(neuXind.indin1|neuXind.indin3) & neuXind.Cctx, neuXind.Dctx),2))


figure; hold all
hc = histogram(ICsigC.ICwcfg1_presentations.SP_ICvsRC  );
histogram(ICsigC.ICwcfg1_presentations.SP_ICvsRC(ICsigC.ICwcfg1_presentations.ICencoder), 'BinEdges', hc.BinEdges )
ranksum(ICsigC.ICwcfg1_presentations.SP_ICvsRC(ICsigC.ICwcfg1_presentations.ICencoder), ICsigC.ICwcfg1_presentations.SP_ICvsRC)


