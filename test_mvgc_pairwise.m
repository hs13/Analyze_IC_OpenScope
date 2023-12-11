%% TODO: GRANGER CAUSALITY BETWEEN NEURONS' PSTH% AND SVM PSTH

% test whether mvgc is calculated pairwise...

X_Cctx = X(neuXind.Cctx,:,:);
%% Model order estimation (<mvgc_schema.html#3 |A2|>)
% Matrix must be positive definite error with 'chol' 

% Calculate information criteria up to specified maximum model order.

ptic('\n*** tsdata_to_infocrit\n');
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X_Cctx,momax,icregmode);
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
[A_Cctx,SIG_Cctx] = tsdata_to_var(X_Cctx,morder,regmode);
assert(~isbad(A_Cctx),'VAR estimation failed - bailing out');
ptoc;

% Report information on the estimated VAR, and check for errors.
%
% _IMPORTANT:_ We check the VAR model for stability and symmetric
% positive-definite residuals covariance matrix. _THIS CHECK SHOULD ALWAYS BE
% PERFORMED!_ - subsequent routines may fail if there are errors here. If there
% are problems with the data (e.g. non-stationarity, colinearity, etc.) there's
% also a good chance they'll show up at this point - and the diagnostics may
% supply useful information as to what went wrong.

info = var_info(A_Cctx,SIG_Cctx);
assert(~info.error,'VAR error(s) found - bailing out');

%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Calculate time-domain pairwise-conditional causalities from VAR model parameters
% by state-space method [4]. The VAR model is transformed into an equivalent state-
% space model for computation. Also return p-values for specified test (F-test or
% likelihood-ratio test; this is optional - if p-values are not required, then it
% is not necessary to supply time series |X|, regression mode |regmode|, or test
% specification |tstat|).

ptic('*** var_to_pwcgc... ');
[F_Cctx,pval_Cctx] = var_to_pwcgc(A_Cctx,SIG_Cctx,X_Cctx,regmode,tstat);
ptoc;

% Check for failed GC calculation

assert(~isbad(F_Cctx,false),'GC calculation failed - bailing out');

% Significance-test p-values, correcting for multiple hypotheses.

sig_Cctx = significance(pval_Cctx,alpha,mhtc);

% Plot time-domain causal graph, p-values and significance.

figure(2); clf;
sgtitlex('Pairwise-conditional Granger causality - time domain');
subplot(1,3,1);
plot_pw(F_Cctx);
title('Pairwise-conditional GC');
subplot(1,3,2);
plot_pw(pval_Cctx);
title(['p-values (' tstat '-test)']);
subplot(1,3,3);
plot_pw(sig_Cctx);
title(['Significant at \alpha = ' num2str(alpha)]);

%%

figure; clf;
sgtitlex('Pairwise-conditional Granger causality - time domain');
subplot(2,3,1);
plot_pw(F(neuXind.Cctx,neuXind.Cctx));
title('Pairwise-conditional GC');
subplot(2,3,2);
plot_pw(pval(neuXind.Cctx,neuXind.Cctx));
title(['p-values (' tstat '-test)']);
subplot(2,3,3);
plot_pw(sig(neuXind.Cctx,neuXind.Cctx));
title(['Significant at \alpha = ' num2str(alpha)]);

subplot(2,3,4);
plot_pw(F_Cctx);
title('Pairwise-conditional GC');
subplot(2,3,5);
plot_pw(pval_Cctx);
title(['p-values (' tstat '-test)']);
subplot(2,3,6);
plot_pw(sig_Cctx);
title(['Significant at \alpha = ' num2str(alpha)]);


figure
plot(pval_Cctx, pval(neuXind.Cctx,neuXind.Cctx), '.')
