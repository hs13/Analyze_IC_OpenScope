%% response to IC vs RE (compare latency across areas)
%sizeCI: each stim is 0.25s, inter-trial interval is 0.5s, drifting grating
% baseline period: -100 to 0ms (alternatively -200 to 0ms would work)
%szvec = [0, 4, 8, 16, 32, 64];

kerwinhalf = 2; kersigma = 1;
kerwinhalf = 4; kersigma = 2;
% kerwinhalf = 7; kersigma = 3;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));

whichvisblock = 'ICwcfg1_presentations';
neuopt = 1;
switch neuopt
    case 0
        neutitle = 'non-center-RF neurons';
    case 1
        neutitle = 'center-RF neurons';
    case 2
        neutitle = 'IC-encoders';
end

lw = 3; fs = 16;
[~,probehierarchy]=sort(visind);
% figure %('Position', [100 100 1000 400])
% annotation('textbox', [0.05 0.9 0.8 0.1], 'string', [nwbsessions{ises} ' All Neurons'], 'FontSize', 14, 'EdgeColor', 'none')
for pltopt = 1:2 % 1 real edges, 2 illusory contours
    switch pltopt
        case 1
            figtitle = [neutitle ' Real Edge'];
        case 2
            %         figtitle = [neutitle ' REt'];
            %     case 3
            figtitle = [neutitle ' Illusory Contour'];
    end
    figure%('Position', [0 0 400 400])
    hold all
    for a = 1:2%numel(probes)
        iprobe = probehierarchy(a);
        
        switch neuopt
            case 0
                neuoi = sesneuoiagg{iprobe}<=4 & RFCIagg(iprobe).RFindclassic~=1 & RFCIagg(iprobe).pRFclassic<0.05;
            case 1
                neuoi = sesneuoiagg{iprobe}<=4 & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
            case 2
                neuoi = sesneuoiagg{iprobe}<=4 & ICsigagg.(whichvisblock)(iprobe).ICencoder1==1;
        end
        switch pltopt
            case 1
                ttoi = ismember(ICtrialtypes, [506 511]);
            case 2
                %             ttoi = ismember(ICtrialtypes, [1105 1109]);
                %         case 3
                ttoi = ismember(ICtrialtypes, [106 111]);
        end
        
        temp = psthavg_ctxagg.(whichvisblock){iprobe}(:,ttoi,neuoi);
        % plot(psthtli, smooth(squeeze(mean(temp, [2,3])), 5), '-')
        temp = convn(temp, kergauss, 'same');
        plot(psthtli, squeeze(mean(temp, [2,3])), 'LineWidth', lw)
        %     temppsth = squeeze(mean(temp, 2));
        %     shadedErrorBar(psthtli, squeeze(mean(temppsth,2)), squeeze(std(temppsth,0,2)/sqrt(nnz(neuoi))), {'LineWidth', lw}, 1)
    end
    set(gca, 'FontSize', fs)
    legend(visareas{probehierarchy})
    % xlim([0 400])
    % ylim ([0 12])
    xlim([-50 200])
    set(gca, 'XGrid', 'on')
    title(figtitle)
    % title(sprintf('%s %s neurons N=%d/%d\nsizeCI_presentations', whicharea, neudesc, nnz(neuoi), nnz(neuctx)), 'interpreter', 'none')
    xlabel('Time (ms)')
    ylabel('Rate (Hz)')
end

%% IC vs RE in each area
whichvisblock = 'ICwcfg1_presentations';

kerwinhalf = 2; kersigma = 1;
kerwinhalf = 4; kersigma = 2;
% kerwinhalf = 7; kersigma = 3;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));

fs=10; lw = 1;
figure
for neuopt = 1:2
    switch neuopt
        case 1
            neutitle = 'center-RF neurons';
        case 2
            neutitle = 'non-center-RF neurons';
    end
    
    % figure%('Position', [100 0 1800 1200])
    % annotation('textbox', [0.1 0.9 0.9 0.1], 'string', ['Neuropixels: ' neutitle ' ' whichvisblock], 'edgecolor', 'none', 'fontsize', fs, 'interpreter', 'none')
    for a = 1:2%numel(probes)
        iprobe = probehierarchy(a);
        subplot(2,2,2*(a-1)+neuopt)
        switch neuopt
            case 1
                neuoi = sesneuoiagg{iprobe}<=4 & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
            case 2
                neuoi = sesneuoiagg{iprobe}<=4 & RFCIagg(iprobe).RFindclassic~=1 & RFCIagg(iprobe).pRFclassic<0.05;
        end
        
        hold all
        for typi = 1:2
            switch typi
                case 1
                    ttoi = ismember(ICtrialtypes, [506 511]);
                    ttcol = [0 0 1];
                case 2
                    ttoi = ismember(ICtrialtypes, [106 111]);
                    ttcol = [0 0.7 0];
            end
            
            temppsth = convn(squeeze(psthavg_ctxagg.(whichvisblock){iprobe}(:,ttoi,neuoi)), kergauss, 'same');
            plot(psthtli, squeeze(mean(temppsth,[2 3])), '-', 'Color', ttcol, 'LineWidth', lw)
            %shadedErrorBar(psthtli/1000, squeeze(mean(temppsth,2)), squeeze(std(temppsth,0,2)/sqrt(nnz(neuoi))), {'-', 'Color', ttcol, 'LineWidth', 1}, 1)
        end
        xlim([-50 200])
        ylim([2 15])
        %     ylim([0 25])
        set(gca, 'FontSize', fs, 'XGrid', 'on')
        if a==1 && neuopt==1
            legend({'RE', 'IC'})
        end
        xlabel('Time (s)')
        ylabel('Rate (Hz)')
        title(sprintf('%s %s', visareas{iprobe}, neutitle), 'interpreter', 'none')
        %     title(sprintf('%s %s N=%d', visareas{iprobe}, neutitle, nnz(neuoi)), 'interpreter', 'none')
    end
end


%% ctr-CRF vs non-ctr-CRF, IC-encoder vs inducer-encoder
whichvisblock = 'ICwcfg1_presentations';

kerwinhalf = 2; kersigma = 1;
kerwinhalf = 4; kersigma = 2;
% kerwinhalf = 7; kersigma = 3;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));

fs=16; lw = 3;

a = 1;
iprobe = probehierarchy(a);
figure
hold all
for neuopt = 1:2
    switch neuopt
        case 1
            yyaxis left
            neutitle = 'center-RF neurons';
            neuoi = sesneuoiagg{iprobe}<=4 & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
            neucol = [0.5 0 1];
        case 2
            yyaxis right
            neutitle = 'non-center-RF neurons';
            neuoi = sesneuoiagg{iprobe}<=4 & RFCIagg(iprobe).RFindclassic~=1 & RFCIagg(iprobe).pRFclassic<0.05;
            neucol = [1 0.7 0];
    end
    ttoi = ismember(ICtrialtypes, [106 111]);
    temppsth = convn(squeeze(psthavg_ctxagg.(whichvisblock){iprobe}(:,ttoi,neuoi)), kergauss, 'same');
    plot(psthtli, squeeze(mean(temppsth,[2 3])), '-', 'Color', neucol, 'LineWidth', lw)
    %shadedErrorBar(psthtli/1000, squeeze(mean(temppsth,2)), squeeze(std(temppsth,0,2)/sqrt(nnz(neuoi))), {'-', 'Color', ttcol, 'LineWidth', 1}, 1)
    if neuopt==1
        ylabel('Rate (Hz)')
    end
end
xlim([-50 200])
%     ylim([2 15])
%     ylim([0 25])
set(gca, 'FontSize', fs, 'XGrid', 'on')
legend({'center-RF', 'non-center-RF'})
xlabel('Time (s)')
ylabel('Rate (Hz)')
%     title(sprintf('%s %s', visareas{iprobe}, neutitle), 'interpreter', 'none')
%     title(sprintf('%s %s N=%d', visareas{iprobe}, neutitle, nnz(neuoi)), 'interpreter', 'none')


figure
hold all
for neuopt = 1:2
    switch neuopt
        case 1
            %yyaxis left
            neutitle = 'center-RF neurons';
            neuoi = sesneuoiagg{iprobe}<=4 & RFCIagg(iprobe).RFindclassic==1 & RFCIagg(iprobe).pRFclassic<0.05;
            neucol = [0.5 0 1];
        case 2
            %yyaxis right
            neutitle = 'non-center-RF neurons';
            neuoi = sesneuoiagg{iprobe}<=4 & RFCIagg(iprobe).RFindclassic~=1 & RFCIagg(iprobe).pRFclassic<0.05;
            neucol = [1 0.7 0];
    end
    ttoi = ismember(ICtrialtypes, [106 111]);
    temp = convn(squeeze(psthavg_ctxagg.(whichvisblock){iprobe}(:,ttoi,neuoi)), kergauss, 'same');
    tempbase = squeeze(mean(temp(psthtli>=-200 & psthtli<0,:,:), 'all'));
    temppeak = max(squeeze(mean(temp(psthtli>=0 & psthtli<400,:,:), [2,3])));
    plot(psthtli, (squeeze(mean(temp, [2,3]))-tempbase)/(temppeak-tempbase), '-', 'Color', neucol, 'LineWidth', lw)
    %shadedErrorBar(psthtli/1000, squeeze(mean(temppsth,2)), squeeze(std(temppsth,0,2)/sqrt(nnz(neuoi))), {'-', 'Color', ttcol, 'LineWidth', 1}, 1)
    %         if neuopt==1
    %             ylabel('Rate (Hz)')
    %         end
end
xlim([-50 200])
%     ylim([2 15])
%     ylim([0 25])
ylim([-0.5 1])
set(gca, 'FontSize', fs, 'XGrid', 'on')
legend({'center-RF', 'non-center-RF'})
xlabel('Time (s)')
ylabel('Normalized Activity')
%     title(sprintf('%s %s', visareas{iprobe}, neutitle), 'interpreter', 'none')
%     title(sprintf('%s %s N=%d', visareas{iprobe}, neutitle, nnz(neuoi)), 'interpreter', 'none')

%% ctr-CRF vs non-ctr-CRF, IC-encoder vs inducer-encoder
whichvisblock = 'ICwcfg1_presentations';

kerwinhalf = 2; kersigma = 1;
kerwinhalf = 4; kersigma = 2;
% kerwinhalf = 7; kersigma = 3;
kergauss = normpdf( (-kerwinhalf:kerwinhalf)', 0,kersigma);
kergauss = (kergauss/sum(kergauss));

fs=16; lw = 3;

a = 1;
iprobe = probehierarchy(a);
figure
hold all
for neuopt = 1:2
        switch neuopt
            case 1
            yyaxis left
                neutitle = 'IC-encoders';
                neuoi = sesneuoiagg{iprobe}<=4 & ICsigagg.(whichvisblock)(iprobe).ICencoder1==1;
                neucol = [0 0.7 0];
            case 2
            yyaxis right
                neutitle = 'inducer-encoders';
                neuoi = sesneuoiagg{iprobe}<=4 & ICsigagg.(whichvisblock)(iprobe).indenc13==1;
                neucol = [0 0 0];
        end
        ttoi = ismember(ICtrialtypes, [106 111]);
    temppsth = convn(squeeze(psthavg_ctxagg.(whichvisblock){iprobe}(:,ttoi,neuoi)), kergauss, 'same');
    plot(psthtli, squeeze(mean(temppsth,[2 3])), '-', 'Color', neucol, 'LineWidth', lw)
    %shadedErrorBar(psthtli/1000, squeeze(mean(temppsth,2)), squeeze(std(temppsth,0,2)/sqrt(nnz(neuoi))), {'-', 'Color', ttcol, 'LineWidth', 1}, 1)
    if neuopt==1
        ylabel('Rate (Hz)')
    end
end
xlim([-50 200])
%     ylim([2 15])
%     ylim([0 25])
set(gca, 'FontSize', fs, 'XGrid', 'on')
legend({'IC-enc', 'ind-enc'})
xlabel('Time (s)')
ylabel('Rate (Hz)')
%     title(sprintf('%s %s', visareas{iprobe}, neutitle), 'interpreter', 'none')
%     title(sprintf('%s %s N=%d', visareas{iprobe}, neutitle, nnz(neuoi)), 'interpreter', 'none')


figure
hold all
for neuopt = 1:2

        switch neuopt
            case 1
                neutitle = 'IC-encoders';
                neuoi = sesneuoiagg{iprobe}<=4 & ICsigagg.(whichvisblock)(iprobe).ICencoder1==1;
                neucol = [0 0.7 0];
            case 2
                neutitle = 'inducer-encoders';
                neuoi = sesneuoiagg{iprobe}<=4 & ICsigagg.(whichvisblock)(iprobe).indenc13==1;
                neucol = [0 0 0];
        end
        ttoi = ismember(ICtrialtypes, [106 111]);
    temp = convn(squeeze(psthavg_ctxagg.(whichvisblock){iprobe}(:,ttoi,neuoi)), kergauss, 'same');
    tempbase = squeeze(mean(temp(psthtli>=-200 & psthtli<0,:,:), 'all'));
    temppeak = max(squeeze(mean(temp(psthtli>=0 & psthtli<400,:,:), [2,3])));
    plot(psthtli, (squeeze(mean(temp, [2,3]))-tempbase)/(temppeak-tempbase), '-', 'Color', neucol, 'LineWidth', lw)
    %shadedErrorBar(psthtli/1000, squeeze(mean(temppsth,2)), squeeze(std(temppsth,0,2)/sqrt(nnz(neuoi))), {'-', 'Color', ttcol, 'LineWidth', 1}, 1)
    %         if neuopt==1
    %             ylabel('Rate (Hz)')
    %         end
end
xlim([-50 200])
%     ylim([2 15])
%     ylim([0 25])
ylim([-0.5 1])
set(gca, 'FontSize', fs, 'XGrid', 'on')    
legend({'IC-enc', 'ind-enc'})
xlabel('Time (s)')
ylabel('Normalized Activity')
%     title(sprintf('%s %s', visareas{iprobe}, neutitle), 'interpreter', 'none')
%     title(sprintf('%s %s N=%d', visareas{iprobe}, neutitle, nnz(neuoi)), 'interpreter', 'none')
