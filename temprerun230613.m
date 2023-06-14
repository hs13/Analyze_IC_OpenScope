nwbsessions = {'sub_1171903433','sub_1172968426','sub_1172969394','sub_1174569641', ...
    'sub_1175512783','sub_1176214862','sub_1177693342','sub_1181314060', ...
    'sub_1181585608','sub_1182593231','sub_1183369803','sub_1186544726', ...
    'sub_1189891322','sub_1194090570'};
Nsessions = numel(nwbsessions);
datadir = '/Users/hyeyoung/Documents/DATA/OpenScopeData/000248/';

for ises = 1:Nsessions
    sesclk = tic;
    fprintf('%d/%d %s\n', ises, Nsessions, nwbsessions{ises})
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    
    try
        load([pathpp 'trackmouseeye.mat'])
        valideyetracking = true;
    catch
        fprintf('Skipping session%d %s -- EyeTracking error in nwb\n', ises, nwbsessions{ises})
        valideyetracking = false;
        %continue
    end
    gazedistthresh = 20;
    visblocks = {'ICkcfg0_presentations','ICkcfg1_presentations','ICwcfg0_presentations','ICwcfg1_presentations', ...
        'RFCI_presentations','sizeCI_presentations'}; %,'spontaneous_presentations'};
    probes = {'A', 'B', 'C', 'D', 'E', 'F'};
    for iprobe = 1:numel(probes)
        
        ppfn = sprintf('%spostprocessed_probe%s.mat', pathpp, probes{iprobe});
        load(ppfn)
        load(sprintf('%svisresponses_probe%s.mat', pathpp, probes{iprobe}))
        
        ICtrialtypes = [0 101 105 106 107 109 110 111 506 511 1105 1109 1201 1299 ...
            1301 1302 1303 1304 1305 1306 1307 1308];
        ICblocks = 1:4;
        %ICsig = struct();
        for b = ICblocks
            disp(visblocks{b})
            tloi = psthtli>0 & psthtli<=400;
            R = squeeze(1000*mean(psth.(visblocks{b})(tloi,:,:), 1))';
            trialorder = ICtrialtypes( vis.(visblocks{b}).trialorder + 1);
            %ICsig.(visblocks{b}) = analyzeStaticICtxi(R, trialorder);
            
            Nneurons = size(R, 1);
            BICREltt = [0 106 506 111 511]; % Blank, IC and REl trial type
            SP_BICREl = NaN(Nneurons, numel(BICREltt)-1);
            Pmww_BICREl = NaN(Nneurons, numel(BICREltt)-1);
            blanktrials = trialorder==0;
            for typi = 2:numel(BICREltt)
                temptrials = trialorder==BICREltt(typi);
                if nnz(blanktrials)==0 || nnz(temptrials)==0
                    continue
                end
                labels = [zeros(1,nnz(blanktrials)) ones(1,nnz(temptrials))];
                
                for ci = 1:Nneurons
                    scores = [R(ci, blanktrials) R(ci, temptrials)];
                    
                    [X,Y,T,AUC,OPTROCPT] = perfcurve(labels,scores,'1');% , 'NBoot',Nshuf);
                    SP_BICREl(ci,typi-1) = AUC;
                    
                    Pmww_BICREl(ci,typi) = ranksum(R(ci, blanktrials), R(ci, temptrials));
                end
            end
            disp('calculated SP_BICREl')
            if ~isequal(SP_BICREl, ICsig.(visblocks{b}).SP_BICREl)
                error('check SP_BICREl/Pmww_BICREl calculation')
            end
            ICsig.(visblocks{b}).Pmww_BICREl = Pmww_BICREl;
        end
        
        save(sprintf('%svisresponses_probe%s.mat', pathpp, probes{iprobe}), ...
            'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'ICsig', 'RFCI', 'RFCIspin', ...
            'sizeCI', 'oriparams', 'ori4params', '-v7.3')
        
        if valideyetracking
            vrfgfn= sprintf('%svisresponses_fixedgaze%dpix_probe%s.mat', pathpp, gazedistthresh, probes{iprobe});
            load(vrfgfn)
            
            ICtrialtypes = [0 101 105 106 107 109 110 111 506 511 1105 1109 1201 1299 ...
                1301 1302 1303 1304 1305 1306 1307 1308];
            ICblocks = 1:4;
            %ICsig_fixedgaze = struct();
            for b = ICblocks
                disp([visblocks{b} ' fixed gaze'])
                tloi = psthtli>0 & psthtli<=400;
                tempR = squeeze(1000*mean(psth.(visblocks{b})(tloi,:,:), 1))';
                temptrialorder = ICtrialtypes( vis.(visblocks{b}).trialorder + 1);
                temptrialsfixedgaze = trialmaxdistmodecom.(visblocks{b})<gazedistthresh & ~triallikelyblink.(visblocks{b});
                validICfix = all(ismember(ICtrialtypes, unique(temptrialorder(temptrialsfixedgaze))));
                if validICfix
                    %ICsig_fixedgaze.(visblocks{b}) = analyzeStaticICtxi(tempR(:,temptrialsfixedgaze), temptrialorder(temptrialsfixedgaze) );
                    R = tempR(:,temptrialsfixedgaze);
                    trialorder = temptrialorder(temptrialsfixedgaze);
                    
                    Nneurons = size(R, 1);
                    BICREltt = [0 106 506 111 511]; % Blank, IC and REl trial type
                    SP_BICREl = NaN(Nneurons, numel(BICREltt)-1);
                    Pmww_BICREl = NaN(Nneurons, numel(BICREltt)-1);
                    blanktrials = trialorder==0;
                    for typi = 2:numel(BICREltt)
                        temptrials = trialorder==BICREltt(typi);
                        if nnz(blanktrials)==0 || nnz(temptrials)==0
                            continue
                        end
                        labels = [zeros(1,nnz(blanktrials)) ones(1,nnz(temptrials))];
                        
                        for ci = 1:Nneurons
                            scores = [R(ci, blanktrials) R(ci, temptrials)];
                            
                            [X,Y,T,AUC,OPTROCPT] = perfcurve(labels,scores,'1');% , 'NBoot',Nshuf);
                            SP_BICREl(ci,typi-1) = AUC;
                            
                            Pmww_BICREl(ci,typi) = ranksum(R(ci, blanktrials), R(ci, temptrials));
                        end
                    end
                    disp('calculated SP_BICREl fixed gaze')
                    if ~isequal(SP_BICREl, ICsig_fixedgaze.(visblocks{b}).SP_BICREl)
                        error('check SP_BICREl/Pmww_BICREl calculation')
                    end
                    ICsig_fixedgaze.(visblocks{b}).Pmww_BICREl = Pmww_BICREl;
                    
                %else
                %    ICsig_fixedgaze.(visblocks{b}) = struct();
                %    disp('skipped')
                end
            end
            
            save(vrfgfn, 'meanFRvec', 'sponFRvec', 'ICtrialtypes', 'spinmaxdistmodecom', 'spinlikelyblink', ...
                'ICsig_fixedgaze', 'RFCI_fixedgaze', 'RFCIspin_fixedgaze', ...
                'sizeCI_fixedgaze', 'oriparams_fixedgaze', 'ori4params_fixedgaze', '-v7.3')
        end
    end
    toc(sesclk)
end