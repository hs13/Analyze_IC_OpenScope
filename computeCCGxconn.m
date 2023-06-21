% simply return peak & trough values and indices within the designated time period
% TODO: THIS CODE CAN DEFINITELY BE MADE FASTER BY ONLY CHECKING PEAKS WHEN
% MAX IS AT BOUNDARIES
function [excpeak, inhtrough] = computeCCGxconn(CCGtli, CCGin, Tptlim)

% the 2X in the 3rd dimension of each peak/trough matrix are in the order of:
% peak value (Hz), peak ms
Nunits = size(CCGin,1);
excpeak = NaN(Nunits, Nunits, 2);
inhtrough = NaN(Nunits, Nunits, 2);

for ipre = 1:Nunits
    for ipost = ipre+1:Nunits
        tempCCG = squeeze(CCGin(ipre, ipost, :));
                
        pretopostpeak = NaN(2,1);
        pretoposttrough = NaN(2,1); 
        posttoprepeak = NaN(2,1); 
        posttopretrough = NaN(2,1);
        
        [pks, pklocs]=findpeaks(tempCCG);
        
        if nnz(CCGtli(pklocs)>=1 & CCGtli(pklocs)<=Tptlim)
            % maximum local peak value within the temporal constraitns
            [v, mpki]=max(pks(CCGtli(pklocs)>=1 & CCGtli(pklocs)<=Tptlim));
            pklocsqual = pklocs(CCGtli(pklocs)>=1 & CCGtli(pklocs)<=Tptlim);
            mtli = pklocsqual(mpki);
            pretopostpeak(1) = v;
            pretopostpeak(2) = CCGtli(mtli);
        else
            % there is no peak in the designated time period, pretopostpeak and pretopostexc are NaN and false
        end
        
        if nnz(CCGtli(pklocs)<=-1 & CCGtli(pklocs)>=-Tptlim)
            % maximum local peak value within the temporal constraitns
            [v, mpki]=max(pks(CCGtli(pklocs)<=-1 & CCGtli(pklocs)>=-Tptlim));
            pklocsqual = pklocs(CCGtli(pklocs)<=-1 & CCGtli(pklocs)>=-Tptlim);
            mtli = pklocsqual(mpki);
            posttoprepeak(1) = v;
            posttoprepeak(2) = -CCGtli(mtli);
        else
            % there is no peak in the designated time period, pretopostpeak and pretopostexc are NaN and false
        end
        
        
        [negtroughs, troughlocs]=findpeaks(-tempCCG);
        troughs = -negtroughs;
        
        if nnz(CCGtli(troughlocs)>=1 & CCGtli(troughlocs)<=Tptlim)
            % maximum local peak value within the temporal constraitns
            [v, mpki]=min(troughs(CCGtli(troughlocs)>=1 & CCGtli(troughlocs)<=Tptlim));
            troughlocsqual = troughlocs(CCGtli(troughlocs)>=1 & CCGtli(troughlocs)<=Tptlim);
            mtli = troughlocsqual(mpki);
            pretoposttrough(1) = v;
            pretoposttrough(2) = CCGtli(mtli);
        else
            % there is no trough in the designated time period, posttopretrough and posttopreinh are NaN and false
        end
        
        if nnz(CCGtli(troughlocs)<=-1 & CCGtli(troughlocs)>=-Tptlim)
            % maximum local peak value within the temporal constraitns
            [v, mpki]=min(troughs(CCGtli(troughlocs)<=-1 & CCGtli(troughlocs)>=-Tptlim));
            troughlocsqual = troughlocs(CCGtli(troughlocs)<=-1 & CCGtli(troughlocs)>=-Tptlim);
            mtli = troughlocsqual(mpki);
            posttopretrough(1) = v;
            posttopretrough(2) = -CCGtli(mtli);
        else
            % there is no trough in the designated time period, posttopretrough and posttopreinh are NaN and false
        end
        
        
        excpeak(ipre,ipost,:) = pretopostpeak;
        inhtrough(ipre,ipost,:) = pretoposttrough;
        
        if ipost ~= ipre
            excpeak(ipost,ipre,:) = posttoprepeak;
            inhtrough(ipost,ipre,:) = posttopretrough;
        end
    end
end
end