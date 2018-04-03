function conc_mM=DCEFunc_Enh2Conc_SPGR(enhancementPct,T10_s,TR_s,TE_s,FA_deg,r1_permMperS,r2s_permMperS)
% calculate concentration for SPGR sequence
% conc_mM: array of concentration estimates; each column = 1 time series
% enhancementPct: array of enhancements including; each column = 1 time series
% T10_s: row vector of pre-contrast T1 values; each entry corresponds to 1
% time series
% TR_s, TE_s: acquisition parameters
% FA_deg: row vector of flip angles; each value corresponds to a different
% time series
% r1_permMperS, r2_permMperS: relaxivity values

conc_mM=nan(size(enhancementPct));

NTimePoints=size(enhancementPct,1);
NSeries=size(enhancementPct,2);

S0=1; % this value is generally unknown and doesn't affect enhancement for SPGR
T2s0_s = 1; % this value is generally unknown and doesn't affect enhancement for SPGR

%% calculate concentration values to try
maxC=8; %maximum concentration
maxNegC=0.1; %maximum negative concentration
minPctInc=1; %minimum percent increase from previous value
minAbsInc=0.0001; %minimum absolute increase from previous value
CTry_mM=-maxNegC;
while 1
    c = CTry_mM(1,end)+max(minAbsInc,(minPctInc/100)*abs(CTry_mM(1,end))); %start at max negative value and increase by pct or abs increment (whichever is greater)
    if c>maxC; break; else
        CTry_mM = [CTry_mM c];
    end
end

T2s_try_s = ( T2s0_s^-1 + r2s_permMperS * CTry_mM ).^-1; %create row vector of T2s for each test concentration

for iSeries=1:NSeries %loop through different time series (e.g. voxels or ROIs)
    
    if isnan(T10_s(1,iSeries)); continue; end; %skip voxels without T1 value
    
    predictedSignal_pre = DCEFunc_getSPGRSignal(S0,T10_s(1,iSeries),T2s0_s,TR_s,TE_s,FA_deg(1,iSeries)); %calculate pre-contrast signal
    
    T1_try_s = ( T10_s(1,iSeries)^-1 + r1_permMperS * CTry_mM ).^-1; %create row vector of T1 for each test concentration
    
    T1_try_filt_s= T1_try_s(T1_try_s > 0); % select only test concentrations that lead to positive T1 values
    T2s_try_filt_s = T2s_try_s(T1_try_s > 0);
    CTry_filt_mM = CTry_mM(T1_try_s > 0);
    NTry_filt=size(CTry_filt_mM,2);
    
    predictedSignal_post = DCEFunc_getSPGRSignal(S0,T1_try_filt_s,T2s_try_filt_s,TR_s,TE_s,FA_deg(1,iSeries)); %calculate row vector of signals with a value for each test concentration
    predictedEnhancementPct = 100 * ( (predictedSignal_post - predictedSignal_pre) ./ predictedSignal_pre ); %calculed predicted enhancements
    
    diff = abs(repmat(predictedEnhancementPct,[NTimePoints 1]) - repmat(enhancementPct(:,iSeries),[1 NTry_filt])); %absolute difference between predicted and actual enhancement for each test conc and time point
    [y,minIdx]=min(diff.'); % choose test concentration values with min difference for each time point
    conc_mM(:,iSeries)=CTry_filt_mM(minIdx); % set concentration value for each time point
    
end

end


