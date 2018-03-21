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

S0=1;
%CTry_mM=[-0.1:0.0001:0.1 0.101:0.001:1 1.01:0.01:8.00]; %calculate enhancement for all trial conc. values and select conc. corresponding to closest enhancement
CTry_mM=[-0.1:0.00001:0.1 0.1001:0.0001:1 1.001:0.001:8.00]; %calculate enhancement for all trial conc. values and select conc. corresponding to closest enhancement
T2s0_s = 1; % this value is generally unknown and doesn't affect enhancement for SPGR

T2s_try_s = ( T2s0_s^-1 + r2s_permMperS * CTry_mM ).^-1; %create row vector of T2s for each test concentration

for iSeries=1:NSeries %loop through different time series (e.g. voxels or ROIs)
    
    if isnan(T10_s(1,iSeries)); continue; end; %skip voxels without T1 value
    
    predictedSignal_pre = DCEFunc_getSPGRSignal(S0,T10_s(1,iSeries),T2s0_s,TR_s,TE_s,FA_deg(1,iSeries)); %calculate pre-contrast signal
    
    T1_try_s = ( T10_s(1,iSeries)^-1 + r1_permMperS * CTry_mM ).^-1; %create row vector of T1 for each test concentration
    
    T1_try_filt_s= T1_try_s(T1_try_s > 0); % select only test concentrations that lead to positive T1 values
    T2s_try_filt_s = T2s_try_s(T1_try_s > 0);
    CTry_filt_mM = CTry_mM(T1_try_s > 0);
    
        
    predictedSignal_post = DCEFunc_getSPGRSignal(S0,T1_try_filt_s,T2s_try_filt_s,TR_s,TE_s,FA_deg(1,iSeries)); %calculate row vector of signals with a value for each test concentration
    predictedEnhancementPct = 100 * ( (predictedSignal_post - predictedSignal_pre) ./ predictedSignal_pre ); %calculed predicted enhancements
    
    for iTimePoint=1:NTimePoints %loop through time points
        
        diff = abs(predictedEnhancementPct - enhancementPct(iTimePoint,iSeries)); %absolute difference between predicted and actual enhancement for each test conc
        [y,minIdx]=min(diff); % choose test conc with min difference
        conc_mM(iTimePoint,iSeries)=CTry_filt_mM(minIdx); % set concentration value
        
    end
end

end


%% local function to predict SPGR signal
% function SPGRSignal = getSPGRSignal(S0,T1_s,T2s_s,TR_s,TE_s,beta_deg)
% %beta_rad = 2*pi*(beta_deg/360);
% SPGRSignal = S0 .* ...
%     (((1-exp(-TR_s./T1_s))*sin(2*pi*(beta_deg/360))) ./ (1-exp(-TR_s./T1_s)*cos(2*pi*(beta_deg/360))) ) ...
%     .* exp(-TE_s./T2s_s) ;
% end