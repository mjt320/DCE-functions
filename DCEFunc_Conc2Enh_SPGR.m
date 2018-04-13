function enhancementPct=DCEFunc_Conc2Enh_SPGR(conc_mM,T10_s,TR_s,TE_s,FA_deg,r1_permMperS,r2s_permMperS)
% calculate signal enhancement from contrast agent concentration for SPGR sequence
% OUTPUT:
% enhancementPct: array of enhancements including; each column = 1 time series
% INPUT:
% conc_mM: array of concentration estimates; each column = 1 time series
% T10_s: row vector of pre-contrast T1 values; each entry corresponds to 1
% time series
% TR_s, TE_s: acquisition parameters for SPGR/FLASH sequence
% FA_deg: row vector of flip angles; each value corresponds to a time series
% r1_permMperS, r2_permMperS: T1 and T2 relaxivity values. Ignores T2'
% effects.

enhancementPct=nan(size(conc_mM));

NTimePoints=size(conc_mM,1);
NSeries=size(conc_mM,2);

S0=1; % this value doesn't affect enhancement for SPGR
T2s0_s = 1; % this value doesn't affect enhancement for SPGR

for iSeries=1:NSeries %loop through different time series (e.g. voxels or ROIs)
    signal_pre = DCEFunc_getSPGRSignal(S0,T10_s(1,iSeries),T2s0_s,TR_s,TE_s,FA_deg(1,iSeries)); %calculate pre-contrast signal
    for iTimePoint=1:NTimePoints %loop through time points
        T2s_s = ( T2s0_s^-1 + r2s_permMperS * conc_mM(iTimePoint,iSeries) ).^-1; %calculate T2* post-enhancement
        T1_s = ( T10_s(1,iSeries)^-1 + r1_permMperS * conc_mM(iTimePoint,iSeries) ).^-1; %calculate T1 post-enhancement
        
        signal_post = DCEFunc_getSPGRSignal(S0,T1_s,T2s_s,TR_s,TE_s,FA_deg(1,iSeries)); %calculate post-contrast signal
        enhancementPct(iTimePoint,iSeries) = 100 * ( (signal_post - signal_pre) / signal_pre ); %calculate enhancement
    end
end

end
