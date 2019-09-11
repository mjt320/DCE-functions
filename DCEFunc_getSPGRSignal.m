function SPGRSignal=DCEFunc_getSPGRSignal(S0,T1_s,T2s_s,TR_s,TE_s,FA_deg)
%calculate signal for SPGR sequence
%OUTPUT:
%vector of signal values
%INPUT:
%S0=fully-relaxed signal
%T1_s=T1 relaxation time
%T2s_s=T2* relaxation time
%TR_s, TE_s, FA_deg: acquisition parameters

SPGRSignal = S0 .* ...
    (((1-exp(-TR_s./T1_s)).*sin(2*pi*(FA_deg/360))) ./ (1-exp(-TR_s./T1_s).*cos(2*pi*(FA_deg/360))) ) ...
    .* exp(-TE_s./T2s_s) ;

end