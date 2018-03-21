function SPGRSignal=DCEFunc_getSPGRSignal(S0,T1_s,T2s_s,TR_s,TE_s,beta_deg)
% calculate signal for SPGR sequence

%beta_rad = 2*pi*(beta_deg/360);
SPGRSignal = S0 .* ...
    (((1-exp(-TR_s./T1_s))*sin(2*pi*(beta_deg/360))) ./ (1-exp(-TR_s./T1_s)*cos(2*pi*(beta_deg/360))) ) ...
    .* exp(-TE_s./T2s_s) ;

end