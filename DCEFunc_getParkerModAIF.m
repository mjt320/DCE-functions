function [t_s,Cp_AIF_mM]=DCEFunc_getParkerModAIF(t_res_s,t_acq_s,t_start_s,Hct)
%returns modified version of population average plasma AIF (Parker 2006) in mMol and time points t
%modifications are described in Heye et al., Neuroimage (2016)
%t_s: time points are placed at intervals of t_res_s, starting at t_res_s/2
%(assuming linear k-space sampling and generates time relative to start of
%acquisition)
%Cp_AIF_mM: column vector of arterial plasma concentrations
%t_res_s=temporal resolution
%t_acq_s=duration of acquisition
%t_start_s=time of start of injection
%Hct=hematocrit (needed as Parker forumula is for Cb not Cp)

t_s=((t_res_s/2):t_res_s:(t_acq_s-t_start_s)).'; %calculate time points (starting at t_res/2)
t_min=t_s/60; %parker uses minute units, so convert time;

%%set Parker parameters defined in his paper
A1=0.809; A2=0.330;
T1=0.17046; T2=0.365;
sigma1=0.0563; sigma2=0.132;
s=38.078; tau=0.483;
alpha=1.050; beta=0.1685;
alpha=3.1671; beta=1.0165;
alpha2=0.5628; beta2=0.0266;

Cb_mM=(A1/(sigma1*sqrt(2*pi)))*exp(-((t_min-T1).^2)/(2*sigma1^2)) + ...
    (A2/(sigma2*sqrt(2*pi)))*exp(-((t_min-T2).^2)/(2*sigma2^2)) + ...
    (alpha*exp(-beta*t_min) + alpha2*exp(-beta2*t_min))./(1+exp(-s*(t_min-tau)));

Cp_AIF_mM=Cb_mM/(1-Hct); %convert to Cp

t_preContrast_s=fliplr(t_start_s-t_res_s/2:-t_res_s:0).'; %pre-contrast time points (calculate backwards to zero, then reverse, so that time interval is constant)
%t_preContrast_s=(0:t_res_s:t_start_s).'; %pre-contrast time points

t_s=[ t_preContrast_s ; t_s + t_start_s]; %add pre-contrast time points
Cp_AIF_mM=[zeros(size(t_preContrast_s)) ; Cp_AIF_mM]; % add pre-contrast concentrations

end