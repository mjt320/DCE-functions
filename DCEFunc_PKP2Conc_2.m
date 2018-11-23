function [Ct_mM, C_cp_mM, C_e_mM] = H2O_PKP2Conc(tRes_s,Cp_AIF_mM,PKP,model,opts)
% OUTPUT:
% Ct_mM: column vector giving overall tissue concentration in mM
% C_cp_mM: column vector giving local concentration in capillary plasma
% C_e_mM: column vector giving local EES concentration
% INPUT:
% tRes_s = time resolution of data in seconds
% Cp_AIF_mM = column vector giving AIF plasma concentration in mM
% PKP = struct containing PK parameters (vP, vE, PS_perMin, FP_mlPer100gPerMin)
% model = string to specify model:
%   Patlak
%   2CXM = two-compartment exchange model

N=size(Cp_AIF_mM,1);
Ct_mM = nan(N,1);
C_cp_mM = nan(N,1);
C_e_mM = nan(N,1);

switch model
    %% Approach is to calculate the propagators for capillary plasma and the 
    % EES as defined in Sourbron and Buckley (2012) Phys. Med. Biol.
    % C_cp and C_e are then convolutions of these propagators with C_p_AIF
    % Note that C_p_AIF is discrete and the propagators are continuous.
    % Therefore the convolution is calculated exactly for each time step,
    % then summed over all time steps.
    case 'Patlak' % calculates IRF for Patlak model
        IRF = PatlakIRF();
        Ct_mM = conv(Cp_AIF_mM.',IRF,'full').';
        Ct_mM = Ct_mM(1:N);
    case '2CXM'        
        v = PKP.vE + PKP.vP;
        T = v/PKP.FP_mlPer100gPerMin;
        Tc = PKP.vP/(PKP.FP_mlPer100gPerMin);
        Te = PKP.vE/PKP.PS_perMin;
        sig_p = ((T+Te) + sqrt((T + Te)^2 - (4*Tc*Te)))/(2*Tc*Te);
        sig_n = ((T+Te) - sqrt((T + Te)^2 - (4*Tc*Te)))/(2*Tc*Te);
        sig_rat = (sig_p*sig_n)/(sig_p - sig_n);
        sig_p_rec = -1/sig_p;
        sig_n_rec = -1/sig_n;
        [h_cp, h_e] = prop2CXM();
        
        %% Calculate C_cp and C_e by convolution of AIF with propogators
        C_cp_mM = conv(Cp_AIF_mM.',h_cp,'full').';
        C_cp_mM = C_cp_mM(1:N); % remove extra entries so that C_cp is same length as AIF
        C_e_mM = conv(Cp_AIF_mM.',h_e,'full').';
        C_e_mM = C_e_mM(1:N);
        %% Calculate C_t from C_cp and C_e
        Ct_mM = (PKP.vP*C_cp_mM) + (PKP.vE*C_e_mM);
        
    otherwise
        error('Model not recognised.');
end

%% Function to calculate IRF for Patlak
    function PatlakIRF() % function to calculate discrete IRF by taking mean for each time point
        IRF=nan(1,N);
        IRF(1)=PKP.vP + (PKP.PS_perMin/2)*(tRes_s/60); % IRF at time zero
        IRF(2:N)=PKP.PS_perMin * (tRes_s/60); % IRF at time zero+t_res, zero+2*t_res, ...
    end

%% Functions to calculate propogators for 2CXM
    function [h_cp,h_e] = prop2CXM()
        h_cp=nan(1,N);
        h_e=nan(1,N);
        h_cp(1)=hcp_2CXMIntegral(0,(tRes_s/60)/2); % h_cp for first time point
        h_e(1)=he_2CXMIntegral(0,(tRes_s/60)/2); % h_e for first time point
        for iTime=2:N
            h_cp(iTime)=hcp_2CXMIntegral( (iTime-1-0.5)*(tRes_s/60), (iTime-1+0.5)*(tRes_s/60) );
            h_e(iTime)=he_2CXMIntegral( (iTime-1-0.5)*(tRes_s/60), (iTime-1+0.5)*(tRes_s/60) );
        end
    end

%% Functions to calculate exact integral of propogators over ranges for 2CXM
    function prop_hcp_Int = hcp_2CXMIntegral(t1,t2)
        prop_hcp_Int = (sig_rat*((sig_n_rec*(1-(Te*sig_n))*exp(-1*t2*sig_n))+(sig_p_rec*((Te*sig_p)-1)*exp(-1*t2*sig_p))))...
           - (sig_rat*((sig_n_rec*(1-(Te*sig_n))*exp(-1*t1*sig_n))+(sig_p_rec*((Te*sig_p)-1)*exp(-1*t1*sig_p))));
    end

    function prop_he_Int = he_2CXMIntegral(t1,t2)
        prop_he_Int = (sig_rat*((sig_n_rec*exp(-1*t2*sig_n))-(sig_p_rec*exp(-1*t2*sig_p))))...
            - (sig_rat*((sig_n_rec*exp(-1*t1*sig_n))-(sig_p_rec*exp(-1*t1*sig_p))));
    end


end


