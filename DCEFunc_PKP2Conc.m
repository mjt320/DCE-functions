function [Ct_mM, IRF, C_cp_mM, C_e_mM] = DCEFunc_PKP2Conc(tRes_s,Cp_AIF_mM,PKP,model,opts)
% OUTPUT:
% Ct_mM: column vector giving overall tissue concentration in mM
% IRF: Impulse Response Function
% C_cp_mM: column vector giving local concentration in capillary plasma
% C_e_mM: column vector giving local EES concentration
% INPUT:
% tRes_s = time resolution of data in seconds
% Cp_AIF_mM = column vector giving AIF plasma concentration in mM
% PKP = struct containing PK parameters (vP, vE, PS_perMin, FP_mlPer100gPerMin)
% model = string to specify model:
%   Patlak
%   2CXM = two-compartment exchange model

% create empty arrays
N=size(Cp_AIF_mM,1); % size of arrays
Ct_mM = nan(N,1); % tissue concentrations
C_cp_mM = nan(N,1); % capillary concentrations
C_e_mM = nan(N,1); % EES concentrations
IRF =  nan(N,1); % Impulse Response Function
h_cp=nan(1,N); % capillary propogators
h_e=nan(1,N); % EES propogators

switch model
    %% Approach is to calculate the propagators for capillary plasma and the 
    % EES as defined in Sourbron and Buckley (2012) Phys. Med. Biol.
    % C_cp and C_e are then convolutions of these propagators with C_p_AIF
    % Note that C_p_AIF is discrete and the propagators are continuous.
    % Therefore the convolution is calculated exactly for each time step,
    % then summed over all time steps.
    
    case 'Patlak' % calculates IRF for Patlak model
        if ~isfield(PKP,'vE'); PKP.vE=1; end %set nominal value for vE so that general propogator approach works for Patlak model
        h_cp(1)=1 ; % h_cp for first time point
        h_e(1)=(PKP.PS_perMin/(2*PKP.vE))*(tRes_s/60); % h_e for first time point
        for iTime=2:N
            h_cp(iTime)=0;
            h_e(iTime)=(PKP.PS_perMin/PKP.vE) * (tRes_s/60);
        end
        
    case '2CXM'
        % calculate constants
        v = PKP.vE + PKP.vP;
        T = v/((1/100)*PKP.FP_mlPer100gPerMin);
        Tc = PKP.vP/((1/100)*PKP.FP_mlPer100gPerMin);
        Te = PKP.vE/PKP.PS_perMin;
        sig_p = ((T+Te) + sqrt((T + Te)^2 - (4*Tc*Te)))/(2*Tc*Te);
        sig_n = ((T+Te) - sqrt((T + Te)^2 - (4*Tc*Te)))/(2*Tc*Te);
        sig_rat = (sig_p*sig_n)/(sig_p - sig_n);
        sig_p_rec = -1/sig_p;
        sig_n_rec = -1/sig_n;
        % calculate propogators for 2CXM
        h_cp(1) = hcp_2CXMIntegral(0,(tRes_s/60)/2); % h_cp for first time point
        h_e(1) = he_2CXMIntegral(0,(tRes_s/60)/2); % h_e for first time point
        for iTime=2:N
            h_cp(iTime) = hcp_2CXMIntegral( (iTime-1-0.5)*(tRes_s/60), (iTime-1+0.5)*(tRes_s/60) );
            h_e(iTime) = he_2CXMIntegral( (iTime-1-0.5)*(tRes_s/60), (iTime-1+0.5)*(tRes_s/60) );
        end
    otherwise
        error('Model not recognised.');
end

%% Calculate C_cp and C_e by convolution of AIF with propogators
C_cp_mM = conv(Cp_AIF_mM.',h_cp,'full').';
C_e_mM = conv(Cp_AIF_mM.',h_e,'full').';
C_cp_mM = C_cp_mM(1:N); % remove extra entries so that C_cp and C_e are same length as AIF
C_e_mM = C_e_mM(1:N);

%% calculate IRF
IRF = (h_e*PKP.vE) + (h_cp*PKP.vP);

%% Calculate C_t from C_cp and C_e
Ct_mM = (PKP.vP*C_cp_mM) + (PKP.vE*C_e_mM);

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
