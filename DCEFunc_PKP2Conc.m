function [Ct_mM, IRF, c_cp_mM, c_e_mM] = DCEFunc_PKP2Conc(tRes_s,Cp_AIF_mM,PKP,model,opts)
% OUTPUT:
% Ct_mM: array of overall tissue concentrations in mM, 1 col per time series
% IRF: array containing Impulse Response Functions, 1 col per time series
% c_cp_mM: array giving local concentration in capillary plasma, 1 col per time series
% c_e_mM: array giving local EES concentrations, 1 col per time series
% INPUT:
% tRes_s = time resolution of data in seconds
% Cp_AIF_mM = column vector giving AIF plasma concentration in mM
% PKP = struct containing PK parameters (vP, vE, PS_perMin, FP_mlPer100gPerMin)
%   each field is a row vector containing parameter corresponding to each time series
% model = string to specify model:
%   Patlak
%   2CXM = two-compartment exchange model
%
% Authors: MJT, JMB, CM

NFrames=size(Cp_AIF_mM,1); % number of time frames
NSeries=size(PKP.vP,2); % number of time series

% create empty arrays
Ct_mM = nan(NFrames,NSeries); % tissue concentrations
c_cp_mM = nan(NFrames,NSeries); % capillary concentrations
c_e_mM = nan(NFrames,NSeries); % EES concentrations
IRF = nan(NFrames,NSeries); % Impulse Response Functions
h_cp = nan(NFrames,NSeries); % capillary propogators
h_e = nan(NFrames,NSeries); % EES propogators

switch model
    %% Approach is to calculate the propagators for capillary plasma and the
    % EES as defined in Sourbron and Buckley (2012) Phys. Med. Biol.
    % c_cp and c_e are then convolutions of these propagators with C_p_AIF
    % Note that C_p_AIF is discrete and the propagators are continuous.
    % Therefore the convolution is calculated exactly for each time step,
    % then summed over all time steps.
    
    case 'Patlak' % calculates IRF for Patlak model
        if ~isfield(PKP,'vE'); PKP.vE=1*ones(1,NSeries); end %set nominal value for vE so that general propogator approach works for Patlak model
        h_cp(1,:)=1 ; % h_cp for first time point
        h_e(1,:)=(PKP.PS_perMin./(2*PKP.vE))*(tRes_s/60); % h_e for first time point
        h_cp(2:NFrames,:)=0; % h_cp for subsequent time points
        h_e(2:NFrames,:)=(repmat(PKP.PS_perMin,[NFrames-1,1])./repmat(PKP.vE,[NFrames-1,1])) * (tRes_s/60); % h_e for subsequent time points
        
    case '2CXM'
        % calculate constants
        v = PKP.vE + PKP.vP;
        T = v./((1/100).*PKP.FP_mlPer100gPerMin);
        Tc = PKP.vP./((1/100).*PKP.FP_mlPer100gPerMin);
        Te = PKP.vE./PKP.PS_perMin;
        sig_p = ((T+Te) + sqrt((T + Te).^2 - (4*Tc.*Te)))./(2*Tc.*Te);
        sig_n = ((T+Te) - sqrt((T + Te).^2 - (4*Tc.*Te)))./(2*Tc.*Te);
        sig_rat = (sig_p.*sig_n)./(sig_p - sig_n);
        sig_p_rec = -1./sig_p;
        sig_n_rec = -1./sig_n;
        % calculate propogators for 2CXM
        h_cp(1,:) = hcp_2CXMIntegral(0,(tRes_s/60)/2); % h_cp for first time point
        h_e(1,:) = he_2CXMIntegral(0,(tRes_s/60)/2); % h_e for first time point
        iTimes=(2:NFrames).'; % calculate h_cp and h_e for all subsequent time points
        h_cp(iTimes,:) = hcp_2CXMIntegral( (iTimes-1-0.5)*(tRes_s/60), (iTimes-1+0.5)*(tRes_s/60) );
        h_e(iTimes,:) = he_2CXMIntegral( (iTimes-1-0.5)*(tRes_s/60), (iTimes-1+0.5)*(tRes_s/60) );

    otherwise
        error('Model not recognised.');
end

%% Calculate c_cp and c_e by convolution of AIF with propogators
c_cp_mM = conv2(1,Cp_AIF_mM.',h_cp.','full').';
c_e_mM = conv2(1,Cp_AIF_mM.',h_e.','full').';
c_cp_mM = c_cp_mM(1:NFrames,:); % remove extra entries so that c_cp and c_e are same length as AIF
c_e_mM = c_e_mM(1:NFrames,:);

%% calculate IRF
IRF = (h_e.*repmat(PKP.vE,[NFrames,1])) + (h_cp.*repmat(PKP.vP,[NFrames,1]));

%% Calculate C_t from c_cp and c_e
Ct_mM = (repmat(PKP.vP,[NFrames,1]).*c_cp_mM) + (repmat(PKP.vE,[NFrames,1]).*c_e_mM);

%% Functions to calculate exact integrals over propogators for 2CXM
    function prop_hcp_Int = hcp_2CXMIntegral(t1,t2)
        prop_hcp_Int = (sig_rat.*((sig_n_rec.*(1-(Te.*sig_n)).*exp(-1*t2*sig_n))+(sig_p_rec.*((Te.*sig_p)-1).*exp(-1*t2*sig_p))))...
            - (sig_rat.*((sig_n_rec.*(1-(Te.*sig_n)).*exp(-1*t1*sig_n))+(sig_p_rec.*((Te.*sig_p)-1).*exp(-1*t1*sig_p))));
    end

    function prop_he_Int = he_2CXMIntegral(t1,t2)
        prop_he_Int = (sig_rat.*((sig_n_rec.*exp(-1*t2*sig_n))-(sig_p_rec.*exp(-1*t2*sig_p))))...
            - (sig_rat.*((sig_n_rec.*exp(-1*t1*sig_n))-(sig_p_rec.*exp(-1*t1*sig_p))));
    end


end
