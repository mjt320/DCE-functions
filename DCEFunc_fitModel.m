function [PKP, CtModelFit_mM]=DCEFunc_fitModel(tRes_s,Ct_mM,Cp_AIF_mM,model,opts)
% Fit concentration curve using model to generate pharmacokinetic
% parameters. At present, only the Patlak model is supported but other
% models may be added on request.
% Output:
% PKP = struct containing pharmacokinetic parameters (vP, vE, PS_perMin,
% FP_mlPer100gPerMin); each parameter is a row vector, with each entry
% corresponding to a time series
% CtModelFit_mM = best-fit concentration curve; each column corresponds to
% a time series.
% Input:
% tRes_s = time resolution in seconds for input data
% Ct_mM: array of tissue concentrations in mM; each column corresponds to a
% time series
% Cp_AIF_mM = column vector giving AIF plasma (not blood!) concentration in mM
% model = string to specify model:
%   'Patlak': non-linear implementation of Patlak.
%   'PatlakLinear': classic linear Patlak plot approach. Can lead to
%       unstable results for low values of Cp
%   'PatlakFast': multiple linear regression implementation of Patlak model
%       (recommended)
% opts = struct containing options:
%   NIgnore = number of points to ignore in cost function. For Patlak,
%       early data points may be excluded to reduce CBF effects
%   PatlakFastRegMode = regression mode for multiple linear regression Patlak fitting:
%       'linear' for standard regression (default),
%       'robust' for robust regression

if ~isfield(opts,'PatlakFastRegMode'); opts.PatlakFastRegMode='linear'; end %default to linear regression


N=size(Ct_mM,2); %number of time series
NTime=size(Ct_mM,1); %number of time points

CtModelFit_mM=nan(NTime,N);

switch model %process data using specified model/implementation
    case 'Patlak' %non-linear implementation
        
        PKP.vP=nan(1,N);
        PKP.PS_perMin=nan(1,N);
        
        fitParams0=[opts.init_vP opts.init_PS_perMin];
        
        for iSeries=1:N %loop through different time series (e.g. voxels or ROIs)
            if sum(isnan(Ct_mM(:,iSeries))) > 0; continue; end % skip concentration profiles containing one or more NaNs
            %% Fit Patlak model to estimated concentration curves
            [fitParams,resnorm,residual,exitflag,output]=lsqcurvefit(@fPatlak,fitParams0,[],Ct_mM(opts.NIgnore+1:end,iSeries)...
                ,[-inf -inf],[inf inf],optimoptions(@lsqcurvefit,'Display','off','TolFun',1e-6,'TolX',1e-6,'TypicalX',fitParams0));
            CtModelFit_mM(:,iSeries)=[nan(opts.NIgnore,1); fPatlak(fitParams,[])]; %get best fit curve
            PKP.vP(1,iSeries)=fitParams(1); %get best fit parameters
            PKP.PS_perMin(1,iSeries)=fitParams(2);
        end
        
        
    case 'PatlakLinear' %classic "Patlak Plot" implmentation
        
        PKP.vP=nan(1,N);
        PKP.PS_perMin=nan(1,N);
        
        PKP.PatlakX=nan(NTime,N);
        PKP.PatlakY=nan(NTime,N);
        
        PKP.PatlakYFit=nan(NTime,N);
        PKP.PatlakR2=nan(1,N);
        
        for n=1:NTime %calculate X and Y variables for linear Patlak model
            PKP.PatlakX(n,:) = ( (tRes_s/60)*(sum(Cp_AIF_mM(1:n-1,1),1) + 0.5*Cp_AIF_mM(n,1)) ) / Cp_AIF_mM(n,1); %integral of Cp / Cp
            PKP.PatlakY(n,:) = Ct_mM(n,:) / Cp_AIF_mM(n,1); % Ct / Cp
        end
        
        for iSeries=1:N %loop through series and do linear regression
            if sum(isnan(Ct_mM(:,iSeries))) > 0; continue; end % skip concentration profiles containing one or more NaNs
            [b,bint,r,rint,s] = regress(PKP.PatlakY(opts.NIgnore+1:end,iSeries),[ones(NTime-opts.NIgnore,1) PKP.PatlakX(opts.NIgnore+1:end,iSeries)]);
            PKP.vP(1,iSeries) = b(1);
            PKP.PS_perMin(1,iSeries) = b(2);
            PKP.PatlakYFit(:,iSeries) = [ones(NTime,1) PKP.PatlakX(:,iSeries)] * b;
            PKP.PatlakR2(1,iSeries) = s(1);
            CtModelFit_mM(:,iSeries)=DCEFunc_PKP2Conc(tRes_s,Cp_AIF_mM,struct('vP',b(1),'PS_perMin',b(2)),'Patlak',[]); %get best fit curve
            CtModelFit_mM(:,iSeries)=[nan(opts.NIgnore,1); CtModelFit_mM(opts.NIgnore+1:end,iSeries)];
        end
        
    case 'PatlakFast' %multiple linear regression implementation (recommended)
        %linear Patlak using Cp and Integral Cp as regressors - vP and K are the coefficients; uses matrix divide for speed
        
        PKP.vP=nan(1,N);
        PKP.PS_perMin=nan(1,N);
        intCp_AIF_mM_min=nan(NTime,1);
        beta=nan(2,N);
        
        for iFrame=1:NTime %calculate Cp integral to use as regressor
            intCp_AIF_mM_min(iFrame,1)=(0.5*Cp_AIF_mM(iFrame) + sum(Cp_AIF_mM(1:iFrame-1,1),1)) * (tRes_s/60);
        end
        
        reg=[Cp_AIF_mM(opts.NIgnore+1:end,:) intCp_AIF_mM_min(opts.NIgnore+1:end,:)]; % put both regressors into a matrix
        
        switch opts.PatlakFastRegMode
            case 'linear'
                beta(:,:) = reg \ Ct_mM(opts.NIgnore+1:end,:); % regression to calculate coefficients
            case 'robust'
                for iSeries=1:N
                    if sum(isnan(Ct_mM(:,iSeries))) > 0; continue; end % skip concentration profiles containing one or more NaNs
                    beta(:,iSeries)= robustfit(reg,Ct_mM(opts.NIgnore+1:end,iSeries),'huber',[],'off'); % robust regression to calculate coefficients
                end
        end
        
        PKP.vP(1,:)=beta(1,:);
        PKP.PS_perMin=beta(2,:);
        
        CtModelFit_mM = [ nan(opts.NIgnore,N) ; reg*beta ]; % calculate best fit concentration
        
    case 'PatlakFastExtra' % EXPERIMENTAL!
        %permits additional non-Patlak regressors in opts.extraRegs
        
        PKP.vP=nan(1,N);
        PKP.PS_perMin=nan(1,N);
        intCp_AIF_mM_min=nan(NTime,1);
        beta=nan(2+size(opts.extraRegs,2),N);
        
        for iFrame=1:NTime %(calculate integral to use as regressor)
            intCp_AIF_mM_min(iFrame,1)=(0.5*Cp_AIF_mM(iFrame) + sum(Cp_AIF_mM(1:iFrame-1,1),1)) * (tRes_s/60);
        end
        
        reg=[Cp_AIF_mM(opts.NIgnore+1:end,:) intCp_AIF_mM_min(opts.NIgnore+1:end,:) opts.extraRegs(opts.NIgnore+1:end,:)];
        beta(:,:) = reg \ Ct_mM(opts.NIgnore+1:end,:);
        
        PKP.vP(1,:)=beta(1,:);
        PKP.PS_perMin=beta(2,:);
        
        CtModelFit_mM = [ nan(opts.NIgnore,N) ; reg*beta ];
end

    function CtModel_mM=fPatlak(p,t) %function to minimise when using non-linear implementation
        testPKP.vP=p(1); testPKP.PS_perMin=p(2);
        CtModel_mM=DCEFunc_PKP2Conc(tRes_s,Cp_AIF_mM,testPKP,model,[]);
        CtModel_mM=CtModel_mM(opts.NIgnore+1:end);
    end

end


