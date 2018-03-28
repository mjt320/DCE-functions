function [PKP, CtModelFit_mM]=DCEFunc_fitModel(tRes_s,Ct_mM,Cp_AIF_mM,model,opts)
% Fit concentration curve using model to generate pharmacokinetic parameters
% Output:
% PKP = struct containing PK parameters (vP, vE, PS_perMin,
% FP_mlPer100gPerMin); each parameter is a row vector, with each value
% corresponding to a column in in Ct_mM
% CtModelFit_mM = best-fit concentration curve
% Input:
% t_s = time in seconds corresponding to each data point
% Ct_mM: array of tissue concentrations in mM; each column corresponds to a
% time series
% Cp_AIF_mM = column vector giving AIF plasma concentration in mM
% model = string to specify model ('Patlak' or '2CXM' or 'PatlakLinear').
% opts = struct containing options (NIgnore = number of points to ignore in
% cost function)


N=size(Ct_mM,2);
NTime=size(Ct_mM,1);

CtModelFit_mM=nan(NTime,N);

switch model
    case 'Patlak'
        
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
        
        
    case 'PatlakLinear'
        
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
        
    case 'PatlakFast'
        %linear Patlak using Cp and Integral Cp as regressors - vP and K are the regressors
        %uses matrix divide so should be fast
        
        PKP.vP=nan(1,N);
        PKP.PS_perMin=nan(1,N);
        intCp_AIF_mM_min=nan(NTime,1);
        beta=nan(2,N);
        
        for iFrame=1:NTime %(calculate integral to use as regressor)
            intCp_AIF_mM_min(iFrame,1)=(0.5*Cp_AIF_mM(iFrame) + sum(Cp_AIF_mM(1:iFrame-1,1),1)) * (tRes_s/60);
        end
        
        reg=[Cp_AIF_mM intCp_AIF_mM_min];
        beta(:,:) = reg \ Ct_mM;
        
        PKP.vP(1,:)=beta(1,:);
        PKP.PS_perMin=beta(2,:);
        
end


    function CtModel_mM=fPatlak(p,t)
        testPKP.vP=p(1); testPKP.PS_perMin=p(2);
        CtModel_mM=DCEFunc_PKP2Conc(tRes_s,Cp_AIF_mM,testPKP,model,[]);
        CtModel_mM=CtModel_mM(opts.NIgnore+1:end);
    end

end


