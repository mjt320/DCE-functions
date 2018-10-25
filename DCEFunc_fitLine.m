function slopeFit=DCEFunc_fitLine(tRes_s,signal,opts)
% WORK IN PROGRESS
% Fit signal curve using model to get slope (for QA use)
% Output:
% slopeFit = struct containing output parameters:
%   signalFit, slope_perS in /s, intercept (value at first measurement),
%   slope_PctPerS normalised to mean signal in %/s,
%   changePct: % change over acquisition (from first to last measurement) relative to mean signal
%   each parameter is a row vector, with each value corresponding to a column in signal
% Input:
% t_s = time resolution in seconds
% signal: array of signals; each column corresponds to a
% time series
% opts = struct containing options:


N=size(signal,2);
NTime=size(signal,1);


t=0:tRes_s:NTime*tRes_s; %defining time in this way means the intercept represents the first measurement

slopeFit.signalFit=nan(NTime,N);
slopeFit.slope_perS=nan(1,N);
slopeFit.intercept=nan(1,N);
slopeFit.slope_PctPerS=nan(1,N);
slopeFit.changePct=nan(1,N);


for n=1:N %loop through different time series (e.g. voxels or ROIs)
    if sum(isnan(signal(:,n))) > 0; continue; end % skip concentration profiles containing one or more NaNs
    
    %% Fit linear model
    
    meanSignal=mean(signal(:,n));
    
    [p,S]=polyfit(t,signal(:,n),1);
    slopeFit.signalFit(:,n)=polyval(p,t);
    slopeFit.slope_perS(1,n)=p(2);
    slopeFit.intercept(1,n)=p(1);
    slopeFit.slope_PctPerS=100*(p(2)/meanSignal);
    slopeFit.changePct=100* ( ( polyval(p,NTime*tRes_s) - polyval(p,0) ) / meanSignal );
    
end


end