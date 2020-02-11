function enhancementPct=DCEFunc_Sig2Enh(signal,baselineScansIdx)
%calculate signal enhancement
%INPUT:
%signal: array of DCE signals including pre-contrast values; each column
%represents a different time series (e.g. voxel, ROI)
%baselineScans: a row vector of indices corresponding to precontrast scans;
%the mean of these scans is used as the baseline when calculating signal
%enhancement
%OUTPUT:
%enhancementPct: array of percentage enhancement values for each time point including pre-contrast measurements; each column
%represents a different time series (e.g. voxel, ROI)

NSignals=size(signal,2); %number of series
enhancementPct=nan(size(signal));

baselineSignal=mean(signal(1:baselineScansIdx,:),1); %calculate baseline signal for each series

for iSignal=1:NSignals %loop through series
    enhancementPct(:,iSignal) = 100*((signal(:,iSignal) - baselineSignal(1,iSignal)) / baselineSignal(1,iSignal)); %calculate enhancement relative to baseline signal
end

end