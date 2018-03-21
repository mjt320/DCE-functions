function enhancementPct=DCEFunc_Sig2Enh(signal,baselineScansIdx)
%signal: array of DCE signals including pre-contrast; each column
%represents a different time series (e.g. voxel, ROI)
%baselineScans: row vector of indices corresponding to precontrast scans

NSignals=size(signal,2);
enhancementPct=nan(size(signal));

baselineSignal=mean(signal(baselineScansIdx,:),1);

for iSignal=1:NSignals
    enhancementPct(:,iSignal) = 100*((signal(:,iSignal) - baselineSignal(1,iSignal)) / baselineSignal(1,iSignal));
end

end