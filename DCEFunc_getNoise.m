function [tNoise,tSNR,tSignal]=DCEFunc_getNoise(signal2D,noiseFrames)
%signal: array of DCE signals including pre-contrast;  each column
%represents a different time series (e.g. voxel, ROI)
%noiseFrames: indices of time frames to use for calculation
%tNoise2D and tSNR2D: row vector of temporal noise and tSNR values
%calculate as std and mean/std following detrending


NSignals=size(signal2D,2);

tNoise=nan(1,NSignals);
tSNR=nan(1,NSignals);
tSignal=nan(1,NSignals);

%%calculate mean signal over noiseFrames
tSignal=mean(signal2D(noiseFrames,:),1);

%%detrend
signalDT2D=detrend(signal2D);

%%calculate noise properties
tNoise=std(signalDT2D(:,:),1);
tSNR=tSignal./tNoise;

end
