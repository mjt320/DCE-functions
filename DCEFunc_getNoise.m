function [tNoise,tSNR,tSignal]=DCEFunc_getNoise(signal2D,noiseFrames)
%WORK IN PROGRESS
%Function to calculate temporal noise level of time series (for QA use)
%Calculated as std and mean/std following linear detrending of signal
%INPUT:
%signal2D: array of signals including baseline signal; each column represents a time series (e.g. voxel, ROI)
%noiseFrames: indices of time frames to use for noise calculation
%OUTPUT:
%tNoise and tSNR: row vectors of temporal noise and tSNR values for each
%time series

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
