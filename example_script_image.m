%%Example script demonstrating how to fit a 4D DCE-MRI image using the Patlak model

clear; close all;

%% Set up paths (USER DEPENDENT!)
procRoot='/ISIS/proc5';
addpath('/usr/local/spm/spm12');
addpath([procRoot '/software/GitHub/InvSVDs-processing']);
addpath([procRoot '/software/GitHub/utilities']);
addpath([procRoot '/software/GitHub/DCE-functions']);

%% Set acquisition parameters
tRes_s=39.62; %temporal resolution
NFrames=32; %number of time frames
FA_deg=15; %SPGR flip angle
TR_s=0.0034; %SPGR TR
TE_s=0.0017; %SPGR TE

%% Set options
Hct=0.45; % hematocrit
r1_permMperS=5.0; % relaxivity (/mMol/s) [gadovist at 3T]
r2s_permMperS=7.1; % relaxivity (/mMol/s) [gadovist at 3T]
DCENFramesBase=3; % time points for calculating baseline signal
NIgnore=7; % time points to exclude in fitting (e.g. to reduce CBF effects)

%% Read 4D DCE-MRI image (which should have been previously realigned)
SI4DHdr=spm_vol(['./rDCE.nii']); %read 4D DCE-MRI image
[SI4D,temp]=spm_read_vols(SI4DHdr);
dims=[size(SI4D,1) size(SI4D,2) size(SI4D,3)]; %get dimensions of image

%% Convert 4D to a 2D array of time series (each column corresponds to a voxel) for ease of processing
SI2D=DCEFunc_reshape(SI4D);

%% Load relative flip angle (k) map and T1 map (which should be co-registered to DCE-MRI space)
[T1Map_s,temp]=spm_read_vols(spm_vol(['./rT1.nii'])); %load T1 map
[kMap,temp]=spm_read_vols(spm_vol(['./rk.nii'])); %load relative flip angle (k) map
FAMap_deg=kMap * FA_deg; %calculate flip angle for each voxel
%if a flip angle map isn't available, the nominal flip angle can be used: FAMap_deg=ones(dims)*FA_deg

%% Calculate enhancement and concentrations in each voxel
enhancement2DPct=DCEFunc_Sig2Enh(SI2D,1:DCENFramesBase); %generates 1 column of enhancements per voxel
conc2DmM=DCEFunc_Enh2Conc_SPGR(enhancement2DPct,T1Map_s(:).',TR_s,TE_s,FAMap_deg(:).',r1_permMperS,r2s_permMperS); %generates 1 column of concentrations per voxel
conc4DmM=DCEFunc_reshape(conc2DmM,dims); %reshape concentration map to 4D

%% Get VIF
AIFMaskData=measure4D(conc4DmM,'./VIF.nii'); %measure time-concentration profile using VIF ROI mask image
Cp_AIF_mM=AIFMaskData.mean/(1-Hct); %correct for Hct to get plasma concentration

%% Fit Patlak model
[PatlakParams, CtModelFit2DFast_mM]=DCEFunc_fitModel(tRes_s,conc2DmM,Cp_AIF_mM,'PatlakFast',struct('NIgnore',NIgnore)); %do the fitting

%% Write output images
SPMWrite4D(SI4DHdr,DCEFunc_reshape(CtModelFit2DFast_mM,dims),'.','PatlakModelFit',16); %write 4D image of model fit 
SPMWrite4D(SI4DHdr(1),DCEFunc_reshape(PatlakParams.vP,dims),'.','vP',16); %write 3D vP image
SPMWrite4D(SI4DHdr(1),DCEFunc_reshape(PatlakParams.PS_perMin,dims),'.','PSperMin',16); %write 3D PS image


