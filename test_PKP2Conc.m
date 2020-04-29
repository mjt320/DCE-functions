clc;clear

NSeries=500;

tRes_s=2.3;
Cp_AIF_mM=[1:10].' + rand(10,1);
vPList=1*rand(1,NSeries);
PSList=1e-2*rand(1,NSeries);
vEList=1*rand(1,NSeries);
FPList=50*rand(1,NSeries);
model='2CXM';
opts=[];

Ct_mM=nan(10,NSeries);
IRF=nan(10,NSeries);
c_cp_mM=nan(10,NSeries);
c_ce_mM=nan(10,NSeries);

for nSeries=1:NSeries
    PKP.vP=vPList(nSeries);
    PKP.PS_perMin=PSList(nSeries);
    PKP.vE=vEList(nSeries);
    PKP.FP_mlPer100gPerMin=FPList(nSeries);
    
    [Ct_mM(:,nSeries), IRF(:,nSeries), c_cp_mM(:,nSeries), c_ce_mM(:,nSeries)] = DCEFunc_PKP2Conc_old(tRes_s,Cp_AIF_mM,PKP,model,opts);
end


PKP2.vP=vPList;
PKP2.PS_perMin=PSList;
PKP2.vE=vEList;
PKP2.FP_mlPer100gPerMin=FPList;
[Ct_mM_2, IRF_2, c_cp_mM_2, c_ce_mM_2] = DCEFunc_PKP2Conc(tRes_s,Cp_AIF_mM,PKP2,model,opts);

[Ct_mM Ct_mM_2];
[IRF IRF_2];
[c_cp_mM c_cp_mM_2];
[c_ce_mM c_ce_mM_2];

max(abs([Ct_mM(:)-Ct_mM_2(:)]))
max(abs([IRF(:)-IRF_2(:)]))
max(abs([c_cp_mM(:)-c_cp_mM_2(:)]))
max(abs([c_ce_mM(:)-c_ce_mM_2(:)]))