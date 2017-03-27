function [ROCmeasure, ConfusMat] = calcROCmeasure(TP,FP,Ptotal,N)
%==========================================================================
% Filename: calcROC.m (function).
% 
% Description:  Calculates the receiver operating characteristic (ROC)
%               curve.
%
% Usage:        [ROCmeasure, ConfusMat] = calcROCmeasure(TP,FP,Ptotal,N)
%
% Input:        TP:     True Positive
%               FP:     False Positive
%               Ptotal: Total number of positive
%               N:      Total number = Ptotal + Ntotal, with Ntotal as
%                       total number of negative.
%
% Output:       ROCmeasure
%                 .TPR = TP./Ptotal:              %True positive rate
%                 .FPR = FP./Ntotal:              %False positive rate
%                 .ACC = (TP+TN)/N:              %Accuracy
%                 .SPC = 1-FPR:                   %Specificity
%                 .PPV = TP./(TP+FP):             %Positive predictive value
%                 .NPV = TN./(TN+FN):             %False predictive value
%                 .FDR = FP./(FP+TP):             %False discovery rate
%                 .MCC = (TP.*TN-FP.*FN)./...     %Matthews Correlation
%                        sqrt(Ptotal.*Ntotal.*... %Coefficient
%                       (TP+FN).*(FP+TN))
%               ConfusMat: Confusion Matrix including TP, FP, FN, and TN.
%
%
% History:
%   - Created:  21/10/2008
%   - Modified:  
%
% Special remarks:
%
% Copyright (C) Carsten Stahlhut, DTU IMM 2008
%==========================================================================

Ntotal = N-Ptotal;
FN = Ptotal-TP;                 %False Negative
TN = Ntotal-FP;                 %True Negative

TPR = TP./Ptotal;              %True positive rate, sensitivity
TPR(isnan(TPR)) = 0;
FPR = FP./Ntotal;              %False positive rate
ACC = (TP+TN)/N;               %Accuracy
SPC = 1-FPR;                   %Specificity
PPV = TP./(TP+FP);             %Positive predictive value
NPV = TN./(TN+FN);             %False predictive value
FDR = FP./(FP+TP);             %False discovery rate
MCC = (TP.*TN-FP.*FN)./...
  sqrt(Ptotal.*Ntotal.*(TP+FN).*(FP+TN)); %Matthews Correlation Coefficient


ROCmeasure.TPR = TPR;
ROCmeasure.FPR = FPR;
ROCmeasure.ACC = ACC;
ROCmeasure.SPC = SPC;
ROCmeasure.PPV = PPV;
ROCmeasure.NPV = NPV;
ROCmeasure.FDR = FDR;
ROCmeasure.MCC = MCC;

ConfusMat.TP = TP;                        %True Positive
ConfusMat.FP = FP;                        %False Positive
ConfusMat.FN = FN;                        %False Negative
ConfusMat.TN = TN;                        %True Negative
