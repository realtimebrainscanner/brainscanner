function AUC = calcAUCROC(TPR,FPR)
%==========================================================================
% Filename: calcAUCROC.m (function).
% 
% Description:  Calculates the Area Under the receiver operating
%               characteristic (ROC) Curve (AUC).
%
% Usage:        AUC = calcAUCROC(TPR,FPR)
%
% Input:        TPR: True Positive Rate
%               FPR: False Positive Rate
%
% Output:       AUC: Area Under the ROC Curve.
%
% History:
%   - Created:  31/10/2008
%   - Modified:  
%
% Special remarks:
%
% Copyright (C) Ricardo Henao and Carsten Stahlhut, DTU Informatics 2009
%==========================================================================

n = length( TPR);
AUC = sum( ( FPR(2:n) - FPR(1:n-1) ).*( TPR(2:n)+TPR(1:n-1) ) )/2;