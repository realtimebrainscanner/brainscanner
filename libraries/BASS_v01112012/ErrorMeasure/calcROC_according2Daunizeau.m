function ROC = calcROC(STrue,S,opts)
%==========================================================================
% Filename: calcROC.m (function).
% 
% Description:  Calculates the receiver operating characteristic (ROC)
%               curve.
%
% Usage:        ROC = calcROC(STrue,S,opts)
%
% Input:        STrue: Simulated dipole activities, Nd x Nt (dipole x time)
%               S:     Estimated dipole activities, Nd x Nt (dipole x time)
%
% Output:       ROC
%                   .tprate: True positive rate.
%                   .fprate: False positive rate
%
% History:
%   - Created:  30/05/2008
%   - Modified:  
%
% Special remarks:
%               [1]: J. Daunizeau and K.J. Friston, "A mesostate-space
%                    model for EEG and MEG", NeuroImage 38, 67-81, 2007.
%
% Copyright (C) Carsten Stahlhut, DTU IMM 2008
%==========================================================================

[Nd, Nt] = size(STrue);

STrue = STrue/max(abs(STrue(:)));
S = S/max(abs(S(:)));

E1 = sum(STrue.^2,2);
E1 = E1/max(E1);

E2 = sum(S.^2,2);
E2 = E2/max(E2);

[Ethres, iEthres] = sort(E2,'descend');


TrueClass = E1 > 0;
P = sum(TrueClass);       %Number of total positive (true active sources)
N = Nd-P;                   %Number of total negative (true inactive sources)


TP = 0;
FP = 0;
tprate = nan*ones(Nd,1);
fprate = tprate;

for i=1:length(Ethres)
    if TrueClass(iEthres(i))
        TP = TP+1;              %True positive
    else
        FP = FP+1;              %False positive
    end
    tprate(i) = TP/P;           %True positive rate
    fprate(i) = FP/N;           %False positive rate
end

ROC.tprate = tprate;
ROC.fprate = fprate;
ROC.Ethres = Ethres;
ROC.iEthres = iEthres;

figure
stairs(fprate,tprate)
% xlabel('False positive rate')
% ylabel('True positive rate')
xlabel('1-Specificity')         %Specificity = 1-fprate --> 1-Specificity = fprate
ylabel('Sensitivity')           %Sensitivity = tprate
grid on
title('ROC')
