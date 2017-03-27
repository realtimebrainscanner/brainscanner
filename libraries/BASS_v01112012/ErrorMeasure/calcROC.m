function ROC = calcROC(X,Xref,opts)
%==========================================================================
% Filename: calcROC.m (function).
% 
% Description:  Calculates the receiver operating characteristic (ROC)
%               curve.
%
% Usage:        ROC = calcROC(Xref,X,kappa,opts)
%
% Input:        X:    Estimated dipole activities, Nd x Nt (dipole x time)
%               Xref: True parameters
%               opts:
%                 .kappa: Vector with different statistical thresholds
%                         kappa. Default is the resolution of X, i.e. X is
%                         sorted in descending order and the
%                 .dataNorm: Normalize data if not already done. Default is
%                 0 (no normalization).
%
% Output:       ROC
%                   .tprate: True positive rate. (sensitivity)
%                   .fprate: False positive rate (1-specificity)
%                   .sensitivity: TP/(TP + FN)
%                   .specificity: TN/(TN + FP)
%                   .posCount: TP + FP
%                   .negCount: TN + FN
%
% History:
%   - Created:  30/05/2008
%   - Modified: 13/04/2009: Order of input arguments has been changed to
%                           follow the same order as the rest of the
%                           calc-functions.
%               09/03/2010: Xref was not supported for Nt>1, which is now
%                           supported.
%
% Special remarks:
%
%               [1]: J. Daunizeau and K.J. Friston, "A mesostate-space
%                    model for EEG and MEG", NeuroImage 38, 67-81, 2007.
%               [2]: C. Grova et al., "Evaluation of EEG localization
%                    methods using realistic simulations of interictal
%                    spikes", NeuroImage 26, 734-753, 2006.
%
% Copyright (C) Carsten Stahlhut, DTU IMM 2008
%==========================================================================

try dataNorm = opts.dataNorm; catch dataNorm = 0; end;
try flag1 = opts.show_fig; catch flag1 = 0; end;

[Nd, Nt] = size(Xref);

if dataNorm
    Xref = Xref/max(abs(Xref(:)));
    X = X/max(abs(X(:)));

    E1 = sum(Xref.^2,2);
    E1 = E1/max(E1);

    E2 = sum(X.^2,2);
    E2 = E2/max(E2);
    
    Xref = E1;
    X = E2;
    Eref = E1;
    
else
    Eref = sum(Xref.^2,2);
    E2 = sum(X.^2,2);
end

% try kappa = opts.kappa; catch [kappa, iEthres] = sort(X,'descend'); end


[kappa, ikappa] = sort(E2,'descend');
TrueClass = Eref > 0;
% [kappa, ikappa] = sort(X,'descend');
% TrueClass = Xref > 0;

P = sum(TrueClass);       %Number of total positive (true active sources)
N = Nd-P;                 %Number of total negative (true inactive sources)

TP = 0;
FP = 0;
tprate = nan*ones(Nd,1);
fprate = tprate;

for i=1:length(kappa)
    if TrueClass(ikappa(i))
        TP = TP+1;              %True positive
    else
        FP = FP+1;              %False positive
    end
    tprate(i) = TP/P;           %True positive rate
    fprate(i) = FP/N;           %False positive rate
end

ROC.tprate = tprate;
ROC.fprate = fprate;
ROC.sensitivity = tprate;
ROC.specificity = 1-fprate;
ROC.posCount = P;
ROC.negCount = N;
ROC.kappa = kappa;
ROC.ikappa = ikappa;

if flag1
    figure
    stairs(fprate,tprate)
    % xlabel('False positive rate')
    % ylabel('True positive rate')
    xlabel('1-Specificity')         %Specificity = 1-fprate --> 1-Specificity = fprate
    ylabel('Sensitivity')           %Sensitivity = tprate
    grid on
    title('ROC')
end