function Res = calcMetrics(M,A,S,Sref,opts)
%==========================================================================
% Filename: calcMetrics.m (function).
% 
% Description:  Calculates a number of error measures:
%                   - VE of M
%                   - GOF on M
%                   - LE of S: Need to specify opts.vert and
%                              opts.vert_ref_center
%                   - MSE of S
%                   - DF of S
%                   - ROC: Receiver operating curve
%                   - AUC: Area under the ROC curve
%
% Usage:        Res = calcMetrics(M,A,S,Sref,opts)
%
% Input:        M: Observations
%               A: Dictionary used to estimate S with.
%               Sref: Simulated dipole activities, Nd x Nt (dipole x time)
%               S:     Estimated dipole activities, Nd x Nt (dipole x time)
%               opts:
%                   .vert: vertices corresponding to S (Nd x 3
%                   .vert_ref_center: Vertex point of center of activity (1 x 3)
%
% Output:       Res: A set of metrics according to the description above.
%
% History:
%   - Created:  21/04/2009
%   - Modified: 22/05/2009: Calculation of LE is now included.
%               31/05/2009: If opts.tselect is not specified tselect is
%               dependent on size of M if smaller than time samples
%               tselect=1 else tselect=11.
%               31/05/2009: AUC included
%
% Special remarks:
%               [1]: Grova et al. (2006): Note however that they in this
%               paper only define it for one time sample
%
% Copyright (C) Carsten Stahlhut, DTU Informatics 2009
%==========================================================================



% Res.NameMethod = fnameMat;
%         Q = load(fnameMat);
try Res.LogEv = opts.LogEv; catch Res.LogEv = NaN; end
try 
    tselect = opts.tselect;
catch 
    if size(M,2)>=11
        tselect = 11;
    else
        tselect = 1;
    end
end

% Variance explained at sensor level
Mrecon = A*S;
Res.VE = calcVE(M,Mrecon);

% Goodness of fit on M
Res.GOF_M = calcGOF(M,Mrecon);

% Localization error in mm
try
    vert = opts.vert;
    vert_ref_center = opts.vert_ref_center;
    Res.LE = calcLE(S,vert,vert_ref_center);
catch
%     disp('Localization error not calculated.')
%     vertFname = ''
%     vert_ref_Fname = '';
%     load(vertFname)
%     load(vert_ref_Fname)
end


% Mean square error on estimated sources

Res.MSE = calcMSE(S,Sref);
% Res.MSE_tselect = calcMSE(S,Sref,11);
Res.MSE_tselect = calcMSE(S,Sref,tselect);  %CS 31052009

% Degree of focalization
% disp('if more than one source - DF not calculated correctly')
iSim = single(sum(Sref.^2,2)>0);        %Note if more than one source all a regarded as one single source - this should be changed
Res.DF = calcDF(S,Sref,iSim,tselect);

% Energy
[E Eref] = calcNormEnergy(S,Sref,tselect);
Res.E = E;
Res.Eref = Eref;

%
Res.SSE = calcSSE(S,Sref);

% Receiver operating curves
ROC = calcROC(E,Eref);
Res.ROC = ROC;

% Area under the ROC curve
AUC = calcAUCROC(ROC.tprate,ROC.fprate);
Res.AUC = AUC;

Res = orderfields(Res);      %alphabetic ordering