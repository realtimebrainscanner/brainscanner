function Res = calcObjectiveMetrics(M,A,S,Sref,opts)
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
%                   - EMD: Earth Mover Distance
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
%   - Created:  07/10/2012
%
% Special remarks:
%               [1]: Grova et al. (2006): Note however that they in this
%               paper only define it for one time sample
%               Code based on the previous calcMetrics.m and
%               calcMetrics_UCSF.m files
%
% Copyright (C) Carsten Stahlhut, DTU Informatics 2012
%==========================================================================

if nargin<5
    opts = [];
end

% Res.NameMethod = fnameMat;
%         Q = load(fnameMat);
try Res.LogEv = opts.LogEv; catch Res.LogEv = NaN; end

if isfield(opts,'tselect'), tselect = opts.tselect; else tselect = 1:size(S,2); end

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
try
    Res.MSE = calcMSE(S,Sref);
    % Res.MSE_tselect = calcMSE(S,Sref,11);
    Res.MSE_tselect = calcMSE(S,Sref,tselect);  %CS 31052009
    disp('MSE calculated.')
    
    % Degree of focalization
    % disp('if more than one source - DF not calculated correctly')
    iSim = single(sum(Sref.^2,2)>0);        %Note if more than one source all a regarded as one single source - this should be changed
    Res.DF = calcDF(S,Sref,iSim,tselect);
    disp('DF calculated.')
    
    % Energy
    [E Eref] = calcNormEnergy(S,Sref,tselect);
    Res.E = E;
    Res.Eref = Eref;
    disp('Energy calculated.')
    
    %
    Res.SSE = calcSSE(S,Sref);
    disp('SSE calculated.')
    
    %     [CC list_active] = calcCorrCoef(Sref,S,struct('choice','all'));
    %     Res.CorrCoeff.CC = CC;
    %     Res.CorrCoeff.list_active = list_active;
catch
    
end


try
    vert = opts.vert;
    vert0 = opts.vert_ref_center;
catch
    vert = NaN;
    vert0 = NaN;
end

try
    [Aprime HR FR hits false_alarms] = calcAprime(S,vert,vert0);
    Res.Aprime.Aprime = Aprime;
    Res.Aprime.HR = HR;
    Res.Aprime.FR = FR;
    Res.Aprime.hits = hits;
    Res.Aprime.false_alarms = false_alarms;
    Res.Aprime.vert0 = vert0;
    disp('Aprime calculated.')
catch
    
end

try
    % Receiver operating curves
    ROC = calcROC(E,Eref);
    Res.ROC = ROC;
    
    % Area under the ROC curve
    AUC = calcAUCROC(ROC.tprate,ROC.fprate);
    Res.AUC = AUC;
    disp('AUC calculated.')
catch
    
end



%% Earth Mover Distance Measure

try
    vert = opts.vert;
    vert0 = opts.vert_ref;
    
    % Distance measure - Euclidean distance between dipoles - this is the
    % cost of moving from one to anouther.
    if isfield(opts,'distMatrix')
        load(opts.distMatrix);
        %      distMatrix = opts.distMatrix;
    else
        distMatrix = zeros(Nd,Nd);
        for i=1:Nd
            distMatrix(:,i) = sum(bsxfun(@minus,vert0,vert(i,:)).^2,2);
        end
    end
    
    [Nd Nt] = size(Sref);
    
    EMDtime = zeros(1,Nt);
    for t=1:Nt
        w1 = Sref(:,t)';
        iw1 = find(w1~=0);
        w2 = S(:,t)';
        iw2 = find(w2~=0);
        distMatrix_sub = distMatrix(iw1,:);
        distMatrix_sub = distMatrix_sub(:,iw2);
        [e,Flow]=emd_mex(w1(iw1),w2(iw2),distMatrix_sub);
        EMDtime(t) = e;
    end
    
    Res.EMD.total = sum(EMDtime);
    Res.EMD.EMDtime = EMDtime;
    
    disp('EMD calculated.')
    
catch
    
end

%         w1=[0.4, 0.2, 0.2, 0.1, 0.1];
%         w2=[0.6, 0.2, 0.1];
%         C= [ 3, 5, 2;
%              0, 2, 5;
%              1, 1, 3;
%              8, 4, 3;
%              7, 6, 5 ];
%
%         [e,Flow]=emd_mex(w1,w2,C)
%

Res = orderfields(Res);      %alphabetic ordering