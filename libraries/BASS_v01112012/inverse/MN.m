function [X alpha beta logEv] = MN(Y,A,opts)
%==========================================================================
% Filename: MN.m (function).
% 
% Description: Minimum norm method - Bayesian formulation with EM-updates
%              of the hyperparameters of the noise variance and the source
%              variance.
%
% Input:        Y: Observations (size: Nc x Nt)
%               A: Mixing matrix (forward model)
%               opts:
%                   .alpha: Init hyperparameter of inverse source variance
%                   .beta: Init hyperparameter of the inverse noise
%                          variance.
%                   .maxIter: Maximum number of iterations.
%
% Output:       X: Regressions coefficients (sources)
%               alpha: Estimated hyperparameter of the inverse source
%                      variance.
%               beta: Estimated hyperparameter of the inverse noise
%                     variance.
%               logEv: Log-evidence value
%
% Special remarks:
%               Nc: Number of channels
%               Nd: Number of coefficients/sources
%               Nt: Number of samples
%
% Example:      [X alpha beta logEv] = MN(Y,A)
%               [X alpha beta logEv] = MN(Y,A,opts)
%
% History:
%   - Created:  02/05/2011
%   - Modified:
%
% Copyright (C) Carsten Stahlhut, DTU Informatics 2011
%==========================================================================

if nargin<3
    opts = [];
end

% Check for optional specification of parameters
if isfield(opts,'alpha') alpha = opts.alpha; else alpha = 1; end
if isfield(opts,'beta') beta = opts.beta; else beta = 1/10; end
if isfield(opts,'maxIter') maxIter = opts.maxIter; else maxIter = 1000; end

%Pre-compute
[Nc,Nd] = size(A);
Nt = size(Y,2);
Ic = speye(Nc);
AAt = A*A';
logEv = zeros(1,maxIter);
YYt = Y*Y';
SIGMAfulldiag = zeros(Nd,Nt);

% Now start iterative estimation
for ite=1:maxIter
    %----------------
    %Update X
    %----------------
    SigmaXinvalpha = pinv(AAt + alpha/beta*Ic);
    X = A'*SigmaXinvalpha*Y;
    
    %----------------
    %Update alpha
    %----------------
    for i=1:Nd
        SIGMAfulldiag(i,1) = 1/alpha-1/alpha*...
                A(:,i)'*SigmaXinvalpha*A(:,i);
    end
    SIGMAfulldiag = repmat(SIGMAfulldiag(:,1),[1 Nt]);

    v1 = (SIGMAfulldiag + X.^2)/2;
    alpha = Nt*Nd/(2*sum(v1(:)));
%     invalpha = sum(2*sum(v1,2))/(Nt*Nd);    %Changed 25/2-09
%     alpha = 1./invalpha;

    %----------------
    %Update beta
    %----------------
    Yhat = A*X;
    AinvDAt = AAt/alpha;
    v2 = 1/2*diag( YYt - 2*Yhat*Y' +...
        Nt*(AinvDAt - AinvDAt * alpha*SigmaXinvalpha * AinvDAt) + Yhat*Yhat');
    beta = Nt*Nc/(2*sum(v2));           %Expected value of beta.
    
    %Calculate log-evidence
%     logEv_temp = logEvMN(X,SigmaXinvalpha);
%     logEv(ite) = logEv_temp;
    logEv(ite) = -Nt*Nc/2*log(2*pi) +Nt/2*sum(log(svd(SigmaXinvalpha*alpha))) -...
        1/2*trace(SigmaXinvalpha*YYt);
end
ite
%EOF
end
