function [A, beta, M, llh] = MARDv2(alphas, beta, Phi, T)
% ARD style M-SBL
tolerance = 1e-4;
maxIterations = 300;
alphaLowerBound = 1e-6;
alphaUpperBound = 1e6;

llh = zeros(1,maxIterations);
modelSize = size(Phi,2);

steps = size(T,2);

activeSet = 1:modelSize;

% C=zeros(N,N);
for k=2:maxIterations
    %% Compute diagonal of posterior covariance
    PhiAInv = bsxfun(@rdivide, Phi, alphas');
    C = (1/beta*eye(size(Phi,1)) + PhiAInv*Phi');
    
    L = chol(C,'lower');
    R = L\PhiAInv;
    diagSigma = (1./alphas) - sum(R.^2,1)';
    
    
    %% Compute posterior mean : mN = beta * (Sigma*(Phi'*t))
    PhiTT = (Phi'*T);
    AInvPhiTT = bsxfun(@rdivide, PhiTT, alphas);
    
    M = beta*(AInvPhiTT - (R'*(R*PhiTT)));
    M2 = M.^2;

    Mi2Normalized = sum(M2, 2)/steps;    
    
    % Make sure alpha has correct dimensions
    if size(alphas,2) > size(alphas,1)
        alphas = alphas';
    end
    
    gamma = 1 - alphas.*diagSigma;    
    alphas = gamma./(Mi2Normalized);
    
    %% Determine current active set
    activeIdx = (alphas < alphaUpperBound);  % & (alphaLowerBound < alphas);
    
    Phi(:,~activeIdx) = [];
    alphas(~activeIdx) = [];

    activeSet = activeSet(activeIdx);
    
    if isempty(alphas)
%         disp('Pruned all weights');
        break;
    end   

    %% Noise variance update
    %||T - Phi*M||_F^2
%     diagSigma(~activeIdx) = [];    
%     frobSquared=trace((T-Phi*M)'*(T-Phi*M));
%     noiseVar = ((1/steps)*frobSquared)/(N - modelSize + sum(alphas.*diagSigma));
%     beta=1/noiseVar;                          %% This will cause problems
%     beta = max(1e-10, min(1e10, 1/noiseVar));   %% This gives bad results
   
    
    logdetC = 2*sum(log(diag(L)));
    LT = L\T;
    TCInvT = LT(:)'*LT(:);     % sum(LT.^2);
    
    llh(k) = -0.5*(steps*logdetC + TCInvT);   %7.85
    
    if abs(llh(k)-llh(k-1)) < tolerance*abs(llh(k-1));
%         SigmaInv = A + beta * (Phi'*Phi);
%         mN = beta * (SigmaInv\(Phi'*t));
%         disp('Converged at');
%         k
        break;
    end
     
    
end
M(~activeIdx,:) = [];

Mtemp = zeros(modelSize, steps);
A = 1e6*ones(1,modelSize);
Mtemp(activeSet,:) = M;
A(activeSet) = alphas;

M=Mtemp;

% beta = 1/noiseVar;
llh = llh(1:k);

end