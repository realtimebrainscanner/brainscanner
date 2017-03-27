function [A, beta, M, llh] = MARDv1(alphas, beta, Phi, T)
% ARD style M-SBL
tolerance = 1e-4;
maxIterations = 300;
alphaLowerBound = 1e-6;
alphaUpperBound = 1e6;

llh = zeros(1,maxIterations);
modelSize = size(Phi,2);
N = size(T,1);
steps = size(T,2);

PhiT = Phi';
PhiTPhi =  PhiT*Phi;


% %%% Gamma: equivalent to A^-1
% A = zeros(modelSize);
% for j=1:modelSize
%     A(j,j) = 1/A(j,j);
% end

activeSet = 1:modelSize;
zeroIndexes = zeros(1,modelSize);
% C=zeros(N,N);
for k=2:maxIterations
    %% Compute diagonal of posterior covariance
%     A = diag(alphas);
    
    SigmaInv = diag(alphas) + beta * PhiTPhi;
    SigmaInvU = chol(SigmaInv);
%     SigmaU = inv(SigmaInvU);
    SigmaU = SigmaInvU\eye(size(SigmaInvU,1));
%     Sigma = SigmaU*SigmaU';  %A^-1 = L^-1'*L^-1 = U^-1 * U^-1'
    diagSigma = sum(SigmaU.^2, 2); % MRA: We only need the diagonal of the covariance during iterations
    
    
    %% Compute posterior mean : mN = beta * (Sigma*(Phi'*t))
    M = beta * (SigmaInvU\(SigmaInvU'\(Phi'*T))); % MRA: We prefer to use cholesky decomp rather than Sigma directly.
    M2 = M.^2;
    
    Mi2Normalized = sum(M2, 2)/steps;
    
    % Make sure alpha has correct dimensions
    if size(alphas,2) > size(alphas,1)
        alphas = alphas';
    end
    
    gamma = 1 - alphas.*diagSigma;    
    
    alphas = gamma./(Mi2Normalized);
    
    %% Determine current active set
    activeIdx = alphas < alphaUpperBound;
    
    Phi(:,~activeIdx) = [];
    PhiTPhi(~activeIdx,:) = [];
    PhiTPhi(:,~activeIdx) = [];
    alphas(~activeIdx) = [];
    M(~activeIdx,:) = [];

    
    activeSet = activeSet(activeIdx);
    
    if isempty(M)
        disp('Pruned all weights');
        break;
    end   

    %% Noise variance update
    %||T - Phi*M||_F^2
%     diagSigma(~activeIdx) = [];    
%     frobSquared=trace((T-Phi*M)'*(T-Phi*M));
%     noiseVar = ((1/steps)*frobSquared)/(N - modelSize + sum(alphas.*diagSigma));
%     beta=1/noiseVar;                          %% This will cause problems
%     beta = max(1e-10, min(1e10, 1/noiseVar));   %% This gives bad results

%     TCInvT = 0;
%     for j=1:steps
%         b=L'\T(:,j);
%         TCInvT = TCInvT + b'*b;
%     end    
%     TCInvT_alt = sum(diag(LT'*LT)); % Same as above
    
    %%% What is this????
    PhiAInv = Phi*diag(sqrt(1./alphas));
    C = (1/beta)*eye(N) + PhiAInv*PhiAInv';   % Used in order to avoid numerical instability
    
%     C = (1/beta)*eye(N) + Phi*diag(1./alphas)*Phi';
    
    L=chol(C);
    logdetC = 2*sum(log(diag(L)));
    LT = L'\T;
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
Mtemp = zeros(modelSize, steps);
A = 1e6*ones(1,modelSize);
Mtemp(activeSet,:) = M;
A(activeSet) = alphas;

M=Mtemp;

% beta = 1/noiseVar;
llh = llh(1:k);

end