function [V,m,X,FreeEnergy,f,M,beta,k] = teVGGD(A,Y,gamma,opts)
%
% Description:      Solves the augmented inverse problem: Y = A * s * X
%                   using the time-expanded Variational Garrote (teVG).
%
% Input:            Y:  Data matrix of size KxT, where T is the number of
%                       measurements.
%                   A:  Design matrix/forward model of size KxN
%                   opts:   see 'Settings'
%                   gamma:  Sparsity level
%
% Output:           X:  Feature matrix of size NxT
%                   m:  Variational mean/expectation of the state 's'.
%                       Can be thresholded at 0.5 to only keep sources with
%                       high evidence.
%                   V:  The solution matrix; V = X.*repmat(m,1,T);
%                   FreeEnergy: The solution's variational free energy
%--------------------------References--------------------------------------
% The Variational Garrote was originally presented in
% Kappen, H. (2011). The Variational Garrote. arXiv Preprint
% arXiv:1109.0486. Retrieved from http://arxiv.org/abs/1109.0486
% and
% Kappen, H.J., & Gómez, V. (2014). The Variational Garrote. Machine
% Learning, 96(3), 269–294. doi:10.1007/s10994-013-5427-7
%
% Time-expanded Variational Garrote reference
% Hansen, S.T., Stahlhut, C. & Hansen, L.K. (2013). Expansion
% of the Variational Garrote to a Multiple Measurement Vectors Model. In
% Twelfth Scandinavian Conference on Artificial Intelligence. ed. / M.
% Jaeger. IOS Press, (pp. 105–114).
%
% teVG with gradient descent
% Hansen, S.T., & Hansen, L.K. (2015). EEG source reconstruction
% performance as a function of skull conductance contrast. In 2015 IEEE
% International Conference on Acoustics, Speech and Signal Processing
% (ICASSP) (pp. 827–831). IEEE. doi:10.1109/ICASSP.2015.7178085
%-----------------------------Author---------------------------------------
% Sofie Therese Hansen, DTU Compute
% March 2016
% -------------------------------------------------------------------------

% Settings:
try max_iter = opts.max_iter; catch; max_iter = 100; end; % Maximum number of iterations
try tol = opts.tol; catch; tol = 1e-5; end; % Convergence criterium
try eta0 = opts.eta0; catch; eta0=1e-3; end; % Learning rate for gradient descent
try updEta = opts.updEta; catch; updEta = 1; end; % Update of learning rate
try etastart = opts.etastart; catch; etastart=200; end; % When to reset learning rate
try run_prune = opts.run_prune; catch; run_prune=1; end; % When to to prune M
try prune = opts.prune; catch; prune=1e-2; end; % When to to prune M
try pn = opts.pnorm; catch; pn=0.8; end; % When to to prune M
try m0 = opts.m0; catch; m0 = zeros(size(A,2),1); end; % initializa m

% Initialize
[K,N] = size(A);
T = size(Y,2);
%m0 = zeros(N,1);
m0 = max(m0,sqrt(eps));
m0 = min(m0,1-sqrt(eps));
m = m0;
A = center(A);
Y = center(Y);
chi_nn =  1/K*sum(A.^2)';
k = 0;term = 0;f = zeros(max_iter,1);
eta = eta0;
M=NaN(N,100);
keep_list = 1:N;
Nold = N;
% Iterate solution until convergence
while k <max_iter && term==0
     k = k+1;

    if k>4 && (min(mprun) < prune ) && run_prune == 1;
        %keyboard
        index = find(mprun > prune);
        A = A(:,index);            % corresponding columns in Phi
        keep_list = keep_list(index);
        m = m(index);
        N = size(m,1);
        chi_nn = chi_nn(index);
        if isempty(keep_list);X=X(keep_list,:);break;end
    end;
    C = eye(K)+1/K*A*spdiags(m./(1-m)./chi_nn,0,N,N)*A';
    yhat = C\Y;
    yhaty = yhat.*Y;
    beta = T*K/sum(yhaty(:));
    lambda = beta*yhat;
    X = bsxfun(@times, 1./((K*beta)*(1-m).*chi_nn), (A'*lambda));
    z = Y-1/beta*lambda;
    X2 = X.^2;
    m = max(m,sqrt(eps));
    m = min(m,1-sqrt(eps));
    dFdm = beta*K/2*sum(repmat((1-2*m).*chi_nn,1,T).*X2,2)-gamma+log(m./(1-m))-sum(X.*(A'*lambda),2);
    m = m -eta*dFdm;
    m = max(m,sqrt(eps));
    m = min(m,1-sqrt(eps));
    % *** Update hyperparameters ***

    f(k) = -T*K/2*log(beta/(2*pi))+beta/2*sum(sum((z-Y).^2))...
        +K*beta/2*sum(sum(X2.*repmat(chi_nn.*m.*(1-m),1,T)))...
        -gamma*sum(m)+N*log(1+exp(gamma))...
        +sum(m.*log(m)+(1-m).*log(1-m))...
        +sum(sum(lambda.*(z-A*(repmat(m,1,T).*X)))); % Free energy
    
    %Check for convergence
    if k>5
        if f(k)==f(k-5) && f(k)==f(k-1);
            term = 1;
        elseif k>10
            if abs(f(k)-f(k-5))<tol && abs(f(k)-f(k-1))<tol
                term = 1;
            end
        end
    end
    % Update learning rate
    if term==0 && updEta==1 && k>2
        if f(k)>f(k-1)
            eta = eta/2;
        elseif f(k)<f(k-1)
            eta = eta*1.1;
        end
    end
    if k>etastart && k<etastart+350
        eta = 1e-3;
    end
    m2 = m.^2;
    mprun = (m2).^(1-pn/2);
    %Xprun = (sum(X2,2)/T).^(1-pn/2);
   %M(:,k)=m; 
end
mold=m;
m=zeros(Nold,1);
m(keep_list)=mold;
Xold = X;
X=zeros(Nold,T);
X(keep_list,:)= Xold;

FreeEnergy = f(k);
V = X.*repmat(m,1,T);