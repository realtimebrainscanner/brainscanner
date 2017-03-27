function [S W invgamma invalpha ite] = SBLwSTE(Y,L,Psi,Phi,Sigma_E,opts)
%==========================================================================
% Filename: SBLwSTE.m (function).
%
% Description: Sparse-Bayesian Learning with space-time events structure.
%              This model is an extension of the the Aquavit algorithm [1].
%              The SBLwSTE function can both handle spatial and temporal
%              basis functions. It is an implementation of the general
%              model
%
%                    Y = L * (Psi*W*Phi + S) + E   Sensor noise E
%
% Input:
%            Y: Observations                    [Nc x Nt]
%            A: Mixing matrix (forward model)   [Nc x Nd]
%          PSi: Spatial design matrix           [Nd x Nm]
%          Phi: Temporal design matrix          [Nk x Nt]
%      Sigma_E: Noise covariance matrix         [Nc x Nc]
%         opts:
%             .maxIter: Maximum number of interations
%             .init.invgamma: 1 (default)
%             .init.invalpha: 1 (default)
%             .init.beta: 1 (default)
%             .group_alpha: 1 (default)
%             .tol_Sigma_E: Relative tolerance on the eigenvalues in Sigma_E
%                           1e-6 (default)
%             .flag_upRule: 1 - EM
%                           2 - Mackay (default)
%                           3 - Wipf (not implemented yet)
%
% Output:      S: Current source 'noise'
%              W: Spatio-temporal maps (regression coefficients)
%          invgamma: Inverse precision parameters on sources.
%          invalpha: Inverse precision parameters on spatio-temporal events.
%            ite: Number of iterations performed
%
% Special remarks:
%   The spatio-temporal maps W is transposed compared to Aquavit [1]!
%
%           Nc: Number of channels
%           Nt: Number of samples
%           Nd: Number of sources
%           Nm: Number of spatial basis functions
%           Nk: Number of temporal basis functions
%
% References
%   [1]: Stahlhut, C.; Attias, H. T.; Wipf, D.; Hansen, L. K. & Nagarajan,
%        S. S.; "Sparse Spatio-Temporal Inference of Electromagnetic Brain
%        Sources", Machine Learning in Medical Imaging (MLMI) 2010,
%        Springer-Verlag Berlin Heidelberg, 2010, 157-164
%
%   [2]: Stahlhut, C.; Attias, H. T.; Wipf, D.; Hansen, L. K. & Nagarajan,
%        S. S.; "Probabilistic M/EEG source imaging from sparse spatio-
%        temporal evetn structure", Machine Learning and Interpretation in
%        Neuroimaging (MLINI'12), 2012, 1-6.
%
%  Author: Carsten Stahlhut
%
%  Copyright (c) DTU Informatics, 2012.
%==========================================================================


MIN_DMU         = 1e-12;

[Nc, Nt] = size(Y);
Nd = size(L,2);
Nk = size(Phi,1);
Nm = size(Psi,2);


if nargin<6
    opts = [];
end
%% Extract optional parameter specification
if isfield(opts,'maxIter'), maxIter = opts.maxIter; else maxIter = 100; end;
if isfield(opts,'flag_upRule'), flag_upRule = opts.flag_upRule; else flag_upRule = 2; end
if isfield(opts,'tol_Sigma_E'), tol_Sigma_E = opts.tol_Sigma_E; else tol_Sigma_E = 1e-6; end
if isfield(opts,'tol_svd'), tol_svd = opts.tol_svd; else tol_svd = 1e-6; end
if isfield(opts,'flag_OrthoPhis')
    flag_OrthoPhis = opts.flag_OrthoPhis;
else
    flag_OrthoPhis = false;
    %     Rphi = Phi*Phi';
    %     detRphi = det(Rphi);
    %     clear Rphi
    %     if abs(1-detRphi)<detTol, flag_OrthoPhis = true; else falg_OrthoPhis = false; end
    
end


% Get user specified init of hyperparameters
if ~isfield(opts,'init'), opts.init = []; end
if isfield(opts.init,'invgamma'), invgamma = opts.init.invgamma; else invgamma = ones(Nd,1); end
if isfield(opts.init,'invalpha'), invalpha = opts.init.invalpha; else invalpha = ones(Nm,Nk); end
if isfield(opts.init,'beta'), beta = opts.init.beta; else beta = 1; end
if isfield(opts.init,'W'), W = opts.init.W; else W = zeros(Nm,Nk);; end
%if isfield(opts,'group_alpha'), group_alpha = opts.group_alpha; else group_alpha = 1; end

% try flag_calcLogEv = opts.flag_calcLogEv; catch flag_calcLogEv = 0; end


% ----------- Extract optional parameter specification done ---------------



% if flag_OrthoPhis
%     [S W gamma ite] = SBLwSTE_ortho(Y,L,Phi,Sigma_E,opts);
% else
%     [S W gamma ite] = SBLwSTE_nonortho(Y,L,Phi,Sigma_E,opts);
% end


%% Reserved space for parameters

if size(invalpha,1)<Nm
    index_k = ones(Nm,1);                    %Do not use space on copy Z times D, since D_k = D.
else
    index_k = (1:Nm)';
end

SigmaWiAll_diag = W;
OmegaWiAll_diag = W;
SigmaS_diag = zeros(Nd,1);
S = zeros(Nd,Nt);
Res_w = Y;     %Residual
Res_s = Y;     %Residual

%Create waitbar
[handle_waitbar,timer_t0] = waitbar_create(mfilename);

logEv = zeros(1,maxIter);
Nactive = zeros(1,maxIter);

%% Robust estimation of inverse Sigma_E
[U,SV,V] = svd(Sigma_E);
sv = diag(SV);
sv = sv( sv > sv(1)*tol_Sigma_E);
Nu = length(sv);
invSigma_e2 = U(:,1:Nu)*diag(1./sqrt(sv))*U(:,1:Nu)';
invSigma_e = U(:,1:Nu)*diag(1./sv)*U(:,1:Nu)';
logdet_invSigma_e = 2*log(det(U)) - sum(log(sv));
clear U SV V sv


% invSigma_e = inv(Sigma_E);

%Pre-compute
AX_all = kron(Phi',L*Psi);
C = kron(speye(Nt),invSigma_e);
% YPhit = Y*Phi';
% Karr = true(1,Nk);          %Array with indices of the active basis

%Rphi = Phi*Phi';


invgammaSqR = invgamma.^(1/2);

% YX = Y*X;                       %Pre-computation  / Barreto
% XtX = X'*X;                     %Pre-computation  / Barreto
Ic = speye(Nc);                 %Pre-computation

ite = 0;
dmu = inf;

while (ite< maxIter && dmu>MIN_DMU)
    ite = ite+1;
    
    %---------------------------------
    %% Check for active set - pruning
    %---------------------------------
    [iac_gamma, invgammaSqR, invgammaSqR_active, Nac_gam] = active_set(invgammaSqR,1e-3);
    [iac_alpha, invalpha, invalpha_active, Nac_alpha] = active_set(invalpha,1e-5);
    
    
    %----------------------------------------
    %% Update of spatio-temporal weights - W
    %----------------------------------------
    OmegaWiAll_diag(:) = 0; %19072010
    
    AX = AX_all(:,iac_alpha);         %Active part of X_all associated with active set of W.
        
    % invOmega = invOmega_L * invOmega_L' )
    % Omega = (invOmega_L)^(-1)' * invOmega_L^(-1);
    
%     % Inverse with cholesky factorization
%     disp('chol')
%     tic
%     Dtilde_idiag = (1:Nac_alpha+1:Nac_alpha^2)';
%     temp = array_add_diagElements( invalpha_active(:).^(-1) , AX' * C * AX  , Dtilde_idiag);
%     invOmega_L = chol( temp ,'lower' );
%     Omega_L = pinv(invOmega_L);
%     diagOmega = sum(Omega_L.^2,1)';                     %OK
%     w = Omega_L' * (Omega_L * (AX'*C) * Res_s(:) );     %OK
%     toc
% %     % Debug of usage of chol
% %     w2 = (Omega_L' * Omega_L) * ( (AX'*C) * Res_s(:) );
% %     compare_variables(w,w2,'Check calc w:')
% %     diagOmega2 = diag(Omega_L' * Omega_L);      %
% %     compare_variables(diagOmega,diagOmega2,'Check diagOmega1 and diagOmega2:')
    
    %Robust estimate Omega by SVD and check of tolerance on singular values
    % SVD is approximately 2 times as fast as the chol implementation
    Dtilde = spdiags(invalpha_active(:).^(-1),0,Nac_alpha,Nac_alpha);
    [U,Sv,V] = svd(full(AX' * C * AX + Dtilde),'econ');
    sv = diag(Sv);
    sv = sv( sv > sv(1)*tol_svd);
    Nu = length(sv);
    sv = sv(1:Nu); U = U(:,1:Nu);
    Omega = U * spdiags(1./sv,0,Nu,Nu) * U';
    w = Omega * ( AX'*C * Res_s(:) );
    diagOmega = diag(Omega);
%     diagOmega3 = diag(Omega);
%     compare_variables(diagOmega3,diagOmega2,'Check diagOmega3 and diagOmega2:')

%     logdetOmega = 2*log(det(U)) - sum( log(sv) );
    
    Yw = AX*w;                      
    Res_w(:) = Y(:) - Yw;           %Data residual based on W-term

    
    
    %---------------------------------
    %% Update of sources S
    %---------------------------------
    A = L(:,iac_gamma);
    AinvSigmaE = A'*beta*invSigma_e;
    
    
%     if Nac_gam>Nc
        AinvGamma = bsxfun(@times,A,invgammaSqR_active'.^2);

        %SVD formulation
        [U Sv V] = svd( AinvGamma*A' + Sigma_E / beta ) ;
        sv = diag(Sv);
        sv = sv( sv > sv(1)*tol_svd);
        Nu = length(sv);
        sv = sv(1:Nu); U = V(:,1:Nu);
        SigmaY = U * diag(1./sv) * U';

        diagSigma_S = invgammaSqR_active.^2 - diag(AinvGamma' * SigmaY * AinvGamma);
        
%         logdet_invSigmaY = 2*log(U) + sum(log(sv));
%         logdetNoiseCov = -( Nc*log(beta) + logdet_invSigma_e );
%         logdetZ = logdet_invSigmaY - logdetNoiseCov;
%         logdet_invSigmaS = Nt*(logdetZ -sum(2*log(invgammaSqR_active)) );

    %     % Cholesky formulation
    %     invSigmaY_low = chol( AinvGamma*A' + Sigma_E / beta ,'lower' ) ;
    %     SigmaY_low = pinv(invSigmaY_low);
    %     SigmaY = SigmaY_low' * SigmaY_low;

        S(iac_gamma,:) = AinvGamma' * SigmaY * Res_w;
        
%     else     %number of active dipoles smaller than number of channels - don't use inversion lemma in this situation.
%         %         Sigma_S = pinv(AinvSigmaE*A + Gamma);
%         
%         disp('Fever source terms V than Nc')
%         
%         keyboard
%         
%         Gamma = spdiags(invgammaSqR_active.^(-2),0,Nac_gam,Nac_gam);
%         [U Sv V] = svd( AinvSigmaE*A + Gamma ) ;
%         sv = diag(Sv);
%         sv = sv( sv > sv(1)*tol_svd);
%         Nu = length(sv);
%         sv = sv(1:Nu); U = V(:,1:Nu);
%         SigmaS = U * diag(1./sv) * U';
% 
%         S(iac_gamma,:) = SigmaS * AinvSigmaE * Res_w;
%         
%         diagSigma_S = diag(SigmaS);
%     end    
    
    Res_s = Y - A*S(iac_gamma,:);
        
    %----------------------------------------
    %% Update S hyperparameters -  gamma
    %----------------------------------------
%     r = mean( reshape( S(:) -  X*W(iac_alpha) ,[Nd, Nt]).^2 , 2 );
    r = mean( S(iac_gamma,:).^2 , 2 );
    %     H_ii = diag(AinvSigmaE * A * Sigma_S );
    H_ii = 1 - invgammaSqR_active.^(-2) .* diagSigma_S;
    invgammaSqR(iac_gamma) = sqrt(r./H_ii);
    
    if sum(H_ii<0)
        disp('Elements of H_ii are negative')
        keyboard
    end
    
    %----------------------------------------
    %% Update W hyperparameters - alpha
    %---------------------------------------
    G_jj = 1 - invalpha_active(:).^(-1) .* diagOmega;      % NB tjek at dette er korrekt!
    %     G_jj = 1 - diag( bsxfun(@times, invalpha_active(:).^(-1) , Omega_L') * Omega_L );
    %     G_jj = 1 - invalpha_active(:).^(-1).*diag(Omega_L'*Omega_L);
    
    if sum(G_jj<0)
        disp('Elements of G_jj are negative')
        figure, plot(G_jj), title('G_{jj}')
        figure, plot(1./invalpha_active), title('alpha_{active}')
        keyboard
    end
    
    invalpha(iac_alpha) = w.^2 ./ G_jj;
    
        if mod(ite,20)
            waitbar_update(handle_waitbar,timer_t0,maxIter,ite);
        end
        
        
%     like = -logdet_invSigmaS + logdetOmega;
    
end

W(:) = 0;
W(iac_alpha) = w;
invgamma = invgammaSqR.^2;

close(handle_waitbar)

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Internal functions


function M2 = array_add_diagElements(v1,M1,M1_idiag)

M2 = M1;
M2(M1_idiag) = M2(M1_idiag) + v1;

end

function [iac, g, g_active, Nac] = active_set(g,thresh)

% iac = false(size(g));
iac = g > max(g(:))*thresh;      %Remove very large alpha values
g(~iac) = 0;
Nac = sum(iac(:));       %Number of active elements in invgamma
g_active = g(iac);

% if Nac==0
%     disp('Nac is 0 - no active hyperparameters')
%     keyboard
% end
end


function B = sparse_diag(X)

N = size(X(:),1);
B = sparse(1:N,1:N,X(:),N,N);

end

function B = sparse_kron_diag(N,X)
% B = kron( I_N , X )

M = size(X(:),1);
NM = M*N;
B = sparse(1:NM,1:NM,repmat(X(:),[N 1]),NM,NM);

end

function compare_variables(a,b,disp_text)

disp(disp_text)
[max_val imax] = max(abs(a-b));
max_val/b(imax)

end



