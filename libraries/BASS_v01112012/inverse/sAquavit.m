function [S W invgamma invalpha LB] = sAquavit(Y,L,Psi,Phi,Sigma_E,opts)
%==========================================================================
% Filename: sAquavit.m (function).
%
% Description: Sparse Bayesian Learning with spatio-temporal structure.
%              The model is an extension of the Aquavit algorithm [1] with
%              spatial basis function, i.e. Aquavit space (sAquavit) [2].
%              The sAquavit function can both handle spatial and temporal
%              basis functions. It is an implementation of the general
%              model
%
%                    Y = L * (Psi*W*Phi + S) + E   Sensor noise E
%                    p(Y|W,S) = prod N(yn| L * (Psi*W*Phi_n + S_n) , Sigma_E
%                    p(W) = prod N(wk| 0, invalpha )
%                    p(S) = prod N(sn | 0, invgamma)
%
%              This function factorizes the variational posterior
%              distribution q(W) = prod q(w_k), i.e. in columns of W
%              (across temporal functions).
%              This function can both handle orthonormal and
%              non-orthonormal temporal basis functions.
%
% Input:
%            Y: Observations                    [Nc x Nt]
%            A: Mixing matrix (forward model)   [Nc x Nd]
%          Psi: Spatial design matrix           [Nd x Nm]
%          Phi: Temporal design matrix          [Nk x Nt]
%      Sigma_E: Noise covariance matrix         [Nc x Nc]
%
%         opts:   Optional parameters
%             .maxIter: Maximum number of interations
%             .init.flag_MVAB: Initialize W,S,invgamma,invalpha based on
%                             minimum-variance beamforming with half of the
%                             empirical data covariance splitted to each
%                             term S or W
%             .init.invgamma: 1 (default)
%             .init.invalpha: 1 (default)
%             .init.beta: 1 (default)      -  Noise learning currently not possible
%             .tol_Sigma_E: Relative tolerance on the eigenvalues in Sigma_E
%                           1e-6 (default)
%             .flag_upRule: 1 - EM  (Not implemented yet)
%                           2 - Mackay
%                           3 - Wipf (default)
%             .flag_calcVBbound: calc lower bound on/off (default: true)
%             .flag_OrthoPhis: True if Phi is orthonormal (default: false).
%
% Output:      S: Current source 'noise'
%              W: Spatio-temporal maps (regression coefficients)
%          invgamma: Inverse precision parameters on sources S.
%          invalpha: Inverse precision parameters on spatio-temporal maps.
%            ite: Number of iterations performed
%
% Special remarks:
%   Note the factorization of q(W) = prod q(w_k) is a natural consequence
%   when orthonormal temporal basis functions are used, i.e. Phi*Phi' = I.
%
%           Nc: Number of channels
%           Nt: Number of samples
%           Nd: Number of dipoles/sources/voxels
%           Nm: Number of spatial basis functions
%           Nk: Number of temporal basis functions
%
% References
%   [1]: Stahlhut, C.; Attias, H. T.; Wipf, D.; Hansen, L. K. & Nagarajan,
%        S. S.; "Sparse Spatio-Temporal Inference of Electromagnetic Brain
%        Sources", Machine Learning in Medical Imaging (MLMI) 2010,
%        Springer-Verlag Berlin Heidelberg, 2010, 157-164.
%
%   [2]: Stahlhut, C.; Attias, H. T.; Wipf, D.; Hansen, L. K. & Nagarajan,
%        S. S.; "Probabilistic M/EEG source imaging from sparse spatio-
%        temporal evetn structure", Machine Learning and Interpretation in
%        Neuroimaging (MLINI'12), 2012, 1-6.
%
%  Author: Carsten Stahlhut, DTU Informatics
%==========================================================================

if nargin<6
    opts = [];
end

if isfield(opts,'tol_Sigma_E'), tol_Sigma_E = opts.tol_Sigma_E; else tol_Sigma_E = 1e-6; end

if isfield(opts,'flag_OrthoPhis')
    flag_OrthoPhis = opts.flag_OrthoPhis;
else
    flag_OrthoPhis = false;
end



[U,SV,V] = svd(Sigma_E);
sv = diag(SV);
sv = sv( sv > sv(1)*tol_Sigma_E);
Nu = length(sv);
invSigma_e2 = U(:,1:Nu)*diag(1./sqrt(sv))*U(:,1:Nu)';

A = invSigma_e2*L;
Yp = invSigma_e2*Y;
    

if flag_OrthoPhis
    [S W invgamma invalpha LB] = sAquavit_OrthoNormal(Yp,A,Psi,Phi,1,opts);
else
    [S W invgamma invalpha LB] = sAquavit_NonOrthoNormal(Yp,A,Psi,Phi,1,opts);
end


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  NON-ORTHONORMAL TEMPORAL BASIS FUNCTION FORMULATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S W invgamma invalpha LB] = sAquavit_NonOrthoNormal(Y,L,Psi,Phi,Sigma_E,opts)


MIN_DMU         = 1e-12;

[Nc, Nt] = size(Y);
Nd = size(L,2);
Nk = size(Phi,1);
Nm = size(Psi,2);
%Pre-compute
F = L*Psi;

%% Extract optional parameter specification
if isfield(opts,'maxIter'), maxIter = opts.maxIter; else maxIter = 100; end;
if isfield(opts,'flag_upRule'), flag_updateRule = opts.flag_upRule; else flag_updateRule = 2; end
if isfield(opts,'flag_calcVBbound'), flag_calcVBbound = opts.flag_calcVBbound; else flag_calcVBbound = true; end
if isfield(opts,'tol_Sigma_E'), tol_Sigma_E = opts.tol_Sigma_E; else tol_Sigma_E = 1e-6; end
if isfield(opts,'tol_svd'), tol_svd = opts.tol_svd; else tol_svd = 1e-6; end


% Get user specified init of hyperparameters
if ~isfield(opts,'init'), opts.init = []; end
if isfield(opts.init,'flag_MVAB'), flag_MVAB = opts.init.flag_MVAB; else flag_MVAB = false; end
if isfield(opts.init,'invgamma'), invgamma2 = sqrt( opts.init.invgamma ); else invgamma2 = ones(Nd,1); end
if isfield(opts.init,'invalpha'), invalpha2 = sqrt( opts.init.invalpha ); else invalpha2 = ones(Nm,1); end
if isfield(opts.init,'beta'), beta = opts.init.beta; else beta = 1; end
if isfield(opts.init,'W'), W = opts.init.W; else W = zeros(Nm,Nk); end
if isfield(opts.init,'S'), S = opts.init.S; else S = zeros(Nd,Nt); end


% Initialize with Minimum Variance Beamformer solution with hCyy=0.5*Cyy
% (half of the empirical data covariance) associated with S term and
% the other half 0.5*Cyy associated with the W-term.    
if flag_MVAB    
    if ~isfield(opts.init,'invgamma')
        S = mvab(L,Y,0.5,tol_Sigma_E);      %Internal MVAB method
        invgamma2 = sqrt( mean( S.^2,2 ) );
    end
    
    if ~isfield(opts.init,'invalpha')
        W = mvab(F,Y*Phi',0.5,tol_Sigma_E); %Internal MVAB method
        invalpha2 = sqrt( mean( W.^2,2 ) );        
    end
    disp('Initialization to MVAB solution with half of the empirical data covariance to each term.')
end
% ----------- Extract optional parameter specification done ---------------



%% Reserved space for parameters


%Create waitbar
[handle_waitbar,timer_t0] = waitbar_create(mfilename);

logEv = zeros(1,maxIter);
Nactive = zeros(1,maxIter);
Ic = speye(Nc);                 %Pre-computation

%% Robust estimation of inverse Sigma_E
if size(Sigma_E,1)==1
    Sigma_E = Ic;
    invSigma_e2 = Sigma_E;
    invSigma_e = Sigma_E;
    logdet_invSigma_e = 0;
else
    [U,SV,V] = svd(Sigma_E);
    sv = diag(SV);
    sv = sv( sv > sv(1)*tol_Sigma_E);
    Nu = length(sv);
%     invSigma_e2 = U(:,1:Nu)*diag(1./sqrt(sv))*U(:,1:Nu)';
    invSigma_e = U(:,1:Nu)*diag(1./sv)*U(:,1:Nu)';
    invSigma_e2 = (invSigma_e2)^(1/2);
    logdet_invSigma_e = 2*log(det(U)) - sum(log(sv));
    clear U SV V sv
end


% Temporal basis functions
% Ensure numerical stability
Rphi = Phi*Phi';
scalePhi = 1 / sqrt(max(diag(Rphi)));
Phi = scalePhi * Phi;
Rphi = Phi*Phi';


eps1 = 1e-8;
invalpha2 = max(invalpha2,eps1);
invgamma2 = max(invgamma2,eps1);

Yw = F*W*Phi;       % Data described by W-term
Ys = L*S;           % Data described by S-term

ite = 0;
dmu = inf;

Nac_alpha = Nm;
Nac_gamma = Nd;

LB = -inf*ones(1,maxIter);
hfig = figure;

Csq_W = cell(1,Nk);
diagOmega_W = cell(1,Nk);
logdet_wAll = zeros(1,Nk);
while (ite< maxIter && dmu>MIN_DMU)
    ite = ite+1;
    
    %----------------------------
    %% Update of W
    %----------------------------
    invDsq = spdiags(invalpha2,0,Nac_alpha, Nac_alpha);
    iw = randperm(Nk);      %Random order update of w_k
    
    X = F*invDsq;
    for j=1:Nk
        k=iw(j);
        
        k_ = (1:Nk)~=k;
        W_ = W(:, k_ );
        Yw_ = F*W_*Rphi(k_,k);
        
        Yres = (Y-Ys)*Phi(k,:)' - Yw_;
        
        [W(:,k) diagOmega_wk Csq_Wk logdet_wAll(k)] = gausPost_invlemma( X ,invDsq,Sigma_E/Rphi(k,k), Yres/Rphi(k,k) , tol_svd, flag_calcVBbound);
        Csq_W{k} = Csq_Wk;
        diagOmega_W{k} = diagOmega_wk;
    end
    logdet_w = sum(logdet_wAll);
    Yw = F*W*Phi;    
    
    
    %---------------------------------
    %% Update of sources S
    %---------------------------------
    invGamma_sq = spdiags(invgamma2,0,Nac_gamma, Nac_gamma);
    [S diagSigma_S Csq_s logdet_s] = gausPost_invlemma(L*invGamma_sq,invGamma_sq,Sigma_E, Y-Yw , tol_svd, flag_calcVBbound);
    logdet_s = Nt * logdet_s;
    Ys = L*S;
    
    
    % Update hyperparameter for W: invalpha2
    % Update hyperparameter for S: invgamma2
    %----------------------------------------
    switch flag_updateRule
        case 0 % EM updates
            error('EM-updates not implemented yet')
            disp('EM update rule')

        case 1  % Mackay update rule
            invalpha2 = mackay_update(W,invalpha2,diagOmega, 2, eps1);
            invgamma2 = mackay_update(S,invgamma2,diagSigma_S, 2, eps1);
            
        case 2  % Wipf update rule
            invalpha2 = wipf_update(W,Csq_W,F,2,eps1);  % W * Nk as wipf_update calculates the mean over columns and it is the sum that should be used.
            invgamma2 = wipf_update(S,{Csq_s},L,2,eps1);

        otherwise
            error('Specified update rule not acceptabel - type ´help sAquavit´ to see choices')
    end
    
    
    %---------------------------------
    %% Update lower bound
    %---------------------------------
    if flag_calcVBbound
        
        Ey = (Y-Ys-Yw);
        Ey2 = trace(invSigma_e * (Ey*Ey') );
                
        % Es2 = sum_n || s_n ||^2_Gamma
        SgammaSq = bsxfun(@times,S,1./invgamma2);
        Es2 = sum(sum(SgammaSq.^2));
        
        % Ew2 = sum_k || w_k ||^2_alpha
        WalphaSq = bsxfun(@times,W,1./invalpha2);
        Ew2 = sum(sum(WalphaSq.^2));

        % Complete lower bound
        % logdet_s = log(| Gamma |) - log( | invSigma_S | )
        LB(ite) = 0.5*(Nt*logdet_invSigma_e -Ey2 +...
                  logdet_s -Es2 + logdet_w  -Ew2);             
    end
    
    
    %% Plot
    figure(hfig)
    subplot(2,2,1);plot((1:ite),LB(1:ite));
        title(['Lower bound: ' int2str(ite) ' / ' int2str(maxIter)]);
        xlabel('iteration');
        set(gca(),'XLim',[0 ite]);
    subplot(2,2,2);plot((1:Nd),invgamma2.^2);
        title(['S: Voxel power: ' num2str(length(invgamma2)) ' / ' num2str(Nd)]);
        xlabel('voxel index');
        set(gca(),'XLim',[1 Nd]);
    subplot(2,2,3);plot((1:Nm),invalpha2.^2);
        title(['W: Region power: ' num2str(length(invalpha2)) ' / ' num2str(Nm)]);
        xlabel('Region index');
        set(gca(),'XLim',[1 Nm]);
    drawnow

                   
    
    if mod(ite,20)
        waitbar_update(handle_waitbar,timer_t0,maxIter,ite);
    end
    
end

invalpha = invalpha2.^2;
invgamma = invgamma2.^2;

close(handle_waitbar)

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ORTHONORMAL TEMPORAL BASIS FUNCTION FORMULATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S W invgamma invalpha LB] = sAquavit_OrthoNormal(Y,L,Psi,Phi,Sigma_E,opts)


MIN_DMU         = 1e-12;

[Nc, Nt] = size(Y);
Nd = size(L,2);
Nk = size(Phi,1);
Nm = size(Psi,2);

%Pre-compute
F = L*Psi;

%% Extract optional parameter specification
if isfield(opts,'maxIter'), maxIter = opts.maxIter; else maxIter = 100; end;
if isfield(opts,'flag_upRule'), flag_updateRule = opts.flag_upRule; else flag_updateRule = 2; end
if isfield(opts,'flag_calcVBbound'), flag_calcVBbound = opts.flag_calcVBbound; else flag_calcVBbound = true; end
if isfield(opts,'tol_Sigma_E'), tol_Sigma_E = opts.tol_Sigma_E; else tol_Sigma_E = 1e-6; end
if isfield(opts,'tol_svd'), tol_svd = opts.tol_svd; else tol_svd = 1e-6; end


% Get user specified init of hyperparameters
if ~isfield(opts,'init'), opts.init = []; end
if isfield(opts.init,'flag_MVAB'), flag_MVAB = opts.init.flag_MVAB; else flag_MVAB = false; end
if isfield(opts.init,'invgamma'), invgamma2 = sqrt( opts.init.invgamma ); else invgamma2 = ones(Nd,1); end
if isfield(opts.init,'invalpha'), invalpha2 = sqrt( opts.init.invalpha ); else invalpha2 = ones(Nm,1); end
if isfield(opts.init,'beta'), beta = opts.init.beta; else beta = 1; end
if isfield(opts.init,'W'), W = opts.init.W; else W = zeros(Nm,Nk); end
if isfield(opts.init,'S'), S = opts.init.S; else S = zeros(Nd,Nt); end
%if isfield(opts,'group_alpha'), group_alpha = opts.group_alpha; else group_alpha = 1; end

% Initialize with Minimum Variance Beamformer solution with hCyy=0.5*Cyy
% (half of the empirical data covariance) associated with S term and
% the other half 0.5*Cyy associated with the W-term.    
if flag_MVAB    
    if ~isfield(opts.init,'invgamma')
        S = mvab(L,Y,0.5,tol_Sigma_E);      %Internal MVAB method
        invgamma2 = sqrt( mean( S.^2,2 ) );
    end
    
    if ~isfield(opts.init,'invalpha')
        W = mvab(F,Y*Phi',0.5,tol_Sigma_E); %Internal MVAB method
        invalpha2 = sqrt( mean( W.^2,2 ) );        
    end
    disp('Initialization to MVAB solution with half of the empirical data covariance to each term.')
end

% ----------- Extract optional parameter specification done ---------------



%% Reserved space for parameters
%Create waitbar
[handle_waitbar,timer_t0] = waitbar_create(mfilename);

logEv = zeros(1,maxIter);
Nactive = zeros(1,maxIter);
Ic = speye(Nc);                 %Pre-computation

%% Robust estimation of inverse Sigma_E
if size(Sigma_E,1)==1
    Sigma_E = Ic;
    invSigma_e2 = Sigma_E;
    invSigma_e = Sigma_E;
    logdet_invSigma_e = 0;
else
    [U,SV,V] = svd(Sigma_E);
    sv = diag(SV);
    sv = sv( sv > sv(1)*tol_Sigma_E);
    Nu = length(sv);
%     invSigma_e2 = U(:,1:Nu)*diag(1./sqrt(sv))*U(:,1:Nu)';
    invSigma_e = U(:,1:Nu)*diag(1./sv)*U(:,1:Nu)';
    invSigma_e2 = (invSigma_e2)^(1/2);
    logdet_invSigma_e = 2*log(det(U)) - sum(log(sv));
    clear U SV V sv
end


eps1 = 1e-8;
invalpha2 = max(invalpha2,eps1);
invgamma2 = max(invgamma2,eps1);

Yw = F*W*Phi;       % Data described by W-term
Ys = L*S;           % Data described by S-term


ite = 0;
dmu = inf;

Nac_alpha = Nm;
Nac_gamma = Nd;

LB = -inf*ones(1,maxIter);
hfig = figure;

while (ite< maxIter && dmu>MIN_DMU)
    ite = ite+1;
    
    %----------------------------
    %% Update of W
    %----------------------------
    invDsq = spdiags(invalpha2,0,Nac_alpha, Nac_alpha);
    [W diagOmega Csq_w logdet_w] = gausPost_invlemma(F*invDsq,invDsq,Sigma_E, (Y-Ys)*Phi' , tol_svd, flag_calcVBbound);    
    logdet_w = logdet_w * Nk;
    Yw = F*W*Phi;    
    
    %---------------------------------
    %% Update of sources S
    %---------------------------------
    invGamma_sq = spdiags(invgamma2,0,Nac_gamma, Nac_gamma);
    [S diagSigma_S Csq_s logdet_s] = gausPost_invlemma(L*invGamma_sq,invGamma_sq,Sigma_E, Y-Yw , tol_svd, flag_calcVBbound);
    logdet_s = logdet_s*Nt;
    Ys = L*S;
    
    
    % Update hyperparameter for W: invalpha2
    % Update hyperparameter for S: invgamma2
    %----------------------------------------
    switch flag_updateRule
        case 0      % EM updates
            error('EM-updates not implemented yet')
            disp('EM update rule')

        case 1      % Mackay update rule
            invalpha2 = mackay_update(W,invalpha2,diagOmega, 2, eps1);
            invgamma2 = mackay_update(S,invgamma2,diagSigma_S, 2, eps1);
            
        case 2      % Wipf update
            invalpha2 = wipf_update(W,{Csq_w},F,2,eps1);
            invgamma2 = wipf_update(S,{Csq_s},L,2,eps1);

        otherwise
            error('Specified update rule not acceptabel - type ´help sAquavit´ to see choices')
    end
    
    
    %---------------------------------
    %% Update lower bound
    %---------------------------------
    if flag_calcVBbound
        
        Ey = (Y-Ys-Yw);
        Ey2 = trace(invSigma_e * (Ey*Ey') );
                
        % Es2 = sum_n || s_n ||^2_Gamma
        SgammaSq = bsxfun(@times,S,1./invgamma2);
        Es2 = sum(sum(SgammaSq.^2));
        
        % Ew2 = sum_k || w_k ||^2_alpha
        WalphaSq = bsxfun(@times,W,1./invalpha2);
        Ew2 = sum(sum(WalphaSq.^2));

        % Complete lower bound
        % logdet_s = log(| Gamma |) - log( | invSigma_S | )
        LB(ite) = 0.5*(Nt*logdet_invSigma_e -Ey2 +...
                  logdet_s -Es2 + logdet_w  -Ew2);
    end
    
    
    %% Plot
    figure(hfig)
    subplot(2,2,1);plot((1:ite),LB(1:ite));
        title(['Lower bound: ' int2str(ite) ' / ' int2str(maxIter)]);
        xlabel('iteration');
        set(gca(),'XLim',[0 ite]);
    subplot(2,2,2);plot((1:Nd),invgamma2.^2);
        title(['S: Voxel power: ' num2str(length(invgamma2)) ' / ' num2str(Nd)]);
        xlabel('voxel index');
        set(gca(),'XLim',[1 Nd]);
    subplot(2,2,3);plot((1:Nm),invalpha2.^2);
        title(['W: Region power: ' num2str(length(invalpha2)) ' / ' num2str(Nm)]);
        xlabel('Region index');
        set(gca(),'XLim',[1 Nm]);
    drawnow
  
                
    if mod(ite,20)
        waitbar_update(handle_waitbar,timer_t0,maxIter,ite);
    end
    
end

invalpha = invalpha2.^2;
invgamma = invgamma2.^2;

close(handle_waitbar)

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Internal functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mu diagCov Csq logdetZ] = gausPost_invlemma(X,Psq,R,Y,tol_svd,flag_calcVBbound)
% X = L * Psq , with L being the forward model and Psq the
% inv-precision^(1/2)

Csq = inverse_svd( full(X * X' + R) , tol_svd);
Mu = Psq * X' * (Csq'*Csq) * Y;

% Use matrix inversion lemma to get the diagonal of cov-matrix
diagXXii = sum( ( Csq * X ).^2,1)';   % diag( X' * Csq' * Csq * X )
diagCov = diag(Psq).^2 .* (1 - diagXXii );

if ~isempty(find(diagCov<0))
    warning('Diagonal elements negative in Covariance matrix')
    keyboard
end


if flag_calcVBbound
    logdetZ = log( det(R)) - log( det(X * X' + R));
else
    logdetZ = NaN;
end

end


function Mu = gausPost_invlemma_wCsq(X,Psq,Csq,Y,tol_svd)
% X = L * Psq , with L being the forward model and Psq the
% inv-precision^(1/2)

Mu = Psq * X' * (Csq'*Csq) * Y;


end





function [Csq U sv] = inverse_svd(invC,tol_svd)

[U Sv V] = svd( invC ) ;
sv = diag(Sv);
sv = sv( sv > sv(1)*tol_svd);
Nu = length(sv);
sv = sv(1:Nu); U = U(:,1:Nu);
% C = U * diag(1./sv) * U';
Csq = diag(1./sqrt(sv) + 0.5*eps/sqrt(sv(1)) ) * U';

end


% function [invhyp2_new] = wipf_update(X,Csq,A,dim_mean,eps1)
% 
% 
% X2mean = mean(X.^2,dim_mean);
% CsqA = Csq*A;
% negGrad = sum(CsqA.^2,1)';
% invhyp_new = ( (1./negGrad) .* X2mean ).^(1/2);
% invhyp2_new = max( invhyp_new.^(1/2) ,eps1);
% 
% if max(abs(imag(invhyp2_new(:))))
%     warning('complex number')
%     keyboard
% end
% 
% end


function [invhyp2_new] = wipf_update(X,Csq,A,dim_mean,eps1)
% This formulation can handle cell-arrays of Csq

X2mean = mean(X.^2,dim_mean);
CsqA = cellfun(@(Csq) Csq * A , Csq,'UniformOutput', false);
aiCai = cellfun(@(CsqA) sum(CsqA.^2,1)', CsqA, 'UniformOutput',false);   
negGrad = mean( cell2mat(aiCai), 2);

invhyp_new = ( (1./negGrad) .* X2mean ).^(1/2);
invhyp2_new = max( invhyp_new.^(1/2) ,eps1);


% % debug for numerical problems
% if max(abs(imag(invhyp2_new(:))))
%     warning('complex number')
%     keyboard
% end

end



function [invhypSqR_new] = mackay_update(X,invhypSqR,diagPostCovX,dim_mean,eps1)
% Mackay updates of they hyperparameters associated X, diagPostCovX is the
% diagonal of the posterior covariance matrix associated X.
% dim_mean specifies the dimension to be averaged across.
% eps1 serves a variable agains numerical instability.

r = mean(X.^2 , dim_mean);
H_ii = 1 - invhypSqR.^(-2) .* diagPostCovX;

invhypSqR_new = sqrt(r./H_ii);
invhypSqR_new = max(invhypSqR_new,eps1);


% % debug for numerical problems
% if sum(isinf(invhypSqR_new))
%     disp('Elements inf')
% end
% 
% if sum(H_ii<=0)
%     disp('Elements of H_ii are negative or zero')
% end


end


function X = mvab(L,Y,CovRatio,tol_svd)
% L: lead field
% Y: data
% CovRatio: The ratio of the data covariance of what X should assign to.
%           For normal mvab CovRatio=1 - X is based of the total
%           data-covariance.
% X: output source solution
%------------------------------
hCyy = CovRatio*(Y*Y')/size(Y,2);
[U Sv V] = svd(hCyy);
sv = diag(Sv);
sv = sv( sv > sv(1)*tol_svd);
Nu = length(sv);
invhCyy = U(:,1:Nu)*diag(1./sv)*U(:,1:Nu)';
invhCyy2 = (invhCyy)^(1/2);

Z = invhCyy2*L;
invZZii = 1./sum(Z.^2,1)';
X = bsxfun(@times, Z'*invhCyy2*Y , invZZii);
end