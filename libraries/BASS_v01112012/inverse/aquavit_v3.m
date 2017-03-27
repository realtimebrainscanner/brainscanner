function [J W gamma alpha ite] = aquavit_v3(Y,L,Phi,Sigma_E,opts)
%==========================================================================
% Filename: aquavit_v3.m - version 3 (function).
% 
% Description: VB statio-temporal reconstruction
%              Implementation of [1]. The model is given by
%                       Y = A*S + E   Sensor noise E
%                       S = W'*Phi + V;
%
%
% Input:    
%            Y: Observations             [Nc x Nt]
%            L: Lead-field matrix        [Nc x Nd]
%          Phi: Temporal design matrix   [Nk x Nt]
%      Sigma_E: Noise covariance matrix  [Nc x Nc]
%         opts:
%             .maxIter: Maximum number of interations
%             .init.gamma: 1 (default)
%             .init.alpha: 1 (default)
%             .init.beta: 1 (default)
%             .group_alpha: 1 (default)
%             .index_groups: Nk x Nd matrix specifying which cluster
%                            alpha_ki belongs to. This allows spatio-
%                            temporal sub clusters. Default: Nk*Nd groups.
%             .tol_Sigma_E: Relative tolerance on the eigenvalues in Sigma_E
%                           1e-6 (default)
%             .flag_upRule: 1 - EM
%                           2 - Mackay (default)
%                           3 - Wipf
%             .thresh.alpha_rel:
%             .thresh.alpha_abs:
%             .thresh.gamma_rel:
%             .thresh.gamma_abs:
%
%
% Output:      S: Current sources
%              W: Spatio-temporal maps (regression coefficients)
%          gamma: Precision parameters on sources and basis.
%            ite: Number of iterations performed
%
% Special remarks:
%   [1]: Stahlhut, C.; Attias, H. T.; Wipf, D.; Hansen, L. K. & Nagarajan,
%        S. S.; "Sparse Spatio-Temporal Inference of Electromagnetic Brain
%        Sources", Machine Learning in Medical Imaging (MLMI) 2010,
%        Springer-Verlag Berlin Heidelberg, 2010, 157-164
%
%           Nc: Number of channels
%           Nt: Number of samples
%           Nd: Number of sources
%           Nk: Number of basis functions
%
%           The function was previously called vbst2 --> aquavit1.m
%   
%   Updates: - SVD-applied on Sigma_E when calculating inv(Sigma_E) with a
%              specified tol on the eigenvalues --> more stabile.
%            - Faster inversion of Sigma_J - actually not calculating this
%              anymore. Calcs a modified A.
%            - 12-02-2011: Wipf updates included. More flexible           
%
% Author:
%   Carsten Stahlhut, DTU Informatics, Technical University of Denmark.
%==========================================================================

MIN_DMU         = 1e-12;

[Nc, Nt] = size(Y);
Nd = size(L,2);
Nk = size(Phi,1);


%% Extract optional parameter specification
try maxIter = opts.maxIter; catch maxIter = 100; end;
% Get user specified init of hyperparameters
try gamma = opts.init.gamma; catch gamma = ones(Nd,1); end
try alpha = opts.init.alpha; catch alpha = ones(Nk,Nd); end
try beta = opts.init.beta; catch beta = 1; end
try group_alpha = opts.group_alpha; catch group_alpha = 0; end
try tol_Sigma_E = opts.tol_Sigma_E; catch tol_Sigma_E = 1e-6; end

% try flag_calcLogEv = opts.flag_calcLogEv; catch flag_calcLogEv = 0; end
try flag_upRule = opts.flag_upRule; catch flag_upRule = 3; end

try alpha_rel = opts.alpha_rel; catch alpha_rel = 1e6; end   %Relative threshold
try alpha_abs = opts.alpha_abs; catch alpha_abs = 1e6; end   %Absolute threshold
try gamma_rel = opts.gamma_rel; catch gamma_rel = 1e6; end   %Relative threshold
try gamma_abs = opts.gamma_abs; catch gamma_abs = 1e6; end   %Absolute threshold



W = zeros(Nk,Nd);
if size(alpha,2)<Nd
    index_k = ones(Nd,1);                    %Do not use space on copy Z times D, since D_k = D.
else
    index_k = (1:Nd)';
end
if group_alpha
    if isfield(opts,'index_groups'), Igroup = opts.index_groups; else Igroup = reshape(1:Nk*Nd,[Nk Nd]); end
else
    Igroup = 1:Nk*Nd;
end

% Ni_group = size(Igroup,1);
% if Ni_group>1, Ni_group = sum(logical(Igroup),1); end

SigmaWiAll_diag = W;
OmegaWiAll_diag = W;
SigmaJ_diag = zeros(Nd,1);
J = zeros(Nd,Nt);
Jhat = J;
E_Ji = zeros(size(gamma));

%Create waitbar
[handle_waitbar,timer_t0] = waitbar_create(mfilename);

logEv = zeros(1,maxIter);
Nactive = zeros(1,maxIter);

%% Robust estimation of inverse Sigma_E
[U,S,V] = svd(Sigma_E);
s = diag(S);
s = s( s > s(1)*tol_Sigma_E);
Nu = length(s);
invSigma_e2 = U(:,1:Nu)*diag(1./sqrt(s))*U(:,1:Nu)';
invSigma_e = U(:,1:Nu)*diag(1./s)*U(:,1:Nu)';
invSigmae2L = invSigma_e2*L;
clear U S V s
% invSigma_e = inv(Sigma_E);
%Pre-compute
% YPhit = Y*Phi';
% Karr = true(1,Nk);          %Array with indices of the active basis
Rphi = Phi*Phi';

iac = true(size(alpha));
iac_gamma = true(size(gamma));

% YX = Y*X;                       %Pre-computation  / Barreto
% XtX = X'*X;                     %Pre-computation  / Barreto
Ic = speye(Nc);                 %Pre-computation


ite = 0;
dmu = inf;

while (ite< maxIter && dmu>MIN_DMU)
    ite = ite+1;
%     K_active = find(Karr);

    %---------------------------------
    % Check for active set - pruning
    %---------------------------------    
    %Active set
%     iac( alpha > min(alpha(:))*1e6  ) = false;      %Remove very large alpha values
    iac( alpha > min(alpha(:))*alpha_rel  ) = false;      %Remove very large alpha values
    alpha(~iac) = inf;

%     iac_gamma( gamma > min(gamma(:))*1e6  ) = false;      %Remove very large gamma values
    iac_gamma( gamma > min(gamma(:))*gamma_rel  ) = false;      %Remove very large gamma values
    gamma(~iac_gamma) = inf;
    Nac_gam = sum(iac_gamma(:));
    
    
%     Ngold = Ng;
    Ng = sum(iac(:));
    Nactive(ite) = Ng;
%     gamma_t = gamma(iac);
    gamma_t = gamma(iac_gamma);
    
%     A = L(:,iac);
%     invGamma = spdiags(1./gamma_t(:),0,Ng,Ng);

    active_i = find(iac);
%     alpha(:,active_i) = inf;        %Force to inf, since J is zero here
%     iac_alpha( alpha > min(alpha(:))*1e6  ) = false;      %Remove very large gamma values
%     alpha(~iac_alpha) = inf;

    
    %---------------------------------
    % Update of J
    %--------------------------------- 
    Jold = J;
    
    Gamma = spdiags(gamma_t,0,Nac_gam,Nac_gam);
    invGammaSq = spdiags( gamma_t.^(-1/2),0,Nac_gam,Nac_gam);
    Afast = invSigmae2L(:,iac_gamma)*invGammaSq;
    Id = speye(Nac_gam);

    A = L(:,iac_gamma);
    AinvSigmaE = A'*beta*invSigma_e;
    
    if Nac_gam>Nc
        [U,Lam,V] = svd(Afast,'econ');
        V = Lam*V';
        lam = diag(Lam);
        Nl = length(lam);
        
        P = V'*diag(1./(1 + lam.*lam))*V;
        %             P = V'*pinv(Ic + Lam*Lam')*V;
        Sigma_J = invGammaSq*(Id - P)*invGammaSq;
    else
        Sigma_J = pinv(AinvSigmaE*A + Gamma);
        P = diag(1-gamma_t.*diag(Sigma_J));
    end
    
    J(iac_gamma,:) = Sigma_J*(AinvSigmaE*Y + Gamma*W(:,iac_gamma)'*Phi);
    J(~iac_gamma,:) = W(:,~iac_gamma)'*Phi;
    
    %         keyboard
    
    %     Gamma = spdiags(gamma,0,Nd,Nd);
    %     AinvSigmaE = A'*beta*invSigma_e;
    %     Sigma_J = pinv(AinvSigmaE*A + Gamma);
    %     J = Sigma_J*(AinvSigmaE*Y + Gamma*W'*Phi);
    
    SigmaJ_diag = diag(Sigma_J);
    %     SigmaJ_diag(iac) = invgam - sum(AinvGamma.*(invSigma_m*AinvGamma),1)';

    %---------------------------------
    % Update of W
    %---------------------------------
%     SigmaWiAll_diag(:) = 0;
    OmegaWiAll_diag(:) = 0; %19072010
    W(:) = 0;
    %Active set
    
    for i=1:Nd
%         active_i
%         i = active_i(j);
        
        Di = alpha(:,index_k(i));
        iac_i = iac(:,index_k(i));
        
        Rphi_ac = Rphi(iac_i,:);
        Rphi_ac = Rphi_ac(:,iac_i);
        
%         Sigma_wi = pinv(Rphi_ac + diag( Di(iac_i) ))/gamma(i);
%         wi = Sigma_wi*gamma(i)*Phi(iac_i,:)*J(i,:)';
%         SigmaWiAll_diag(iac_i,i) = diag(Sigma_wi);
        Omega_wi = pinv(Rphi_ac + diag( Di(iac_i) ));
        wi = Omega_wi*Phi(iac_i,:)*J(i,:)';
        W(iac_i,i) = wi;
        OmegaWiAll_diag(iac_i,i) = diag(Omega_wi);      
        
        Jhat(i,:) = wi'*Phi(iac_i,:);
    end

    
    switch flag_upRule
        case 1  %EM
            %---------------------------------
            % Update of gamma
            %---------------------------------       
%             gamma = (Nt+1)./( E_wtw +  );          
%             E_Ji = Nt*SigmaJ_diag + sum(J.^2,2);
            E_Ji(iac_gamma) = Nt*SigmaJ_diag + sum(J(iac_gamma,:).^2,2); %Change 15072010
            E_Ji(~iac_gamma) = sum(J(~iac_gamma,:).^2,2);
            
%             invgamma = ( Nk./gamma + E_Ji - J*Phi'*W )/(Nt+1);
            invgamma = ( Nk./gamma + E_Ji - sum(J.*(W'*Phi),2) )/(Nt+1);
            gamma = 1./(invgamma + max(invgamma)*1e-10);
            
            %---------------------------------
            % Update of alpha
            %---------------------------------        
        %     alpha = 1./bsxfun(@times,(SigmaWiAll_diag + W.^2),gamma);
            GamTemp = repmat(gamma',[Nk 1]);
%             alpha(iac) = 1./((SigmaWiAll_diag(iac) + W(iac).^2).*GamTemp(iac) );
            alpha(iac) = 1./(OmegaWiAll_diag(iac) + (W(iac).^2).*GamTemp(iac) );
            
            %Group basis functions
            if group_alpha
                invalpha_mean = ones(Nk,1)*mean(1./alpha,1);
                alpha(iac) = 1./invalpha_mean(iac);
            end
            
        case 2  %MacKay
            %---------------------------------
            % Update of gamma
            %---------------------------------       
%             Z = AinvSigmaE*A*Sigma_J;
%             z = diag(Z);       
            try
                z = 1 - gamma_t.*SigmaJ_diag;
            catch
                 if isempty(gamma_t),   z = 1;   end
            end
            sumV2 = sum((J-Jhat).^2,2);     % source noise contribution
            penW = zeros(Nk,Nd);            % 12042011: to be checked 
            penW(iac) = alpha(iac).*(W(iac).^2);
            r = sumV2 + sum(penW,1)';
%             r = sum((J-Jhat).^2,2) + sum(alpha(iac).*(W(iac).^2),1)';
% %             r = sum((J-W'*Phi).^2,2) + sum(alpha.*(W.^2),1)';
            gamma_temp = Nt*z./r(iac_gamma);
            if sum(isnan(gamma_temp(:)))
                keyboard
            end
            gamma(iac_gamma) = gamma_temp;
          
            %---------------------------------
            % Update of alpha
            %---------------------------------                    
%             Z = 1 - alpha(iac).*OmegaWiAll_diag(iac);     %NB husk active set
%             GamTemp = repmat(gamma',[Nk 1]);
%             R = (W(iac).^2).*GamTemp(iac);
%             alpha(iac) = Z./R;
% 
%             %Group basis functions
%             if group_alpha
%                 invalpha_mean = ones(Nk,1)*mean(1./alpha,1);
%                 alpha(iac) = 1./invalpha_mean(iac);
%             end
            GamTemp = repmat(gamma',[Nk 1]);            
            Z = inf*ones(size(iac));
            Z(iac) = 1 - alpha(iac).*OmegaWiAll_diag(iac);     %Eq.(15) multi with D_i^{-1} in Wipf's note
            R = (W.^2).*GamTemp;
            invalpha_subgroups = (sum_subgroups( R,Igroup ) ./...
                                 sum_subgroups(Z,Igroup) );
            alpha = copy_subgroup_val(1./invalpha_subgroups,Igroup,Nk,Nd);
            alpha(~iac) = inf;
            
                                    
        case 3  %Wipf
            z = gamma_t.*diag(P);         %Eq. (13) in Wipf's note
%             z = diag(Gamma - Gamma*Sigma_J*Gamma);
            sumV2 = sum((J-Jhat).^2,2);     % source noise contribution
            penW = zeros(Nk,Nd);  penW(iac) = alpha(iac).*(W(iac).^2);
            r = sumV2 + sum(penW,1)';
%             r = sum((J-Jhat).^2,2) + sum(alpha(iac).*(W(iac).^2),1)';
            gamma(iac_gamma) = sqrt( Nt*z./r(iac_gamma) );    %Eq. (19) in Wipf's note
            
            %---------------------------------
            % Update of alpha
            %---------------------------------                    
            GamTemp = repmat(gamma',[Nk 1]);
%             if group_alpha
                R = (W.^2).*GamTemp;
                Z = inf*ones(size(iac));
                Z(iac) = 1 - alpha(iac).*OmegaWiAll_diag(iac);     %Eq.(15) multi with D_i^{-1} in Wipf's note
                alpha_subgroups = sqrt(...
                    sum_subgroups( alpha.*Z,Igroup ) ./...
                    sum_subgroups(R,Igroup) );
                alpha = copy_subgroup_val(alpha_subgroups,Igroup,Nk,Nd);
                alpha(~iac) = inf;
%             else
%                 Z = 1 - alpha(iac).*OmegaWiAll_diag(iac);     %Eq.(15) multi with D_i^{-1} in Wipf's note
%                 R = (W(iac).^2).*GamTemp(iac);
%                 alpha(iac) = sqrt( (alpha(iac).*Z ) ./ R);    %Eq.(22) in Wipf's note
%             end
            
            
        otherwise
            error('Unknow update rule')
    end
    
%     %Group basis functions
%     if group_alpha
%         invalpha_mean = ones(Nk,1)*mean(1./alpha,1);
%         alpha(iac) = 1./invalpha_mean(iac);
%     end
    

    %---------------------------------
    % Update of beta
    %---------------------------------        

    %Change in mu
    [dmu idmu] = max(abs( J(:)-Jold(:) ) );
    dmu = dmu/max(abs(J(idmu)));
    
    %Update Waitbar
    waitbar_update(handle_waitbar,timer_t0,maxIter,ite)

end

% gamma(~iac) = inf;
% J(~iac,:) = 0;
% mu_x(~iac) = 0;
% G = reshape(mu_x,Nd,Nk);
% clear mu_x
% J = G*Phi;

logEv = logEv(1:ite);       %Delete un-used space
Nactive = Nactive(1:ite);           %Delete un-used space

close(handle_waitbar)

end

function sumx_subgroup = sum_subgroups(X,I)
% I is a J x K matrix. The k'th column consists of the indices of the
% elements that is within the k'th subgroup. Note that the number of
% elements within the k'th group and the g'th group don't have to be equal.
% The remaining spots in I just needs to be filled with zeros.
Nrows = size(X,1);
Y = [zeros(Nrows,1) X];
J = I + Nrows;      %Changing the indices numbers so it corresponds to same
                    %elements in Y.
sumx_subgroup = sum( Y(J), 1);

end

function X = copy_subgroup_val(y,I,Nrows,Ncols)
% Function intended to copy values from subgroups to all the elements
% within each of the subgroups.
X = zeros(Nrows,Ncols+1);
J = I + Nrows;
[N K] = size(J);        % N: #elements in a group, K:#subgroups

X( J ) = ones(N,1)*y(:)';
X(:,1) = [];

% for k=1:K
%     X(J(:,k )) = y(k);
% end

end


% function sumx_subgroup = sum_subgroups_cell(X,I)
% % I keeps the indices for each subgroup. I consists of a cell for each
% % group, since the same number of elements within a subgroup is not
% % necessarily equal. I.e. I consists of K subgroups with the length of the
% % k'th cell J_k long - J_k elements in this group.
% 
% K = length(I);
% sumx_subgroup = zeros(1,K);
% for k=1:K
%     sumx_subgroup(k) = sum( X(I{k}) );
% end
% 
% end
