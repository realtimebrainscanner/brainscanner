function [S W invAlpha invBeta,opts] = BayesMNandLORETA_svd(Y,A,opts)
%==========================================================================
% Filename: BayesMNandLORETA.m
%
% Description:  Bayesian formulation of the minimum norm (MN) and low
%               resolution brain electromagnetic tomography (LORETA).
%               Hyperparameters are updated with expectation-maximization
%               steps.
%
% Input:        Y: Observations [ Nc x Ns ]
%               A: Lead field matrix / Forward model [ Nc x Nd ]
%               opts:
%                 .maxIter: maximum number of EM-steps
%                 .K: Spatial coherence matrix - prior. Default identity
%                     corresponding to MN.
%                 .Sigma_e: Fixed noise covariance matrix [ Nc x Nc ]
%                 .invAlpha: Hyperparameter for the sources [ scalar ]
%                 .invBeta: Hyperparameter for the noise variane [ scalar ]
%                 .W: Init weights used to obtain sources S [Nd x Nc]
%                 .S: Init sources [Nd x Ns]
%                 .thresh: Threshold on singular values when inverting a
%
% Output:       S: New source estimates based on updated invAlpha, invBeta,
%                  and W.
%               W: Updated weights.
%               invAlpha: Hyperparameter for the sources [ scalar ];
%               invBeta: Hyperparameter for the noise variance [ scalar ];
%
% Last revision: 2012/03/21
%
% Special remarks:
%               opts.flag_update.hyp_trace should only be set to true when
%               it is the MN-method that is used.
%
% Authors: Carsten Stahlhut
%
% Copyright (c) DTU Informatics, Technical University of Denmark
%==========================================================================

[Nc Nd] = size(A);

%% Check for optional variables
if isfield(opts,'maxIter'), maxIter = opts.maxIter; else maxIter = 1; end
if isfield(opts,'K'), K = opts.K; else K = speye(Nd); end
if isfield(opts,'invK'), invK = opts.invK; else invK = pinv(full(K)); opts.invK = invK; end
if isfield(opts,'Sigma_e'), Sigma_e = opts.Sigma_e; else Sigma_e = eye(Nc); end
if isfield(opts,'invAlpha'), invAlpha = opts.invAlpha; else invAlpha = 1; end
if isfield(opts,'invBeta'), invBeta = opts.invBeta; else invBeta = 1; end
if isfield(opts,'AKAt'), AKAt = opts.AKAt; else AKAt = A*K*A'; opts.AKAt = AKAt; end
if isfield(opts,'KAt'), KAt = opts.KAt; else KAt = K*A'; opts.KAt = KAt; end
if isfield(opts,'thresh'), thresh = opts.thresh; else thresh = 1e-6; end
if isfield(opts,'flag_update')
    flag_update = opts.flag_update;
else
    flag_update.S = true;
    flag_update.alpha = true;
    flag_update.beta = true;
    flag_update.hyp_trace = false;      % false: works for MN and LORETA, true: only for MN
end
if ~isfield(flag_update,'S'), flag_update.S = true; end
if ~isfield(flag_update,'alpha'), flag_update.alpha = true; end
if ~isfield(flag_update,'beta'), flag_update.beta = true; end

if (isfield(opts,'W') && isfield(opts,'S'))
    W = opts.W;
    S = opts.S;
else
    %% Calculate weights W and sources S
    inputMatrix = invAlpha * AKAt + invBeta*Sigma_e;
    [UInputMat singInputMat]=svd(inputMatrix);
    sing = diag(singInputMat);
    sing = sing(sing>sing(1)*thresh);
    UInputMat = UInputMat(:,1:length(sing));
    W = invAlpha * KAt * UInputMat * spdiags(1./sing,0,length(sing),length(sing)) * UInputMat';
    S = W * Y;
end


if isfield(opts,'flag_debug'), flag_debug = opts.flag_debug; else flag_debug = false; end
if flag_debug
    invAlpha_arr = zeros(1,maxIter);
    invBeta_arr =  zeros(1,maxIter);
else
    invAlpha_arr = [];
    invBeta_arr = [];
end

invSigma_e = pinv(Sigma_e);

% S = W * Y;      %DISP new!!
% Y
% keyboard

for ite=1:maxIter
    ite
    
    %% Update hyperparameter alpha and beta
    [invAlpha_ invBeta_] = updateHyperparameters(Y,S,A,K,invK,W,KAt,invSigma_e,invAlpha,flag_update);
    if flag_update.alpha
        invAlpha = invAlpha_;
    end
    if flag_update.beta
        invBeta = invBeta_;
    end
    
    if flag_debug, invAlpha_arr(ite) = invAlpha; invBeta_arr(ite) = invBeta; end
    
    %% Calculate weights W and sources S
    inputMatrix = invAlpha * AKAt + invBeta*Sigma_e;
    [UInputMat singInputMat]=svd(inputMatrix);
    sing = diag(singInputMat);
    sing = sing(sing>sing(1)*thresh);
    UInputMat = UInputMat(:,1:length(sing));
    W = invAlpha * KAt * UInputMat * spdiags(1./sing,0,length(sing),length(sing))  * UInputMat';
%     inputMatrix = invAlpha * AKAt + invBeta*Sigma_e;
%     W = invAlpha * KAt * pinv(inputMatrix);
    if flag_update.S
        S = W * Y;
    end
%     disp('Now reading for printing W')
%     keyboard
end

if flag_debug
    figure
    plot(1:maxIter,invAlpha_arr), xlabel('ite'),ylabel('\alpha^{-1}')
    grid on
    
    figure
    plot(1:maxIter,invBeta_arr), xlabel('ite'),ylabel('\beta^{-1}')
    grid on
end

end


function [invAlpha_ invBeta_ ] = updateHyperparameters(Y,S,A,K,invK,W,KAt,invSigma_e,invAlpha,flag_update)
%==========================================================================
% Filename: updateHyperparameters.m
%
% Description:  Hyperparameters are updated according to EM-steps.
%
% Input:        A: Lead field matrix / Forward model [ Nc x Nd ]
%               opts:
%
% Output:       invAlpha: Hyperparameter for the sources [ scalar ];
%               invBeta: Hyperparameter for the noise variance [ scalar ];
%==========================================================================
% flag_update.hyp_trace = false;

[Nc Ns] = size(Y);
Nd = size(S,1);

if flag_update.hyp_trace

    %---------------------
    % Solution 1

    %%%%%%%%%%%%%%%%%%%%%%%
    % Update invAlpha
    %%%%%%%%%%%%%%%%%%%%%%%
    SigmaS_SSt = Ns*(invAlpha * K - W*KAt'*invAlpha) + S*S';
    invAlpha_ = trace(invK*SigmaS_SSt)/(Nd*Ns);

    % Sigma_S = invAlpha * K - W*KAt'*invAlpha;
    % invAlpha_2 = trace(K*(Ns*Sigma_S + S*S'))/(Nd*Ns);


    %%%%%%%%%%%%%%%%%%%%%%%
    % update invBeta
    %%%%%%%%%%%%%%%%%%%%%%%
    B = A'*invSigma_e*A;
    Ex_stXs = trace(B*SigmaS_SSt);
    % Ex_stXs = trace(B*(Ns*Sigma_S + S*S'));
    Ey2 = trace(invSigma_e*(Y*Y'-2*A*S*Y')) + Ex_stXs;
    invBeta_ = Ey2/(Nc*Ns);

else

    
    %---------------------
    % Solution 2

    %%%%%%%%%%%%%%%%%%%%%%%
    % Update invAlpha2
    %%%%%%%%%%%%%%%%%%%%%%%
    Sigma_S = invAlpha * K - W*A*K*invAlpha;
    tempSmean = 0;
    for t=1:Ns
        tempSmean = tempSmean + S(:,t)'*invK*S(:,t);
    end
    Es = Ns*trace(invK*Sigma_S) + tempSmean;
    invAlpha_ = Es/(Ns*Nd);


    % Sigma_S = invAlpha * K - W*KAt'*invAlpha;
    % % invAlpha_2 = trace(K*(Ns*Sigma_S + S*S'))/(Nd*Ns);
    % invAlpha_2 = trace(invK*(Ns*Sigma_S + S*S'))/(Nd*Ns);
    %%%%%%%%%%%%%%%%%%%%%%%
    % update invBeta
    %%%%%%%%%%%%%%%%%%%%%%%
%     invSigma_e = pinv(Sigma_e);
    B = A'*invSigma_e*A;
    Ey = 0;
    for t=1:Ns
        Ey = Ey + Y(:,t)'*invSigma_e*Y(:,t) -...
            2*Y(:,t)'*invSigma_e*A*S(:,t) + S(:,t)'*B*S(:,t);
    end
    Ey = Ey + Ns * trace(B*Sigma_S);
    invBeta_ = Ey/(Nc*Ns);

end


end