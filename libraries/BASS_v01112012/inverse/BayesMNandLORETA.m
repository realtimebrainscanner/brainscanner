function [S W invAlpha invBeta,opts] = BayesMNandLORETA(Y,A,opts)
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
%
% Output:       S: New source estimates based on updated invAlpha, invBeta,
%                  and W.
%               W: Updated weights.
%               invAlpha: Hyperparameter for the sources [ scalar ];
%               invBeta: Hyperparameter for the noise variance [ scalar ];
%
% Last revision: 2011/09/30
%
% Special remarks:
%
% Authors: Carsten Stahlhut
%         Arkadiusz Stopczynski
%
% Copyright (c) DTU Informatics, Technical University of Denmark
%==========================================================================

[Nc Nd] = size(A);

%% Check for optional variables
if isfield(opts,'maxIter'), maxIter = opts.maxIter; else maxIter = 1; end
if isfield(opts,'K'), K = opts.K; else K = speye(Nd); end
if isfield(opts,'invK'), invK = opts.invK; else invK = pinv(full(K)); opts.invK = invK; end
if isfield(opts,'Sigma_e'), Sigma_e = opts.Sigma_e; else Sigma_e = speye(Nc); end
if isfield(opts,'invAlpha'), invAlpha = opts.invAlpha; else invAlpha = 1; end
if isfield(opts,'invBeta'), invBeta = opts.invBeta; else invBeta = 1; end
if isfield(opts,'AKAt'), AKAt = opts.AKAt; else AKAt = A*K*A'; opts.AKAt = AKAt; end
if isfield(opts,'KAt'), KAt = opts.KAt; else KAt = K*A'; opts.KAt = KAt; end
if (isfield(opts,'W') && isfield(opts,'S'))
    W = opts.W;
    S = opts.S;
else
    inputMatrix = invAlpha * AKAt + invBeta*Sigma_e;
    W = invAlpha * KAt * pinv(inputMatrix);
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


S = W * Y;      %DISP new!!
% Y
% keyboard

for ite=1:maxIter
    %% Update hyperparameter alpha and beta
    [invAlpha invBeta] = updateHyperparameters(Y,S,A,K,invK,W,KAt,Sigma_e,invAlpha);
    if flag_debug, invAlpha_arr(ite) = invAlpha; invBeta_arr(ite) = invBeta; end
    
    %% Calculate weights W and sources S
    inputMatrix = invAlpha * AKAt + invBeta*Sigma_e;
    W = invAlpha * KAt * pinv(inputMatrix);
    S = W * Y;
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


function [invAlpha_ invBeta_ ] = updateHyperparameters(Y,S,A,K,invK,W,KAt,Sigma_e,invAlpha)
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

[Nc Ns] = size(Y);
Nd = size(S,1);

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
B = A'*Sigma_e*A;
Ex_stXs = trace(B*SigmaS_SSt);
% Ex_stXs = trace(B*(Ns*Sigma_S + S*S'));
Ey2 = trace(Sigma_e*(Y*Y'-2*A*S*Y')) + Ex_stXs;
invBeta_ = Ey2/(Nc*Ns);


% %---------------------
% % Solution 2
%
%% update invAlpha:
% % temp = W*A;
% % temp2 = temp*K;
% % temp3 = temp2*invAlpha;
% % temp4 = invAlpha*full(K);
% % temp5 = temp4-temp3;
% %
%%%%%%%%%%%%%%%%%%%%%%%
% Update invAlpha
%%%%%%%%%%%%%%%%%%%%%%%
% Sigma_S = invAlpha * K - W*A*K*invAlpha;
% tempSmean = 0;
% for t=1:Ns
%     tempSmean = tempSmean + S(:,t)'*K*S(:,t);
% end

% tempSmean2 = 0;
% for j=1:Nd
%     for i=1:Nd
%         for t=1:Ns
%             tempSmean2 = tempSmean2 + S(i,t)*K(i,j)*S(j,t);
%         end
% %         keyboard
%     end
% end

% % Solution 3 - MN updates K = eye(Nd)
% % SIGMAfulldiag = repmat(SIGMAfulldiag(:,1),[1 Nt]);
% SIGMAfulldiag = repmat(diag(Sigma_S),[1 Ns]);
% v1 = (SIGMAfulldiag + S.^2)/2;
% invAlpha_3 = sum(2*sum(v1,2))/(Ns*Nd);    %Changed 25/2-09
% [invAlpha_ invAlpha_3]

% invAlpha_new = (Ns*trace(K*Sigma_S) + tempSmean)/(Nd*Ns)
%
%
% %%%%%%%%%%%%%%%%%%%%%%%
% % update invBeta
% %%%%%%%%%%%%%%%%%%%%%%%
% Ey= Ns * trace(B*Sigma_S);
% for t=1:Ns
%     Ey = Ey + y(:,t)'*Sigma_e*y(:,t) -...
%         2*y(:,t)'*Sigma_e*A*S(:,t) + S(:,t)'*B*S(:,t);
%
% end
% invBeta_new = Ey/(Nc*Ns)

end