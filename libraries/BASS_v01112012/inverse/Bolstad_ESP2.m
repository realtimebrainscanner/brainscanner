function [Theta X lambda CardA pwiseApprox lambda_arr lmax ThetaAll] = Bolstad_ESP2(Y,H,opts)
%==========================================================================
% Filename: Bolstad_ESP2.m (function).
% 
% Usage: [Theta X lambda CardA pwiseApprox lambda_arr lmax ThetaAll] =
%           Bolstad_ESP2(Y,H,opts)
%
% Description: This method is based on Bolstad et al.'s paper about event
%              sparse penalty (ESP) regularization. The ESP seeks a
%              solution that is composed by a small number of
%              "space-time-events" (STEs). "Each STE is a spatio-temporal
%              signal defined in terms of a group of basis functions" [1].
%              In contrast to the mfile Bolstad_ESP this version is vector
%              based when updating Theta's eq.(13) in [1].
%
%                   min_Theta || Y - H*S*Theta*T' || + lambda sum_ij ( ||
%                   vec Theta ||_p )
%
% Input:    
%           Y: Observations             [M x T]
%           H: Lead-field matrix        [M x G]
%           opts:
%               .maxIter: Maximum number of interations
%               .Theta: Init valu
%               .spatial.S: Spatial bases - default indentity matrix (cell
%                           array)
%               .spatial.I: Number of spatial submatrices;
%               .temporal.T: T_j temporal submatrix with an orthogonal basis for a
%                            temporal subspace.
%               .temporal.J: Number of J temporal submatrices.
%               .e: default 0.01
%               .lambda: default according to eq.(15) in [1].  
%
% Output:   Theta: Spatio-temporal dictionary at the given lambda value.
%               X: S*Theta*T' (Source solution)
%          lambda: Optimal lambda value according to heuristic approach
%           CardA:
%     pwiseApprox: Piecewise approximation graphs at the different lambda
%                  values in lambda_arr. Each row is associated with a
%                  lambda value in lambda_arr.
%      lambda_arr: Penalty array
%            lmax: Maximum lambda value used
%        ThetaAll: Cell-array with All spatio-temporal dictionaries
%                  associated with lambda_arr.
%
% Special remarks:
%
%           M: Number of channels
%           N: Number of samples
%           D: Number of sources
%           I: Number of spatial submatrices
%           J: Number of temporal submatrices
%
%           Note when S and T are set to default this methods is actually
%           performing ME estimates.
%
%           Depending on how the temporal submatrices are specified we
%           can have choose to have sparsity in temporal basis functions
%           (i.e. J = #temporal basis functions and submatrix T_j is then a
%           vector with a single temporal basis function) or put all
%           temporal basis functions into a single submatrix (i.e. J=1 -
%           T_j submatrix with all the temporal basis functions)
%
% References:
% [1] A. Bolstad, B. V. Veen, and R. Nowark, "Space-time event sparse
%     penalization for magneto-/electroencephalography", NeuroImage
%     46, pp.1066-1081, 2009.
%
% Author:
%   Carsten Stahlhut
%   Technical University of Denmark, DTU Informatics (c) 2010
%==========================================================================

% *** Control parameters ***
MIN_DMU         = 1e-12;  

%Create waitbar
[handle_waitbar,timer_t0] = waitbar_create(mfilename);

[M N] = size(Y);
D = size(H,2);

try maxIter = opts.maxIter; catch maxIter = 10000; end;
try
    S = opts.spatial.S;
    I = opts.spatial.I;
catch
    S = speye(D);
    I = D;
end;

try
    T = opts.temporal.T;
    J = opts.temporal.J;
catch
    T = speye(N);
    J = 1;
end;
try e = opts.e; catch e = 0.01; end;

try sigma2 = opts.sigma2; catch sigma2 = 1; end;

%% Start EM iterations

normTT = norm(full(T*T'));
HS = H*S;
normHSSH = norm(HS*HS');
c = 1/(normTT*normHSSH);
alpha2 = sigma2*c * 0.99;   %ensure that alpha2 <= sigma2*c

NcolS = size(S,2);
NcolT = size(T,2);
Q = NcolS/I;
P = NcolT/J;
W = Q*P;        %Number of elements in submatrix Theta_ij
K = I*J;        %Number of submatrices in Theta (#spatial basis times #temporal basis)

R = 1:K;
R = reshape(R,[I,J]);
R = kron(R,ones(Q,P));
[dummy iR] = sort(R(:),'ascend');

try Theta = opts.Theta; catch Theta = zeros(NcolS,NcolT); end;
theta = Theta(iR); theta = reshape(theta,[W K]);

try
    lambda = opts.lambda;
    lmax = max(lambda);
catch
    temp_lam = zeros(I,J);
    ir = 0;
    for i=1:I
        for j=1:J
            ir = ir+1;
            iiS = (i-1)*Q+1:i*Q;
            iiT = (j-1)*P+1:j*P;
            temp_lam(i,j) = norm(S(:,iiS)'*H'*Y*T(:,iiT),'fro');
%             temp_lam(i,j) = norm(S(iR==ir)'*H'*Y*T(iR==ir),'fro');
        end
    end
    lmax = 2*max(temp_lam(:));
    reg_array = 0.02:e:1-e;
    reg_array = 0:e:1-e;
    lambda = fliplr(reg_array)*lmax;
end
lambda_arr = lambda;
Nlam = length(lambda_arr);

% rho = zeros(1,maxIter);
CardA = NaN*ones(1,Nlam);       %The cardinal of STEs
CardA0 = rank(H)/W;             %W: #elements in Theta_ij
% CardA0 = rank(H)/K;             %K: #submatrices
% CardA0 = rank(H)/(NcolS*NcolT); %NcolS*NcolT = #elements in Theta

ThetaAll = cell(1,Nlam);
[ThetaAll{:}] = deal(zeros(size(Theta)));

%% Now iteratively update Theta and Z, while decreasing the penalty weight
%% parameter lambda
for ilam=1:Nlam
    lambda = lambda_arr(ilam);
    
    for ite=1:maxIter

        %--------------------------------------------
        % Update Z according to eq.(14) in [1]
        %--------------------------------------------
        res = Y - H*S*Theta*T';
        Z = Theta + c*S'*H'*res*T;
%         rho(ite) = sqrt(res(:)'*res(:));

        %--------------------------------------------
        % Update Theta according to eq.(13) in [1]
        %--------------------------------------------
        g = alpha2*lambda/(2*sigma2);

        Zvec = Z(iR);
        Zvec = reshape(Zvec,[W K]);
        norm_z = sqrt(sum(Zvec.^2,1));
        index = find(norm_z > g);

        Theta_old = Theta;
        theta(:) = 0;
        theta(:,index) = repmat((1 - g./norm_z(index)),[W 1]).*Zvec(:,index);
        Theta(iR) = theta(:);


        dmu = max(max( abs(Theta_old(:) - Theta(:) ) ));
        if (dmu < MIN_DMU)  break;  end;

        %Update Waitbar
        waitbar_update(handle_waitbar,timer_t0,maxIter*Nlam,ite+(ilam-1)*maxIter)
    end
    
    %Piecewise approximation for penalty weight parameter selection
    CardA(ilam) = length(index);
%     aQ = (CardA(ilam)-CardA0)/lambda^2;
%     pwiseApprox(ilam,ilam:Nlam) = aQ*lambda_arr(ilam:Nlam).^2 + CardA0;
    
    ThetaAll{ilam} = Theta;
end

%% Penalty weight parameter selection based on linear/quadratic approx
% Make piecewise approximation
pwiseApprox = zeros(Nlam,Nlam);
for ilam=1:Nlam
    %Linear term
    aL = CardA(ilam)/(lambda_arr(ilam)-lmax);
    pwiseApprox(ilam,1:ilam) = aL*lambda_arr(1:ilam) - aL*lmax;

    %Quadratic term
    aQ = (CardA(ilam)-K)/(lambda_arr(ilam)^2-0);
    bias = K;
    pwiseApprox(ilam,ilam:Nlam) = aQ*(lambda_arr(ilam:Nlam).^2) + bias;    
end

%% Find lambda value based on the piecewise approximation that leads to the
%% smallets squared error relative to the CardA.
E = repmat(CardA,[Nlam 1]) - pwiseApprox;
E = sum(E.^2,2);
[dummy ilam_op] = min(E);
lambda = lambda_arr(ilam_op);
Theta = ThetaAll{ilam_op};

X = S*Theta*T';
% res = Y - H*X;
% rho(ite+1) = sqrt(res(:)'*res(:));

close(handle_waitbar)

figure
plot(lambda_arr/lmax,CardA)
hold on
plot(lambda_arr/lmax,pwiseApprox(ilam_op,:),'--r')
xlabel('\lambda/\lambda_{max}')
ylabel('#STEs in estimate')
grid on
set(gca,'XDir','reverse')

% keyboard
end
