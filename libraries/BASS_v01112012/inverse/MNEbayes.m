%==========================================================================
% Filename: MNEbayes.m (function).
% 
% Description:  This algorithm uses expectation maximization (EM) to
%               estimate the distribution of the sources of the
%               EEG-signals. Assumes certain system matrix.
%               Model: - M: Gaus likelihood with mean AS and covariace
%                           with different Beta parameters for each channel
%                           and the variance at a given channel is assumed
%                           time-invariant.
%                      - S: Gaus prior with different Alpha parameters for
%                                 each source and time.
%                      - beta: Control the noise variance for the model.
%                              Delta functions (i.e. ML estimates).
%                      - Alpha: A hyperparameter on each source at each
%                               time index to control the relance of that
%                               given source (ARD).
%                               Delta functions (i.e. ML estimates).
%
% Input:        M:      Target values (Recorded EEG signal).
%               A:      Forward field.
%               param:  Structure with user specified parameters.
%
% Output:       S
%
% History:
%   - Created:  30/12/2007
%   - Modified: 12/05/2008: Calculation of relative A_error moved to out of
%                           this function. The A_error in this function is
%                           now the absolute error.
%               27/05/2008: Calculation of LB - problems with
%               log(det(LAMBDA)), changed to sum(log(svd(LAMBDA)))
%               10/03/2011: Renamed and part of BASS toolbox
%
% Special remarks:
%               Original filename: algorithm_EM3full_MNE
%               To be done: change global variable to locally and let param
%               be optional variable, i.e. use try/catch or isfield
%
% Copyright (C) Carsten Stahlhut (cs), DTU Informatics 2011
%==========================================================================
function varout = MNEbayes(M, A0, param, varargin)
global v2 Nt Nc B idiagB beta

slash_op = '/';
if isfield(param,'convcrit'), convcrit=param.convcrit; else convcrit=1e-8; end

if nargin<5
    gamma = 0.1;
    if nargin<2
        error('Not enough input arguments - see help algorithmIV for more information.')
    end
end

% A = A0;
[Nc, Nd] = size(A0);
Nt = size(M,2);
iActive = true(1,Nd);                   %Active sources.
NdActive = sum(iActive);                %Number of active sources.

scaling = 10^ceil( abs( log10( var(M(:)) ) )/2);
scaling = 1;
M = M*scaling;

%==========================================================================
% alpha, beta, og gamma boer vel ikke vaelges som parameter som skal
% initialiseres hvis der anvendes hyperprior, idet de herved vaelges ud fra
% parametrene, dvs.:
%==========================================================================
if param.hyp.init_alpha_true == 1
%     alpha = param.hyp.alpha;              %Hyperparameter.
    alpha = param.hyp.alpha(:,1);
else
    range = 10^ceil( abs( log10( var(M(:)) ) )/2);
    
    alpha = 1;                  %Random choice.
    
    ny_alpha = param.hyp.ny_alpha;          %Only for init.
    zeta_alpha = param.hyp.zeta_alpha;      %Only for init.
    alpha = ny_alpha/zeta_alpha;
    alpha = alpha/((range/scaling)^2);        %CS added 03/03/2009
end

if param.hyp.init_beta_true == 1
    beta = param.hyp.beta;              %Parameter.
else
%     beta = 100/var(M(:));                %According to Tipping.
    
    %CS changed 25/02/2009
    beta = 1/var(M(:));                           %Random choice
end

% gamma = param.hyp.gamma;              %Hyperparameter.
    

if (size(alpha,1) < Nd)
%     keyboard
    alpha = repmat(alpha,[Nd 1]);
end
invDt = sparse(Nd,Nd);
idiag = 1:Nd+1:Nd^2;
invDt(idiag) = 1./alpha;

% if (size(alpha,2) < Nt)
% %     keyboard
%     alpha = repmat(alpha,[1 Nt]);
% end

if (length(beta) < Nc)
%     keyboard
    beta = repmat(beta,[Nc 1]);
end
B = sparse(Nc,Nc);
idiagB = 1:Nc+1:Nc^2;
B(idiagB) = beta;

if param.hyp.BetaCommon
    BetaUpMethod = @betaStepCommon;
else
    BetaUpMethod = @betaStepDifferent;
end

% if (length(gamma) < Nc)
% %     keyboard
%     gamma = repmat(gamma,[Nc 1]);
% end


% % Ic = eye(Nc);
% A_error = NaN*ones(param.Nite_A+1,1);

% A0_En = sum(A0(:).^2);
% A_error(1) = sum((A0(:)-param.ATrue(:)).^2)/A0_En;

% % A_error(1) = sum((A0(:)-param.ATrue(:)).^2);
% A_error(:) = sum((A0(:)-param.ATrue(:)).^2);

%--------------------------------------------------------------------------
% Start algorithm.
%--------------------------------------------------------------------------
c = clock;
disp(sprintf(['==============================================================\n',...
     '%s started %02.0f-%02.0f-%.0f at the time %02.0f:%02.0f:%02.0f \n'],...
    func2str(param.method),c(3),c(2),c(1),c(4),c(5),c(6)))
clear c


% %---------------------
% % Subspace
% %---------------------
% %Reduced space
% [U LAMB V] = svd(A0,'econ');
% Nk = size(V,2);
% VT = V';
% A0r = A0*V;            %Reduce A0 to subspace ( A0 = U*LAMB*V') 
% % A0r = U*LAMB            %Is the same as A0*V, since U*LAMB*V'*V = U*LAMB
MY = A0;
A = A0;
Nk = Nd;

% PSI = eye(Nk);  %In reduced space PSI = V'*(Id/gammaj)*V = Ik/gammaj
% PSI = diag(1./gamma);
% THETA = eye(Nk)/gamma;

% % D = diag(alpha);
% % Id = eye(Nd);
% ExAtBA = 0;
% for j=1:Nc
%     ExAAtj = (PSI + MY(j,:)'*MY(j,:));
%     ExAtBA = ExAtBA + beta(j)*ExAAtj;
% end

% MY = A0';
% keyboard
% L = zeros(param.Nite_A+1,1);
% % L = zeros(param.Nite_A*2+1,1);
% Lmat = zeros(param.Nite_A*2+1,5);
s_error = inf*ones(param.Nite_A,1);
varout.add.ln_evid = inf*ones(param.Nite_A,1);              %Changed 14/12
varout.add.m_error = inf*ones(param.Nite_A+1,1);

% while dE_A < last_error
% while change in lower bound is larger than a specified max-change.
% ETA = zeros(Nk,Nt);
S = zeros(Nd,Nt);
% ExStAAtS = zeros(Nc,1);
% ExAtA = ExStAAtS;
SIGMAfulldiag = zeros(NdActive,Nt);

% varout.add.alpha_ite = zeros(Nd,Nt,param.Nite_A+1);
varout.add.beta_ite = inf*ones(Nc,param.Nite_A+1);
varout.add.gamma_ite = inf*ones(Nc,param.Nite_A+1);

% varout.add.alpha_ite(:,:,1) = alpha;
varout.add.beta_ite(:,1) = beta;
% varout.add.gamma_ite(:,1) = gamma;

h = waitbar(0,['Time left: Unknown.  -  0% done.']);
set(h,'Name',['Processing ' func2str(param.method)])
t0 = clock;
% 
nrun = 100;
% v2array = zeros(param.Nite_A);
% for n=1:param.Nite_A                                %Iterative approach.
maxiter = param.Nite_A;
s_error = inf*ones(maxiter+1,1);
RMSE_S  = nan*ones(maxiter+1,1);

varout.add.LBlong = inf*ones(5,maxiter+1);
Fcost = 1e100;
dFcost = inf;
n = 0;
% for n=1:param.Nite_A                            %Iterative approach.
nmin = 100;
nmin = maxiter

optsRMSE.win = 'full';
while (abs(dFcost/Fcost)>convcrit & n<maxiter) | n<nmin
    n = n+1;

%     keyboard
%     if ~mod(n,100)
%         n
%     end
    %---------------------------------------
    % s-step (E-step in EM algorithm)
    %---------------------------------------
%     temp = ExAtBA +...
%             VT*(repmat(alpha(iActive),[1 Nk]).*V);
    LAMBDA = inv(B)+A*invDt*A';
    invLAMBDA = inv(LAMBDA);
       
%     SIGMA = invDt - invDt*A'*(invLAMBDA*(A*invDt));    
%     ETA = (SIGMA*A')*(B*M);

%     S = invDt*(A'*B*M) - invDt*A'*((invLAMBDA*A*invDt*A')*(B*M));
        S = invDt*A'*invLAMBDA*M;       %CS added 25/02/2009

% %     S(iActive,:) = St(iActive,:);
% %     S(iActive,:) = V*ETA;          %Transformation back to the originally space.
% 
%     s_error(n) = sum((S(:)-param.STrue(:)).^2);
% %      tempRMSE = calcRMSE(param.STrue,S,optsRMSE);
% %     RMSE_S(n) = tempRMSE;
% 
    Mre = A*S;
    varout.add.m_error(n) = sum((Mre(:)-M(:)).^2);
%     clear Mre

    %---------------------------------------------
    %No A-step
    %---------------------------------------------
%     A_error(n+1) = sum((A0(:)-param.ATrue(:)).^2);
%     A_error(n+1) = sum(sum(abs(A-param.ATrue)));

   
%     for t=1:Nt
        for i=1:NdActive
            SIGMAfulldiag(i,1) = invDt(i,i)-invDt(i,i)*...
                A(:,i)'*invLAMBDA*A(:,i)*invDt(i,i);
        end
        SIGMAfulldiag = repmat(SIGMAfulldiag(:,1),[1 Nt]);
%     end

    v1 = (SIGMAfulldiag + S.^2)/2;

%     ExStAAtS = zeros(Nc,1);
% %     sumExsst = (ETA*ETA' + SIGMAsum);     %sum over E[s*s'] for all t.
%     sumExsst = (S*S' + Nt*SIGMA);     %sum over E[s*s'] for all t.
%     for j=1:Nc
%         ExStAAtS(j) = MY(j,:)*sumExsst*(MY(j,:)');
%     end
%     v2old = sum((M.^2)/2 - (A0*S).*M,2) + ExStAAtS/2; %Note A0 is used since A=A0

    AinvDAt = A*invDt*A';
    v2 = 1/2*diag( M*M' - 2*(A*S)*M' +...
        Nt*(AinvDAt - AinvDAt * invLAMBDA * AinvDAt) + Mre*Mre');
%     keyboard
%     v2array(n) = sum(abs(v2old-v2));
%     v3 = ExAtA/2 - sum(MY.*A0r,2) + sum(A0r.*A0r,2)/2;

   

    %----------------------------------------------------------------------
    % Estimation of parameters. (M-step in EM-algorithm)
    %----------------------------------------------------------------------
    if param.hyp.hyp_alpha          %Update when ~0.
%         alpha = 1./(2*v1);          %Expected value of alpha.
        alpha = Nt./(2*sum(v1,2));          %Expected value of alpha.
                                    %v1 is multiplied by 2 since v1
                                    %includes divided by 2 in

%         invDt(idiag) = 1./alpha;
%         invalpha = (2*sum(v1,2))/Nt;    %Changed 27/5-08
        invalpha = sum(2*sum(v1,2))/(Nt*Nd);    %Changed 25/2-09
%         minInvAlpha = max(invalpha(:))*1e-6;
%         idummy = invalpha < minInvAlpha;
%         invalpha(idummy) = minInvAlpha;
        invDt(idiag) = invalpha;
        alpha = 1./invalpha;
        clear idummy invalpha
    end
    if param.hyp.hyp_beta          %Update when ~0.            
            feval(BetaUpMethod)
    end
    
    varout.add.beta_ite(:,n+1) = beta;
%     varout.add.gamma_ite(:,n+1) = gamma;
    
    %Changed 27/5-08
%     logLike = -Nt*Nc/2*log(2*pi) - Nt/2*log(det(LAMBDA)) -...
%         1/2*trace(invLAMBDA*(M*M'));
    logLike = -Nt*Nc/2*log(2*pi) - Nt/2*sum(log(svd(LAMBDA))) -...
        1/2*trace(invLAMBDA*(M*M'));

%     E = beta(1)*sum(v2(:))*2 + alpha(1)*sum(v1(:))*2;
%     logEv = 0.5*(Nt*( Nd*(log(alpha/(2*pi))) + sum(log(beta/(2*pi))) ) -...
%         ( E ) );
    
    LB = logLike;

%     L(n) = LB;
    varout.add.ln_evid(n) = LB;
    
    dFcost = LB-Fcost;      %Change in lower bound (cost)
    Fcost = LB;             %Update cost

    
    %---------------------------------------------
    % Update waitbar with estimated time left.
    %---------------------------------------------
    Totalsec = etime(clock,t0)/n*(param.Nite_A-n);
    TimeLeftHour = floor(Totalsec/3600);
    TimeLeftMin = floor(mod(Totalsec,3600)/60);
    TimeLeftSec = Totalsec - (TimeLeftHour*3600+TimeLeftMin*60);
    waitbar(n/param.Nite_A,h,...
     sprintf('Time left: %3.0f h %02.0f min %02.0f sec.  -  %2.0f %c done.',...
     TimeLeftHour,TimeLeftMin,TimeLeftSec,floor(100*n/param.Nite_A),'%'))
 
%       if (~mod(n,500))
%         fname_mat = sprintf(['Mat_files' slash_op 'Temp' param.A_choice '_',...
%                             param.fmethod '_Aite%.0f'],param.Nite_A)
%         save(fname_mat)
%         fname_txt = sprintf(['Mat_files' slash_op 'Txt_' param.A_choice '_',...
%                             param.fmethod '_Aite%.0f.txt'],param.Nite_A)
%         fid = fopen(fname_txt,'w');
% %         fprintf(fid,'Iterations: %8.0f\n',n);
%         fprintf(fid,'Iteration number: %8.0f\nTime left: %3.0f h %02.0f min %02.0f sec.  -  %2.0f %c done.',...
%      n,TimeLeftHour,TimeLeftMin,TimeLeftSec,floor(100*n/param.Nite_A),'%');
%         fclose(fid);
%       end

end

%---------------------------------------------------------
% Update s so it corresponds to the newest estimate of A.
%---------------------------------------------------------
    LAMBDA = inv(B)+A0*invDt*A0';
    invLAMBDA = inv(LAMBDA);
       
%     SIGMA = invDt - invDt*A'*(invLAMBDA*(A*invDt));    
%     ETA = (SIGMA*A')*(B*M);
%     S = invDt*(A'*B*M) - invDt*A'*((invLAMBDA*A*invDt*A')*(B*M));
    S = invDt*A'*invLAMBDA*M;       %CS added 25/02/2009

%     S(iActive,:) = V*ETA;          %Transformation back to the originally space.

% s_error(n+1) = sum((S(:)-param.STrue(:)).^2);
% % tempRMSE = calcRMSE(param.STrue,S,optsRMSE);
% % RMSE_S(n+1) = tempRMSE;

Mre = A*S;
varout.add.m_error(n+1) = sum((Mre(:)-M(:)).^2);
clear Mre

% figure,
% plot(1:param.Nite_A,v2array)

% EMcalcLB
% L(n+1) = LB;

varout.Amap = A0;
varout.smap = S/scaling;
varout.add.iter = n+1;
varout.alpha = repmat(alpha,[1 Nt]);
varout.beta = beta;
% varout.gamma = gamma;

% varout.add.s_error = s_error;
% varout.add.RMSE_S = RMSE_S;
% varout.add.A_error = A_error/A_error(1);
% varout.add.A_error = A_error;       %Relative error calculated in main.
% % varout.add.ln_evid = L;     %lower bound
% % varout.add.L_A = L(2:2:end);
% % varout.add.L_A = L;
% % varout.add.Lmat = Lmat;


close(h)
c = clock;
disp(sprintf(['\n%s finished %02.0f-%02.0f-%.0f at the time %02.0f:%02.0f:%02.0f \n',...
    '==============================================================\n'],...
    func2str(param.method),c(3),c(2),c(1),c(4),c(5),round(c(6))))

clear c
% 
% fname_mat = sprintf(['Mat_files' slash_op 'Finished' param.A_choice '_',...
%                             param.fmethod '_Aite%.0f'],param.Nite_A);
% save(fname_mat)

if isfield(param,'CombineMethods')
    figure
        stem(1:length(alpha),alpha,'.r')
        xlabel('Source number, i')
        ylabel('{\alpha}_i')
        title('Estimated values of the hyperparameter {\bf \alpha}')
        grid
%         legend('True','Estimated')
        fname = 'Alpha_valuesEstimated';
        if param.fig.print
            print(gcf,param.fig.dformat,...
                [param.fig.fsubdir fname param.fig.format])
        end

    
    c = clock;
    disp(sprintf(['Combination of ',...
     '%s with %s \nstarted %02.0f-%02.0f-%.0f at the time %02.0f:%02.0f:%02.0f \n'],...
    func2str(param.method),param.CombineMethods{1}(:,:),c(3),c(2),c(1),c(4),c(5),c(6)))

%     LBarray = varout.add.ln_evid;
%     orgfMethod = param.fmethod;
    CombineMethods = param.CombineMethods;
    N_CombMethods = length(CombineMethods);
    
    KeepFields.fmethod = param.fmethod;
    KeepFields.A_choice = param.A_choice;
    KeepFields.fig = param.fig;
    KeepFields.STrue = param.STrue;
    KeepFields.ATrue = param.ATrue;
    KeepFields.MTrue = param.MTrue;
    KeepFields.Nt = param.Nt;


    for iCM = 1:N_CombMethods
        fmethod = CombineMethods{iCM}(:,:);
        param = feval(str2func(fmethod));
        param.fmethod = ['Comb' fmethod];
        
        param.hyp.init_alpha_true = 1;          %Use the estimated alpha and
        param.hyp.init_beta_true = 1;           %beta as initialization in the
        param.hyp.alpha = varout.alpha;                %next algorithm.
        param.hyp.beta = varout.beta;
        param.hyp.init_gamma_true = 1;
%         param.hyp.gamma = varout.gamma;
        A0 = varout.Amap;
        
        param.A_choice = KeepFields.A_choice;
        param.fig = KeepFields.fig;
        param.STrue = KeepFields.STrue;
        param.ATrue = KeepFields.ATrue;
        param.MTrue = KeepFields.MTrue;
        param.Nt = KeepFields.Nt;
        
        param.fig.fsubdir = ['Figures' slash_op param.A_choice,...
            slash_op 'Comb' func2str(param.method) slash_op fmethod slash_op];
    
%         keyboard
        keep KeepFields N_CombMethods CombineMethods param M A0 varargin
    
%         keyboard
        varout = feval(param.method,M,A0,param,varargin);
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBSFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================================================================
function [] = betaStepCommon(varargin)
global v2 Nt Nc B idiagB beta

beta_mean = Nt*Nc/(2*sum(v2));           %Expected value of beta.
beta(:) = beta_mean;
B(idiagB) = beta;
%==========================================================================


                                            
%==========================================================================
function [] = betaStepDifferent(varargin)
global v2 Nt Nc B idiagB beta

beta = Nt./(2*v2);          %Multiply v2 with 2 due to v2 was divided
B(idiagB) = beta;           %with 2 when calculated.

%==========================================================================
