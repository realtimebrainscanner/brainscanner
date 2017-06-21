 function [gamma_mean,gamma_median,error_val] = teVGGD_wcross(A,Y,opts)
%
% Description:      Performs cross-validation to estimate the sparsity
%                   parameter, gamma, in the time-expanded Variational
%                   Garrote (teVG) which solves the augmented inverse
%                   problem: Y = A * s * X.
%
% Input:            Y:  Data matrix of size KxT, where T is the number of
%                       measurements.
%                   A:  Design matrix/forward model of size KxN
%                   opts: see 'Settings'
%
% Output:           gamma_mean: Mean optimal gamma.
 %                  gamma_median: Median optimal gamma.
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
try min_gamma = opts.min_gamma; catch; min_gamma = -100; end; % Minimum gamma value
try max_gamma = opts.max_gamma; catch; max_gamma = -10; end; % Maximum gamma value
try n_gamma = opts.n_gamma; catch; n_gamma = 10; end; % Number of gamma values
try Cf = opts.Cf; catch;  Cf=4; end; % Number of folds in the cross-validation
try opts.max_iter = opts.max_iter; catch; opts.max_iter = 100; end; % Maximum number of iterations

gamma_all = linspace(min_gamma,max_gamma,n_gamma); % Gamma sequence
K = size(A,1);

% Cross-validation indices
rng(5)% Seed for reproducibility
indices = crossvalind('Kfold', K, Cf);

% Preallocation
error_val = zeros(Cf,n_gamma);

%---------Cross-validation loop---------------
for cf = 1:Cf
    %fprintf('Running cross fold %d of %d\n',cf,Cf)
    A_val  = center(A(indices==cf,:));
    Y_val = center(Y(indices==cf,:));
    A_train = center(A(indices~=cf,:));
    Y_train = center(Y(indices~=cf,:));
    for i = 1:n_gamma;
        V = teVGGD(A_train,Y_train,gamma_all(i),opts);
        % Calculate normalized mean squared error on the validation set
        error_val(cf,i) = mean(mean((A_val*(V)-Y_val).^2))/mean(mean(Y_val.^2));
    end
    
end
% Find minimum
[~, imin] = min(error_val,[],2);
gamma_median = median(gamma_all(imin));
gamma_mean = mean(gamma_all(imin));
