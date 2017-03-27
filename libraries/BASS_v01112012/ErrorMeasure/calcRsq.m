function Rsq = calcRsq(STrue,S,opts)
%==========================================================================
% Filename: calcRsq.m (function).
% 
% Description:  Calculates the coefficient of determination also known as
%               the R^2-value. Measure used to perform goodness of fit
%               evaluation.
%
%               This function may be used in two different ways;
%                   1. calculating the temporal accuracy of each dipole,
%                      i.e. a R^2-value is calculated for each dipole.
%                   2. calculating the temporal accuracy between the N-most
%                      active dipoles, i.e. a N x N matrix with R^2-values
%                      is calculated. By combining this measure with MaxLE
%                      (maximum localization error) one obtains both
%                      spatial and temporal accuracy of the N-most active
%                      dipoles.
%
% Usage:        Rsq = calcRsq(STrue,S)
%
% Input:        STrue:  Simulated dipole activities, (Nd x Nt).
%               S:      Estimated dipole activities, (Nd x Nt).
%               opts
%                   .RsqMethod: 'all' (compare each simulated dipole with
%                                      the corresponding estimated, based
%                                      on dipole number)
%                               'Nmost' (N-most active dipole, default).
%                   .N: Number of most active dipoles to be compared (only
%                       relevant if .RsqMethod = 'Nmost'. Default is .N=1.
%
% Output:       Rsq:    R^2 value (coefficient of determination), (Nd x 1)
%                       corresponding to a R^2-value for each dipole. It
%                       thes
%
% History:
%   - Created:  31/05/2008
%   - Modified:  
%
% Special remarks:
%               Nd: Number of dipoles.
%               Nt: Number of time points.
%
% Copyright (C) Carsten Stahlhut, DTU IMM 2008
%==========================================================================


try, RsqMethod = opts.RsqMethod; catch RsqMethod = 'Nmost'; end
try, N = opts.N; catch N = 1; end

% NB er ikke faerdig!!!
% benyt evt corrcoef til at beregne correlation coefficient.

ypred = polyval(P2,var1); 
dev = var2 - mean(2); 		
SST = sum(dev.^2); 		
resid = var2 - ypred; 		
SSE = sum(resid.^2);
normr = sqrt(SSE);	% residual norm
Rsq = 1 - SSE/SST; 