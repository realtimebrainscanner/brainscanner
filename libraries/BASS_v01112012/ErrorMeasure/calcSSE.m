function SSE = calcSSE(S,Sref)
%==========================================================================
% Filename: calcSSE.m (function).
% 
% Description:  Calculates the log of the sum of squared errors on the
%               estimated dipole activities normalized with respect to the
%               simulated dipole activity:
%
% Usage:        SSE = calcSSE(S,Sref)
%
% Input:        S:      Estimated dipole activities.
%               Sref:  Simulated dipole activities.
%
% Output:       SSE:    Sum of squared error score.
%
% History:
%   - Created:  30/05/2008
%   - Modified: 15/04/2009: No longer the log
%
% Special remarks:
%       [1] Daunizeau et al. 2007
%
% Copyright (C) Carsten Stahlhut, DTU IMM 2008
%==========================================================================

% SSE = log10( sum((S(:)-Sref(:)).^2) / (Sref(:)'*Sref(:)) );

SSE = sum((S(:)-Sref(:)).^2) / sqrt((S(:)'*S(:)) * (Sref(:)'*Sref(:)) );