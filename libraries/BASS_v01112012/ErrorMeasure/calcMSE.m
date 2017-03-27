%==========================================================================
% Filename: calcMSE.m (function).
% 
% Description:  Calculates the mean squared error
%
% Usage:        MSE = calcMSE(Strue,S,opts)
%
% Input:        S:     Estimated dipole activities, Nd x Nt (dipole x time)
%               Sref:  Simulated dipole activities, Nd x Nt (dipole x time)
%               opts
%                   .win: Select window positioned around maximum absolute
%                         amplitude in original data. Choices: 'full' or
%                         'select'. Default is 'select'
%                   .twin: Specify number of time point around the
%                          maximum. This field only relevant if
%                          .win='select'. Default is .twin = 9.
%
% Output:       MSE:   Root mean square error normalized such that the
%                       zero-solution (S=0) leads to RMSE = 1.
%
% History:
%   - Created:  12/03/2009
%   - Modified:  
%
% Special remarks:
%               [1]: Grova et al. (2006): Note however that they in this
%               paper only define it for one time sample
%
% Copyright (C) Carsten Stahlhut, DTU Informatics 2009
%==========================================================================
function MSE = calcMSE(S,Sref,tselect)

try dummy = tselect; catch tselect = 1:size(S,2); end

MSE = sum(sum((S(:,tselect)-Sref(:,tselect)).^2))/...
    sum(sum(Sref(:,tselect).^2));