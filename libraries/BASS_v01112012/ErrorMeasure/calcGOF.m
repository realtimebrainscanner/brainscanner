function GOF = calcGOF(M,Mrecon)
%==========================================================================
% Filename: calcGOF.m (function).
% 
% Description:  Calculates the Goodness of fit
%
% Usage:        GOF = calcGOF(M,Mrecon)
%
% Input:        M:     
%               Mrecon:
%
% Output:       GOF: Goodness of fit
%
% History:
%   - Created:  31/05/2009
%   - Modified:  
%
% Special remarks:
%
% Copyright (C) Carsten Stahlhut, DTU Informatics 2009
%==========================================================================

GOF = 1-sum((Mrecon(:)-M(:)).^2)/...
        (sum(M(:).^2));    %Goodness of fit