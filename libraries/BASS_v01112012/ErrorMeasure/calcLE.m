function [LE iSmax] = calcLE(S,vert,vert_ref,distMeas)
%==========================================================================
% Filename: calcLE.m (function).
% 
% Description:  Calculates the localization error between the maximum
%               dipole in S and the dist.
%
% Usage:        LE = calcLE(S,vert,vert_ref,distMeas)
%
% Input:        S: The estimated activity
%               vert: Indicies on dipoles the localization error (LE) should
%                   be measured with respect to. vi should keep the most
%                   active dipoles.
%               vert_ref: Type of distance measure: 'Euclidean' or
%                         'Geodesic'. Currently only 'Euclidean' is
%                         available.
%
% Output:       LE:  Localization error in [mm].
%
% History:
%   - Created:  15/04/2009
%   - Modified: 
%
% Special remarks:
%               distMeas not used at the moment. A new version of this
%               function should also include the Geodesic distance.
%
% Copyright (C) Carsten Stahlhut, DTU IMM 2008
%==========================================================================

dummy = max(S,[],2);
[dummy iSmax] = max(dummy,[],1);

%%Euclidean distance
LE = sqrt( sum( (vert(iSmax,:)-vert_ref).^2  ) );