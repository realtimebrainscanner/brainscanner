function [] = waitbar_update(h,t0,maxIter,ite)
%==========================================================================
% Filename: waitbar_update.m (function).
% 
% Description: Update waitbar. Note to create a waitbar use the function
%              waitbar_create.m
%
% Input:    
%             h: Handle to waitbar (obtained from waitbar_create.m)
%            t0: Timer object (obtained from waitbar_create.m)
%       maxIter: Maximum number of iteration
%           ite: Current iteration
%
% Output:      
%
% Authors:
%   Carsten Stahlhut
%   Technical University of Denmark, DTU Informatics (c) 2010
%==========================================================================

Totalsec = etime(clock,t0)/ite*(maxIter-ite);
TimeLeftHour = floor(Totalsec/3600);
TimeLeftMin = floor(mod(Totalsec,3600)/60);
TimeLeftSec = Totalsec - (TimeLeftHour*3600+TimeLeftMin*60);

waitbar(ite/maxIter,h,...
sprintf('Time left: %3.0f h %02.0f min %02.0f sec.  -  %2.0f %c done.',...
TimeLeftHour,TimeLeftMin,TimeLeftSec,floor(100*ite/maxIter),'%'))
