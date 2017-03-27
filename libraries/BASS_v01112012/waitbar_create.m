function [h,t0] = waitbar_create(function_name)
%==========================================================================
% Filename: waitbar_create.m (function).
% 
% Description: Create waitbar. To update the waitbar use the function
%              waitbar_update.m
%
% Input:   function_name: A string with the name of the file that is
%                         beeing preprocessed.
% Output:      
%             h: Handle to waitbar (obtained from waitbar_create.m)
%            t0: Timer object (obtained from waitbar_create.m)
%
% Authors:
%   Carsten Stahlhut
%   Technical University of Denmark, DTU Informatics (c) 2010
%==========================================================================

h = waitbar(0,['Time left: Unknown.  -  0% done.']);
% set(h,'Name',['Processing ' func2str(function_name)])
set(h,'Name',['Processing ' function_name])
t0 = clock;