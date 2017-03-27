function RMSE = calcRMSE(STrue,S,opts)
%==========================================================================
% Filename: calcRMSE.m (function).
% 
% Description:  Calculates the root mean squared error, which measures the
%               discrepancy bwtween the original and reconstructed source
%               distributions without any global scaling effect, [1].
%               Measure:
%
%               RMSE = || Strue/s1 - S/s2 ||_fro / || Sture/s1 - 1 ||_fro
%
%               where s1 represents the amplitude of the largest (absoulute
%               value) dipole in STrue and similar for s2 (just in S).
%               source amplitudes. Only nine time points around the maximum
%               of the original data are considered, according to [1].
%
% Usage:        RMSE = calcRMSE(Strue,S,opts)
%
% Input:        STrue: Simulated dipole activities, Nd x Nt (dipole x time)
%               S:     Estimated dipole activities, Nd x Nt (dipole x time)
%               opts
%                   .win: Select window positioned around maximum absolute
%                         amplitude in original data. Choices: 'full' or
%                         'select'. Default is 'select'
%                   .twin: Specify number of time point around the
%                          maximum. This field only relevant if
%                          .win='select'. Default is .twin = 9.
%
% Output:       RMSE:   Root mean square error normalized such that the
%                       zero-solution (S=0) leads to RMSE = 1.
%
% History:
%   - Created:  30/05/2008
%   - Modified:  
%
% Special remarks:
%               [1]: C. Phillips, J. Mattout, M.D. Rugg, P. Maquet, and
%                    K.J. Friston, "An empirical Bayesian solution to the
%                    source reconstruction problem in EEG", NeuroImage 24,
%                    997-1011, 2005.
%
% Copyright (C) Carsten Stahlhut, DTU IMM 2008
%==========================================================================

try, win = opts.win; catch win = 'select'; end
try, twin = opts.twin; catch twin = 9; end

[Nd, Nt] = size(STrue);


if strcmp(win,'select')
    [s1 is1] = max(abs(STrue(:)));
    tmax = ceil(is1/Nd);

    t_win = 9;
    if (tmax> twin)
        tstart = tmax-twin;
    else
        tstart = 1;
    end

    if (tmax > Nt-twin)
        tstop = Nt;
    else
        tstop = tmax+twin;
    end

    % Select a window positioned around tmax.
    STrue = STrue(:,tstart:tstop);
    S = S(:,tstart:tstop);
end

s1 = max(abs(STrue(:)));
s2 = max(abs(S(:)));
if (s2>eps)                 % Check for zero-solution.
    S = STrue/s1 - S/s2;    % S changed to keep difference between true and est.
else
    S = STrue/s1 - 1;       % Zero-solution, S=0 => S* = eps/eps =1
end


%Normalize such that zero-solution S=0 results in RMSE=1. If S=0 =>
%S*=eps/eps = 1, i.e. all elements in S* are 1.
Stemp = STrue/s1 -1;            %Normalize with zero-solution eps/eps=1
s3 = Stemp(:)'*Stemp(:);
RMSE = sqrt( S(:)'*S(:) /s3 );  %Note S keeps the difference.
