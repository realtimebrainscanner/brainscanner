function data_out = contrasts(data_in,C)
%==========================================================================
% Filename: contrasts.m (function).
% 
% Description:  
%
% Input:        data_in: input data
%                     C: Specify contrast
%
% Output:       data_out: output data with contrast applied
%
% History:
%   - Created:  20/6/2010
%   - Modified: 
%
% Special remarks: 
%
% Copyright (C) Carsten Stahlhut, DTU Informatics 2010
%==========================================================================

% ic = setdiff(D.meegchannels,D.badchannels);
% data_sens = D(ic,:,:);
Ncont = size(C,1);

if iscell(data_in)
    Ntrl = length(data_in);
    [Nrows,Ncol] = size(data_in{1});
    data_out = cell(1,Ncont);
    [data_out{:}] = deal(sparse(Nrows,Ncol));
    for i=1:Ncont
        for itrl=1:Ntrl
            data_out{i} = data_out{i} + data_in{itrl}*C(i,itrl);
        end
    end
else
    [Nc,Ns,Ntrl] = size(data_in);
    data_out = zeros(Nc,Ns,Ncont);
    for i=1:Ncont
        for itrl=1:Ntrl
            data_out(:,:,i) = data_out(:,:,i) + data_in(:,:,itrl)*C(i,itrl);
        end
    end
end
