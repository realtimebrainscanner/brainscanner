%==========================================================================
% Filename: rescaling.m
%
% Description: Performs a rescaling of a signal/image to a specified
%              dynamic range.
%           
% Function call:
%           out = rescaling(in,range_out,range_in,param)
%
% Input:    in:         Signal/image to be compressed.
%           range_out:  Dynamic range, which the signal/image should be
%                       compressed to. range_out is a 1x2 vector with the
%                       minimum and the maximum output value.
%                       range_out = [min(out(:)) max(out(:))].
%           range_in:   The dynamic range of the signal/image, which should
%                       be rescaled to range_out. Is a 1x2 vector with the
%                       lower and the higher bounds of the dynamic range.
%
% Output:   out:    The rescaled signal/image.
%
% Last updated: 21/10/2009.
%
% Special remarks:
%
% Copyright (C) Carsten Stahlhut, DTU Informatics 2009
%==========================================================================
function out = rescaling(in,range_out,range_in)

if nargin<3
    range_in(1) = min(in(:));
    range_in(2) = max(in(:));
end

dif_in = diff(range_in);
dif_out = diff(range_out);

out = (in-range_in(1))*(dif_out/dif_in)+range_out(1);

index_min = out < range_out(1);         %Suppression of values out of
index_max = out > range_out(2);         %the specified dynamic range.
out(index_min) = range_out(1);
out(index_max) = range_out(2);