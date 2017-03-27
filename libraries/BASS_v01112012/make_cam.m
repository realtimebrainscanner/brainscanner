function [CAM, MAP, thresh, levels] = make_cam(im_back,AM,param)
%==========================================================================
% Filename: make_cam.m (function).
% 
% Description:  Generates a thresholded color activity map.
%
% Input:        im_back: Background image to be illustrated, e.g. im_back
%                        could be the curvature when showing activity on a
%                        inflated brain.
%                    AM: Activity (the sources)
%                 param
%                   .thresh: Threshold for activity (default: 0.5 ~ all
%                            activity with amplitudes less than 50% of
%                            the maximum amplitude will be thresholded)
%                     .bMAP: Color specification for background image
%                            (default: 2 gray values)
%                     .aMAP: Color specification for activity map
%                            (default: jet(128) )
%               
%
% Output:       CAM: Thresholded color activity map
%               MAP: Colormap - [bMAP; aMAP]
%            thresh: Applied threshold for activity
%            levels: Number of bMAP levels and aMAP levels [size(bMAP,1)
%                    size(aMAP,1)]
%
% Last updated:  14/12/2010
%
% Author: Carsten Stahlhut
%         Technical University of Denmark, DTU Informatics
%==========================================================================

try
    MAP = param.MAP;
catch
    try 
        bMAP = param.bMAP;
    catch 
        bMAP = gray(10);
        bMAP = [bMAP(2,:); bMAP(7,:)]; 
    end
    try aMAP = param.aMAP; catch aMAP = jet(128); end

    MAP = [bMAP; aMAP];
%     MAP = [bMAP;...     %Make colormap corresponding to
%            jet(param.color_levels)];
end

try thresh = param.thresh; catch thresh = 0.5; end

levels = [size(bMAP,1) size(aMAP,1)];
range_back = [1 size(bMAP,1)];          %Rescale image to the user
CAM = rescaling(im_back,range_back); %specified bmode graylevels.

Nt = size(AM,2);
CAM = repmat(CAM,[1 Nt]);

index.pos = abs(AM)>=max(abs(AM(:)))*thresh;
range_cMAP = size(bMAP,1)+[1 size(aMAP,1)];          %Rescale image to the user
CAM_temp = rescaling(AM,range_cMAP); %specified bmode graylevels.
CAM(index.pos) = CAM_temp(index.pos);

