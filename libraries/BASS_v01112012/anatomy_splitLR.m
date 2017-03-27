function [regionLR_vert regionLR_list] = anatomy_splitLR(MNIvert,region_vert,region_list)
%==========================================================================
% Filename: anatomy_splitLR.m (function).
%
% Description:  Function to split region into left and right hemisphere
%               regions.
%
% Input:        MNIvert: Vertices coordinates [Nv x 3] in MNI space.
%               region_vert: Logical array keeping the information of which
%                            vertices being in a given region. Rows:
%                            region number, columns: vertex number.
%               region_list: list with regions specifying which vertices
%                            are included in the given regions
%
% Output:       regionLR_vert: Logical array keeping the information of which
%                            vertices being in a given region. Rows:
%                            region number, columns: vertex number.
%               regionLR_list: list with regions specifying which vertices
%                            are included in the given regions
%
% History:
%   - Created:  27/09/2012
%   - Modified:
%
% Author: Carsten Stahlhut
% Copyright (C) DTU Informatics 2012
%==========================================================================



x = MNIvert(:,1);
iR = x>0;
iL = ~iR;

figure
hold on
hL = plot3(MNIvert(iL,1),MNIvert(iL,2),MNIvert(iL,3),'.b');
hR = plot3(MNIvert(iR,1),MNIvert(iR,2),MNIvert(iR,3),'.r');
grid on
axis equal


Nregions = length(region_list);


regionL_vert = bsxfun(@times,region_vert,iL');  %Left region
regionR_vert = bsxfun(@times,region_vert,iR');  %Right region
regionLR_vert = [regionL_vert; regionR_vert];


% Left
for ir=1:Nregions
    regionLR_list(ir).name = [region_list(ir).name '_L'];
    regionLR_list(ir).ivert = find(regionL_vert(ir,:));
end

% Right
for ir=1:Nregions
    regionLR_list(ir+Nregions).name = [region_list(ir).name '_R'];
    regionLR_list(ir+Nregions).ivert = find(regionR_vert(ir,:));
end
