function region = anatomy_brodmann()
%==========================================================================
% Filename: anatomy_lookup.m (function).
%
% Description:  Creates a list with names on different Brodmann areas
%               corresponding to the desired ones.
%
% Input:        area_no: Brodmann area number (Numbers should be between
%                        1-47)
%
% Output:       region: List with brodmann region names.
%
% History:
%   - Created:  19/07/2011
%   - Modified: 
%
% Author: Carsten Stahlhut
% Copyright (C) DTU Informatics 2011
%==========================================================================

%Nr = length(area_no);
%region = cell(Nr,1);

% for i=1:Nr
%     region{i} = sprintf('brodmann_area_%.0f',area_no(i));
% end

region = {
    'brodmann_area_1'
'brodmann_area_2'
'brodmann_area_3'
'brodmann_area_4'
'brodmann_area_5'
'brodmann_area_6'
'brodmann_area_7'
'brodmann_area_8'
'brodmann_area_9'
'brodmann_area_10'
'brodmann_area_11'
'brodmann_area_12'
'brodmann_area_13'
'brodmann_area_14'
'brodmann_area_15'
'brodmann_area_16'
'brodmann_area_17'
'brodmann_area_18'
'brodmann_area_19'
'brodmann_area_20'
'brodmann_area_21'
'brodmann_area_22'
'brodmann_area_23'
'brodmann_area_24'
'brodmann_area_25'
'brodmann_area_26'
'brodmann_area_27'
'brodmann_area_28'
'brodmann_area_29'
'brodmann_area_30'
'brodmann_area_31'
'brodmann_area_32'
'brodmann_area_33'
'brodmann_area_34'
'brodmann_area_35'
'brodmann_area_36'
'brodmann_area_37'
'brodmann_area_38'
'brodmann_area_39'
'brodmann_area_40'
'brodmann_area_41'
'brodmann_area_42'
'brodmann_area_43'
'brodmann_area_44'
'brodmann_area_45'
'brodmann_area_46'
'brodmann_area_47'
}
