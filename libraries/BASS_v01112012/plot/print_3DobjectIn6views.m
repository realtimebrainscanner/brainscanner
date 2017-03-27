function fig_names = print_3DobjectIn6views(fname,fig_path_doc,opts)
%==========================================================================
% Filename: print_3DobjectIn6vies.m (function).
% 
% Description:  
%
% Input:        
%
% Output:       fig_names: cell-array with the printed figures
%
% Example:      handles = plot_3Dbrain(vert,face,data)
%               handles = plot_3Dbrain(vert,face,data,opts)
%
% History:
%   - Created:  31/03/2012
%
% Copyright (C) Carsten Stahlhut, DTU Informatics 2012
%==========================================================================
if nargin<3
    opts.hfig = gcf;
    
    if nargin<2
        fig_path_doc = pwd;
        
        if nargin<1
            fname = 'brain';
        end
    end
end
if isfield(opts,'hfig'), figHandle = opts.hfig; else figHandle = gcf; end
if isfield(opts,'figFormat'), figFormat = opts.figFormat; else figFormat = '.eps'; end
if isfield(opts,'figDFormat'), figDFormat = opts.figDFormat; else figDFormat = '-depsc2'; end

figure(figHandle)

fig_names = cell(6,1);
view(-90,0)
fig_name = fullfile(fig_path_doc,[fname '_view1']);
saveas(gcf,[fig_name '.fig'])
print(gcf,figDFormat,[fig_name figFormat])
fig_names{1} = fig_name;

view(-90,90)
fig_name = fullfile(fig_path_doc,[fname '_view2']);
saveas(gcf,[fig_name '.fig'])
print(gcf,figDFormat,[fig_name figFormat])
fig_names{2} = fig_name;

%     view(1,-90)
view(-90,-90)
fig_name = fullfile(fig_path_doc,[fname '_view3']);
saveas(gcf,[fig_name '.fig'])
print(gcf,figDFormat,[fig_name figFormat])
fig_names{3} = fig_name;

view(90,0)
fig_name = fullfile(fig_path_doc,[fname '_view4']);
saveas(gcf,[fig_name '.fig'])
print(gcf,figDFormat,[fig_name figFormat])
fig_names{4} = fig_name;


view(90,90)
fig_name = fullfile(fig_path_doc,[fname '_view5']);
saveas(gcf,[fig_name '.fig'])
print(gcf,figDFormat,[fig_name figFormat])
fig_names{5} = fig_name;

view(90,-90)
fig_name = fullfile(fig_path_doc,[fname '_view6']);
saveas(gcf,[fig_name '.fig'])
print(gcf,figDFormat,[fig_name figFormat])
fig_names{6} = fig_name;