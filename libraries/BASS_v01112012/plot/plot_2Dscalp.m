function hfig = plot_2Dscalp(data,Cpos,Cnames,opts)
%==========================================================================
% Filename: plot_2Dscalp.m (function).
% 
% Description:  
%
% Input:        data: Scalp data ( channels x 1 )
%               Cpos: 2D coordinates of the channels ( 2 x channels)
%               Cnames: Names of the channels (cell-array)
%               opts:
%                   .hfig: Figure handle if a previous figure should be
%                          used.
%                   .flag_Cnames: Show channels names on plot.
%                   .flag_contour: Show contour of a scalp with a nose.
%
% Output:       hfig: Handles to figure
%
% Copyright (C) Carsten Stahlhut, DTU Informatics 2010
%==========================================================================
if nargin<4, opts.hfig = figure; end
if isfield(opts,'hfig'), hfig = opts.hfig; else hfig = figure; end
if isfield(opts,'Cpos_scale'), Cpos_scale = opts.Cpos_scale; else Cpos_scale = 1; end
if isfield(opts,'flag_Cnames'), flag_Cnames = opts.flag_Cnames; else flag_Cnames = true; end
if isfield(opts,'flag_contour'), flag_contour = opts.flag_contour; else flag_contour = true; end
if isfield(opts,'resolution'), res = opts.resolution; else res = 0.001; end
if isfield(opts,'flag_colorbar'), flag_colorbar = opts.flag_colorbar; else flag_colorbar = true; end
if isfield(opts,'flag_electrodes'), flag_electrodes = opts.flag_electrodes; else flag_electrodes = true; end
if isfield(opts,'flag_restrict_circle'), flag_restrict_circle = opts.flag_restrict_circle; else flag_restrict_circle = false; end

figure(hfig);
hold on

Cpos = (Cpos - 0.5)*Cpos_scale;
Cpos = Cpos + 0.5;

Cheigh = repmat(max(data),[1 length(data)]);

x = min(Cpos(1,:)):res:max(Cpos(1,:));
y = min(Cpos(2,:)):res:max(Cpos(2,:));

[x1,y1] = meshgrid(x,y);
xp = Cpos(1,:)';
yp = Cpos(2,:)';

z = griddata(xp, yp, data, x1, y1);

if flag_restrict_circle
    [THETA,RHO] = cart2pol(x1(:)-0.5,y1(:)-0.5);
    z(RHO>0.5) = NaN;
end

%% Check for outside scalp

%Headfigure:
if flag_contour
    t = 0:2*pi/1000:2*pi;
    x_head = cos(t)/2+0.5;
    y_head = sin(t)/2+0.5;
    plot3(x_head,y_head,repmat(max(data),[1 length(x_head)]),'k','LineWidth',1)
    x_nose = [-0.2 0 0.2];
    y_nose = sin(acos(x_nose));
    y_nose(2)= y_nose(2)+0.1;
    plot3(x_nose/2+0.5,y_nose/2+0.5,repmat(max(data),[1 length(x_nose)]),'k','LineWidth',1)
end
% 
% figure
surface(x,y,z);
shading('interp')
hold on
if flag_electrodes
    plot3(Cpos(1,:), Cpos(2,:), Cheigh, '.k');
end
if flag_Cnames
    text(Cpos(1,:),Cpos(2,:),Cheigh,Cnames,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','Bottom','fontsize',10)
end

axis equal
% title('2D scalp map')
if flag_colorbar
    h3 = colorbar;
    set(get(h3,'Xlabel'),'String','[\muV]')
    set(gca,'Xticklabel',[],...
        'Yticklabel',[]);
    axis([-0.005 1.005 0 1.05])
end

set(gca,'Visible','off')