function handles = setup3DBrain(vert,face,data,opts)
%==========================================================================
% Filename: setup_3Dbrain.m (function).
% 
% Description:  
%
% Input:        vert: Vertices
%               face: Faces
%               data: data (activity) that will be shown on brain
%               opts:
%                   .cMAP: colormap
%                   .taxis: Time axis in [s].
%                   .flag_interp: Do shading interpolation - default true.
%                   .flag_interactive: Allows interactive visualization
%                                      of the data using a sliderbar. Note
%                                      opts.taxis should also be specified.
%                                      Default: false
%                   .flag_colorbar: Show colorbar - default true.
%                   .flag_relval: Rescale to relative values with
%                                 max(abs(data)) as 1 - default false.
%                   .FaceAlpha: Transparency - default 1 no transparency
%
% Output:       handles: Handles to patch
%
% Example:      handles = plot_3Dbrain(vert,face,data)
%               handles = plot_3Dbrain(vert,face,data,opts)
%
% History:
%   - Created:  28/11/2008
%   - Modified: 27/12/2010: Now with slider to toogle through the source
%                           acitivity
%               22/01/2011: It is now possible to use the subplot function
%                           by using subplot outside this function.
%
% Copyright (C) Carsten Stahlhut, DTU Informatics 2011
%==========================================================================

if (nargin<3 || isempty(data))
    data = zeros(size(vert,1),1);
end

if nargin<4
    opts.dummy = [];
end

% try cMAP = opts.cMAP; catch cMAP = hsv(256); end;
try cMAP = opts.cMAP; catch; cMAP = jet(256); end;

try crange = opts.crange; catch; crange = [min(data(:)) max(data(:))]; end
try cBrain = opts.cBrain; catch; cBrain = [0.5 0.5 0.5]; end
% try fs = opts.fs; catch; fs = 50; end;
try thresh = opts.thresh; catch; thresh = 0.5; end
% try taxis = opts.taxis; catch; taxis = []; end;
try flag_interp = opts.flag_interp; catch; flag_interp = true; end
% try flag_interactive = opts.flag_interactive; catch; flag_interactive = false; end

try flag_colorbar = opts.flag_colorbar; catch; flag_colorbar = true; end
try hfig = opts.hfig; catch; hfig = figure; end
try ax = opts.axes; catch; ax = axes('Parent',hfig,'Position',[.13 .15 .78 .75]); end;
try flag_relval = opts.flag_relval; catch; flag_relval = false; end
try data_unit = opts.data_unit; catch; data_unit = []; end
try FaceAlpha = opts.FaceAlpha; catch; FaceAlpha = 1; end

if isfield(opts,'crangeSym'), crangeSym = opts.crangeSym; else crangeSym = false; end
if crangeSym
    crange = [-max(abs(crange)) max(abs(crange))];
end

[Nd, Nt] = size(data);


if ~isfield(opts,'cdata')
    
if flag_relval
    data_max = max(abs(data(:)));
    data = data/data_max;
end

[data_max, t_max] = max( max( abs(data),[],1 ) );
% data = data(:,t_max);
ithresh = find(abs(data(:)) <= data_max*thresh);
clear t_max

%% Color look up
if diff(crange)~=0
    %Rescale data according to crange
    scale_out = length(cMAP);
    % range_out = [0 scale_out-1];
    range_out = [1 scale_out];      %Can ensure that activity should be zero to get the gray color
    cdata = rescaling(data,range_out,crange);
    cdata = cMAP(round(cdata(:)),:);
%     keyboard
else
    cdata = zeros(Nd,1,3);
end
%cBrain
cdata(ithresh,:) = repmat(cBrain,[length(ithresh) 1]);
cdata = reshape(cdata,[Nd Nt 3]);
cdata = permute(cdata,[1 3 2]);

else
    cdata = opts.cdata;     % Nd x 3 x Nt
end

handles.vert = vert;
handles.face = face;
handles.cdata = cdata;

handles.hfig = figure(hfig);
handles.axes = ax;

set(gcf,'Renderer','OpenGL')
set(gca,'Visible','off')


handles.patch = patch('vertices',handles.vert,'faces',handles.face,'FaceVertexCData',handles.cdata(:,:,1));

colormap(cMAP)

if diff(crange)~=0, caxis(crange), else caxis([0 1]), end

% set(handles.patch,'FaceColor',cBrain,'EdgeColor','none','FaceAlpha',FaceAlpha);
handles.patch.FaceColor = cBrain;
handles.patch.EdgeColor = 'none';
handles.patch.FaceAlpha = FaceAlpha;


if flag_interp, shading interp, else shading flat, end

lighting gouraud

light = findobj(handles.hfig,'Type','Light');
if isempty(light)
    camlight;
    lightangle(0,270);lightangle(90,0);lightangle(0,0);%lightangle(90,0);
else
    camlight(light(1));
    lightangle(light(2), 0,270); lightangle(light(3), 90,0); lightangle(light(4), 0,0);%lightangle(90,0);
end;

%  camlight
%  zoom off
%  lightangle(0,270);lightangle(90,0),lightangle(0,0),%lightangle(90,0);
%material([.01 .2 .4 .7 .8]);%
%material([.1 .1 .4 .5 .4]);

%view(-168,12)
view(0,0)
axis image
% hold on
% if flag_colorbar
%     hcbar = colorbar;
%     if ~isempty(data_unit)
%         set(get(hcbar,'XLabel'),'String',data_unit)             %Colorbar label
%     end
% end
% set(gcf,'Toolbar','figure')

end