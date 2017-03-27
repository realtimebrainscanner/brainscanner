function [handles,hpatch] = plot_multipleLayers(data,hLayer,opts)
%==========================================================================
% Filename: show_3Dbrain.m (function).
% 
% Description:  Fast sampling from multivariate normal distribution with
%
% Input:        data
%               hLayer: Cell-array with the different layers to be
%               illustrated. Each cell consists of a layer in gii (gifti)
%               format.
%               opts:
%
% Output:       S: Variable sampled from multivariate normal distribution,
%                  size (Nd x Ns).
%               cov_S:
%
% Example:      hLayer = {mesh.tess_scalp, mesh.tess_oskull,...
%                         mesh.tess_iskull, mesh.tess_ctx};
%
% History:
%   - Created:  28/11/2008
%   - Modified: 
%
% Special remarks:
%
% Copyright (C) Carsten Stahlhut, DTU Informatics 2008
%==========================================================================

% Original file: plot_eeg_brain - however changed the in
try title_text = opts.title_text; catch title_text = ''; end
try fs = opts.fs; catch fs = 1; end
try dataType = opts.dataType; catch dataType = 'real'; end
try crange = opts.crange; catch crange = []; end
try data_unit = opts.data_unit; catch data_unit = []; end
try cMAP = opts.cMAP; catch cMAP = hot(256); end;
try alpha_arr = opts.alpha; catch alpha_arr = 0.2*ones(1,length(hLayer)); alpha_arr(end)=1; end;
try FaceColor_arr = opts.FaceColor_arr; catch FaceColor_arr = []; end;

switch dataType
    case 'abs'
        data = full(abs(data));
        if isempty(data_unit)
            data_unit = 'Am';
%             data_unit = 'nAm';
        end

    case 'absNorm'                  %Normalised
        data = full(abs(data));
        data = data/max(data(:));

    case 'real'
        data = full(data);
        if isempty(data_unit)
            data_unit = 'Am';
%             data_unit = 'nAm';
        end

    case 'realNorm'                 %Normalised
        data = full(data);
        data = data/max(abs(data(:)));
        
    case 'energy'
        data = sum(full(data.^2),2);
        if isempty(data_unit)
            data_unit = '(Am)^2';
%             data_unit = 'nAm';
        end

    case 'energyNorm'               %Normalised
        data = sum(full(data.^2),2);
        data = data/max(data(:));
        
    case 'Zscore'
        try sig2 = opts.data_variance; catch error('No variance specified'); end
        %Ensure data and sig2 have the same size
        if size(data,2)~=size(sig2,2)
            sig2 = repmat(sig2,[1 size(data,2)]);
        end
        if size(data,1)~=size(sig2,1)
            sig2 = repmat(sig2,[size(data,1) 1]);
        end
        data = data./sqrt(sig2);
        if isempty(data_unit)
            data_unit = 'Zscore';
        end
        
    case {'PPV'}
        param.bmode_levels = 1;
        param.flow_levels = 210;
        param.cmax = crange(2);                 %Maximum
        param.cmin = crange(1);                 %Minimum
        param.Nlabels = 5;              %Number of labels on colorbar.
        
        cMAP = hot(256);
%         %Option 1 - write is assigned the value -1
        cMAP = cMAP(1:param.flow_levels,:);
        cMAP = [1 1 1; cMAP];

%         %Option 2 - black is assigned the value -1
%         cMAP = cMAP(end-param.flow_levels+1:end,:);
%         cMAP = cMAP(size(cMAP,1):-1:1,:);       %Reverse order
%         cMAP = [0 0 0; cMAP];
                
        range_in = [param.cmin param.cmax];
        range_out = [param.bmode_levels param.flow_levels];
        idata = find(data >= 0);
        data(idata) = rescaling(data(idata),range_out,range_in);
        data(data<0) = 0;
        crange = [0 param.flow_levels];         %New rescaling

    case {'maxDist'}
        param.bmode_levels = 1;
        param.flow_levels = 210;
        param.cmax = crange(2);                 %Maximum
        param.cmin = crange(1);                 %Minimum
        param.Nlabels = 5;              %Number of labels on colorbar.
        
        cMAP = hot(256);
        %Option 1 - write is assigned the value -1
        cMAP = cMAP(1:param.flow_levels,:);
        cMAP = cMAP(size(cMAP,1):-1:1,:);       %Reverse order
        cMAP = [1 1 1; cMAP];

%         %Option 2 - black is assigned the value -1
%         cMAP = cMAP(end-param.flow_levels+1:end,:);
%         cMAP = cMAP(size(cMAP,1):-1:1,:);       %Reverse order
%         cMAP = [0 0 0; cMAP];
                
        range_in = [param.cmin param.cmax];
        range_out = [param.bmode_levels param.flow_levels];
        idata = find(data >= 0);
        data(idata) = rescaling(data(idata),range_out,range_in);
        data(data<0) = 0;
        crange = [0 param.flow_levels];         %New rescaling
        
    case {'SNR'}
        cMAP = hot(256);
%         cMAP = cMAP(10:240,:);

    case {'cost'}
        cMAP = hot(128);
%         cMAP = cMAP(1:210,:);
       
    otherwise
        error('This option is not supported')
end

if isempty(crange)
    crange = [min(data(:)) max(data(:))];
end
if isempty(data_unit)
    data_unit = '';
end

% crange = [min(srcs_disp(:)) max(srcs_disp(:))];
crange;
NLayer = length(hLayer);

% handles.vert = vert;
% handles.face = face;
handles.data = data;
handles.dataType = dataType;
handles.data_unit = data_unit;

%------------------------------
% Make figure
%------------------------------
figure
h = axes('Visible','off');
hpatch = NaN*ones(1,NLayer);

for iL=1:NLayer
    htemp = gifti(hLayer(iL));
    handles.vert = htemp.vertices;
    handles.face = htemp.faces;
    clear htemp
    
%     if iL<=NLayer
    if iL<NLayer
%         handles.fig1 = patch('vertices',handles.vert,...
%             'faces',handles.face,'FaceVertexCData',(NLayer+1-iL)/NLayer*zeros(size(handles.vert,1),1));
        handles.fig1 = patch('vertices',handles.vert,...
            'faces',handles.face,'FaceVertexCData',-(NLayer+1-iL)/NLayer*ones(size(handles.vert,1),1));
%         set(handles.fig1,'FaceAlpha',0.3)
        set(handles.fig1,'FaceAlpha',alpha_arr(iL))
        set(handles.fig1,'FaceColor',[.5 .5 .5],'EdgeColor','none')
        cMap_layer = copper(20);
        colormap(cMap_layer)
    else
        handles.fig1 = patch('vertices',handles.vert,...
            'faces',handles.face,'FaceVertexCData',data(:,1));
        set(handles.fig1,'FaceAlpha',1)
        set(handles.fig1,'FaceColor',[.5 .5 .5],'EdgeColor','none');
        colormap(cMAP)
    end
    hpatch(iL) = handles.fig1;

    
%     colormap(cMAP)
    caxis(crange);
%     set(handles.fig1,'FaceColor',[.5 .5 .5],'EdgeColor','none');
    shading interp
    hold off
    
%     keyboard
end

    colormap(cMAP)
    try caxis(crange);, end
    set(handles.fig1,'FaceColor',[.5 .5 .5],'EdgeColor','none');
    shading interp

lighting gouraud
camlight
zoom off
lightangle(0,270);lightangle(270,0),lightangle(0,0),lightangle(90,0);
material([.1 .1 .4 .5 .4]);
% view(-45,50)
% view(-149,30)
% view(-162,22)
view(170,8)
view(-168,12)
axis image

% handles.colorbar = colorbar;
% dummy = get(handles.colorbar,'OuterPosition');
% dummy(2) = 0.1;
% % dummy(2) = 0.15;
% % bar_ratio = 0.65/dummy(4);
% dummy(4) = 0.8;
% % dummy(4) = 0.75;
% set(handles.colorbar,'OuterPosition',dummy)
% set(h,'Position',[-0.01 0.05 0.95 0.95])
% % text(55,40,50,title_text,'HorizontalAlignment',...
% %     'center','VerticalAlignment',...
% %     'middle','fontsize',14)
% 
% if ~strcmp(dataType,'cost')
%     dummy = get(handles.colorbar,'Position');
%     dummy(3) = 0.05;
%     set(handles.colorbar,'Position',dummy);
% else
%     dummy = get(handles.colorbar,'Position');
%     dummy(3) = 0.07;
%     set(handles.colorbar,'Position',dummy);
% end
% 
% if (strcmp(dataType,'PPV') | strcmp(dataType,'maxDist'))
%     h2 = handles.colorbar;
% %     set(h2,'Position',get(h,'position'))
%     set(h2,'Ylim',[param.bmode_levels param.flow_levels])
%     set(h2,'YTick',linspace(param.bmode_levels,...
%         param.flow_levels,param.Nlabels))
%     set(h2,'YTickLabel',num2str((round(100*...
%         linspace(param.cmin,param.cmax,param.Nlabels))/100)'))
%     set(get(h2,'XLabel'),'String',data_unit)             %Colorbar label
% %     set(get(h2,'Title'),'String','   Hastighed')
%     handles.colorbar = h2;
% else    
%     set(get(handles.colorbar,'XLabel'),'String',data_unit)             %Colorbar label
% end



% % title(title_text)
% 
% % dummy = get(handles.colorbar,'OuterPosition');
% % dummy(2) = 0.1;
% % % dummy(2) = 0.15;
% % % bar_ratio = 0.65/dummy(4);
% % dummy(4) = 0.8;
% % % dummy(4) = 0.75;
% % set(handles.colorbar,'OuterPosition',dummy)
% % 
% % set(h,'Position',[-0.01 0.05 0.95 0.95])
% % 
% % text(55,40,50,title_text,'HorizontalAlignment',...
% %     'center','VerticalAlignment',...
% %     'middle','fontsize',14)
% 

% for ti=1:size(data,2)
%     set(handles.fig1,'FaceVertexCdata',data(:,ti));
%     drawnow
%     pause(1/(fs*4))
% end
% 
% if nargin>2
%     show_3Dbrain(data,handles.vert,handles.face,opts)
% else
%     show_3Dbrain(data,handles.vert,handles.face)
% end