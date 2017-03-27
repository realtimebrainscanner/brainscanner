function [handles, dv, view_info] = plot_3DbrainInflated(vert,face,data,opts)
%==========================================================================
% Filename: plot_3Dbrain.m (function).
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
%               07/10/2012: Order of input 
%
% Copyright (C) Carsten Stahlhut, DTU Informatics 2012
%==========================================================================
dv = [];

if (nargin<3 || isempty(data))
    data = zeros(size(vert,1),1);
end

if nargin<4
    opts.dummy = [];
end

% try cMAP = opts.cMAP; catch cMAP = hsv(256); end;
try cMAP = opts.cMAP; catch cMAP = jet(256); end;
try cBrain = opts.cBrain; catch cBrain = [0.7 0.7 0.7]; end
try fs = opts.fs; catch fs = 50; end;
try thresh = opts.thresh; catch thresh = 0.5; end
try taxis = opts.taxis; catch taxis = []; end;
try flag_interp = opts.flag_interp; catch flag_interp = true; end
try flag_interactive = opts.flag_interactive; catch flag_interactive = false; end
try flag_colorbar = opts.flag_colorbar; catch flag_colorbar = true; end
try hfig = opts.hfig; catch hfig = figure; end
try flag_relval = opts.flag_relval; catch flag_relval = false; end
try data_unit = opts.data_unit; catch data_unit = []; end
try FaceAlpha = opts.FaceAlpha; catch FaceAlpha = 1; end

try iter_inflate = opts.iter_inflate ; catch iter_inflate  = inf; end

if isfield(opts,'crangeSym'), crangeSym = opts.crangeSym; else crangeSym = false; end

try flag_save_movie = opts.flag_save_movie; catch flag_save_movie = false; end
try flag_rot = opts.flag_rot; catch flag_rot = false; end
try fname_movie = opts.fname_movie; catch fname_movie = 'movie3DbrainInflated'; end
try view_coor = opts.view_coor; catch view_coor = [168 12]; end
try view_drot = opts.view_drot; catch view_drot = [0 0]; end

if flag_relval
    data_max = max(abs(data(:)));
    data = data/data_max;
end
try crange = opts.crange; catch crange = [min(data(:)) max(data(:))]; end

if crangeSym
    crange = [-max(abs(crange)) max(abs(crange))];
end

[Nd, Nt] = size(data);

%------------------------------
% Make figure
%------------------------------
figure
patch_structure = patch('Faces',face,'Vertices',vert,'FaceVertexCData',max(data,[],2),...
      'FaceColor','flat');
M = spm_mesh_inflate(patch_structure,iter_inflate,0);
C = spm_mesh_curvature_modified(patch_structure);

Nlevels_background = 20;
cMAP_back = gray(Nlevels_background);
% cMAP_back = cMAP_back([2 7],:);

    scale_out_back = size(cMAP_back,1);
    range_out_back = [1 scale_out_back];      %Can ensure that activity should be zero to get the gray color
    cdata_back = rescaling(C,range_out_back,[min(C(:)) max(C(:))]);
    cdata_back = cMAP_back(round(cdata_back(:)),:);


if ~isfield(opts,'cdata')


[data_max, t_max] = max( max( abs(data),[],1 ) );
% data = data(:,t_max);
ithresh = find(abs(data(:)) <= data_max*thresh);
clear t_max


% ithresh = abs(data) <= data_max*thresh;
cdata_back = repmat(cdata_back,[Nt 1]);



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
% cdata(ithresh,:) = repmat(cBrain,[length(ithresh) 1]);
cdata(ithresh,:) = cdata_back(ithresh,:);

cdata = reshape(cdata,[Nd Nt 3]);
cdata = permute(cdata,[1 3 2]);

else
    cdata = opts.cdata;     % Nd x 3 x Nt
end






temp = get(M);
handles.vert = temp.Vertices;
handles.face = temp.Faces;
handles.cdata = cdata;

close


handles.hfig = figure(hfig);
set(gcf,'Renderer','OpenGL')
set(gca,'Visible','off')
handles.patch = patch('vertices',handles.vert,...
    'faces',handles.face,'FaceVertexCData',handles.cdata(:,:,1));
colormap(cMAP)
if diff(crange)~=0, caxis(crange), else caxis([0 1]), end
set(handles.patch,'FaceColor',cBrain,'EdgeColor','none','FaceAlpha',FaceAlpha);
if flag_interp, shading interp, else shading flat, end
lighting gouraud
camlight
zoom off
lightangle(0,270);lightangle(270,0),lightangle(0,0),lightangle(90,0);
material([.1 .1 .3 .5 .4]);
view(view_coor(1),view_coor(2))

axis image
hold on
if flag_colorbar
    hcbar = colorbar;
    if ~isempty(data_unit)
        set(get(hcbar,'XLabel'),'String',data_unit)             %Colorbar label
    end

end
set(gcf,'Toolbar','figure')
% set(handles.fig1,'facecolor','w');



%---------------------------------------------------------------------
% Multiple time points show movie or toogle through the data using
% interactive plotting
%---------------------------------------------------------------------
if flag_interactive
    if ~isempty(taxis)
        minorStep = 1/(length(taxis)-1);
        majorStep = minorStep*5;
        handles.taxis = taxis;

        handles.slideText = uicontrol('Style','text',...
                'Position',[440 20 75 20],...
                'String',[num2str(round(taxis(1)*1e3)) ' ms']);

        handles.slidebar = uicontrol('Style', 'slider',...
            'Min',1,'Max',length(taxis),'Value',1,...
            'Position', [125 20 300 20],...
            'SliderStep',[minorStep majorStep],...
            'Callback', {@sliderUpdate,handles});        
    else
        error('Chosing interactive plot requires the field taxis to be specified in opts.')
    end
else
    
    if flag_save_movie
        nFrames = size(data,2);
        % Preallocate movie structure.
        mov(1:nFrames) = struct('cdata', [],...
            'colormap', []);
    end
    
    if ~isempty(taxis)
        % Create textbox
        ha1 = annotation(handles.hfig,'textbox',...
            [0.74 0.017 0.21 0.067],...
            'String',{sprintf('Time: %3.3f sec',taxis(1) ) },...
            'FitBoxToText','off');
        
        for ti=1:size(data,2)
            set(ha1,'String',{sprintf('Time: %3.3f sec',taxis(ti) ) } )
            set(handles.patch,'FaceVertexCdata',handles.cdata(:,:,ti));
            drawnow
            pause(1/(fs))
            
            
            if flag_save_movie
                mov(ti) = getframe;
            end
            
            if flag_rot
                view(view_coor(1),view_coor(2))
                view_coor = view_coor + view_drot;
                if abs(view_coor(2))>=90
                    view_drot(2) = -view_drot(2);
                    view_coor(2) = view_coor(2) + 2*view_drot(2);
                end
            end
                        
        end
        
        
        if flag_save_movie
            save(fname_movie,'mov')
        end
        

    else
        if size(data,2)>1
            ha1 = annotation(handles.hfig,'textbox',...
                [0.74 0.017 0.21 0.067],...
                'String',{sprintf('Time: %3.3f sec',1 ) },...
                'FitBoxToText','off');

            for ti=1:size(data,2)
                set(ha1,'String',{sprintf('Sample: %3.0f',ti ) } )
                set(handles.patch,'FaceVertexCdata',handles.cdata(:,:,ti));
                drawnow
                pause(1/(fs))
                
                
                if flag_save_movie
                    mov(ti) = getframe(ha1);
                end 
            end
            
            disp('Now playing movie by movie-fnc')
            if flag_save_movie
                save(fname_movie,'mov')
            end

            
        end
    end
       
end


view_info.coor = view_coor;
view_info.drot = view_drot;

end


function sliderUpdate(hObj,event,info)
%==========================================================================
% Filename: sliderUpdate.m (function).
% 
% Description:  Called to set FaceVertexCdata of patch in figure axes
%               when user moves the slider control
%
% Input:        hObj: Handle to the slider object.
%               event: Faces
%               info: handles with info
%
% Output:       None
%
% Copyright (C) Carsten Stahlhut, DTU Informatics 2011
%==========================================================================
    ti = round(get(hObj,'Value'));
    set(info.patch,'FaceVertexCdata',info.cdata(:,:,ti));
    set(info.slideText,'String',[num2str(round(info.taxis(ti)*1e3)) ' ms'])
end


