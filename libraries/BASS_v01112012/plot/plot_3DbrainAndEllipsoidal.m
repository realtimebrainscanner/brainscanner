function [handles iVar DipColor] = plot_3DbrainAndEllipsoidal(vert,face,index,Var,opts)
%==========================================================================
% Filename: plot_3DbrainAndEllipsoidal.m (function).
% 
% Description:  
%
% Input:        vert: Vertices
%               face: Faces
%               index: Indices of vertices where ellipsoidal to be shown -
%                      dipoles.
%               opts:
%                   .flag_legend: Show legend if 1. Default off, i.e. = 0.
%                   .ellipScale: Scale ellipses. Default 10.
%                   .cMAP: Color map. Default is hot(256).
%                   .cBrain: Color of the brain. Default = [0.7 0.7 0.7].
%                   .DipColor: Color of ellipses
%                   .Ndip: #dipoles to be shown. Default is 5. Note that
%                          this field is only used if input 'index' is
%                          empty.
%                   .src_data: Show source activity. Default is zero.
%                   .src_thresh: Threshold for source data. Default 50% of
%                                maximum response in the source data.
%                   .flag_colorbar: false/true (off/on). Default off: false.
%
% Output:       handles: Handles to patch
%               iVar: Indices of dipoles (shown as ellipses) plotted.
%               DipColor: Color of the dipoles associated with iVar.
%
% History:
%   - Created:  01/06/2010
%   - Modified: 13/12/2010: Now possibly to both show ellipses and source
%                           data.
%
% Special remarks: 
%               When Var is approximately 50 the size seems good.
%               
%               Work to be done - use of colorbar does not corresponds to
%               the patch. The colorbar corresponds to the 
%
% Copyright (C) Carsten Stahlhut, DTU Informatics 2010
%==========================================================================
Nd = size(vert,1);
Ndip = length(index);

if Ndip==0
    try Ndip = opts.Ndip; catch Ndip = 5; end;
    [nVar iVar] = sort(Var,'descend');
    iVar = iVar(1:Ndip);
    Var = nVar(1:Ndip);
else
    iVar = index;
end

Ndip = sum(Var~=0);
iVar = iVar(1:Ndip);

try flag_legend = opts.flag_legend; catch flag_legend = false; end

try 
    cMAP = opts.cMAP; 
catch
    clevels = 256;
    cMAP = jet(clevels);
%     cMAP = hsv(256);
%     cMAP=cMAP(1:round(clevels-0.1*clevels),:);    %if cMAP based on hsv or hot
end

try cBrain = opts.cBrain; catch cBrain = [0.7 0.7 0.7]; end
try ellipScale = opts.ellipScale; catch ellipScale = 10; end
try flag_colorbar = opts.flag_colorbar; catch flag_colorbar = 0; end
try hfig = opts.hfig; catch hfig = figure; end

try 
    DipColor = opts.DipColor;
    dummy = DipColor(Ndip,:);     %Check that DipColor at least consists of Ndip-fields
catch
%     DipColor = jet(Ndip);
    DipColor = lines(Ndip);
%     col = ['b','g','r','c','m','y','k','w'];
%     tmp = ceil(Ndip./numel(col));
%     DipColor = repmat(col,1,tmp);
end

try 
    ellipSize = opts.ellipSize;
catch
    ellipSize = 'fixed';
    Var(:) = 6;
end

try data = opts.src_data; catch data = zeros(Nd,1); end
try thresh = opts.src_thresh; catch thresh = 0.5; end

[data_max, t_max] = max( max( abs(data),[],1 ) );
data = data(:,t_max);
% data_max = max( abs(data(:,t_max)) );
try crange = opts.crange; catch crange = [min(data(:)) max(data(:))]; end
ithresh = find(abs(data(:)) <= data_max*thresh);

Var = ellipScale*Var;
if size(Var,2)==1
    Var = repmat(Var,[1 3]);
end


%% Color look up
if diff(crange)~=0
    %Rescale data according to crange
    scale_out = length(cMAP);
    % range_out = [0 scale_out-1];
    range_out = [1 scale_out];      %Can ensure that activity should be zero to get the gray color
    cdata = rescaling(data,range_out,crange);
    cdata = cMAP(round(cdata),:);
else
    cdata = zeros(Nd,3);
end
cdata(ithresh,:) = repmat(cBrain,[length(ithresh) 1]);

handles.vert = double(vert);        %Ensures that vert is in double format
handles.face = double(face);
handles.cdata = cdata;


%% Show source activity on brain
handles.hfig  = figure(hfig);
set(gca,'Visible','off')
handles.fig1 = patch('vertices',handles.vert,...
    'faces',handles.face,'FaceVertexCData',handles.cdata);
colormap(cMAP)
caxis(crange)
set(handles.fig1,'FaceColor',cBrain,'EdgeColor','none');
shading interp
lighting gouraud
camlight
zoom off
lightangle(0,270);lightangle(270,0),lightangle(0,0),lightangle(90,0);
material([.1 .1 .4 .5 .4]);
view(-168,12)
axis image
hold on

% set(handles.fig1,'FaceColor',[.7 .7 .7],'EdgeColor','none');
% set(handles.fig1,'facecolor','w');
% set(handles.fig1,'facecolor','w','EdgeColor','b');


%% Show ellipses on brain. Ellipses indicate activity with most variance
for i =1:Ndip
    j = iVar(i);
    [x,y,z]= ellipsoid(vert(j,1),vert(j,2),vert(j,3),...
    1.*sqrt(Var(i,1)),1.*sqrt(Var(i,2)),1.*sqrt(Var(i,3)),20);
    he(i) = surfl(x, y, z);
%     colormap(cMAP)
    set(he(i),'edgecolor','none',...
        'facecolor',DipColor(i,:),...
        'facealpha',0.9)
    text_legend{i} = ['#' num2str(i)];
end

if flag_legend
    legend(he,text_legend)
end


if flag_colorbar
    Nlabels = 9;       %number of labels on the colorbar
%     ha2 = axes('position',[0.70 0.12 0.04 0.79]);
%     h2 = colorbar('Location','South');
    h2 = colorbar;
%     set(h2,'Position',get(ha2,'position'))
    set(h2,'Ylim',[0 1])
    set(h2,'YTick',linspace(0,1,Nlabels))
    if crange(2)<10        
        set(h2,'YTickLabel',num2str((round(100*...
            linspace(crange(1),crange(2),Nlabels))/100)'))
    else
        set(h2,'YTickLabel',num2str((round(...
            linspace(crange(1),crange(2),Nlabels)))'))
    end
%     set(get(h2,'XLabel'),'String','[ unit? ]' )             %Colorbar label
%     set(get(h2,'Title'),'String','   Hastighed')
%     set(ha2,'visible','off')
end

