function simple_gui_test
% SIMPLE_GUI2 Select a data set from the pop-up menu, then
% click one of the plot-type push buttons. Clicking the button
% plots the selected data in the axes.

   %  Create and then hide the UI as it is being constructed.
   f = figure('Visible','off','Position',[360,500,450,285]);
   
   %  Construct the components.
   hsurf = uicontrol('Style','pushbutton','String','Surf',...
           'Position',[315,220,70,25]);
   hmesh = uicontrol('Style','pushbutton','String','Mesh',...
           'Position',[315,180,70,25]);
   hcontour = uicontrol('Style','pushbutton',...
           'String','Contour',...
           'Position',[315,135,70,25]); 
   htext = uicontrol('Style','text','String','Select Data',...
           'Position',[325,90,60,15]);
   hpopup = uicontrol('Style','popupmenu',...
           'String',{'Peaks','Membrane','Sinc'},...
           'Position',[300,50,100,25]);
   ha = axes('Units','Pixels','Position',[50,60,200,185]); 
   align([hsurf,hmesh,hcontour,htext,hpopup],'Center','None');
   
   % Make the UI visible.
   f.Visible = 'on';

end