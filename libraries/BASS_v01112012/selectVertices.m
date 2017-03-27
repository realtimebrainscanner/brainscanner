function [viArray handles] = selectVertices(vert,face,data,opts)
%
% Example:
%   [viArray handles] = selectVertices( vert,face,zeros(size(vert,1),1) );
%
try dummy = opts.thresh; catch; opts.thresh = 0; end  % threshold in percent of maximum signal
try dummy = opts.patchSize; catch; opts.patchSize = 3; end% number of incremental surounding spheres for ROI
try dummy = opts.fs; catch; opts.fs = 100; end % speed of movie : 1 = 1 second per sample; default : 50
try dummy = opts.cBrain; catch; opts.cBrain = [0.7 0.7 0.7]; end% only for full brain (not modifiable for inflated brain)
try dummy = opts.cMAP; catch; opts.cMAP = jet(256); end
try dummy = opts.flag_relval; catch; opts.flag_relval = false; end
try dummy = opts.flag_interactive; catch; opts.flag_interactive = false; end % true for slidebar; false for movie
try dummy = opts.hfig; catch; opts.hfig = figure; end
try dummy = opts.flag_colorbar; catch; opts.flag_colorbar = false; end
try dummy = opts.taxis; catch; opts.taxis = 1; end

handles = plot_3Dbrain(vert,face,data,opts);
path_work = pwd;
qinput1 = 's';
qinput2 = 's';
iselect = 0;
viArray = [];
WOIdata = zeros(size(data));
while ~strcmp(qinput1,'e')
    cd(fullfile(spm('Dir'),'external\fieldtrip\plotting\private'))
    qinput2 = 's';
    while ~strcmp(qinput2,'q')
        %             select3dtool(handles.hfig)
        select3dtool(handles.hfig)
        waitforbuttonpress
        qinput2 = input('Push key q for quit when desired position is found:','s');
        qinput1 = input('Push key e to exit the selection tool or any other keys to continue:','s');
    end
    iselect = iselect+1;
    [p v vi face_select facei_select] = select3d;
    viArray(iselect) = vi;
    WOIdata(vi) = 1;
    
    %         close('Select 3-D Tool')
    
    cd(path_work)
    close
    opts.hfig=2;
    handles = plot_3Dbrain(vert,face,WOIdata,opts);
    %         set(handles.patch,'FaceVertexCdata',data);
    drawnow
end