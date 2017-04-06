% Brain Scanner Application
function varargout = BrainScanner(varargin)


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_OpeningFcn, ...
                   'gui_OutputFcn',  @main_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before main is made visible.
function main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main (see VARARGIN)

% Choose default command line output for main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%% Load libraries
if ~exist('MARDv2','file')
    addpath(genpath('functions')); end

if ~exist('rescaling', 'file');
    addpath('libraries/BASS_v01112012'); end

if ~exist('lsl_loadlib','file')
    addpath(genpath('libraries/liblsl-Matlab')); end
try 
    lib = lsl_loadlib(env_translatepath('libraries/liblsl-Matlab:/bin'));
catch
    lib = lsl_loadlib();
end

%% Setup stream and settings 

opts.eegStreamName = 'openvibeSignal'; %'EEG';
opts.bufferrange = 10; % (in seconds)
opts.channelrange = 1:24;

opts.samplingRate = 250;
opts.blockSize = 64;
opts.refreshrate = 1/(2*opts.blockSize/opts.samplingRate); % Maybe change this!

% Pre processing
[opts.filterB, opts.filterA] = butter(4,[1 45]/(opts.samplingRate/2));
opts.artifactRemoval = 0;
opts.filter = 0;
opts.zeromean = 0; 
opts.reref = 0;
opts.standardize = 0; 


% Model setup
basis = load('model/basisFunctions.mat');
basisFunctions = basis.IDX2;
opts.basisFunctions = basisFunctions;

easyCapModelFull = importdata('model/easyCapModel.mat');
opts.forwardModel = easyCapModelFull(:, basisFunctions);
opts.numSources = numel(basisFunctions);

vertface = load('model/vertface');
opts.faces = vertface.face2;
opts.verts = vertface.vert2;
opts.QG = basis.QG;

% basisFunctions = sort(randperm(size(easyCapModelFull,2), numSources));
% forwardModel = easyCapModelFull(:, basisFunctions);
% opts.numSources = 400;

% Other
opts.recoveryMethod = 'MARD';
opts.saveData = 0;
opts.log = 0;

handles.eeg = EEGStream(lib, opts, handles);


initialize_gui(hObject, handles, false);

% UIWAIT makes main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%% Stimuli experiment settings
function runExperiment(handles)

handles.eeg.options.saveData = 1;
handles.eeg.options.log = 1;




% --- Outputs from this function are returned to the command line.
function varargout = main_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.


% Update handles structure
guidata(handles.figure1, handles);



% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if hObject.Value
    if handles.eeg.isConnected
        handles.text13.String = 'Live';
    else
        handles.text13.String = 'Replay';
    end
    
    handles.pushbutton5.Enable = 'off';
    handles.togglebutton6.Enable = 'off';
    handles.eeg.start();
else
    handles.text13.String = 'Off';
    handles.eeg.stop();
    
    handles.togglebutton10.Value = 0;
    handles.togglebutton11.Value = 0;
    handles.togglebutton14.Value = 0;
    
    handles.pushbutton5.Enable = 'on';
    handles.togglebutton6.Enable = 'on';
end


% --- Executes on button press in togglebutton2.
function togglebutton2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if hObject.Value
    handles.text14.String = 'On';
else
    handles.text14.String = 'Off';
end
handles.eeg.options.filter = hObject.Value;


% --- Executes on button press in togglebutton3.
function togglebutton3_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if hObject.Value
    handles.text15.String = 'On';
else
    handles.text15.String = 'Off';
end
handles.eeg.options.reref = hObject.Value;

% --- Executes on button press in togglebutton4.
function togglebutton4_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if hObject.Value
    handles.text16.String = 'On';
else
    handles.text16.String = 'Off';
end
handles.eeg.options.zeromean = hObject.Value;


% --- Executes on button press in togglebutton6.
function togglebutton6_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if hObject.Value
    if handles.eeg.connect()
        handles.text17.String = 'Connected';
        handles.togglebutton1.Enable = 'on';
        handles.text23.String = 'No file';
    else 
        hObject.Value = 0;
    end
else
    if handles.eeg.isConnected
        handles.eeg.disconnect(); end;
    handles.togglebutton1.Enable = 'off';
    handles.text17.String = 'Off';
end


% --- Executes on button press in togglebutton14.
function togglebutton14_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~hObject.Value
    handles.eeg.closeFigure(handles.eeg.TimeFigure); end;

handles.eeg.showTiming = hObject.Value;


% --- Executes on button press in togglebutton10.
function togglebutton10_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~hObject.Value
    handles.eeg.closeFigure(handles.eeg.DataFigure); end;

handles.eeg.showData = hObject.Value;


% --- Executes on button press in togglebutton11.
function togglebutton11_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~hObject.Value
    handles.eeg.closeFigure(handles.eeg.BrainFigure); end;

handles.eeg.showBrain = hObject.Value;


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.eeg.stop();
delete(hObject);


% --- Executes on button press in pushbutton5. 
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fileName, path] = uigetfile('*.csv','Select the saved EEG data file (.csv)');
if fileName
    handles.text23.String = 'File loaded';
    handles.togglebutton1.Enable = 'on';
    handles.eeg.replayFileName = strcat(path,fileName);
else
    disp('No file selected');
end


% --- Executes on button press in togglebutton16.
function togglebutton16_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton16 (see GCBO)
if hObject.Value
    handles.text24.String = 'On';
else
    handles.text24.String = 'Off';
end
handles.eeg.options.saveData = hObject.Value;
handles.eeg.options.log = hObject.Value;


% --- Executes on button press in togglebutton17.
function togglebutton17_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton17 (see GCBO)

if hObject.Value
    handles.text25.String = 'On';
else
    handles.text25.String = 'Off';
end
% runExperiment(handles);
handles.eeg.options.saveData = hObject.Value;
handles.eeg.options.log = hObject.Value;
handles.eeg.showExperiment = hObject.Value;
