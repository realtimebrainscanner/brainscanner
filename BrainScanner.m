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
%if ~exist('MARDv2','file')
%    addpath(genpath('functions')); end

addpath(genpath('functions'))

if ~exist('rescaling', 'file');
    addpath('libraries/BASS_v01112012'); end

if ~exist('lsl_loadlib','file')
    addpath(genpath('libraries/liblsl-Matlab')); end

if ~exist('asr_process','file')
    addpath(genpath('libraries/BCILAB')); end

if ~exist('mtspecgramc','file')
    addpath(genpath('libraries/chronux_2_11')); end

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
opts.blockSize = 32;
opts.refreshrate = 1/(2*opts.blockSize/opts.samplingRate); % Maybe change this!

% Pre processing for viewer
[opts.filterB, opts.filterA] = butter(4,[1 45]/(opts.samplingRate/2));
opts.artefactRemoval = 0;
opts.filter = 0;
opts.reref = 0;
opts.bad_chans=[];%

% Options for experiment
[opts_exp.filterB, opts_exp.filterA] = butter(4,[5 15]/(opts.samplingRate/2));
opts_exp.artefactRemoval = 0;
opts_exp.filter = 0;
opts_exp.reref = 0;
opts_exp.bad_chans=opts.bad_chans;% e.g. [5 21];
opts.experiment = opts_exp;


% artifact removal
opts.artefactRemoval = 0;

% plotting
opts.numSamplesToPlot = 1000;
opts.rangeChannelPlot = 100;

% Model setup
basis = load('model/basis_functions_24ch.mat','Qg');
basisFunctions = basis.Qg;    % 569 basis functions uni-lateral
opts.basisFunctions = basisFunctions;

easyCapModelFull = load('model/Gain_mbraintrain_24ch.mat','Gain');

easyCapModelFull=easyCapModelFull.Gain;
easyCapModelFull(opts.bad_chans,:)=[];
no_chan=size(easyCapModelFull,1);
easyCapModelFull=(eye(no_chan)-1/no_chan)*easyCapModelFull; % set gain to have average reference
opts.forwardModel = easyCapModelFull*basisFunctions';
opts.numSources = size(basisFunctions,1);
channels=load('model/Gain_mbraintrain_24ch.mat','Channel');
opts.channames={channels.Channel(:).Name};clear channels

vertface = load('model/vertface_24ch_ICBMtemp');
opts.faces = vertface.face;
opts.verts = vertface.vert;
opts.numChannels = 24-numel(opts.bad_chans);

% prediction
opts.print_predicted_label = 1;

%opts.QG = basis.QG;

% % Old Model setup
% basis = load('model/basisFunctions.mat');
% basisFunctions = basis.IDX2;    % 776 basis functions uni-lateral
% opts.basisFunctions = basisFunctions;
%
% easyCapModelFull = importdata('model/easyCapModel.mat');
% opts.forwardModel = easyCapModelFull(:, basisFunctions);
% opts.numSources = numel(basisFunctions);
%
% vertface = load('model/vertface');
% opts.faces = vertface.face2;
% opts.verts = vertface.vert2;
% opts.QG = basis.QG;



% basisFunctions = sort(randperm(size(easyCapModelFull,2), numSources));
% forwardModel = easyCapModelFull(:, basisFunctions);
% opts.numSources = 400;

% Other
opts.recoveryMethod = 'teVG';%'MARD';
if strcmp(opts.recoveryMethod,'teVG');
    try 
        load('gamma.mat', 'gmedian');
        opts.gamma = nanmedian(gmedian);
        fprintf('Loaded gamma value = %3.2f from file gamma.mat\n', opts.gamma);
    catch
        opts.gamma=-75;
        fprintf('File gamma.mat not found. Using default value for gamma = %3.2f.\n', opts.gamma);
    end;
end
opts.saveData = 0;
opts.log = 0;
opts.verbose = 1;
opts.trainVG = 0;
opts.numSamplesCalibrationVGData=4*opts.blockSize;

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
        handles.togglebutton18.Enable = 'off';

    else
        hObject.Value = 0;
        handles.togglebutton18.Enable = 'on';

    end
else
    if handles.eeg.isConnected
        handles.eeg.disconnect(); end;
    handles.togglebutton1.Enable = 'off';
    handles.text17.String = 'Off';
            handles.togglebutton18.Enable = 'on';

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
    handles.eeg.replayDataFile=[];
else
    disp('No file selected');
end


% --- Executes on button press in togglebutton16.
function togglebutton16_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton16 (see GCBO)
if hObject.Value
    handles.text24.String = 'On';
    handles.eeg.setup();
else
    handles.text24.String = 'Off';
end
handles.eeg.options.saveData = hObject.Value;
handles.eeg.options.log = hObject.Value;


% --- Executes on button press in togglebutton18.
function togglebutton18_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.eeg.options.trainVG = hObject.Value;

% Hint: get(hObject,'Value') returns toggle state of togglebutton18

% --- Executes on button press in togglebutton19.
function togglebutton19_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if hObject.Value
    handles.text27.String = 'On';
    %[options.filterB, options.filterA] = butter(4,[5 15]/(handles.eeg.options.samplingRate/2));
    %options.samplingRate=handles.eeg.options.samplingRate;
    
  
    handles.eeg.ClassificationModel=trainModel(handles.eeg.options);
else
    handles.text27.String = 'Off';
end
handles.eeg.options.trainModel = hObject.Value;



% --- Executes on button press in apply_asr_button.
function apply_asr_button_Callback(hObject, eventdata, handles)
% hObject    handle to apply_asr_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of apply_asr_button

if hObject.Value
    handles.text_asr.String = 'On';
else
    handles.text_asr.String = 'Off';
end

handles.eeg.options.artefactRemoval = hObject.Value;

% --- Executes on button press in load_asr_button.
function load_asr_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_asr_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of load_asr_button
fprintf('Choose calibration data file:\n');
[filename_data, PathName] = uigetfile('*.csv');
fprintf('Loading data file %s\n',filename_data)
calib_data = csvread([PathName,filename_data], 1, 0);%% columns 0
calib_data = calib_data(:, 1:24);
handles.eeg.asr_state = asr_calibrate(preprocess(calib_data', handles.eeg.options), handles.eeg.options.samplingRate);    
handles.text_viewer_asr_loaded.String = 'Loaded';

% --- Executes on button press in button_filter_exp.
function button_filter_exp_Callback(hObject, eventdata, handles)
% hObject    handle to button_filter_exp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button_filter_exp
if hObject.Value
    handles.text30.String = 'On';
else
    handles.text30.String = 'Off';
end

handles.eeg.options.experiment.filter = hObject.Value;


% --- Executes on button press in button_reref_exp.
function button_reref_exp_Callback(hObject, eventdata, handles)
% hObject    handle to button_reref_exp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button_reref_exp
if hObject.Value
    handles.text31.String = 'On';
else
    handles.text31.String = 'Off';
end
handles.eeg.options.experiment.reref = hObject.Value;

% --- Executes on button press in togglebutton25.
function togglebutton25_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton25


% --- Executes on button press in button_asr_exp.
function button_asr_exp_Callback(hObject, eventdata, handles)
% hObject    handle to button_asr_exp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button_asr_exp

if hObject.Value
    handles.text33.String = 'On';
else
    handles.text33.String = 'Off';
end

handles.eeg.options.experiment.artefactRemoval = hObject.Value;

% --- Executes on button press in button_load_asr_exp.
function button_load_asr_exp_Callback(hObject, eventdata, handles)
% hObject    handle to button_load_asr_exp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button_load_asr_exp
fprintf('Choose calibration data file:\n');
[filename_data, PathName] = uigetfile('*.csv');
fprintf('Loading data file %s\n',filename_data)
calib_data = csvread([PathName,filename_data], 1, 0);%% columns 0
calib_data = calib_data(:, 1:24);
handles.eeg.options.experiment.asr_state = asr_calibrate(preprocess(calib_data', handles.eeg.options), handles.eeg.options.samplingRate);    

handles.text35.String = 'Loaded';
