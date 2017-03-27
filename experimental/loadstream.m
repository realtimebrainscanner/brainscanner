% Load
% opts = arg({'streamname','StreamName'},'EEG',{'EEG','Other'},'LSL stream that should be displayed. The name of the stream that you would like to display.');


% if ~exist('arg_define','file')
%     addpath(genpath(fileparts(mfilename('fullpath')))); end

if ~exist('lsl_loadlib','file')
    addpath(genpath('liblsl-Matlab')); end


try 
    lib = lsl_loadlib(env_translatepath('liblsl-Matlab:/bin'));
catch
    lib = lsl_loadlib();
end



eegStreamName = 'openvibeSignal'; %'EEG';

opts.streamname = eegStreamName;
opts.arg_direct = 1;
opts.bufferrange = 10;
opts.timerange = 1;
opts.datascale = 150;
opts.channelrange = [1:24];
opts.samplingrate = 250;
opts.blockSize = 32;
opts.refreshrate = 1/(2*opts.blockSize/opts.samplingrate); %1;  %50; % 0.18; % Probably change this!
opts.freqfilter = 0;
opts.reref = 0;
opts.standardize = 1; % 0 or 1?
opts.zeromean = 1;
% dummy stuff
opts.parent_fig = [];
opts.parent_ax = [];
opts.pageoffset = 0;
opts.position = [];


opts.lpFilter = designfilt('lowpassiir','FilterOrder',6,'PassbandFrequency',45,'PassbandRipple',0.2, 'SampleRate',opts.samplingrate);
opts.hpFilter = designfilt('highpassiir','FilterOrder',6,'PassbandFrequency',1,'PassbandRipple',0.2, 'SampleRate',opts.samplingrate);
opts.bpFilt = designfilt('bandpassfir','FilterOrder',20, 'CutoffFrequency1',0.2,'CutoffFrequency2',45, 'SampleRate',opts.samplingrate);

% create stream inlet, figure and stream buffer
inlet = EEGStream.CreateInlet(lib,opts);
stream = EEGStream.CreateStreambuffer(opts,inlet.info());    
% [fig,ax,lines] = create_figure(opts,@on_key,@on_close);
% optionally design a frequency filter
if length(opts.freqfilter) == 4
    B = design_bandpass(opts.freqfilter,stream.srate,20,true);
elseif isscalar(opts.freqfilter)
    B = ones(opts.freqfilter,1)/max(1,opts.freqfilter);
else
    error('The FIR filter must be given as 4 frequencies in Hz [raise-start,raise-stop,fall-start,fall-stop] or moving-average length in samples.');
end
disp(inlet.info());

% start a timer that reads from LSL and updates the display
th = timer('TimerFcn',{@EEGStream.OnTimer, opts, inlet, stream, lib},'Period',1.0/opts.refreshrate,'ExecutionMode','fixedRate');
start(th);



