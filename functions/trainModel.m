function out=trainModel(options)
out.options=options;
fs=options.samplingRate;
num_channels = options.numChannels;% number of channels
run_source=input('Run source localization (0/1): ');
out.run_source=run_source;
uf=2;% upsample factor
out.uf=uf;
run_calcGamma=1; % Set to 0 and specify a gamma value if you have it predefined.
gamma_predefined=-65; % Set sparsity parameter
eeg_channel_list=9:10'; % Used if run_source=0;
out.feature_idx = 5:15; %Hz, frequencies of interest
%% Load data

num_trials = 1;

% load data & stim
data = cell(num_trials, 1);
stims = cell(num_trials, 1);

for idx_trial = 1:num_trials
    fprintf('Choose stim file\n');
    [filename_stim,PathName] = uigetfile('*.mat');
    fprintf('Loading stim file %s\n',filename_stim)
    
    stim_data = load([PathName,filename_stim], 'stim');
    stim_data.stim=repmat(stim_data.stim,[uf,1]);
    stims{idx_trial} = stim_data.stim(:);
    %fprintf('done.\n');
    
    fprintf('Choose data file\n');
    [filename_data,PathName] = uigetfile('*.csv');
    fprintf('Loading data file %s\n',filename_data)
    subject_data = csvread([PathName,filename_data], 1, 0);%% columns 0
    
    % remove first 5s.
    prompt='Type number of seconds to delay EEG: ';
    EEGdelay = input(prompt);
    subject_data =  subject_data((EEGdelay*fs)+1:end, 1:24);
    
        % how many samples do have labels for?
    if isfield(options.experiment,'asr_state')
        try
            num_labelled_samples = length(stim_data.stim)*fs+ceil(1/uf*250);
        catch
            num_labelled_samples = length(stim_data.stim)*fs;
        end
        else
        num_labelled_samples = length(stim_data.stim)*fs;
    end
    data{idx_trial} = subject_data(1:num_labelled_samples,:)';
    
    
end;

%% Apply preprocessing block-wise
wl=1/uf*fs;
data_preprocessed = zeros(24-length(options.experiment.bad_chans),size(data{1},2));
for d=1:floor(length(data{1})/wl)
    [data_preprocessed(:,(d-1)*wl+1:wl*d)] = preprocess(data{idx_trial}(:,(d-1)*wl+1:wl*d), options.experiment);
end

%% Artifact removal using ASR
if options.experiment.artefactRemoval
    if isfield(options.experiment, 'asr_state')
        [data_preprocessed, out.asr_state] = asr_process(data_preprocessed , fs, options.experiment.asr_state);
    else
        fprintf('trainModel: Load calibration data before applying ASR!\n');
    end
end


if isfield(options.experiment,'asr_state')
    [data_preprocessed , out.asr_state] = asr_process(data_preprocessed , fs, options.experiment.asr_state);
    data_preprocessed=data_preprocessed(:,floor(0.25*fs):end);
    no_win=floor(length(data_preprocessed)/wl);
    stims{idx_trial}=stims{idx_trial}(1:no_win);
else
    no_win=floor(length(data_preprocessed)/wl);
end
%% Keep structure
data_preprocessed = {data_preprocessed};

%%
train_trial_idx = 1;
%% Source localize
if run_source==1
    basis = load('model/basis_functions_24ch.mat','Qg');
    basisFunctions = basis.Qg;    % 569 basis functions uni-lateral
    opts.basisFunctions = basisFunctions;
    easyCapModelFull = load('model/Gain_mbraintrain_24ch.mat','Gain');
    easyCapModelFull=easyCapModelFull.Gain;
    easyCapModelFull(options.bad_chans,:)=[];
    no_chan=size(easyCapModelFull,1);
    easyCapModelFull=(eye(no_chan)-1/no_chan)*easyCapModelFull; % set gain to have average reference
    forwardModel = easyCapModelFull*basisFunctions';
    %
    if run_calcGamma==1 % Estimate the sparsity parameter, gamma
        wl=1/uf*fs;
        disp('Training source localization parameter')
        opts.run_prune=1;opts.prune=1e-2;opts.pnorm = 1;% Use active set
        opts.min_gamma=-100;
        for d=1:floor(length(data_preprocessed{idx_trial})/wl)
            if ~rem(d,50)
                fprintf('window %d out of %d windows\n',d,floor(length(data_preprocessed{idx_trial})/wl))
            end
            [gamma_mean(d), gamma_median(d)]= teVGGD_wcross(forwardModel,data_preprocessed{train_trial_idx}(:,(d-1)*wl+1:wl*d), opts);
        end
    else
        gamma_median=gamma_predefined; % Use predefined gamma
    end
end
%%
if run_source==1
    gamma_median=median(gamma_median);%64.47;%-24.2;
    out.gamma=gamma_median;
    fprintf('Gamma is set to %d\n',gamma_median)
    opts.run_prune=1;opts.prune=1e-2;opts.pnorm = 1; % Use active set
    wl=1/uf*fs;
    Sources = cell(num_trials, 1);
    opts.m0=zeros(569,1);
    for idx_trial = 1:num_trials
        fprintf('Calculating sources\n');
        for d=1:floor(length(data_preprocessed{idx_trial})/wl);
            [sources,~,~,~] = teVGGD(forwardModel,data_preprocessed{idx_trial}(:,(d-1)*wl+1:wl*d),gamma_median,opts);
            Sources{idx_trial}(:,(d-1)*wl+1:wl*d)=sources;
        end
    end
end

%% If source localization is chosen then find the sources (basisfunctions) that are most correlated with the stimulation

if run_source==1
    stim1=repmat(stims{train_trial_idx},[1,fs/uf])'; % repeat stimuli to match the EEG windows
    stim1=stim1(:);
    clear co
    for c=1:size(basisFunctions,1)%
        [up,lo]=envelope(Sources{train_trial_idx}(c,:));%AllS{1}%data_preprocessed
        tmp=corrcoef(up,stim1');
        co(c)=tmp(2);
    end
    figure,plot(abs(co))
    co(isnan(co))=0;
    [i,j]=sort(abs(co),'descend');
    legend('Absolute value of correlation between sources (basis functions) and stimuli')
    xlabel('Basis function index');
    figure,
    tr=1;
    plot(sum(Sources{tr}(j(1),:).^2,1)),hold on,
    stim=repmat(stims{tr},[1,125])';
    stim=stim(:);
    hold on,plot(stim);
    legend('Source with maximum correlation to stimulation','Stimulation')
    xlabel('Time samples');
    out.BFcorr=j;
end
%% Run source localization again where the support is initialiazed to be active in the two most correlated sources.
if run_source==1
    gamma_median=median(gamma_median);%64.47;%-24.2;
    out.gamma=gamma_median;
    fprintf('Gamma is set to %d\n',gamma_median)
    opts.run_prune=1;opts.prune=1e-2;opts.pnorm = 1;% Use active set
    wl=1/uf*fs;
    Sources = cell(num_trials, 1);
    opts.m0=zeros(569,1);opts.m0(j(1:2),1)=1;
    for idx_trial = 1:num_trials
        fprintf('Calculating sources again\n');
        for d=1:floor(length(data_preprocessed{idx_trial})/wl);
            [sources,~,~,~] = teVGGD(forwardModel,data_preprocessed{idx_trial}(:,(d-1)*wl+1:wl*d),gamma_median,opts);
            Sources{idx_trial}(:,(d-1)*wl+1:wl*d)=sources;
        end
    end
end

%% Extract features
if run_source==1
    for tr=1:num_trials
        y{tr}=((Sources{tr}(j(1:50),:)).^2);%
        out.channel_list = 1:2;
        
    end
else
    y=data_preprocessed;
    out.channel_list=eeg_channel_list; % Change
end
features = cell(num_trials, 1);
figure;
for idx_trial = 1:num_trials
    fprintf('Extracting features for trial %d...', idx_trial);
    feat = [];   
    for idx_channel = out.channel_list
        % compute multitaper spectrogram
        params = struct('Fs', fs, 'tapers', [2, 9]); % timebandwidth product = 2, num tapers = 9
        movingwin = [1/uf, 1/uf]; % window size and window step
        for d=1:no_win
            if sum(y{idx_trial }(idx_channel, (d-1)*wl+1:wl*d))==0
                disp('zero y')
                y{idx_trial}(idx_channel, (d-1)*wl+1:wl*d)=randn(size(y{idx_trial }(idx_channel, (d-1)*wl+1:wl*d)))*1e-16;
            end
            [up((d-1)*wl+1:wl*d),lo]=envelope(y{idx_trial}(idx_channel, (d-1)*wl+1:wl*d));
        end
        [S,t, f] = mtspecgramc(up, movingwin, params );
        feat = [feat log10(S(:, out.feature_idx))];
    end
    
    % extract
    features{idx_trial} = feat;
    
    fprintf('done.\n');
    
    imagesc(t, f, log10(S)'); axis xy;
    title(sprintf('Trial %d', idx_trial));
    ylim([1, 15]);
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    % output
    
    
end;
out.features=features;
fprintf('\n\n')

%% standardize before fitting model
m = nanmean(features{train_trial_idx}, 1);
s = nanstd(features{train_trial_idx}, 1);
for idx_trial = 1:num_trials
    features{idx_trial} = bsxfun(@rdivide, bsxfun(@minus, features{idx_trial}, m), s);
end;
out.m=m;
out.s=s;
%% Train classifier
fprintf('Training classifier...')

train_feat=features{train_trial_idx};
train_stim=stims{train_trial_idx};
excl=1*uf; % number of trials to exclude on both sides of condition transition;
trans=find(diff(train_stim));
trans=repmat(trans,[1,excl*2])+repmat(-excl:excl-1,[length(trans),1]);
%train_feat(trans,:)=[];
%train_stim(trans)=[];

% train
model = fitglm(features{train_trial_idx}, stims{train_trial_idx},'linear','Distribution','binomial', 'Intercept', true);
%model = fitglm(train_feat, train_stim,'linear','Distribution','binomial',
%'Intercept', true); % Use to discard data around task switches
out.model=model;

fprintf('done.\n')

% predict
fprintf('Predicting...')
predicted_stim = cell(num_trials, 1);
for idx_trial = 1:num_trials
    predicted_stim{idx_trial} = predict(model, features{idx_trial});
end
fprintf('done.\n\n')

% map from probabilities to labels
binarize = @(x) 1.0*(x > 0.5);

% compute metrics
classification_acc = zeros(num_trials, 1);
for idx_trial = 1:num_trials
    classification_acc(idx_trial) = mean(stims{idx_trial}(:) == binarize(predicted_stim{idx_trial}(:)));
    if idx_trial == train_trial_idx; traintest = 'training data'; else;  traintest = 'test data'; end   
    % output
    fprintf('Trial %d: Acc = %4.3f (%s)\n', idx_trial, classification_acc(idx_trial), traintest)
end
out.trainingError=classification_acc(idx_trial);
%% Plot

figure;
for idx_trial = 1:num_trials    
    hold on;  
    plot(stims{idx_trial}, 'k-', 'linewidth', 3)
    plot(predicted_stim{idx_trial}, 'r--')
    legend('Stimuli', 'Prediction')
    ylim([-0.1, 1.1]);
    title(sprintf('Trial %d', idx_trial));
    grid on; 
    xlabel('Time');
end