function out=trainModel(options)
out.options=options;
%addpath C:\Users\sofha\Documents\realtime\git\brainscanner\functions
%addpath(genpath('C:\Users\sofha\Documents\realtime\chronux_2_11'))
% sampling rate
fs=options.samplingRate;
% number of channels
num_channels = 24;
run_source=input('Run source localization (0/1): ');
out.run_source=run_source;
uf=2;% upsample factor
out.uf=uf;
run_calcGamma=1;
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
    subject_data =  subject_data((EEGdelay*fs)+1:end, 1:num_channels);
    
    % how many samples do have labels for?
    if isfield(options,'asr_state')
    try
    num_labelled_samples = length(stim_data.stim)*fs+ceil(1/uf*250);
    catch
    num_labelled_samples = length(stim_data.stim)*fs;
    end
    else
    num_labelled_samples = length(stim_data.stim)*fs;
    end
    data{idx_trial} = subject_data(1:num_labelled_samples, 1:num_channels)';
    %data{idx_trial} = subject_data(:, 1:24)';
    
    %fprintf('done.\n\n');
    
    
end;
%% Preprocess

% filter and rereference
data_preprocessed = cell(num_trials, 1);

wl=1/uf*fs;
for idx_trial = 1:num_trials
    for idx_channel = 1:num_channels
        for d=1:floor(length(data{idx_trial})/wl);
            %if d==1;
            %   [data_pre,zf] = filter(filterB, filterA, [data{idx_trial}(idx_channel,(d-1)*wl+1:wl*d),zeros(1,500)]);
            [data_preprocessed{idx_trial}(idx_channel,(d-1)*wl+1:wl*d)] = filtfilt(options.filterB,options.filterA, [data{idx_trial}(idx_channel,(d-1)*wl+1:wl*d)]);
            %else
            %  [data_preprocessed{idx_trial}(idx_channel,(d-1)*wl+1:wl*d),zf] = filter(filterB, filterA, data{idx_trial}(idx_channel,(d-1)*wl+1:wl*d),zf);
            %end
        end
    end
end;
for idx_trial = 1:num_trials
    W=eye(num_channels)-1/num_channels;
    data_preprocessed{idx_trial} = W*data_preprocessed{idx_trial};
end;
if isfield(options,'asr_state')
    [data_preprocessed{idx_trial} , out.asr_state] = asr_process(data_preprocessed{idx_trial} , fs, options.asr_state);
    data_preprocessed{idx_trial}=data_preprocessed{idx_trial}(:,floor(0.25*fs):end);
    no_win=floor(length(data_preprocessed{idx_trial})/wl);
    stims{idx_trial}=stims{idx_trial}(1:no_win);
else
    no_win=floor(length(data_preprocessed{idx_trial})/wl);
end
    
   % stims{idx_trial}=stims{idx_trial}(1:no_win);
%%
train_trial_idx = 1;
%% Source localize
if run_source==1
    opts.bad_chans=[];
    basis = load('model/basis_functions_24ch.mat','Qg');
    basisFunctions = basis.Qg;    % 569 basis functions uni-lateral
    opts.basisFunctions = basisFunctions;
    easyCapModelFull = load('model/Gain_mbraintrain_24ch.mat','Gain');
    easyCapModelFull=easyCapModelFull.Gain;
    easyCapModelFull(opts.bad_chans,:)=[];
    no_chan=size(easyCapModelFull,1);
    easyCapModelFull=(eye(no_chan)-1/no_chan)*easyCapModelFull; % set gain to have average reference
    forwardModel = easyCapModelFull*basisFunctions';
    %
    if run_calcGamma==1
        wl=1/uf*fs;
        disp('Training source localization parameter')
        opts.run_prune=1;opts.prune=1e-2;opts.pnorm = 1;opts.min_gamma=-100;
        for d=1:no_win;%floor(length(data_preprocessed{idx_trial})/wl)
            if ~rem(d,50)
                fprintf('window %d out of %d windows\n',d,floor(length(data_preprocessed{idx_trial})/wl))
            end
            [gamma_mean(d), gamma_median(d)]= teVGGD_wcross(forwardModel,data_preprocessed{train_trial_idx}(:,(d-1)*wl+1:wl*d), opts);
        end
    else gamma_median=-66;
    end
end
%%
if run_source==1
    gamma_median=median(gamma_median);%64.47;%-24.2;
    out.gamma=gamma_median;
    fprintf('Gamma is set to %d\n',gamma_median)
    opts.run_prune=1;opts.prune=1e-2;opts.pnorm = 1;
    wl=1/uf*fs;
    Sources = cell(num_trials, 1);AllS = cell(num_trials, 1);
    opts.m0=zeros(569,1);
    for idx_trial = 1:num_trials
        fprintf('Calculating sources\n');
        for d=1:no_win;%floor(length(data_preprocessed{idx_trial})/wl);
            [sources,msources,~,~] = teVGGD(forwardModel,data_preprocessed{idx_trial}(:,(d-1)*wl+1:wl*d),gamma_median,opts);
            %AllS{idx_trial}(:,(d-1)*wl+1:wl*d)=basisFunctions'*sources;
            Sources{idx_trial}(:,(d-1)*wl+1:wl*d)=sources;
        end
    end
end

%%
if run_source==1
    stim1=repmat(stims{train_trial_idx},[1,125])';
    stim1=stim1(:);
    % figure,plot(stim1);
    % n=16;%6 9;
    % act=find(sum(Sources{train_trial_idx},2));
    % hold on,plot(Sources{1}(act(n),:)/max(abs(Sources{1}(act(n),:))))
    clear co
    for c=1:569%
        [up,lo]=envelope(Sources{train_trial_idx}(c,:));%AllS{1}%data_preprocessed
        %[up,lo]=envelope(AllS{1}(c,:));%%data_preprocessed
        
        tmp=corrcoef(up,stim1');
        co(c)=tmp(2);
    end
    figure,plot(abs(co))
    co(isnan(co))=0;
    [i,j]=sort(abs(co),'descend');
    figure,
    tr=1;
    plot(sum(Sources{tr}(j(1),:).^2,1)),hold on,
    stim=repmat(stims{tr},[1,125])';
    stim=stim(:);
    hold on,plot(stim);
    legend('Source with maximum correlation to stimulation','Stimulation')
    out.BFcorr=j;
end
%%
if run_source==1
    gamma_median=median(gamma_median);%64.47;%-24.2;
    out.gamma=gamma_median;
    fprintf('Gamma is set to %d\n',gamma_median)
    opts.run_prune=1;opts.prune=1e-2;opts.pnorm = 1;
    wl=1/uf*fs;
    Sources = cell(num_trials, 1);AllS = cell(num_trials, 1);
    opts.m0=zeros(569,1);opts.m0(j(1:2),1)=1;
    for idx_trial = 1:num_trials
        fprintf('Calculating sources again\n');
        for d=1:no_win;%floor(length(data_preprocessed{idx_trial})/wl);
            [sources,msources,~,~] = teVGGD(forwardModel,data_preprocessed{idx_trial}(:,(d-1)*wl+1:wl*d),gamma_median,opts);
            %AllS{idx_trial}(:,(d-1)*wl+1:wl*d)=basisFunctions'*sources;
            Sources{idx_trial}(:,(d-1)*wl+1:wl*d)=sources;
        end
    end
end

%% Extract features
if run_source==1
    for tr=1:num_trials
        %[up,lu]=envelope(AllS{tr}(j(1:100),:).^2);
        y{tr}=((Sources{tr}(j(1:50),:)).^2);%
        %y{tr}=sum(AllS{tr}(occ,:).^2);%
        out.channel_list = 1:2;%:10;%10;%9:10';%;occ(361);%1;jh(1)';%find(sum(Sources{train_trial_idx},2))';%[1:45];
        
    end
else
    y=data_preprocessed;
    out.channel_list=9:10';
end
features = cell(num_trials, 1);

% which channels to include

%channel_list = 1:5;%10;%9:10';%;occ(361);%1;jh(1)';%find(sum(Sources{train_trial_idx},2))';%[1:45];

% which frequencies to include?
out.feature_idx = 5:15;

figure;

for idx_trial = 1:num_trials
    
    fprintf('Extracting features for trial %d...', idx_trial);
    
    feat = [];
    
    for idx_channel = out.channel_list
        idx_channel
        % compute multitaper spectrogram
        params = struct('Fs', fs, 'tapers', [2, 9]); % timebandwidth product = 2, num tapers = 9
        movingwin = [1/uf, 1/uf]; % window size and window step
        for d=1:no_win;%floor(length(data_preprocessed{idx_trial})/wl);
            if sum(y{idx_trial }(idx_channel, (d-1)*wl+1:wl*d))==0
                disp('zero y')
                y{idx_trial}(idx_channel, (d-1)*wl+1:wl*d)=randn(size(y{idx_trial }(idx_channel, (d-1)*wl+1:wl*d)))*1e-16;
            end
            [up((d-1)*wl+1:wl*d),lo]=envelope(y{idx_trial}(idx_channel, (d-1)*wl+1:wl*d));
        end
        [S,t, f] = mtspecgramc(up, movingwin, params );
        %S(S==-Inf)=[];%S(1);
        feat = [feat log10(S(:, out.feature_idx))];
        %feat(feat==-Inf)=feat(1);
    end;
    
    % extract
    features{idx_trial} = feat;
    
    fprintf('done.\n');
    
    subplot(2, 2, idx_trial);
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
train_feat(trans,:)=[];
train_stim(trans)=[];

% train
%model = fitcdiscr(features{train_trial_idx}, stims{train_trial_idx}, 'DiscrimType','Linear');
%model = fitcdiscr(features{train_trial_idx}, stims{train_trial_idx}, 'DiscrimType','Quadratic')
model = fitglm(features{train_trial_idx}, stims{train_trial_idx},'linear','Distribution','binomial', 'Intercept', true);
%model = fitglm(train_feat, train_stim,'linear','Distribution','binomial', 'Intercept', true);
out.model=model;

fprintf('done.\n')

% predict
fprintf('Predicting...')
predicted_stim = cell(num_trials, 1);
for idx_trial = 1:num_trials
    predicted_stim{idx_trial} = predict(model, features{idx_trial});
end;
fprintf('done.\n\n')

% map from probabilities to labels
binarize = @(x) 1.0*(x > 0.5);

% compute metrics
classification_acc = zeros(num_trials, 1);
for idx_trial = 1:num_trials
    classification_acc(idx_trial) = mean(stims{idx_trial}(:) == binarize(predicted_stim{idx_trial}(:)));
    if idx_trial == train_trial_idx; traintest = 'training data'; else;  traintest = 'test data'; end;
    
    % output
    fprintf('Trial %d: Acc = %4.3f (%s)\n', idx_trial, classification_acc(idx_trial), traintest)
end;
out.trainingError=classification_acc(idx_trial);
%% Plot

figure;
for idx_trial = 1:num_trials
    
    if idx_trial == train_trial_idx; traintest = 'training data'; else;  traintest = 'test data'; end;
    
    hold on;
    
    plot(stims{idx_trial}, 'k-', 'linewidth', 3)
    plot(predicted_stim{idx_trial}, 'r--')
    legend('Stimuli', 'Prediction')
    ylim([-0.1, 1.1]);
    title(sprintf('Trial %d', idx_trial));
    grid on;
    
end;