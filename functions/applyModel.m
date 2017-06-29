function [predicted_stim, out]=applyModel(self,out,data)

%addpath C:\Users\sofha\Documents\realtime\git\brainscanner\functions
% sampling rate
fs=self.options.samplingRate;
% number of channels
run_source=out.run_source;
uf=out.uf;% upsample factor
num_channels=self.options.numChannels;
%%

% filter and rereference
% data_preprocessed =NaN(size(data));
% 
% 
% for idx_channel = 1:num_channels
%     %if d==1;
%     %   [data_pre,zf] = filter(filterB, filterA, [data{idx_trial}(idx_channel,(d-1)*wl+1:wl*d),zeros(1,500)]);
%     [data_preprocessed(idx_channel,:)] = filtfilt(out.options.filterB, out.options.filterA, [data(idx_channel,:)]);
%     %else
%     %  [data_preprocessed{idx_trial}(idx_channel,(d-1)*wl+1:wl*d),zf] = filter(filterB, filterA, data{idx_trial}(idx_channel,(d-1)*wl+1:wl*d),zf);
% end;
% 
% W=eye(num_channels)-1/num_channels;
% data_preprocessed = W*data_preprocessed;

data_preprocessed = preprocess(data, self.options.experiment);

if isfield(self.options.experiment,'asr_state')
    [data_preprocessed , out.asr_state] = asr_process(data_preprocessed, fs, self.options.experiment.asr_state);
end
%% Source localize
if run_source==1
    
    basis = load('model/basis_functions_24ch.mat','Qg');
    basisFunctions = basis.Qg;    % 569 basis functions uni-lateral
    opts.basisFunctions = basisFunctions;
    easyCapModelFull = load('model/Gain_mbraintrain_24ch.mat','Gain');
    easyCapModelFull=easyCapModelFull.Gain;
    easyCapModelFull(self.options.bad_chans,:)=[];
    no_chan=size(easyCapModelFull,1);
    easyCapModelFull=(eye(no_chan)-1/no_chan)*easyCapModelFull; % set gain to have average reference
    forwardModel = easyCapModelFull*basisFunctions';
    %
    
    
    gamma_median=out.gamma;%64.47;%-24.2;
    opts.run_prune=1;opts.prune=1e-2;opts.pnorm = 1;
    opts.m0=zeros(569,1);opts.m0(out.BFcorr(1:2),1)=1;
    [Sources,~,~,~] = teVGGD(forwardModel,data_preprocessed,gamma_median,opts);
    %AllS{idx_trial}(:,(d-1)*wl+1:wl*d)=basisFunctions'*sources;
end
%% Extract features
if run_source==1
    %[up,lu]=envelope(AllS{tr}(j(1:100),:).^2);
    %BFcorrExist=find(sum(Sources(out.BFcorr,:),2));
    y=((Sources(out.BFcorr(1:50),:)).^2);%
    %y{tr}=sum(AllS{tr}(occ,:).^2);%
    out.channel_list = 1:2';%:10;%10;%9:10';%;occ(361);%1;jh(1)';%find(sum(Sources{train_trial_idx},2))';%[1:45];
else
    y=data_preprocessed;
    out.channel_list=9:10';
end

% which channels to include

%channel_list = 1:5;%10;%9:10';%;occ(361);%1;jh(1)';%find(sum(Sources{train_trial_idx},2))';%[1:45];

% which frequencies to include?
feature_idx = out.feature_idx;%5:15;
feat = [];

for idx_channel = out.channel_list
    % compute multitaper spectrogram
    params = struct('Fs', fs, 'tapers', [2, 9]); % timebandwidth product = 2, num tapers = 9
    movingwin = [1/uf, 1/uf]; % window size and window step
    if sum(y(idx_channel, :))==0
        y(idx_channel, :)=randn(size(y(idx_channel, :)))*1e-16;
    end
    [up,lo]=envelope(y(idx_channel, :));
    [S,t, f] = mtspecgramc(up, movingwin, params );
    %S(S==-Inf)=[];%S(1);
    feat = [feat log10(S(:, feature_idx))];
    %feat(feat==-Inf)=feat(1);
end;

% extract
features = feat;
%features((features==-Inf))=NaN;
out.features=features;

%fprintf('\n\n')

%% standardize before fitting model
%features(isnan(features))=0;
features = bsxfun(@rdivide, bsxfun(@minus, features, out.m), out.s);

%% Train classifier

% predict
%fprintf('Predicting...')
predicted_stim = predict(out.model, features);
out.predicted_stim=predicted_stim;

if self.options.print_predicted_label
    if predicted_stim>0.5
        fprintf('Closed eyes, prob=%3.2f\n\n',predicted_stim)
    else
        fprintf('Open eyes, prob=%3.2f\n\n',predicted_stim)
    end
end;


