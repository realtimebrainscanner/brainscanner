% EEG data visualization / inspection
clear;
if ~exist('rescaling', 'file');
    addpath('libraries/BASS_v01112012'); end

% Load data
path = 'EEGData/';
fileName = 'raw_22-03-2017 12-59-20.csv'; % First half 30 Hz sine and 40 Hz sine last half
tableData = readtable([path fileName]);

colNames = tableData.Properties.VariableNames;
disp(colNames);



%% Timing
% plot(tableData{:, 'TimeStamp'});

timeStamps = tableData{:, 'TimeStamp'};
timeDiff = [];
i = 2;

for i=2:numel(timeStamps)
    timeDiff(i-1) = (timeStamps(i) - timeStamps(i-1))*1000;
end
hist(timeDiff);

fprintf('Run-time: %f \n', timeStamps(end)-timeStamps(1));

%% Data
rawData = tableData{:, colNames(1:22)}';
data = rawData;

% surf(rawData(:,1:1000))



% Rereference
% data = bsxfun(@minus, rawData, mean(rawData));


% Filtering
[b, a] = butter(4,[1 45]/(250/2));
for i=1:22
    data(i,:) = filtfilt(b, a, data(i,:));
end

% Zero mean
data = bsxfun(@minus, data, mean(data, 2));





%% Show frequency plots
% Welch's method

% channel = 2;

Fs = 250;

winSize = 512; % = ~1 second
% size(data,2);

[B, A] = butter(4,[1 45]/(Fs/2));

figure(22);
for i=0:floor(size(data,2)/winSize)-1
    dataBlock = rawData(:,1+(i*winSize):((i+1)*winSize));
    
    for j=1:22
        dataBlock(j,:) = filtfilt(B, A, dataBlock(j,:));
    end
    
    [pxx, f] = pwelch(dataBlock');
    freq = 0:Fs/(2*size(f,1)-1):Fs/2;
    plot(freq,10*log10(pxx(:,1)));
%     h=surf(10*log10(pxx')); view(0,90);
%     set(h,'XData',freq);
    pause(0.5);
end


%%
X = fft(rawData')';
N = size(X,2);

X = 20*log10(abs(X)/N);
X = X(:,1:N/2+1);

freq = 0:Fs/(2*size(X,2)-1):Fs/2;

% dT = (tableData{:, 'TimeStamp'}(end)-tableData{:, 'TimeStamp'}(1))/numel(tableData{:, 'TimeStamp'});
plot(freq, X(1,:));

%% Welch
[pxx, f] = pwelch(rawData');
freq = 0:Fs/(2*size(f,1)-1):Fs/2;

% h=surf(10*log10(pxx'));
% set(h,'XData',freq);
figure;
plot(freq,10*log10(pxx(:,1)));




%% Source localization and 3D brain plotting

easyCapModelFull = importdata('easyCapModel.mat');
basis = load('basisFunctions.mat');
basisFunctions = basis.IDX2;

numSources = numel(basisFunctions);

% basisFunctions = sort(randperm(size(easyCapModelFull,2), numSources));
forwardModel = easyCapModelFull(:, basisFunctions);
sources = [];

blockSize = 64;

[B, A] = butter(4,[1 45]/(250/2));

for i=0:floor(size(rawData,2)/blockSize)-1
    
    dataBlock = rawData(:,1+(i*blockSize):((i+1)*blockSize));
    
    for j=1:size(rawData,1)
        dataBlock(j,:) = filtfilt(B, A, dataBlock(j,:));
    end
    
    init_alphas=1*ones(numSources, 1);
%     [~, ~, Mtemp, ~] = MARDv1(init_alphas, 1, forwardModel, dataBlock);
    [~, ~, M, ~] = MARDv2(init_alphas, 1, forwardModel, dataBlock);
        
%     lambda = 1e5;
%     PhiTPhiReg = forwardModel'*forwardModel + lambda*eye(numSources);
%     M = PhiTPhiReg\(forwardModel' * dataBlock);
    
    sources = [sources M];
%     surf(M);
%     pause(0.5);
    
    if mod(i, 50) == 0
       fprintf('%1.2f %% \n', (i/(floor(size(rawData,2)/blockSize)-1)*100))
    end
end
disp('done');


%%
QG = basis.QG;

vertface = load('vertface');
face = vertface.face2;
vert = vertface.vert2;
% vert = vert(idx,:);

sources = QG(:,basisFunctions) * sources;

%%
close(figure(1));
brainOpts.hfig = figure(1);
brainOpts.axes = axes('Parent',brainOpts.hfig,'Position',[.13 .15 .78 .75]); %self.BrainAxis;
brainHandles = setup3DBrain(vert, face, zeros(size(vert,1),1), brainOpts);

frames = [];

numFrames = floor(size(data,2)/blockSize)-1;

for i=0:numFrames
    
    sourceBlock = sources(:,1+(i*blockSize):((i+1)*blockSize));
    sourceBlock = sum(sourceBlock,2);

    brainHandles = plot3DBrain(brainHandles, sourceBlock);
    brainHandles.axes.View = [i*(360/numFrames) 0];
    
    F = getframe(gca);
    frames = [frames F];
%     drawnow;
%     pause
    pause(0.05);
end
%%

% h= imshow(frames(1).cdata);
f=figure;
movie(f, frames, 1, 10, [30 30 30 30])


%%

% figure(1)
filename = 'brainSpin2.gif';
for i=1:numel(frames)
%       y = x.^n;
%       plot(x,y)
%       drawnow
%       frame = getframe(1);
      im = frame2im(frames(i));
      [imind,cm] = rgb2ind(im,256);
      delay=0.15;
      if i == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',delay);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',delay);
      end
end
disp('done');





%%

%%
f1 = figure(1);
%
opts.hfig = f1;
plot_3Dbrain(vert, face, sum(sources(:,1:128),2), opts);

%%

for i=0:floor(size(data,2)/blockSize)-1
    sourceBlock = sources(:,1+(i*blockSize):((i+1)*blockSize));
    surf(sourceBlock);
%     plot_3Dbrain(vert, face, sourceBlock);
    pause(0.4);
end


