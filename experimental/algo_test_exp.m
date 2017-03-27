% EEG data visualization / inspection
clear;

% Load data
path = 'EEGData/';
fileName = 'raw_22-03-2017 12-59-20.csv'; % First half 30 Hz sine and 40 Hz sine last half
tableData = readtable([path fileName]);

colNames = tableData.Properties.VariableNames;

%% Data
rawData = tableData{:, colNames(1:22)}';
data = rawData;

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
    sources = [sources M];
    
    surf(M);
    pause(0.5);
    
    if mod(i, 50) == 0
       fprintf('%1.2f %% \n', (i/(floor(size(rawData,2)/blockSize)-1)*100))
    end
end