% EEG data visualization / inspection

% Load data
path = 'EEGData/';
% fileName = 'log17-03-2017 11-41-04.csv'; % 90%
% fileName = 'log17-03-2017 12-53-35.csv'; % flawless
% fileName = 'log17-03-2017 13-09-47.csv'; % 91%
% fileName = 'log17-03-2017 13-22-53.csv'; % 
fileName = 'log20-03-2017 10-49-09.csv'; % 

tableData = readtable([path fileName]);

colNames = tableData.Properties.VariableNames;
disp(colNames);

disp(sprintf('%d blocks', numel(tableData{:,'blockSize'})));


%% Timing
plot(tableData{:, 'bufferTimeStamp_end_'}, 'o');

timeStamps = tableData{:, 'timeStamp_end_'};
timeDiff = [];
i = 2;

for i=2:numel(timeStamps)
    timeDiff(i-1) = (timeStamps(i) - timeStamps(i-1))*1000;
end
hist(timeDiff);


%% Duration of loop

subplot(2,1,1); plot(tableData{:,'updateDuration'});
set(gca, 'YScale','log');
subplot(2,1,2); plot(tableData{:,'blockSize'});


%%

blockSizes = tableData{:,'blockSize'};
%bufferBlockSizes = tableData{:,'bufferBlockSize'};

[counts, centers] = hist(blockSizes, [32 64 96]);
bar(centers, counts/sum(counts));
title(sprintf('Histogram of %d blocks', numel(tableData{:,'blockSize'})));
xlabel('Block size');



