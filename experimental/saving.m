
dataHeader = 'Fp1,Fp2,F3,F4,C3,C4,P3,P4,O1,O2,F7,F8,T7,T8,P7,P8,Fz,Cz,Pz,AFz,Cpz,POz,TimeStamp,Event,EventTime,EventFile';

data=randn(22, 64);
timeStamps=randn(1, 64);

% {};
experimentData={};
experimentData{1} = 'hej';
experimentData{2} = 'med';
experimentData{3} = 'dig';

experimentData = [num2cell(NaN(64-1, 3)); experimentData];
%%


dataTemp = [data' timeStamps']';
some=num2cell(dataTemp);

dataToWrite = [some' experimentData]';

%%
self.dataFileFormat = '';
for i=1:22
    self.dataFileFormat = [self.dataFileFormat '%1.4f,'];
end
self.dataFileFormat = [self.dataFileFormat '%1.6f,%s,%s,%1.4f\n'];
% self.dataFileFormat = [self.dataFileFormat '%s,'];


% dlmwrite('testtt.csv', dataToWrite, '-write') ;


fileID = fopen('testt.csv', 'w');

fprintf(fileID, self.dataFileFormat, dataToWrite{:,:});
fclose(fileID);