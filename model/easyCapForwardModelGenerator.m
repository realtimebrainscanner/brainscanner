
data = importdata('labels.mat');

easyCapSensors = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2','F7','F8','T7','T8','P7','P8','Fz','Cz','Pz','AFz','Cpz','POz'};
SPM8Sensors = {'C29','C16','C32','C10','D19','B22','A7','B4','A15','A28','D7','C7','D23','B26','D31','B11','C21','A1','A19','C19','A3','A21'};

indices = [];
for i=1:numel(SPM8Sensors)
    SPM8Sensor = SPM8Sensors(i);
    [found, index] = find(strcmp(SPM8Sensor, data) == 1);
    if found
        disp(sprintf('%s found at index: %i', SPM8Sensor{:}, index));
        indices = [indices index];
    end    
end

fullModel = importdata('forwardmodel_full.mat');

easyCapModel = fullModel(indices,:);
save('easyCapModel','easyCapModel');

%%

