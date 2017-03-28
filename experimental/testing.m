
data = importdata('model/labels.mat');

sensors = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2','F7','F8','T7','T8','P7','P8','Fz','Cz','Pz','AFz','Cpz','POz'};
mappingSensors = {'C29','C16','C32','C10','D19','B22','A7','B4','A15','A28','D7','C7','D23','B26','D31','B11','C21','A1','A19','C19','A3','A21'};

indices = [];
for i=1:numel(mappingSensors)
    mappingSensor = mappingSensors(i);
    [found, index] = find(strcmp(mappingSensor, data) == 1);
    if found
        disp(sprintf('%s found at index: %i', mappingSensor{:}, index));
        indices = [indices index];
    end    
end

fullModel = importdata('model/forwardmodel_full.mat');

easyCapModel = fullModel(indices,:);


%%

% find(newLead(1,:) == 
% newLeadNorm = normc(newLead);

[tf, index] = ismember(A, newLead, 'rows');
sum(tf) + sum(index)


%%

correlations=zeros(128,1);
for i=1:128
    R=corrcoef(newLead(1,:), A(i,:));
    correlations(i)=R(2,1);
end
