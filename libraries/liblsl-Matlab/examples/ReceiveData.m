%% instantiate the library
disp('Loading the library...');
lib = lsl_loadlib();

% resolve a stream...
disp('Resolving an EEG stream...');
result = {};
while isempty(result)
    result = lsl_resolve_byprop(lib,'name','openvibeSignal'); end

% create a new inlet
disp('Opening an inlet...');
inlet = lsl_inlet(result{1});

disp('Now receiving data...');

timeStamps = [];
timeDiff = [];
i = 2;
timeStamps(1) = 0;

while true && i < 10002
    % get data from the inlet
    [vec,ts] = inlet.pull_sample();
    % and display it
%     fprintf('%.2f\t',vec);

%     fprintf('%.5f\n',ts);
    
    timeStamps(i) = ts;
    timeDiff(i) = (ts - timeStamps(i-1))*1000;
    
    i = i+1;
end
inlet.close_stream();
inlet.delete();

timeStamps(1) = [];
timeDiff(1:2) = [];
%%

% plot(timeStamps);

plot(timeDiff);