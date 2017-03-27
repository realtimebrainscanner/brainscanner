% instantiate the library
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

disp('Now receiving chunked data.....');



timeStamps = [];
timeDiff = [];
i = 2;
timeStamps(1) = 0;
data = [];

t0 = tic;
while true && i < 900
    % get chunk from the inlet
    [chunk,stamps] = inlet.pull_chunk();
    for s=1:length(stamps)
        % and display it
%         fprintf('%.2f\t',chunk(:,s));
%         fprintf('%.5f\n',stamps(s));
    end
    
    if ~isempty(stamps) 
        toc(t0)
        t0 = tic;
        
        timeStamps(i) = stamps(end);
        timeDiff(i) = (stamps(end) - timeStamps(i-1))*1000;
        dataIdx = size(chunk,2);
        if size(chunk,2) ~= 32 
            disp('Size issue');
            disp(size(chunk,2));
            i
        end
        data = [data chunk];
        i = i+1;
    else
%         disp('what');
    end
    
    
    pause(0.05);
end


inlet.close_stream();
inlet.delete();

timeDiff(1:2) = [];
timeStamps(1) = [];
%%

plot(timeDiff);
figure;
plot(timeStamps);

%%
dSize = 32;

equalRows = [];

for i=1:(size(data,2)/32)-2
    idx = i;
    % idx = 61;
    
    % figure(9); surf(data(:,dSize*(idx):dSize*(idx+1)));
    % figure(10); surf(data(:,dSize*(idx+1)+1:dSize*(idx+2)+1));
    
    chunk1 = data(:,dSize*(idx):dSize*(idx+1)-1);
    chunk2 = data(:,dSize*(idx+1):dSize*(idx+2)-1);
    
    % numel(find((chunk1 == chunk2) > 0));
    % [C, ia, ib] = intersect(chunk1', chunk2', 'rows');
    % idRows = numel(intersect(chunk1', chunk2', 'rows'));
    
%     equalRows(i) = numel(find((chunk1-chunk2) == 0));
    equalRows(i) = 0;
    for j=1:dSize
        if numel(intersect(chunk1(:,j)', chunk2', 'rows')) == 24
            equalRows(i) = equalRows(i) +1;
            j
        end
    end
    
    % if numel(C) > 0
    %     numel(C)
    % % C
    %     idx
    % end
    % (chunk1 == chunk2)
    % pause
end
plot(equalRows)

