function data = preprocess( data, options )
% Preprocesses data according to the settings specified by options
% Including:
%               Removal of bad channels
%               Filtering
%               Re-refering                

    numChannels = size(data, 1);

    % Remove bad channels
    if ~isempty(options.bad_chans)
        data(options.bad_chans,:)=[];
    end

    if options.filter
        for i=1:numChannels
            data(i,:) = filtfilt(options.filterB, options.filterA, data(i,:));
        end
    end
    
    % Optionally reref the data
    if options.reref
        data = bsxfun(@minus, data, mean(data));
    end


end

