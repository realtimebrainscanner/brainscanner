classdef EEGStream
    
    properties
        timeStamps
    end
    
    methods (Static)
        
        % create an inlet to read from the stream with the given name
        function inlet = CreateInlet(lib,opts)
            % look for the desired device
            result = {};
            disp(['Looking for a stream with name=' opts.streamname ' ...']);
            while isempty(result)
                result = lsl_resolve_byprop(lib,'name',opts.streamname); end
            % create a new inlet
            disp('Opening an inlet...');
            inlet = lsl_inlet(result{1},opts.bufferrange);
        end
        
        % create a new stream buffer to hold our data
        function stream = CreateStreambuffer(opts,info)
            stream.srate = info.nominal_srate();
            stream.chanlocs = struct('labels',EEGStream.DeriveChannelLabels(info));
            stream.buffer = zeros(length(stream.chanlocs),max(max(opts.bufferrange,opts.timerange)*stream.srate,100));
            [stream.nsamples,stream.state] = deal(0,[]);
        end
        
        % derive a list of channel labels for the given stream info
        function channels = DeriveChannelLabels(info)
            channels = {};
            ch = info.desc().child('channels').child('channel');
            while ~ch.empty()
                name = ch.child_value_n('label');
                if name
                    channels{end+1} = name; end %#ok<AGROW>
                ch = ch.next_sibling_n('channel');
            end
            if length(channels) ~= info.channel_count()
                disp('The number of channels in the steam does not match the number of labeled channel records. Using numbered labels.');
                channels = cellfun(@(k)['Ch' num2str(k)],num2cell(1:info.channel_count(),1),'UniformOutput',false);
            end
        end
        
        % update display with new data
        function OnTimer(varargin)
            try
                opts = varargin{3};
                inlet = varargin{4};
                stream = varargin{5};
                lib = varargin{6};
                
                % pull a new chunk from LSL
                [chunk,timestamps] = inlet.pull_chunk();
                if isempty(chunk)
                    return; end
                
                % optionally filter the chunk
                chunk(~isfinite(chunk(:))) = 0;
%                 if ~exist(opts.lpFilter, 'var') && ~exist(opts.hpFilter, 'var')
                if ~exist(opts.bpFilter, 'var')
                    
                    for j=opts.channelsrange
                        chunk = filter(opts.bpFilter,chunk(j,:)); 
                    end
                end
                


                % append it to the stream buffer
                [stream.nsamples,stream.buffer(:,1+mod(stream.nsamples:stream.nsamples+size(chunk,2)-1,size(stream.buffer,2)))] = deal(stream.nsamples + size(chunk,2),chunk);
                
                % extract channels/samples to plot
                samples_to_get = min(size(stream.buffer,2), round(stream.srate*opts.timerange));
                channels_to_get = intersect(opts.channelrange + opts.pageoffset*length(opts.channelrange), 1:size(stream.buffer,1));
                stream.data = stream.buffer(channels_to_get,1+round(mod(stream.nsamples-samples_to_get: stream.srate/opts.samplingrate : stream.nsamples-1,size(stream.buffer,2))));
                [stream.nbchan,stream.pnts,stream.trials] = size(stream.data);
                stream.xmax = max(timestamps) - lsl_local_clock(lib);
                stream.xmin = stream.xmax - (samples_to_get-1)/stream.srate;
                
                
                
%                 h = gcf;
%                 this = get(h, 'Children');
%                 if ~isempty(this)
%                     this = get(this, 'Children');
%                     oldTimeStamps = get(this, 'YData');
%                     updates = get(this, 'XData'); 
%                     if isempty(updates)
%                         updates = 1;
%                     end
%                     oldTimeStamps = [oldTimeStamps (timestamps(end)-updates(end))];
%                     set(this, 'YData', oldTimeStamps, 'XData', [updates updates(end)+1]);
%                 else
%                     plot(1, timestamps(end));
%                 end

%                 disp([timestamps(1) timestamps(end) numel(timestamps)]);


                numSamples = numel(timestamps);
                numChannels = size(chunk,1);
                
                
                % optionally post-process the data
                if opts.reref
                    stream.data = bsxfun(@minus,stream.data,mean(stream.data)); end
                if opts.standardize
                    stream.data = bsxfun(@times,stream.data,1./std(stream.data,[],2)); end
                if opts.zeromean
                    stream.data = bsxfun(@minus, stream.data, mean(stream.data,2)); end
                
                % arrange for plotting
                plotoffsets = (0:stream.nbchan-1)'*opts.datascale;
                plotdata = bsxfun(@plus, stream.data, plotoffsets);
                plottime = linspace(stream.xmin,stream.xmax,stream.pnts);
                
                
                stream.data = bsxfun(@times,chunk,1./std(chunk,[],2));
                stream.data = bsxfun(@minus, stream.data, mean(stream.data,2));
            
                figure(2);
                h = gcf;
                this = get(h, 'Children');
                if ~isempty(this)
                    this = get(this, 'Children');
                    set(this, 'ZData', stream.data);
                else
                    surf(stream.data);
                end
                
                
                
%                 figure(2); surf(stream.data);
%                 view(0,30);
%                 size(plotdata)
                
                
%                 % update graphics
%                 if isempty(lines)
%                     lines = plot(ax,plottime,plotdata);
%                     title(ax,opts.streamname);
%                     xlabel(ax,'Time (sec)','FontSize',12);
%                     ylabel(ax,'Activation','FontSize',12);
%                 else
%                     for k=1:min(length(lines),size(plotdata,1))
%                         set(lines(k),'Xdata',plottime, 'Ydata',plotdata(k,:)); end
%                     for k = size(plotdata,1)+1:length(lines)
%                         set(lines(k),'Ydata',nan(stream.pnts,1)); end
%                 end
%                 
%                 % update the axis limit and tickmarks
%                 axis(ax,[stream.xmin stream.xmax -opts.datascale stream.nbchan*opts.datascale + opts.datascale]);
%                 set(ax, 'YTick',plotoffsets, 'YTickLabel',{stream.chanlocs(channels_to_get).labels});
%                 
%                 drawnow;
            catch e
                % display error message
                fprintf('vis_stream error: %s\noccurred in:\n',e.message);
                for st = e.stack'
                    if ~isdeployed
                        try
                            fprintf('   <a href="matlab:opentoline(''%s'',%i)">%s</a>: %i\n',st.file,st.line,st.name,st.line);
                        catch
                            fprintf('   <a href="matlab:edit %s">%s</a>: %i\n',st.file,st.name,st.line);
                        end
                    else
                        fprintf('   %s: %i\n',st.file,st.line);
                    end
                end
                on_close();
            end
        end
        
        
%             % handle key presses
%     function on_key(key)
%         switch lower(key)
%             case 'enter' 
%                 stop(th);
%         end
%     end
        
        
        
    end %methods
    
end