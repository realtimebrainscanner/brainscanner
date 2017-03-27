classdef simpleSourceLocalization < handle
    
    properties
        Figure                  % Graphics handles
        Axis
        Line
        Surface
        TickerText
        TickerEdit
        
        Figure2
        Surface2
        Axis2
        
        TimerCount
        
        % Source localization
        channels = 22; 
        lambda = 1e5;   % Ridge regularizing factor
        
        numFuncs = 100;
        numActiveFuncs = 10;
        timeSteps = 64;
        trueBeta = 25;
        trueAlpha = 2;
        A
        idx
        idx_freq
        w
        
        % Filter
        lpFilter = designfilt('lowpassiir','FilterOrder',6,'PassbandFrequency',40,'PassbandRipple',0.2, 'SampleRate',500);
        hpFilter = designfilt('highpassiir','FilterOrder',6,'PassbandFrequency',1,'PassbandRipple',0.2, 'SampleRate',500);
%         bpFilter = designfilt('bandstopiir','FilterOrder',2, 'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, 'DesignMethod','butter','SampleRate',500);
        
        % Imported data
        dataPath = ('recordings/');
        fileName = 'rec1.mat';
        data
        
        Timer                   % Timer object to get updated prices
        TimerUpdateRate = 0.5;     % In seconds 
        Some = 'test'
    end
    
    
    
    methods
        
        function app = simpleSourceLocalization
            
            app.Figure = figure('MenuBar','none',...           % Main figure
                'NumberTitle','off','Name','True sources',...
                'Position', [100,100,560,420], ...
                'CloseRequestFcn',@app.closeApp) ;
            app.Axis = axes('Parent',app.Figure,...            % Axis for prices
                'Position',[.13 .15 .78 .75]);            
            app.Surface = surf(zeros(app.numFuncs, app.timeSteps));
            ylabel(app.Axis,'Source index') ;
            xlabel(app.Axis,'Time') ;
%             set(app.Axis,'XTickLabel','') ;
            title(app.Axis,'Filtered signal')
%             set(app.Axis, 'ZLim', [-2, 2]);
            
            app.Figure2 = figure('MenuBar','none',...           % Secondary figure
                'NumberTitle','off','Name','Source localization',...
                'Position', [800,100,560,420], ...
                'CloseRequestFcn',@app.closeApp) ;
            app.Axis2 = axes('Parent',app.Figure2, 'Position',[.13 .15 .78 .75]); 

            app.Surface2 = surf(zeros(app.numFuncs, app.timeSteps));
            ylabel(app.Axis2,'Source index') ;
            xlabel(app.Axis2,'Time') ;
            title(app.Axis2,'Estimated sources')
%             set(app.Axis2, 'ZLim', [-5, 5]);
            
%             figure('Position',[10000,10000,400,300])
            
            initialize(app);
            
            app.Timer = timer;                                 % Create timer
            app.Timer.ExecutionMode = 'fixedRate';
            app.Timer.Period = app.TimerUpdateRate;
            app.Timer.TimerFcn = @app.inputCallback;
            app.Timer.TasksToExecute = size(app.data,2)/app.timeSteps;
            start(app.Timer);
        end
        
        function closeApp(app,hObject,eventdata)
            try
                stop(app.Timer)
                delete(app.Timer)
            end
            delete(app.Figure)
            delete(app.Figure2)
        end      
        
        
        function inputCallback(app,hObject,eventdata)
            % This function runs when the timer updates
            timerUpdate = app.Some;
            try
                measurements = readData(app);
            catch
                errordlg(['Error in ' timerUpdate])
            end
            
            % Filter 
            measurements = filterData(measurements, app);
            
            % Mean signal
%             measurements = zeroMean(measurements, app);
            set(app.Surface,'ZData', zeroMean(measurements, app));
            
            plotFrequencies(measurements, app);
%             M = recover(filteredData, app);
%             disp(measurements);
            
        end
        
        
        function M = recover(measurements, app)
            % Ridge
            M = [];
            for l=1:app.timeSteps
                M(:,l) = (app.A'*app.A + app.lambda*eye(size(app.A, 2)))\(app.A'*measurements(:,l));
            end
            
            % Min norm
            
            % MARD
%             init_alphas=0.1*ones(app.numFuncs,1);
%             [alphas, beta, M, llh] = MARD(init_alphas, 1, app.A, measurements);

            set(app.Surface2,'ZData', M)
        end
        
        
        % Show spectogram of imported data
        function plotFrequencies(data, app)
            winSize = 80;
            fftSignal = fft(data',winSize)'; 
            
            L = size(fftSignal,2);
            P2 = 20*log10(abs(fftSignal/L));
            P1 = P2(:,1:L/2+1);
            P1(:,2:end-1) = 2*P1(:,2:end-1);
            
            f = 500*(0:(L/2))/L;
            
            % Plot
%             set(app.Surface2,'ZData', P1(:,1:10), 'XData', f(:,1:10));
            set(app.Surface2,'ZData', P1, 'XData', f);
            ylabel(app.Axis2,'Sensor index');
            xlabel(app.Axis2,'Frequency');
            zlabel(app.Axis2,'dB');
            title(app.Axis2,'Freq spectrum of sensors');
%             view(0,90);
        end
        
        
        function initialize(app)
            % Syntehtic
%             initializeSynthetic(app);
            % Imported
            initializeImported(app);
            % Real
        end
        
       function initializeImported(app)
            app.TimerCount = 0;
            
            app.data = importdata([app.dataPath app.fileName]);
            app.data(20:21,:) = []; 
            
            easyCapModelFull = importdata('easyCapModel.mat');
            app.A = easyCapModelFull(:, sort(randperm(size(easyCapModelFull,2), app.numFuncs)));
%             app.A = easyCapModelFull(:, 1:app.numFuncs);
        end
        
        function initializeSynthetic(app)
            app.TimerCount = 0;
            app.A = randn(app.channels, app.numFuncs);
            app.w = zeros(app.numFuncs, 1);
            app.idx=sort(randperm(app.numFuncs,app.numActiveFuncs));
            factor = 1;
            app.w(app.idx) = factor*normrnd(0,sqrt(1/app.trueAlpha), [1 size(app.idx)]);
            
            app.idx_freq = zeros(app.numFuncs,1);
            for i=1:app.numFuncs
                app.idx_freq(i) = 0.5*randn;
            end
            
        end
        
        
        % Read data from file
        function [measurements] = readData(app)
            % If synthetic
%             measurements = generateData(app)
            
            % Import
            readRange = [1+app.TimerCount*app.timeSteps:app.TimerCount*app.timeSteps+app.timeSteps];
%             readRange = 
            measurements = app.data(:,readRange);
            
            app.TimerCount = app.TimerCount+1;
        end
        
        
        function [filteredData] = filterData(measurements, app)
            frame = double(measurements);
            for j=1:app.channels
                frame(j,:) = filtfilt(app.lpFilter, frame(j,:));
                frame(j,:) = filtfilt(app.hpFilter, frame(j,:));
                
            end
            filteredData = frame;
        end
        
        
        function [measurements] = zeroMean(measurements, app)
            measurements = bsxfun(@times,measurements,1./std(measurements,[],2));
            measurements = bsxfun(@minus, measurements, mean(measurements,2));
        end
        
        
        % Synthetic data
        function [measurements, x] = generateData(app)
            x=zeros(size(app.A,2), app.timeSteps);
            
            range = [1+app.TimerCount*app.timeSteps:app.TimerCount*app.timeSteps+app.timeSteps];
            for i=1:size(app.A,2)
                x(i,:)=app.w(i)*sin(range*app.idx_freq(i));
            end
            
            y = app.A*x;
            noise = normrnd(0, sqrt(1/app.trueBeta), [app.channels app.timeSteps]);

%             app.TimerCount = app.TimerCount+1;
            measurements = y + noise;
        end
        
    end
end