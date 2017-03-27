

data = importdata('recordings/rec1.mat');
data(20:21,:) = [];
measurements = data(:,1:20);

A = randn(size(data,1),100);
MARD(ones(100,1), 1, A, measurements);


%%
% wo = 50/(500/2);  bw = wo/35;
% [b,a] = iirnotch(wo,bw); 
% fvtool(b,a);

bsFilt = designfilt('bandstopiir','FilterOrder',6, ...
               'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
               'DesignMethod','butter','SampleRate',500);
           
lpFilt = designfilt('lowpassiir','FilterOrder',6, ...
                'PassbandFrequency',40,'PassbandRipple',0.2, 'SampleRate',500);
            
hpFilt = designfilt('highpassiir','FilterOrder',6, ...
         'PassbandFrequency',1,'PassbandRipple',0.2, ...
         'SampleRate',500);

% bpFilt = designfilt('bandpassiir','FilterOrder',6, ...
%          'HalfPowerFrequency1',2,'HalfPowerFrequency2',40, ...
%          'SampleRate',500);
     
bpFilt = designfilt('bandpassfir','FilterOrder',20, ...
         'CutoffFrequency1',0.2,'CutoffFrequency2',45, ...
         'SampleRate',500);     
     
% bpFilt2 = designfilt('bandstopiir','FilterOrder',2, ...
%                'HalfPowerFrequency1',124,'HalfPowerFrequency2',126, ...
%                'DesignMethod','butter','SampleRate',500);     


% testData = double(EEG.data);
testData = double(data);
           
% filteredData = filtfilt(bpFilt,testData);
% filteredData2 = filtfilt(hpFilt, filteredData);

% filteredData3 = filtfilt(bpFilt2, filteredData);

% plot(filteredData); hold on;
% plot(filteredData2);

%%
timeSteps = 80;
channels = 22;
frames = size(testData,2)/timeSteps;

for i = 0 : frames-1
    frameStart = 1+timeSteps*i;
     
    frame = testData(:,frameStart:frameStart+timeSteps-1);
    frame2 = testData(:,frameStart:frameStart+timeSteps-1);
    
    for j=1:channels
        frame(j,:) = filtfilt(lpFilt, frame(j,:));
        frame(j,:) = filtfilt(hpFilt, frame(j,:));
    
        frame2(j,:) = filter(lpFilt, frame2(j,:));
        frame2(j,:) = filter(hpFilt, frame2(j,:));
    end

%     frame2 = testData(:,frameStart:frameStart+timeSteps-1);
%     frame2 = filtfilt(lpFilt,testData(:,frameStart:frameStart+timeSteps-1));
%     frame2 = filtfilt(hpFilt,frame2);
%     frame2 = filter(lpFilt,frame2');


    A = fft(frame')';
%     A_abs = 20*log10(abs(A));
    L = size(A,2);
    P2 = 20*log10(abs(A/L));
    P1 = P2(:,1:L/2+1);
    P1(:,2:end-1) = 2*P1(:,2:end-1);
    
    f = 500*(0:(L/2))/L;
    f1 = figure(1);
%     h=surf(P1); %view(0,90);
    h=surf(P1); view(8,16);
    set(h, 'XData', f);
    set(f1, 'Position', [100,100,560,420]);
    title(sprintf('Frame: %i', i+1));
    
    A = fft(frame2')';
%     fftData = 20*log10(abs(A/L));
%     fftData = fftData(:,1:size(fftData,2)/2+1);
    
    
    P2 = 20*log10(abs(A/L));
    P1 = P2(:,1:L/2+1);
    P1(:,2:end-1) = 2*P1(:,2:end-1);
    
    f = 500*(0:(L/2))/L;
    f2 = figure(2);
%     h=surf(P1); %view(0,90);
    h=surf(P1); view(8,16);
    set(h, 'XData', f);
    title(sprintf('Frame: %i', i+1));
    set(f2, 'Position', [800,100,560,420]);
    
    
    
%     pause(0.01);
    pause
end


%%

   %% Time specifications:
%    Fs = 100;                      % samples per second
%    dt = 1/Fs;                     % seconds per sample
%    StopTime = 1;                  % seconds
%    t = (0:dt:StopTime-dt)';
%    N = size(t,1);
   
   % Fourier Transform:
%    x = frame(1,:); 

   Fs = 500; 

   frameStart = 1079;
   frame = testData(:,frameStart:frameStart+79);
   
   
   
   for j=1:channels
%        frame(j,:) = filtfilt(lpFilt, frame(j,:));
%        frame(j,:) = filtfilt(hpFilt, frame(j,:));
       frame(j,:) = filter(bpFilt, frame(j,:));
   end
   
%    frame = bsxfun(@times, frame, 1./std(frame,[],2));
   frame = bsxfun(@minus, frame, mean(frame,2)); 
   x = frame;%(1,:);
   
   X = (fft(x')');
   
   N = size(X,2);
   % Frequency specifications:
   dF = Fs/N;                      % hertz
%    f = -Fs/2:dF:Fs/2-dF;           % hertz
   
   f = 0:dF:(Fs/2);
   
   % Plot the spectrum:
   figure;
   fftData = 20*log10(abs(X)/N);
   fftData = fftData(:,1:size(fftData,2)/2+1);
   
   if size(fftData,1) > 1
       h = surf(frame);
%        set(h, 'XData', f);
%        set(gca, 'ZScale','log');
   else
       h = plot(f,fftData);
       set(gca, 'YScale','log');
   end

   xlabel('Frequency (in hertz)');
   title('Magnitude Response');

%%
%
range=1:80;
rawData = testData(12,range);
filtData = filtfilt(lpFilt,rawData); % filteredData2(1,range);
filtData = filtfilt(hpFilt,filtData); % filteredData2(1,range);

L = numel(filtData);
Y = fft(rawData);
Yfilt = fft(filtData);
P2raw = 20*log10(abs(Y/L));
P2filt = 20*log10(abs(Yfilt/L));

P1raw = P2raw(1:L/2+1);
P1raw(2:end-1) = 2*P1raw(2:end-1);

P1filt = P2filt(1:L/2+1);
P1filt(2:end-1) = 2*P1filt(2:end-1);

f = 500*(0:(L/2))/L;
figure(99);

plot(f,P1raw); hold on;
plot(f,P1filt, '--'); hold off;
set(gca, 'YScale','log');

%%
% load openloop60hertz, openLoop = openLoopVoltage;

% openLoop = testData(1,:);
openLoop = double(EEG.data(1,:));

% rawData = testData(1,:);
% openLoop = rawData;
rawData = openLoop;
filtData = filtfilt(lpFilt,rawData); % filteredData2(1,range);
filtData = filtfilt(hpFilt,filtData); % filteredData2(1,range);

Fs = 1000;
t = (0:length(openLoop)-1)/Fs;

plot(t,openLoop)
ylabel('Voltage (V)')
xlabel('Time (s)')
title('Open-Loop Voltage with 60 Hz Noise')
grid

bsFilt = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',Fs);

buttLoop = filtData; % filtfilt(bsFilt,openLoop);

plot(t,openLoop,t,buttLoop)
ylabel('Voltage (V)')
xlabel('Time (s)')
title('Open-Loop Voltage')
legend('Unfiltered','Filtered')
grid


fftData = openLoop;
L = numel(fftData);
Y = fft(fftData);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

fftData2 = buttLoop;
Y2 = fft(fftData2);
P22 = abs(Y2/L);
P12 = P22(1:L/2+1);
P12(2:end-1) = 2*P12(2:end-1);

f = 500*(0:(L/2))/L;
figure(99);
plot(f,P1); hold on;
plot(f,P12); hold off;
set(gca, 'YScale','log');
