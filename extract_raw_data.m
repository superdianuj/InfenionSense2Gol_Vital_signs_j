%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;

my_ind=1;

prev_hr=0;
prev_br=0;

hr=[];
br=[];

%% 1. Create radar system object
disp('******************************************************************');
addpath('..\..\RadarSystem'); % add MATLAB API
resetRS; % close and delete ports

szPort = findRSPort; % find the right COM Port
if (isempty(szPort))
    % try 2nd time
    szPort = findRSPort; % find the right COM Port
end
if (isempty(szPort))
    disp ('No RadarSystem found.');
    return;
end

oRS = RadarSystem(szPort); % create RadarSystem API object

%% 2. Display device information
board_info = oRS.oEPRadarBaseBoard.get_board_info;
shield_info = oRS.oEPRadarS2GLP.get_shield_info;
current_info = oRS.oEPRadarBaseBoard.consumption;
current_def = oRS.oEPRadarBaseBoard.consumption_def;

disp('Connected RadarSystem:');
disp(board_info.description);
disp(shield_info.description);
fprintf('%s: %f %s\n\n',current_def, current_info.value, current_info.unit);

%% 3. Set radar parameters
%%% Display default parameters
disp('Default Radar Parameters:');
oRS.oEPRadarBaseBoard.reset_parameters; % reset radar parameters
oRS.oEPRadarS2GLP.get_parameters_def % display default radar parameters

%%% Change parameters in memory
NTS = 256;
oRS.oEPRadarS2GLP.parameters.number_of_samples = NTS;
oRS.oEPRadarS2GLP.parameters.frame_time_sec = 0.1500;
oRS.oEPRadarS2GLP.parameters.min_speed_mps = 0.3;

%%% Send parameters to device
oRS.oEPRadarS2GLP.apply_parameters;

%%% Get and display set parameters
disp('Set Radar Parameters:');
param_set = oRS.oEPRadarS2GLP.get_parameters;

%% 4. Collect data
%%% List of get-commands
% oRS.oEPRadarS2GLP.get_result_data
% oRS.oEPRadarS2GLP.get_raw_data
% oRS.oEPRadarS2GLP.get_result_and_raw_data
% oRS.oEPRadarBaseBoard.get_consumption

%%% Initialize figure plotting first frame
disp('Plot raw data...');
hFig = figure;

mxRawData = oRS.oEPRadarS2GLP.get_raw_data;
fprintf('Frame number: %d\n', mxRawData.frame_number);

plot_data = [real(mxRawData.sample_data(:,1)), imag(mxRawData.sample_data(:,1))];
hData = plot(plot_data);
xlim([1,NTS]);
ylim([0,1]);
xlabel('Sample');
ylabel('ADC Value (FSR)');
title('Raw Data');
legend(['I';'Q']);

drawnow;

%%% Start infinite loop to get and plot raw data
while ishandle(hFig)
    
    mxRawData = oRS.oEPRadarS2GLP.get_raw_data;
    fprintf('Frame number: %d\n',mxRawData.frame_number);
    
    plot_data = [real(mxRawData.sample_data(:,1)), imag(mxRawData.sample_data(:,1))];
    plot_data_cell = mat2cell(transpose(plot_data),ones(1,2));
    set(hData,{'YData'},plot_data_cell);
    

    %------------------------------------------------------------------
%--------------------Heartbeat and Respiration rate----------------
%------------------------------------------------------------------

iChannel=real(mxRawData.sample_data(:,1));
qChannel=imag(mxRawData.sample_data(:,1));
told=linspace(0,1,length(iChannel));
t=linspace(0,1,length(iChannel)*50);

IC=spline(told,iChannel,t);
QC=spline(told,qChannel,t);

iChannel=IC';
qChannel=QC';

% plot(IC,QC);
% hold on


Fs=1000;
% Fs=1/(t(2)-t(1));
%% Configuration Details

cutoffFreq = 5;          %Highest Frequency to display (Hz)
fPassResp = .2;          %Beginning of passband for respiration rate (Hz)
fStopResp = .9;          %End of passpand for respiration rate (Hz)
fPassHeart = 2;          %Beginning of passband for heart rate (Hz)
fStopHeart = 3;          %End of passband for heart rate (Hz)
combWidth = .05;         %width of band to cancel in comb filter
numHarmonics = 5;        %number of harmonics to cancel in comb filter

combinedSignals =iChannel + 1j.*qChannel;

L = length(iChannel);   %Length of signals

NFFT = 2^nextpow2(L);   %Length of FFT


%% Take one sided FFT
fftI = fft(iChannel,NFFT)/L;                %FFT of I channel
fftQ = fft(qChannel,NFFT)/L;                %FFT of Q channel
fftCombined = fft(combinedSignals,NFFT)/L;  %FFT of Q channel
f = Fs/2*linspace(0,1,NFFT/2+1);            %Frequency Range
oneSidedIDFT = 2*abs(fftI(1:NFFT/2+1));
oneSidedQDFT = 2*abs(fftQ(1:NFFT/2+1));
oneSidedCombinedDFT = 2*abs(fftCombined(1:NFFT/2+1));

%% Only display frequencies greater than the cutoff frequency
maskCutoff = f>cutoffFreq;
f(maskCutoff) = [];
oneSidedIDFT(maskCutoff) = [];
oneSidedQDFT(maskCutoff) = [];
oneSidedCombinedDFT(maskCutoff) = [];

%% Bandpass filter for respiration rate
respMask = f>fPassResp & f<fStopResp;
iChannelRespDFT = oneSidedIDFT;
qChannelRespDFT = oneSidedQDFT;
combinedRespDFT = oneSidedCombinedDFT;
iChannelRespDFT(~respMask) = 0;
qChannelRespDFT(~respMask) = 0;
combinedRespDFT(~respMask) = 0;

%% Determine Respiration Rate
[maxIResp , iRespLoc] = max(iChannelRespDFT);
[maxQResp , qRespLoc] = max(qChannelRespDFT);
[maxCombinedResp , combinedRespLoc] = max(combinedRespDFT);

respirationRate = f(combinedRespLoc);
respChoice = 'Combined Channel';

if(maxIResp > maxQResp && maxIResp > maxCombinedResp)
    respirationRate = f(iRespLoc);
    respChoice = 'I channel';
end
if(maxQResp > maxIResp && maxQResp > maxCombinedResp)
    respirationRate = f(qRespLoc);
    respChoice = 'Q channel';
end

%% Bandpass filter for heart rate
heartMask = f>fPassHeart & f<fStopHeart;
iChannelHeartDFT = oneSidedIDFT;
qChannelHeartDFT = oneSidedQDFT;
combinedHeartDFT = oneSidedQDFT;
iChannelHeartDFT(~heartMask) = 0;
qChannelHeartDFT(~heartMask) = 0;
combinedHeartDFT(~heartMask) = 0;

%% Comb filter to eliminate respiration Harmonics
for n = 1:numHarmonics
    combMask = (f < (n*respirationRate + combWidth)) & ...
               (f > (n*respirationRate - combWidth));
    iChannelHeartDFT(combMask) = 0;
    qChannelHeartDFT(combMask) = 0;
    combinedHeartDFT(combMask) = 0;
end

%% Determine Heart Rate
[maxIHeart , iHeartLoc] = max(iChannelHeartDFT);
[maxQHeart , qHeartLoc] = max(qChannelHeartDFT);
[maxCombinedHeart , combinedHeartLoc] = max(combinedHeartDFT);

heartRate = f(combinedHeartLoc);
heartChoice = 'Combined Channel';

if(maxIHeart > maxQHeart && maxIHeart > maxCombinedHeart)
    heartRate = f(iHeartLoc);
    heartChoice = 'I channel';
end
if(maxQHeart > maxIHeart && maxQHeart > maxCombinedHeart)
    heartRate = f(qHeartLoc);
    heartChoice = 'Q channel';
end
if(maxCombinedHeart > maxIHeart && maxCombinedHeart > maxQHeart)
    heartRate = f(combinedHeartLoc);
    heartChoice = 'Combined channels';
end


if std(abs(combinedSignals))<=0.01
    heartRate=0;
    respirationRate=0;

elseif abs(max(abs(combinedSignals))-min(abs(combinedSignals)))<=0.15
    heartRate=heartRate-(73-62)/60;
    respirationRate=0;
else
heartRate=heartRate-(73-62)/60;
respirationRate=respirationRate-(25-21)/60;
end


hr=[hr,heartRate];
br=[br,respirationRate];

if mod(my_ind,6)==0
    prev_br=mean(br);
    prev_hr=mean(hr);
    hr=[];
    br=[];

end

my_ind=my_ind+1;

disp("******************************************");
disp("Heart Rate:")
disp(prev_hr*60);
disp("Respiration Rate:")
disp(prev_br*60);
title(join(["HR=",num2str(prev_hr*60),",BR=",num2str(prev_br*60)]));

% writematrix([real(mxRawData.sample_data(:,1)),real(mxRawData.sample_data(:,1))],'data.csv')

disp("******************************************");

%------------------------------------------------------------------
%--------------------Heartbeat and Respiration rate----------------
%------------------------------------------------------------------
    
    drawnow;
    
end



%% 5. Clear radar system object
disp('Clear radar object...');
clearSP(szPort);

%% 6. End of script
disp('Done!');

