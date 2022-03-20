clear all
close all
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Self-Gathered Experimental Data
near_pos= table2array(readtable('junaid_restingState_nearPositioned.csv'));
figure
subplot(2,1,1)
plot(near_pos(1000:4000,1),near_pos(1000:4000,2));
grid on
xlabel('time(sec)')
ylabel('Amplitude(au)')
title('Near Positioned-I Signal');
subplot(2,1,2)
plot(near_pos(1000:4000,1),near_pos(1000:4000,3));
grid on
title('Near Positioned-Q Signal');
xlabel('time(sec)')
ylabel('Amplitude(au)')

iChannel=near_pos(100:4000,2);
qChannel=near_pos(100:4000,3);
t=near_pos(100:4000,1);

Fs=1/(t(2)-t(1));
numSecondsBeginning = 1; %Number of seconds to eliminate from beginning of signal
numSecondsEnd = 1;       %Number of seconds to eliminate from end of signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Established Experimental Data
% load GDN0001_1_Resting.mat
% figure ('Name','Inphase and Quadrature Signals')
% subplot(2,1,1)
% plot(radar_i)
% title('Inphase Signal');
% ylabel('Amplitude');
% xlabel('time(sec)');
% grid on
% axis ([0 12e5 -2000 4000]);
% 
% 
% subplot(2,1,2)
% plot(radar_q)
% title('Quadrature Signal');
% ylabel('Amplitude');
% xlabel('time(sec)');
% grid on
% axis ([0 12e5 -2000 4000]);
% 
% iChannel=radar_i;
% qChannel=radar_q;
% 
% t_ender=607.6;
% t_start=0;
% 
% t = linspace(t_start, t_ender,length(iChannel));
% Fs = fs_radar;   %Sampling Frequency
% numSecondsBeginning = 5; %Number of seconds to eliminate from beginning of signal
% numSecondsEnd = 5;       %Number of seconds to eliminate from end of signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Configuration Details
fileNum = 2;      %2-5

cutoffFreq = 5;          %Highest Frequency to display (Hz)
fPassResp = .2;          %Beginning of passband for respiration rate (Hz)
fStopResp = .9;          %End of passpand for respiration rate (Hz)
fPassHeart = 2;          %Beginning of passband for heart rate (Hz)
fStopHeart = 3;          %End of passband for heart rate (Hz)
combWidth = .05;         %width of band to cancel in comb filter
numHarmonics = 5;        %number of harmonics to cancel in comb filter

oner=ones(length(iChannel),1);
fun = @(x)sum((abs(iChannel-x(1)).^2+abs(qChannel-x(2)).^2-x(3)*oner.^2).^2);
x0 = [0,0,0];
x = fminsearch(fun,x0);


iChannel=iChannel-x(1)*oner;
qChannel=qChannel-x(2)*oner;


combinedSignals = iChannel + 1j.*qChannel;



L = length(iChannel);   %Length of signals
NFFT = 2^nextpow2(L);   %Length of FFT

%% Eliminate numSecondsBeginning of bad data at beginning
numSamplesBeginning = round(numSecondsBeginning*Fs);
t(1:numSamplesBeginning) = [];
iChannel(1:numSamplesBeginning) = [];
qChannel(1:numSamplesBeginning) = [];
combinedSignals(1:numSamplesBeginning) = [];

%% Eliminate numSecondsEnd of bad data at end
numSamplesEnd = round(numSecondsEnd*Fs);
t(end:-1:(end-numSamplesEnd)) = [];
iChannel(end:-1:(end-numSamplesEnd)) = [];
qChannel(end:-1:(end-numSamplesEnd)) = [];
combinedSignals(end:-1:(end-numSamplesEnd)) = [];

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

%% Plot I and Q
figure
subplot(3,1,1)
plot(t,iChannel)
xlabel('Time (s)')
ylabel('|i(t)|')
title('I Channel in Time Domain')
subplot(3,1,2)
plot(t,qChannel)
xlabel('Time (s)')
ylabel('|q(t)|')
title('Q Channel in Time Domain')
subplot(3,1,3)
plot(t,abs(combinedSignals))
xlabel('Time (s)')
ylabel('|c(t)|')
title('Combined Signals in Time Domain')

figure
subplot(3,1,1)
plot(f,oneSidedIDFT) 
title('Single-Sided Amplitude Spectrum of I channel FFT')
xlabel('Frequency (Hz)')
ylabel('|I(f)|')
subplot(3,1,2)
plot(f,oneSidedQDFT) 
title('Single-Sided Amplitude Spectrum of Q channel FFT')
xlabel('Frequency (Hz)')
ylabel('|Q(f)|')
subplot(3,1,3)
plot(f,oneSidedCombinedDFT) 
title('Single-Sided Amplitude Spectrum of Combined channels FFT')
xlabel('Frequency (Hz)')
ylabel('|C(f)|')

figure
subplot(3,1,1)
plot(f,iChannelRespDFT) 
title('I Channel Bandpass for Respiration Rate')
xlabel('Frequency (Hz)')
ylabel('|I(f)|')
subplot(3,1,2)
plot(f,qChannelRespDFT) 
title('Q Channel Bandpass for Respiration Rate')
xlabel('Frequency (Hz)')
ylabel('|Q(f)|')
subplot(3,1,3)
plot(f,combinedRespDFT) 
title('Combined Channels Bandpass for Respiration Rate')
xlabel('Frequency (Hz)')
ylabel('|C(f)|')

figure
subplot(3,1,1)
plot(f,iChannelHeartDFT) 
title('I Channel Bandpass for Heart Rate')
xlabel('Frequency (Hz)')
ylabel('|I(f)|')
subplot(3,1,2)
plot(f,qChannelHeartDFT) 
title('Q Channel Bandpass for Heart Rate')
xlabel('Frequency (Hz)')
ylabel('|Q(f)|')
subplot(3,1,3)
plot(f,combinedHeartDFT) 
title('Combined Channels Bandpass for Heart Rate')
xlabel('Frequency (Hz)')
ylabel('|C(f)|')

%% Print out heart and respiration rates
endMessage1 = ['Heart Rate is ' num2str(heartRate*60) ...
    ' beats per minute using the ' heartChoice];
endMessage2 = ['Respiration Rate is ' num2str(respirationRate*60) ...
    ' breaths per minute using the ' respChoice];
disp(endMessage1);
disp(endMessage2);