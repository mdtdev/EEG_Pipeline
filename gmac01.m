% Script for extracting power in each frequency range from EEG data block
% with rows being the time course and columns being the electrode axis.
% Have not checked to see if this works yet.
%
% From: https://emotiv.com/forum/forum4/topic18214/messages/
%   Posted by: GMAC Administrator @ emotiv who is the author.
%   Not listed as having any copyright, but note sure.
%
% Assumes: 128 Hz Sampling Freq; matrix has EEG data only; no data
% corrections!
%
% Transcribed with minor edits by MDT on 2016.01.18

fftlength  = 256;                      % make the window for sample length fftlength
hanning    = [1:fftlength]';
hanning_in = 2*pi*(hanning - (fftlength+1)/2)/(fftlength+1); % rescaled x-axis to match sample length
hanning    = (sin(hanning_in)./hanning_in).^2;               % sinc^2
hanning    = repmat(hanning, 1, 14);   % match to number of channels

f = [128/fftlength:128/fftlength:128]; % frequency index for the spectral array

thetaIndex    = find(f>=4 & f<8);
alphaIndex    = find(f>=8 & f<12);
lowBetaIndex  = find(f>=12 & f<16);
highBetaIndex = find(f>=16 & f<25);
gammaIndex    = find(f>=25 & f<40);
totIndex      = find(f>=4 & f<=40);
outdata       = [];

% med     = median(eeg.raw,2);             % remove median of each sample
% eeg.raw = eeg.raw - repmat(med, 1, 14);
% 
% for j=2:size(eeg.raw,1)                  % limit slew rate
%     del = eeg.raw(j,:) - eeg.raw(j-1,:);
%     del = min(del,  ones(1,14)*15);
%     del = max(del, -ones(1,14)*15);
%     eeg.raw(j,:) = eeg.raw(j-1,:) + del;
% end
% 
% % High pass filter section
% 
% a = 0.0078125;  % HPF filter coefs
% b = 0.9921875;
% 
% preVal   = zeros(1,14);
% eeg.filt = zeros(size(eeg.raw));
% 
% for j=2:size(eeg.raw,1)
%     preVal = a * eeg.raw(j,:) + b * preVal;
%     eeg.filt(j,:) = eeg.raw(j,:) - preVal;
% end                                          % end HPF

% Add more vectors to the "eeg" data structure!

% Create eeg.filt

eeg.theta    = [];
eeg.alpha    = [];
eeg.lowBeta  = [];
eeg.highBeta = [];
eeg.gamma    = [];
eeg.tot      = [];
eeg.totmed   = [];

% Main loop!

for k = fftlength:32:size(eeg.filt,1)                            % step through every quarter second starting at first possible sample
    spectrum     = fft(eeg.filt(k-fftlength+1:k,:) .* hanning);  % apply window to HP filtered data
    spectrum     = sqrt(spectrum .* conj(spectrum));             % get magnitude
    eeg.theta    = [eeg.theta; k sum(spectrum(thetaIndex,:))];   % append total spectral power in band, including sample index k
    eeg.alpha    = [eeg.alpha; k sum(spectrum(alphaIndex,:))];
    eeg.lowBeta  = [eeg.lowBeta; k sum(spectrum(lowBetaIndex,:))];
    eeg.highBeta = [eeg.highBeta; k sum(spectrum(highBetaIndex,:))];
    eeg.gamma    = [eeg.gamma; k sum(spectrum(gammaIndex,:))];
    eeg.tot      = [eeg.tot; k sum(spectrum(totIndex,:))];
end

% PLOT

timeRange = 30:2500;
channel   = 8;

figure;
plot(eeg.theta(timeRange,1), eeg.theta(timeRange,channel), 'b');
hold on
plot(eeg.alpha(timeRange,1), eeg.alpha(timeRange,channel), 'r');
title('AF3')   % Hopefully AF3 is 2 -- need to check
legend('theta','alpha');
xlabel('Sample Number (Start of FFT; 32 sample steps)')
ylabel('Amplitude (squared?)');
%ylim([0 2500]);

a2tr = eeg.alpha(:,channel)./eeg.theta(:,channel);

figure;
plot(eeg.alpha(timeRange,1), a2tr(timeRange), 'k');
legend('alpha to theta ratio');
xlabel('Sample Number (Start of FFT; 32 sample steps)')
ylabel('Ratio (Alpha to Theta)');

% Alpha versus Alpha
%
% There are 14 channels, from 2 to 15 (inclusive); channel 1 is not an
% electrode! It is a frame counter moving in steps of 32 samples each step.
% The 256 sample frequency analyses are overlapping.

timeRange = 240:2304;  % Chopped off the first and last minute of recording
L         = 3;         % L and R electrode numbers
R         = 14;
alphaL    = eeg.alpha(timeRange, L);
alphaR    = eeg.alpha(timeRange, R);
alphaDiff = alphaR - alphaL;

plot(eeg.alpha(timeRange,1), alphaDiff, 'k');
title('Difference alphaR - alphaL, AF4-AF3');
xlabel('Sample Number (Start of FFT; 32 sample steps)')
ylabel('Difference in Frontal Alpha');
refline(0,0);        % Add horizontal line to plot
