function AIS = alphaImbalance(eegDataBlob, duration, logbase)

% AIS = alphaImbalance(eegDataBlob, duration, logbase)
%
% Returns alpha imbalance scores a la Harmon-Jones (2006). Result is a
% matrix of alpha imblance scores. This function does no cleaning of the
% data such as removing eye blinks, etc. Also assumes the EEG channels are
% arranged in EPOC order: channel 1 and 14 are homologous sites, as are 2
% and 13, 3 and 12, and so on. AIS is arranged as number of differences
% (channels/2) by Harmon-Jones imbalance scores (1 column). DEFAULTS:
% duration is 2 (2 seconds) and logbase is set to 10.
%
% Expects the EEG time courses in GMAC order time X channels, not the
% natural order used by eeglab!
%
% Alpha power is INVERSELY related to cortical activity, so the larger the
% AIS score, the lower the power. HIGHER scores represent GREATER relative
% LEFT HEMISPHERIC activity.
%
% MDT
% 2016.02.02
% Version 0.0.1 ALPHA

    % SETUP and Defaults
    
    if nargin < 3
        logbase = 10;
    end

    if nargin < 2
        duration = 2;   % Default number of seconds per epoch; arbitrary
    end
    
    % Construct FFT Windows (Mostly from gmac01.m; see for details!)
    
    fftlength  = eegDataBlob.Fs * duration; % Based on duration variable!
    hanning    = [1:fftlength]';
    hanning_in = 2*pi*(hanning - (fftlength+1)/2)/(fftlength+1);
    hanning    = (sin(hanning_in)./hanning_in).^2;
    hanning    = repmat(hanning, 1, 14);
    f          = [128/fftlength:128/fftlength:128];
    
    % The following code selects the locations, in f, where the named
    % frequencies are located. To see the alpha frequencies, for instance,
    % you would use f(alphaIndex), etc. The use of < and <= is to make sure
    % the spectrum is not ever counted double.
    
    thetaIndex    = find(f >=  4 & f <   8);
    alphaIndex    = find(f >=  8 & f <= 12); % Note new bounds
    betaIndex     = find(f >  12 & f <  25); % No more high and low beta
    gammaIndex    = find(f >= 25 & f <= 40);
    totIndex      = find(f >=  4 & f <= 40);
    outdata       = [];
    
    % Remove the median at each timepoint from each of the timecourses;
    % check out the correlations both with and without doing this to see
    % the difference. It is substantial!
    
    med              = median(eegDataBlob.data, 2);
    eegDataBlob.data = eegDataBlob.data - repmat(med, 1, 14);
    
    % Limit the SLEW RATE and do a HIGH PASS filter on the data
    
    maxDelta = 15;  % Maximum change between two data points
    
    for j = 2:size(eegDataBlob.data,1)
        del = eegDataBlob.data(j,:) - eegDataBlob.data(j-1,:); % Get differences
        del = min(del,  ones(1,14)*maxDelta);                  % Diff or maxDelta (whichever smaller)
        del = max(del, -ones(1,14)*maxDelta);                  % Result of above or -maxDelta (whichever bigger)
        eegDataBlob.data(j,:) = eegDataBlob.data(j-1,:) + del; % Apply limits above to the data
    end
    
    a = 0.0078125; % HPF filter coefs (Magic Numbers)
    b = 0.9921875;
    
    preVal                 = zeros(1,14);
    eegDataBlob.filterData = zeros(size(eegDataBlob.data));
    
    for j = 2:size(eegDataBlob.data, 1)
        preVal = a * eegDataBlob.data(j, :) + b * preVal;
        eegDataBlob.filterData(j, :) = eegDataBlob.data(j, :) - preVal;
    end
    
    % Make locations to store frequency powers
    
    eegDataBlob.theta    = [];
    eegDataBlob.alpha    = [];
    eegDataBlob.beta     = [];
    eegDataBlob.gamma    = [];
    eegDataBlob.tot      = [];
    eegDataBlob.totmed   = [];
    
    % Main power computation
    %
    % Overlap is 25% of duration, based loosely on Harmon-Jones again
    
    overlapAmount = 0.25;
    stepSamples   = ceil((1 - overlapAmount)*fftlength);
        
    for k = fftlength:stepSamples:size(eegDataBlob.filterData,1)
        spectrum  = fft(eegDataBlob.filterData(k - fftlength + 1:k, :) .* hanning);
        spectrum  = sqrt(spectrum .* conj(spectrum)); 
        
        eegDataBlob.theta = [eegDataBlob.theta; k sum(spectrum(thetaIndex, :))];
        eegDataBlob.alpha = [eegDataBlob.alpha; k sum(spectrum(alphaIndex, :))];
        eegDataBlob.beta  = [eegDataBlob.beta;  k sum(spectrum(betaIndex, :))];
        eegDataBlob.gamma = [eegDataBlob.gamma; k sum(spectrum(gammaIndex, :))];
        eegDataBlob.tot   = [eegDataBlob.tot;   k sum(spectrum(totIndex, :))];
    end
    
    averageAlphaPowerByChannel = mean(eegDataBlob.alpha(:,2:15));
    
    % Harmon-Jones: Asymmetry indices defined as log right minus log left
    % alpha power. He does not state the log used, so we go with log10:
    
    leftAlpha  = averageAlphaPowerByChannel(1:7);
    rightAlpha = averageAlphaPowerByChannel(14:-1:8);  % Trick to reverse order of channels!
    
    AIS = log(rightAlpha)./log(logbase) - log(leftAlpha)./log(logbase);
end