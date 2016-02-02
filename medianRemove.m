function eegObject = medianRemove(eegObject)

% eegObject = medianRemove(eegObject)
%
% Computes the median of an EEG dataset (in eeglab format) then removes it
% from all of the channel signals. Standard processing step.

dataclump = eegObject';
med       = median(dataclump, 2);
dataClump = dataclump - repmat(med, 1, 14);
eegObject = dataclump';
end


