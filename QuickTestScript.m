% Some quick playing around!
%
% Uses subject 1001 session 1 eeg data from Fall 2015
%
% Load it into eeglab first!!!


Fs   = EEG.srate;
data = EEG.data';

size(data)

chunkEarly.Fs = Fs;
chunkMid.Fs   = Fs;
chunkLate.Fs  = Fs;

chunkEarly.data = data(1:27219, :);
chunkMid.data   = data(27219:(2*27219), :);
chunkLate.data  = data((2*27219):end, :);
size(chunkLate.data)

AISearly = alphaImbalance(chunkEarly);
AISmid   = alphaImbalance(chunkMid);
AISlate  = alphaImbalance(chunkLate);

size(AISearly)

AIStable = table(AISearly', AISmid', AISlate', 'VariableNames', {'Early', 'Middle', 'Late'}, 'RowNames', {'AF4-AF3', 'F8-F7', 'F4-F3', 'FC6-FC5', 'T8-T7', 'P8-P7', 'O2-O1'})