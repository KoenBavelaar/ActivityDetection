clear all; clc; close all;

dataFolder = uigetdir(pwd, 'Select the folder with .dat files');
files = dir(fullfile(dataFolder, '*.dat'));

outputSpectrograms = {};
fileLabels = {};

X_all = [];

for k = 1:length(files)
    fprintf('Processing %s (%d of %d)\n', files(k).name, k, length(files));
    filepath = fullfile(dataFolder, files(k).name);
    
    % Parse metadata
    [person, activity, repetition] = parseFilename(files(k).name);
    
    % .dat uitpakken 
    fileID = fopen(filepath, 'r');
    dataArray = textscan(fileID, '%f');
    fclose(fileID);
    radarData = dataArray{1};
    fc = radarData(1);
    Tsweep = radarData(2) / 1000;
    NTS = radarData(3);
    Bw = radarData(4);
    Data = radarData(5:end);
    fs = NTS / Tsweep;
    nc = length(Data)/NTS*Tsweep / Tsweep;
    Data_time = reshape(Data, [NTS nc]);
    Data_time = hilbert(real(Data_time));
    record_length = length(Data)/NTS * Tsweep;

    % Split into 1.25-second chunks (5/4)
    if record_length >= 5 && record_length <= 10
        chunkDuration = 5/1;  % seconds
        totalSweepsPerSecond = 1 / Tsweep;
        sweepsPerChunk = round(chunkDuration * totalSweepsPerSecond);
        numChunks = floor(size(Data_time, 2) / sweepsPerChunk);

        for chunkIdx = 1:numChunks
            sweepStart = (chunkIdx - 1) * sweepsPerChunk + 1;
            sweepEnd = chunkIdx * sweepsPerChunk;
            if sweepEnd > size(Data_time, 2)
                break;  % safety check
            end

            thisChunk = Data_time(:, sweepStart:sweepEnd);

            % Range-Time FFT
            win = repmat(kaiser(NTS,10),1,size(thisChunk,2));
            tmp = fftshift(fft(thisChunk.*win),1);
            Data_rangetime = tmp(NTS/2+1:NTS,:) + flipud(tmp(1:NTS/2,:));
            Data_rangetime = Data_rangetime(4:end, :);

            % Find and sum relevant range bins
            row_sums = sum(abs(Data_rangetime), 2);
            [~, max_row] = max(row_sums);
            max_row = max(7, max_row);  
            delta = 7;
            Rbin_start = max_row - delta + 3;
            Rbin_stop = max_row + delta + 3;
            myvector = sum(Data_rangetime(Rbin_start:Rbin_stop, :));

            % Spectrogram
            PRF = 1 / Tsweep;
            WindowLength = 150;
            OverlapPercentage = 0.90;
            NFFTPoints = 2 * WindowLength;

            Data_spectrogram = fftshift(spectrogram(myvector, WindowLength, round(WindowLength*OverlapPercentage), NFFTPoints), 1);
            Data_spectrogram = flipud(abs(Data_spectrogram));
            Data_spectrogram = Data_spectrogram(:, 1:286);

            % Sped up spectrogram (15%)
            OverlapPercentage = 0.885;

            Data_spectrogram_sped = fftshift(spectrogram(myvector, WindowLength, round(WindowLength*OverlapPercentage), NFFTPoints), 1);
            Data_spectrogram_sped = flipud(abs(Data_spectrogram_sped));
            Data_spectrogram_sped = Data_spectrogram_sped(:, 1:286);

            
            % Slowed down spectrogram (15%)
            OverlapPercentage = 0.915;

            Data_spectrogram_slowed = fftshift(spectrogram(myvector, WindowLength, round(WindowLength*OverlapPercentage), NFFTPoints), 1);
            Data_spectrogram_slowed = flipud(abs(Data_spectrogram_slowed));
            Data_spectrogram_slowed = Data_spectrogram_slowed(:, 1:286);

            % Normalize
            specNorm = Data_spectrogram / max(Data_spectrogram(:));
            specSped = Data_spectrogram_sped / max(Data_spectrogram_sped(:));
            specSlowed = Data_spectrogram_slowed / max(Data_spectrogram_slowed(:));

            spec3D = cat(3, specNorm, specSped, specSlowed);  % [freq, time, type]

            size(spec3D)
            outputSpectrograms{end+1} = spec3D;
            fileLabels{end+1} = struct( ...
                'filename', files(k).name, ...
                'person', person, ...
                'activity', activity, ...
                'repetition', repetition, ...
                'chunkIndex', chunkIdx ...
            );
        end
    end
end

% Save results
save(fullfile(dataFolder, 'Augmented_Data.mat'), ...
     'outputSpectrograms', 'fileLabels', '-v7.3');

fprintf('Done!');


% ========== Helper: Parse filename ==========
function [person, activity, repetition] = parseFilename(filename)
    [~, name, ~] = fileparts(filename);
    person = str2double(regexp(name, 'P(\d+)', 'tokens', 'once'));
    activity = str2double(regexp(name, 'A(\d+)', 'tokens', 'once'));
    repetition = str2double(regexp(name, 'R(\d+)', 'tokens', 'once'));
end
