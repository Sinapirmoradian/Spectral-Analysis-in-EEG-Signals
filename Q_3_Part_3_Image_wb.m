% Load EEG data for the face condition
f_wb = load('face_wb.mat'); % non-face without baseline
f_wb = f_wb.face_data_wn;

% Load EEG data for the non-face condition
nf_wb = load('nonface_wb.mat'); % non-face without baseline
nf_wb = nf_wb.nonface_data_wn;

% Set the sampling frequency
fs = 256;

% Define frequency bands of interest
frequency_bands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};
frequency_ranges = {[0.5, 4], [4, 8], [8, 13], [13, 30], [30, 100]}; % Specify the frequency ranges for each band

% Perform spectral analysis for each frequency band
psd_results_f = cell(size(frequency_bands, 2),1);
psd_results_nf = cell(size(frequency_bands, 2),1);

for j = 1:size(frequency_bands, 2) % Loop over frequency band
    for i = 1:size(f_wb, 1) % Loop over channels
        % Calculate the PSD for the face condition using the Multitaper method
        [pxx_f, f_axis] = pmtm(squeeze(f_wb(i, :, :)), 3.5, size(f_wb, 2), fs, 'ConfidenceLevel', 0.95);
        freq_indices = f_axis >= frequency_ranges{j}(1) & f_axis <= frequency_ranges{j}(2);
        psd_results_f{j}(i, :) = mean(pxx_f(freq_indices, :), 2);

        % Calculate the PSD for the non-face condition using the Multitaper method
        [pxx_nf, ~] = pmtm(squeeze(nf_wb(i, :, :)), 3.5, size(nf_wb, 2), fs, 'ConfidenceLevel', 0.95);
        psd_results_nf{j}(i, :) = mean(pxx_nf(freq_indices, :), 2);
    end
    
    % Plotting PSD within frequency bands for the face condition
    subplot(size(frequency_bands, 2), 2, 2 * j - 1);
    imagesc(10 * log10(psd_results_f{j}));
    colorbar;
    xlabel('Channel');
    ylabel('Trial');
    title(sprintf('Face Condition - %s Band', frequency_bands{j}, 'without baseline'));
    
    % Plotting PSD within frequency bands for the non-face condition
    subplot(size(frequency_bands, 2), 2, 2 * j);
    imagesc(10 * log10(psd_results_nf{j}));
    colorbar;
    xlabel('Channel');
    ylabel('Trial');
    title(sprintf('Non-Face Condition - %s Band', frequency_bands{j}, 'without baseline'));
end
