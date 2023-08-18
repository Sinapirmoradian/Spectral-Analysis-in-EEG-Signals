% Load EEG data for the face condition
f_wb = load('face_wb.mat');  % non-face without baseline
f_wb = f_wb.face_data_wn;

% Load EEG data for the non-face condition
nf_wb = load('nonface_wb.mat');  % non-face without baseline
nf_wb = nf_wb.nonface_data_wn;

% Set the sampling frequency
fs = 256;

% Define frequency bands of interest
frequency_bands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};
frequency_ranges = {[0.5, 4], [4, 8], [8, 13], [13, 30], [30, 100]}; % Specify the frequency ranges for each band

% Perform spectral analysis for each frequency band
psd_results_f = cell(size(frequency_bands, 2), 1);
psd_results_nf = cell(size(frequency_bands, 2), 1);
psd_results_f_mean = cell(size(frequency_bands, 2), 1);
psd_results_nf_mean = cell(size(frequency_bands, 2), 1);
psd_results_f_sem = cell(size(frequency_bands, 2), 1);
psd_results_nf_sem = cell(size(frequency_bands, 2), 1);
p_value = cell(size(frequency_bands, 2), 1);

for j = 1:size(frequency_bands, 2) % Loop over frequency bands
    % Initialize PSD results matrices
    psd_results_f{j} = zeros(size(f_wb, 1), size(f_wb, 3));
    psd_results_nf{j} = zeros(size(nf_wb, 1), size(nf_wb, 3));
    
    for i = 1:size(f_wb, 1) % Loop over channels
        % Calculate the PSD for the face condition using the Multitaper method
        [pxx_f, f_axis] = pmtm(squeeze(f_wb(i, :, :)), 3.5, size(f_wb, 2), fs, 'ConfidenceLevel', 0.95);
        freq_indices = f_axis >= frequency_ranges{j}(1) & f_axis <= frequency_ranges{j}(2);
        psd_results_f{i, j} = mean(pxx_f(freq_indices, :), 2);

        % Calculate the PSD for the non-face condition using the Multitaper method
        [pxx_nf, ~] = pmtm(squeeze(nf_wb(i, :, :)), 3.5, size(nf_wb, 2), fs, 'ConfidenceLevel', 0.95);
        psd_results_nf{i, j} = mean(pxx_nf(freq_indices, :), 2);
    end
    
    % Calculate mean and standard error of the mean (SEM) for face condition
    psd_results_f_mean{j} = mean(cat(2, psd_results_f{:, j}), 2);
    psd_results_f_sem{j} = std(cat(2, psd_results_f{:, j}), 0, 2) / sqrt(size(f_wb, 1));
    
    % Calculate mean and standard error of the mean (SEM) for non-face condition
    psd_results_nf_mean{j} = mean(cat(2, psd_results_nf{:, j}), 2);
    psd_results_nf_sem{j} = std(cat(2, psd_results_nf{:, j}), 0, 2) / sqrt(size(nf_wb, 1));
    
    % Plot PSD with error bars for the face condition
    figure;
    errorbar(f_axis(freq_indices), 10 * log10(psd_results_f_mean{j}), 10 * log10(psd_results_f_sem{j}));
    xlim([frequency_ranges{j}(1), frequency_ranges{j}(2)]);
    xlabel('Frequency (Hz)');
    ylabel('Power Spectral Density (dB/Hz)');
    title(sprintf('Face Condition - %s Band', frequency_bands{j}));
    
    % Plot PSD with error bars for the non-face condition
    figure;
    errorbar(f_axis(freq_indices), 10 * log10(psd_results_nf_mean{j}), 10 * log10(psd_results_nf_sem{j}));
    xlim([frequency_ranges{j}(1), frequency_ranges{j}(2)]);
    xlabel('Frequency (Hz)');
    ylabel('Power Spectral Density (dB/Hz)');
    title(sprintf('Non-Face Condition - %s Band', frequency_bands{j}));
    
    % Perform t-test between face and non-face conditions
    [~, p_value{j}] = ttest2(psd_results_f_mean{j}, psd_results_nf_mean{j});
    fprintf('T-Test Results for %s Band:\n', frequency_bands{j}, 'without baseline');
    fprintf('Face vs. Non-Face: p-value = %f\n\n', p_value{j}, 'without baseline');
end
for k = 1:numel(frequency_bands)
    p_values(k) = p_value{k,1};
end
% Plot the p-values for each frequency band
figure;
bar(p_values);
ylabel('P-value');
xticks(1:numel(frequency_bands));
xticklabels(frequency_bands);
title('Wilcoxon Test: Face vs Non-face without baseline');
% Plot a significant line at the 5% level
significant_line = 0.05; % Set the significance level
line([0, numel(frequency_bands)+1], [significant_line, significant_line], 'Color', 'r', 'LineStyle', '--');

hold off; % Disable hold on to finalize the plot