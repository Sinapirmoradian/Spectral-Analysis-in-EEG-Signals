% Load EEG data for the face condition
f = load('face.mat');
f = f.face_data;
 
% Load EEG data for the non-face condition
nf = load('nonface.mat');
nf = nf.nonface_data;

% Set the sampling frequency
fs = 256;

% Multitaper method - Face Condition
pxx_f = [];
for i = 1:size(f, 1)
    for j = 1:size(f, 3)
        pxx_f(i, j, :) = pmtm(f(i, :, j), 3.5, size(f, 2), fs, 'ConfidenceLevel', 0.95); % PSD of 126 channels
    end
end

m_psd_f = squeeze(mean(pxx_f, 2));
err_psd_f = squeeze(std(pxx_f, 0, 2) / sqrt(size(f, 3)));

% Multitaper method - Non-Face Condition
pxx_nf = [];
for i = 1:size(nf, 1)
    for j = 1:size(nf, 3)
        pxx_nf(i, j, :) = pmtm(nf(i, :, j), 3.5, size(nf, 2), fs, 'ConfidenceLevel', 0.95); % PSD of 126 channels
    end
end

m_psd_nf = squeeze(mean(pxx_nf, 2));
err_psd_nf = squeeze(std(pxx_nf, 0, 2) / sqrt(size(nf, 3)));

% Perform t-test
alpha = 0.05; % Significance level

[h, p] = ttest2(m_psd_f, m_psd_nf, 'Alpha', alpha, 'Vartype', 'unequal');

% Display t-test results
fprintf('T-Test Results:\n');
fprintf('---------------------------------------------------\n');
fprintf('Null Hypothesis: The means are equal.\n');
fprintf('Alternative Hypothesis: The means are not equal.\n');
fprintf('Significance Level (alpha): %.2f\n', alpha);
fprintf('P-Value: %.6f\n', p);
fprintf('Hypothesis Test Result: ');
if h
    fprintf('Reject the null hypothesis.\n');
else
    fprintf('Fail to reject the null hypothesis.\n');
end
fprintf('---------------------------------------------------\n');

% Plotting p-value
figure;
bar(p);
xlabel('each Frequency point');
ylabel('P-Value');
title('P-Value for each Frequency point of face  and non-face condition with baseline normalization');
% Add legend indicating significance level
hold on;
sig_level = 0.05;
line([0, numel(p) + 1], [sig_level, sig_level], 'Color', 'red', 'LineStyle', '--');
legend('P-Value', 'Significance Level (0.05)');

% Plotting - Face Condition
figure;
subplot(1, 2, 1);
errorbar(10 * log10(m_psd_f'), 10 * log10(err_psd_f'));
xlabel('Hz');
ylabel('dB');
title('Error bar of Multitaper PSD Estimate for all Frequency points of face condition with baseline normalization');

subplot(1, 2, 2);
plot(10 * log10(m_psd_f'));
xlim([0 140]);
xlabel('Hz');
ylabel('dB');
title('Multitaper PSD Estimate for all Frequency points of face condition with baseline normalization');

% Plotting - Non-Face Condition
figure;
subplot(1, 2, 1);
errorbar(10 * log10(m_psd_nf'), 10 * log10(err_psd_nf'));
xlabel('Hz');
ylabel('dB');
title('Error bar of Multitaper PSD Estimate for all Frequency points of non-face condition with baseline normalization');

subplot(1, 2, 2);
plot(10 * log10(m_psd_nf'));
xlim([0 140]);
xlabel('Hz');
ylabel('dB');
title('Multitaper PSD Estimate for all Frequency points of non-face condition with baseline normalization');

% Comparison - Face vs. Non-Face Conditions
m_psd_f_er = mean(m_psd_f, 1);
sem_f_amp = 1.96 * std(m_psd_f) / sqrt(size(f, 1));

m_psd_nf_er = mean(m_psd_nf, 1);
sem_nf_amp = 1.96 * std(m_psd_nf) / sqrt(size(nf, 1));

figure;
subplot(1, 2, 1);
errorbar(10 * log10(m_psd_f_er), sem_f_amp, 'm');
hold on;
plot(10 * log10(m_psd_f_er), 'b', 'LineWidth', 1.5);
ylim([-100 25]);
xlabel('Hz');
ylabel('dB');
title('PSD of face condition of all Frequency points with baseline normalization');
legend('confidence_interval-%95');

subplot(1, 2, 2);
errorbar(10 * log10(m_psd_nf_er), sem_nf_amp, 'm');
hold on;
plot(10 * log10(m_psd_nf_er), 'k', 'LineWidth', 1.5);
ylim([-100 25]);
xlabel('Hz');
ylabel('dB');
title('PSD of non-face condition of all Frequency points with baseline normalization');
legend('confidence_interval-%95');
