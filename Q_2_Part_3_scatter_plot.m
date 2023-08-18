% Load the face and non-face ERP data
f = load('face.mat');
face_erp = f.face_data;
nf = load('nonface.mat');
n_face_erp = nf.nonface_data;

erp_f = [];
erp_nf = [];

t = linspace(-100,1000,282);
% Initialize arrays for storing N170 timing and amplitude
face_n170_timing = zeros(126, 1);
face_n170_amplitude = zeros(126, 1);
nonface_n170_timing = zeros(126, 1);
nonface_n170_amplitude = zeros(126, 1);

% Define the threshold for the N170 component
n170_threshold = -170;

% Find the N170 timing and amplitude for each channel in the face data
for j = 1:126
    erp_f(j, :) = mean(face_erp(j, :, :), 3);
    [~, face_n170_index] = min(abs(erp_f(j, :) - n170_threshold));
    face_n170_timing(j) = t(face_n170_index);
    face_n170_amplitude(j) = erp_f(j, face_n170_index);
end

% Find the N170 timing and amplitude for each channel in the non-face data
for j = 1:126
    erp_nf(j, :) = mean(n_face_erp(j, :, :), 3);
    [~, nonface_n170_index] = min(abs(erp_nf(j, :) - n170_threshold));
    nonface_n170_timing(j) = t(nonface_n170_index);
    nonface_n170_amplitude(j) = erp_nf(j, nonface_n170_index);
end

% Plot the N170 timing and amplitude for face data
figure;
scatter(face_n170_timing, face_n170_amplitude, 'filled', 'b');
hold on;

% Plot the N170 timing and amplitude for non-face data
scatter(nonface_n170_timing, nonface_n170_amplitude, 'filled', 'r');

% Customize the plot
xlim([-120 1000]);
ylim([-100 100]);
xlabel('Time (ms)');
ylabel('Amplitude');
title('N170 Timing and Amplitude Comparison');
legend('Face Data', 'Non-Face Data');
