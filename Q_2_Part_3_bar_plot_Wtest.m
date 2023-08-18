f = load('face.mat');
face_erp = f.face_data;

nf = load('nonface.mat');
n_face_erp = nf.nonface_data;

epc_f = [];
erp_f = [];
err_f = [];

epc_nf = [];
erp_nf = [];
err_nf = [];

t = linspace(-100, 1000, 282);

% PLOT ERP FOR FACE DATA FOR EACH CHANNEL
for j = 1:126
    erp_f(j, :) = mean(face_erp(j, :, :), 3);
    epc_f(j, :, :) = face_erp(j, :, :); % Uncomment this line to calculate error bars
    
    % Find timing and amplitude of approximately -170 ERP for face data
    [~, idx] = min(abs(erp_f(j, :) - (-170)));
    amp = erp_f(j, idx);
    
    % Store time and amplitude for each channel
    timing_f(j) = t(idx);
    amplitude_f(j) = amp;
    
    % Compute ERP for non-face data
    erp_nf(j, :) = mean(n_face_erp(j, :, :), 3);
    epc_nf(j, :, :) = n_face_erp(j, :, :); % Uncomment this line to calculate error bars
    
    % Find timing and amplitude of approximately -170 ERP for non-face data
    [~, idx] = min(abs(erp_nf(j, :) - (-170)));
    amp = erp_nf(j, idx);
    
    % Store time and amplitude for each channel
    timing_nf(j) = t(idx);
    amplitude_nf(j) = amp;
end

% Calculate confidence intervals for face and non-face conditions
conf_int_f = 1.96 * std(epc_f, 0, 3) / sqrt(size(face_erp, 3));
conf_int_nf = 1.96 * std(epc_nf, 0, 3) / sqrt(size(n_face_erp, 3));

conf_int_f = mean(conf_int_f, 2);
conf_int_nf = mean(conf_int_nf, 2);
% Plot timing and amplitude with confidence interval for face and non-face conditions
figure
errorbar(timing_f, amplitude_f, conf_int_f, 'b', 'LineWidth', 1.5)
hold on
errorbar(timing_nf, amplitude_nf, conf_int_nf, 'r', 'LineWidth', 1.5)
hold off
xlabel('Time (ms)')
ylabel('Amplitude')
title('Timing and Amplitude of Approximately -170 ERP')
legend('Face Condition', 'Non-Face Condition')
%%
figure;
hold on
% Plot the timing comparison with confidence intervals
errorbar(1:126, timing_f, conf_int_f, 'bo');
errorbar(1:126, timing_nf, conf_int_nf, 'ro');

% Customize the plot
xlim([0 127]);
xticks(1:126);
xticklabels(1:126);
xlabel('Channel');
ylabel('Timing (ms)');
title('N170 Timing Comparison');
legend('Face Data', 'Non-Face Data');

% Plot amplitude comparison between face and non-face with confidence intervals
figure;
hold on;

% Plot the amplitude comparison with confidence intervals
errorbar(1:126, amplitude_f, conf_int_f, 'bo');
errorbar(1:126, amplitude_nf, conf_int_nf, 'ro');

% Customize the plot
xlim([0 127]);
xticks(1:126);
xticklabels(1:126);
xlabel('Channel');
ylabel('Amplitude');
title('N170 Amplitude Comparison');
legend('Face Data', 'Non-Face Data');
%% WILCOXON STATISTICAL TEST

% Perform Wilcoxon signed-rank test between face and non-face conditions for amplitude
[a_p_value, a_h] = signrank(amplitude_f(1,:), amplitude_nf(1,:));
alpha = 0.05;
% Perform Wilcoxon signed-rank test between face and non-face conditions for timing
[t_p_value, t_h] = signrank(timing_f(1,:), timing_nf(1,:));

p_timing = t_p_value;
p_amplitude = a_p_value;

h_timing = t_h;
h_amplitude = a_h;

% Plot the results of Wilcoxon test for amplitude
figure;
bar(1, p_amplitude, 'b');
hold on;
bar(2, p_timing, 'r');
ylabel('p-value');
ylim([0 1e-20]);
xticks([1 2]);
xticklabels({'Amplitude', 'Timing'});
title('Wilcoxon Test: Face vs Non-face');
legend('p-value of amplitude', 'p-value of timing');

% Add significance indicator if the null hypothesis is rejected
if h_amplitude == 1
    text(1, p_amplitude, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
if h_timing == 1
    text(2, p_timing, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
%%
%% barplot -confidence interval of amplitude face vs nonface
m_nf_amp=mean(amplitude_nf, 2);% non-face
sem_nf_amp=std(amplitude_nf)/sqrt(126);
CI_nf_amp1=(1.96*sem_nf_amp);
CI_nf_amp2=(-1.96*sem_nf_amp);

m_f_amp=mean(amplitude_f, 2);% non-face
sem_f_amp=std(amplitude_f)/sqrt(126);
CI_f_amp1=(1.96*sem_f_amp);
CI_f_amp2=(-1.96*sem_f_amp);
figure
b={'amplitude-non-face95%CI','amplitude-face/95%CI'};%assigning the names x axis
mean_amp=[m_nf_amp,m_f_amp];
err_high=[CI_nf_amp1,CI_f_amp1];
err_low=[CI_nf_amp2,CI_f_amp2];
%sem_amp=[sem_nf_amp,sem_f_amp];
x=categorical(b);
bar(x,mean_amp);
hold on
er=errorbar(x,mean_amp,err_low,err_high,'k','LineWidth',1);
er.LineStyle='none';
ylim([-7 3]);
ylabel('voltage');
title('differences between amplitudes of face and non-face is significant ')
hold off
%%
%% confidence interval latency
m_nf_delay=mean(timing_nf,2);%non face
sem_nf_delay=std(timing_nf)/sqrt(126);
CI_nf_delay1=(1.96*sem_nf_delay);
CI_nf_delay2=-(1.96*sem_nf_delay);

m_f_delay=mean(timing_f,2);% face
sem_f_delay=std(timing_f)/sqrt(126);
CI_f_delay1=+(1.96*sem_f_delay);
CI_f_delay2=-(1.96*sem_f_delay);
figure
b={'latency-non-face95%CI','latency-face/95%CI'};%assigning the names x axis
mean_del=[m_nf_delay,m_f_delay];
err_high=[CI_nf_delay1,CI_f_delay1];
err_low=[CI_nf_delay2,CI_f_delay2];
%sem_del=[sem_nf_amp,sem_f_amp];
x=categorical(b);
bar(x,mean_del,'r');
hold on
er=errorbar(x,mean_del,err_low,err_high,'k','LineWidth',1);
er.LineStyle='none';
ylim([0 240]);
ylabel('time(ms)');
title('differences between latencies of face and non-face is not significant ')
hold off
