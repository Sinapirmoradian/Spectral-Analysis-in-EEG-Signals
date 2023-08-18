
%% Part1
% new named variable
[ALLEEG, face, SET]    = pop_newset(ALLEEG, EEG, 1,'retrieve',1,'study',0);  
[ALLEEG, nonface, SET] = pop_newset(ALLEEG, EEG, 1,'retrieve',2,'study',0); 
fs = 256;
% lowFreq and highFreq for the delta EEG band it will be 0.5-4 Hz
lowFreq = 0.5; 
highFreq = 4; 

face_delta      = pop_eegfiltnew(face, 'locutoff',lowFreq,'hicutoff',highFreq,'plotfreqz',1);
nonface_delta   = pop_eegfiltnew(nonface, 'locutoff',lowFreq,'hicutoff',highFreq,'plotfreqz',1);

% lowFreq and highFreq for the theta EEG band it will be 4-8 Hz
lowFreq = 4; 
highFreq = 8; 

face_theta      = pop_eegfiltnew(face, 'locutoff',lowFreq,'hicutoff',highFreq,'plotfreqz',1);
nonface_theta   = pop_eegfiltnew(nonface, 'locutoff',lowFreq,'hicutoff',highFreq,'plotfreqz',1);

% lowFreq and highFreq for the alpha EEG band it will be 8-13 Hz
lowFreq = 8; 
highFreq = 13; 

face_alpha      = pop_eegfiltnew(face, 'locutoff',lowFreq,'hicutoff',highFreq,'plotfreqz',1);
nonface_alpha   = pop_eegfiltnew(nonface, 'locutoff',lowFreq,'hicutoff',highFreq,'plotfreqz',1);

% lowFreq and highFreq for the beta EEG band it will be 13-30 Hz
lowFreq = 13; 
highFreq = 30; 

face_beta      = pop_eegfiltnew(face, 'locutoff',lowFreq,'hicutoff',highFreq,'plotfreqz',1);
nonface_beta   = pop_eegfiltnew(nonface, 'locutoff',lowFreq,'hicutoff',highFreq,'plotfreqz',1);

% lowFreq and highFreq for the gamma EEG band it will be 30-100 Hz
lowFreq = 30; 
highFreq = 100; 

face_gamma      = pop_eegfiltnew(face, 'locutoff',lowFreq,'hicutoff',[],'plotfreqz',1);
nonface_gamma   = pop_eegfiltnew(nonface, 'locutoff',lowFreq,'hicutoff',[],'plotfreqz',1);

EEG_size = numel(face_gamma.data(:,1,1)); 

newName = 'delta'; 
for ii=1:EEG_size
    oldName = face.setname; 
    HilbertResults(1).setname = sprintf('%s_%s',oldName, newName);
    HilbertResults(1).hilbert(ii,:,:) = hilbert(face_delta.data(ii,:,:));
end

for ii=1:EEG_size
    oldName = nonface.setname; 
    HilbertResults(2).setname = sprintf('%s_%s',oldName, newName);
    HilbertResults(2).hilbert(ii,:,:) = hilbert(nonface_delta.data(ii,:,:));
end

newName = 'theta'; 
for ii=1:EEG_size
    oldName = face.setname; 
    HilbertResults(3).setname = sprintf('%s_%s',oldName, newName);
    HilbertResults(3).hilbert(ii,:,:) = hilbert(face_theta.data(ii,:,:));
end

for ii=1:EEG_size
    oldName = nonface.setname; 
    HilbertResults(4).setname = sprintf('%s_%s',oldName, newName);
    HilbertResults(4).hilbert(ii,:,:) = hilbert(nonface_theta.data(ii,:,:));
end

newName = 'alpha'; 
for ii=1:EEG_size
    oldName = face.setname; 
    HilbertResults(5).setname = sprintf('%s_%s',oldName, newName);
    HilbertResults(5).hilbert(ii,:,:) = hilbert(face_alpha.data(ii,:,:));
end

for ii=1:EEG_size
    oldName = nonface.setname; 
    HilbertResults(6).setname = sprintf('%s_%s',oldName, newName);
    HilbertResults(6).hilbert(ii,:,:) = hilbert(nonface_alpha.data(ii,:,:));
end

newName = 'beta'; 
for ii=1:EEG_size
    oldName = face.setname; 
    HilbertResults(7).setname = sprintf('%s_%s',oldName, newName);
    HilbertResults(7).hilbert(ii,:,:) = hilbert(face_beta.data(ii,:,:));
end

for ii=1:EEG_size
    oldName = nonface.setname; 
    HilbertResults(8).setname = sprintf('%s_%s',oldName, newName);
    HilbertResults(8).hilbert(ii,:,:) = hilbert(nonface_beta.data(ii,:,:));
end


newName = 'gamma'; 
for ii=1:EEG_size
    oldName = face.setname; 
    HilbertResults(9).setname = sprintf('%s_%s',oldName, newName);
    HilbertResults(9).hilbert(ii,:,:) = hilbert(face_gamma.data(ii,:,:));
end

for ii=1:EEG_size
    oldName = nonface.setname; 
    HilbertResults(10).setname = sprintf('%s_%s',oldName, newName);
    HilbertResults(10).hilbert(ii,:,:) = hilbert(nonface_gamma.data(ii,:,:));
end

% Mean of epochs for each electrode
for k=1:numel(HilbertResults)
    HilbertResults(k).meanEpochs = mean(HilbertResults(k).hilbert,3);
end
clearvars k

% Mean of channels 
for k=1:numel(HilbertResults)
    HilbertResults(k).meanEpochsTime = mean(HilbertResults(k).meanEpochs,1);
end
clearvars k
%%
figure();
% Plotting the subplots
subplot(3, 2, 1)
x = HilbertResults(1).meanEpochsTime;
y = rad2deg(angle(x));
F = linspace(0.5, 4, 282);
plot(F, y);
title('delta band-face');
xlabel('frequency (HZ)');
ylabel('phase (degree)');
grid on

subplot(3, 2, 2)
x = HilbertResults(3).meanEpochsTime;
y = rad2deg(angle(x));
F = linspace(4, 8, 282);
plot(F, y);
title('theta band-face');
xlabel('frequency (HZ)');
ylabel('phase (degree)');
grid on

subplot(3, 2, 3)
x = HilbertResults(5).meanEpochsTime;
y = rad2deg(angle(x));
F = linspace(8, 13, 282);
plot(F, y);
title('alpha band-face');
xlabel('frequency (HZ)');
ylabel('phase (degree)');
grid on

subplot(3, 2, 4)
x = HilbertResults(7).meanEpochsTime;
y = rad2deg(angle(x));
F = linspace(13, 30, 282);
plot(F, y);
title('beta band-face');
xlabel('frequency (HZ)');
ylabel('phase (degree)');
grid on

% Adding an empty subplot at the middle
subplot(3, 2, [5, 6])
axis off

% Plotting the last subplot at the middle
x = HilbertResults(9).meanEpochsTime;
y = rad2deg(angle(x));
F = linspace(30, 110, 282);
plot(F, y);
title('gamma band-face');
xlabel('frequency (HZ)');
ylabel('phase (degree)');
grid on

% Adjusting the figure size to fit the subplots
set(gcf, 'Position', [100, 100, 800, 800])

sgtitle('Phase Plots');
%%
figure();
% Plotting the subplots
subplot(3, 2, 1)
x = HilbertResults(2).meanEpochsTime;
y = rad2deg(angle(x));
F = linspace(0.5, 4, 282);
plot(F, y);
title('delta band-nonface');
xlabel('frequency (HZ)');
ylabel('phase (degree)');
grid on

subplot(3, 2, 2)
x = HilbertResults(4).meanEpochsTime;
y = rad2deg(angle(x));
F = linspace(4, 8, 282);
plot(F, y);
title('theta band-nonface');
xlabel('frequency (HZ)');
ylabel('phase (degree)');
grid on

subplot(3, 2, 3)
x = HilbertResults(6).meanEpochsTime;
y = rad2deg(angle(x));
F = linspace(8, 13, 282);
plot(F, y);
title('alpha band-nonface');
xlabel('frequency (HZ)');
ylabel('phase (degree)');
grid on

subplot(3, 2, 4)
x = HilbertResults(8).meanEpochsTime;
y = rad2deg(angle(x));
F = linspace(13, 30, 282);
plot(F, y);
title('beta band-nonface');
xlabel('frequency (HZ)');
ylabel('phase (degree)');
grid on

% Adding an empty subplot at the middle
subplot(3, 2, [5, 6])
axis off

% Plotting the last subplot at the middle
x = HilbertResults(10).meanEpochsTime;
y = rad2deg(angle(x));
F = linspace(30, 110, 282);
plot(F, y);
title('gamma band-nonface');
xlabel('frequency (HZ)');
ylabel('phase (degree)');
grid on

% Adjusting the figure size to fit the subplots
set(gcf, 'Position', [100, 100, 800, 800])

sgtitle('Phase Plots');
%% 
% Mean epochs for each electrode
for k=1:numel(HilbertResults)
    
    HilbertResults(k).meanChannels = mean(HilbertResults(k).hilbert,1);
    HilbertResults(k).meanChannels = squeeze(HilbertResults(k).meanChannels);
end
%% Part_2
for k=1:numel(HilbertResults)
    [time , ep]=size(HilbertResults(k).meanChannels);
    for i = 1:time
        ITPC(k , i) = round(1000*abs(mean(exp(1i*angle(HilbertResults(k).meanChannels(i,:))))))/1000;
    end
end

clearvars k

figure();

% Plotting the first four subplots
subplot(3, 2, 1)
F = linspace(0.5, 4, 282);
plot(F, ITPC(1, :));
title('delta band-face');
xlabel('frequency (HZ)');
ylabel('amplitude');
grid on

subplot(3, 2, 2)
F = linspace(4, 8, 282);
plot(F, ITPC(3, :));
title('theta band-face');
xlabel('frequency (HZ)');
ylabel('amplitude');
grid on

subplot(3, 2, 3)
F = linspace(8, 13, 282);
plot(F, ITPC(5, :));
title('alpha band-face');
xlabel('frequency (HZ)');
ylabel('amplitude');
grid on

subplot(3, 2, 4)
F = linspace(13, 30, 282);
plot(F, ITPC(7, :));
title('beta band-face');
xlabel('frequency (HZ)');
ylabel('amplitude');
grid on

% Adding an empty subplot at the middle
subplot(3, 2, [5, 6])
axis off

% Plotting the last subplot at the middle
F = linspace(30, 110, 282);
plot(F, ITPC(9, :));
title('gamma band-face');
xlabel('frequency (HZ)');
ylabel('amplitude');
grid on

% Adjusting the figure size to fit the subplots
set(gcf, 'Position', [100, 100, 800, 800])

sgtitle('ITPC Plots');
%%
figure();

% Plotting the first four subplots
subplot(3, 2, 1)
F = linspace(0.5, 4, 282);
plot(F, ITPC(2, :));
title('delta band-non face');
xlabel('frequency (HZ)');
ylabel('amplitude');
grid on

subplot(3, 2, 2)
F = linspace(4, 8, 282);
plot(F, ITPC(4, :));
title('theta band-non face');
xlabel('frequency (HZ)');
ylabel('amplitude');
grid on

subplot(3, 2, 3)
F = linspace(8, 13, 282);
plot(F, ITPC(6, :));
title('alpha band-non face');
xlabel('frequency (HZ)');
ylabel('amplitude');
grid on

subplot(3, 2, 4)
F = linspace(13, 30, 282);
plot(F, ITPC(8, :));
title('beta band-non face');
xlabel('frequency (HZ)');
ylabel('amplitude');
grid on

% Adding an empty subplot at the middle
subplot(3, 2, [5, 6])
axis off

% Plotting the last subplot at the middle
F = linspace(30, 110, 282);
plot(F, ITPC(10, :));
title('gamma band-non face');
xlabel('frequency (HZ)');
ylabel('amplitude');
grid on

% Adjusting the figure size to fit the subplots
set(gcf, 'Position', [100, 100, 800, 800])

sgtitle('ITPC Plots');

%%
[bootstat_delta,bootsam_delta] = bootstrp(1000,@corr,ITPC(1,:),ITPC(2,:));
[bootstat_theta,bootsam_theta] = bootstrp(1000,@corr,ITPC(3,:),ITPC(4,:));
[bootstat_alpha,bootsam_alpha] = bootstrp(1000,@corr,ITPC(5,:),ITPC(6,:));
[bootstat_beta,bootsam_beta] = bootstrp(1000,@corr,ITPC(7,:),ITPC(8,:));
[bootstat_gamma,bootsam_gammma] = bootstrp(1000,@corr,ITPC(9,:),ITPC(10,:));

bootstat_delta_face = bootstrp(200,@mean,ITPC(1,:));
bootstat_delta_nonface = bootstrp(200,@mean,ITPC(2,:));
bootstat_theta_face = bootstrp(200,@mean,ITPC(3,:));
bootstat_theta_nonface = bootstrp(200,@mean,ITPC(4,:));
bootstat_alpha_face = bootstrp(200,@mean,ITPC(5,:));
bootstat_alpha_nonface = bootstrp(200,@mean,ITPC(6,:));
bootstat_beta_face = bootstrp(200,@mean,ITPC(7,:));
bootstat_beta_nonface = bootstrp(200,@mean,ITPC(8,:));
bootstat_gamma_face = bootstrp(200,@mean,ITPC(9,:));
bootstat_gamma_nonface = bootstrp(200,@mean,ITPC(10,:));

figure();

subplot(3, 2, 1)
histogram(bootstat_delta_face, 'FaceColor', 'blue');
hold on
histogram(bootstat_delta_nonface, 'FaceColor', 'red');
title(sprintf(['delta band mean face = ' num2str(mean(bootstat_delta_face)), '\ndelta band mean nonface = ' num2str(mean(bootstat_delta_nonface))]));
grid on

subplot(3, 2, 2)
histogram(bootstat_theta_face, 'FaceColor', 'blue');
hold on
histogram(bootstat_theta_nonface, 'FaceColor', 'red');
title(sprintf(['theta band mean face = ' num2str(mean(bootstat_theta_face)), '\ntheta band mean nonface = ' num2str(mean(bootstat_theta_nonface))]));
grid on

subplot(3, 2, 3)
histogram(bootstat_alpha_face, 'FaceColor', 'blue');
hold on
histogram(bootstat_alpha_nonface, 'FaceColor', 'red');
title(sprintf(['alpha band mean face = ' num2str(mean(bootstat_alpha_face)), '\nalpha band mean nonface = ' num2str(mean(bootstat_alpha_nonface))]));
grid on

subplot(3, 2, 4)
histogram(bootstat_beta_face, 'FaceColor', 'blue');
hold on
histogram(bootstat_beta_nonface, 'FaceColor', 'red');
title(sprintf(['beta band mean face = ' num2str(mean(bootstat_beta_face)), '\nbeta band mean nonface = ' num2str(mean(bootstat_beta_nonface))]));
grid on

% Adding an empty subplot at the middle
subplot(3, 2, [5, 6])
axis off

% Plotting the last subplot at the middle
histogram(bootstat_gamma_face, 'FaceColor', 'blue');
hold on
histogram(bootstat_gamma_nonface, 'FaceColor', 'red');
% title({'gamma band mean face = ' num2str(mean(bootstat_gamma_face)), 'gamma band mean nonface = ' num2str(mean(bootstat_gamma_nonface))});
title(['gamma band mean face = ' num2str(mean(bootstat_gamma_face)), '     gamma band mean nonface = ' num2str(mean(bootstat_gamma_nonface))]);
grid on
% Adjusting the figure size to fit the subplots
set(gcf, 'Position', [100, 100, 800, 800])

sgtitle('Histograms');
%% Q3 
[ALLEEG data SET] = pop_newset(ALLEEG, EEG, 1,'retrieve',3,'study',0); 

%%
% lowFreq and highFreq for the delta EEG band it will be 0.5-4 Hz
lowFreq = 0.5; 
highFreq = 4; 
data_delta  = pop_eegfiltnew(data, 'locutoff',lowFreq,'hicutoff',highFreq,'plotfreqz',1);

% lowFreq and highFreq for the theta EEG band it will be 4-8 Hz
lowFreq = 4;
highFreq = 8;

data_theta  = pop_eegfiltnew(data, 'locutoff',lowFreq,'hicutoff',highFreq,'plotfreqz',1);

% lowFreq and highFreq for the alpha EEG band it will be 8-13 Hz
lowFreq = 8; 
highFreq = 13; 

data_alpha = pop_eegfiltnew(data, 'locutoff',lowFreq,'hicutoff',highFreq,'plotfreqz',1);

% lowFreq and highFreq for the beta EEG band it will be 13-30 Hz
lowFreq = 13; 
highFreq = 30; 

data_beta = pop_eegfiltnew(data, 'locutoff',lowFreq,'hicutoff',highFreq,'plotfreqz',1);

% lowFreq and highFreq for the gamma EEG band it will be 30-100 Hz
lowFreq = 30; 
highFreq = 100; 

data_gamma = pop_eegfiltnew(data, 'locutoff',lowFreq,'hicutoff',[],'plotfreqz',1);
EEG_size = numel(data.data(:,1,1)); 

newName = 'delta'; 
for ii=1:EEG_size
    oldName = data.setname; 
    dataRes(1).setname = sprintf('%s_%s',oldName, newName);
    dataRes(1).hilbert(ii,:,:) = hilbert(data_delta.data(ii,:,:));
end

newName = 'theta'; 
for ii=1:EEG_size
    oldName = data.setname; 
    dataRes(2).setname = sprintf('%s_%s',oldName, newName);
    dataRes(2).hilbert(ii,:,:) = hilbert(data_theta.data(ii,:,:));
end

newName = 'alpha'; 
for ii=1:EEG_size
    oldName = data.setname; 
    dataRes(3).setname = sprintf('%s_%s',oldName, newName);
    dataRes(3).hilbert(ii,:,:) = hilbert(data_alpha.data(ii,:,:));
end

newName = 'beta'; 
for ii=1:EEG_size
    oldName = data.setname; 
    dataRes(4).setname = sprintf('%s_%s',oldName, newName);
    dataRes(4).hilbert(ii,:,:) = hilbert(data_beta.data(ii,:,:));
end

newName = 'gamma'; 
for ii=1:EEG_size
    oldName = data.setname; 
    dataRes(5).setname = sprintf('%s_%s',oldName, newName);
    dataRes(5).hilbert(ii,:,:) = hilbert(data_gamma.data(ii,:,:));
end

%%
% Mean of channels 
for k=1:numel(dataRes)
    dataRes(k).meanChannels = mean(dataRes(k).hilbert,1);
    dataRes(k).meanChannels = squeeze(dataRes(k).meanChannels);
end
clearvars k

%%
for k = 1:100 
    for k=1:numel(dataRes)
        dataRes(k).meanChannels = transpose(dataRes(k).meanChannels);
        dataRes(k).meanChannels = dataRes(k).meanChannels(randperm(size(dataRes(k).meanChannels , 1)) , :);
        dataRes(k).meanChannels = transpose(dataRes(k).meanChannels);
    end
    clearvars k

    for k=1:numel(dataRes)
        [time , ep]=size(dataRes(k).meanChannels);
        for i = 1:time
            ITPC_data(k,k , i) = round(1000*abs(mean(exp(1i*angle(dataRes(k).meanChannels(i,:))))))/1000;
        end
    end
    clearvars k
end
ITPC_data_sh = mean(ITPC_data,1);
ITPC_data_sh = squeeze(ITPC_data_sh);
figure();
% Plotting the subplots
subplot(3, 2, 1)
F = linspace(0.5, 4, 282);
plot(F, ITPC_data_sh(1, :));
title('delta band');
xlabel('frequency (HZ)');
ylabel('amplitude');
grid on

subplot(3, 2, 2)
F = linspace(4, 8, 282);
plot(F, ITPC_data_sh(2, :));
title('theta band');
xlabel('frequency (HZ)');
ylabel('amplitude');
grid on

subplot(3, 2, 3)
F = linspace(8, 13, 282);
plot(F, ITPC_data_sh(3, :));
title('alpha band');
xlabel('frequency (HZ)');
ylabel('amplitude');
grid on

subplot(3, 2, 4)
F = linspace(13, 30, 282);
plot(F, ITPC_data_sh(4, :));
title('beta band');
xlabel('frequency (HZ)');
ylabel('amplitude');
grid on

% Adding an empty subplot at the middle
subplot(3, 2, [5, 6])
axis off

% Plotting the last subplot at the middle
F = linspace(30, 110, 282);
plot(F, ITPC_data_sh(5, :));
title('gamma band');
xlabel('frequency (HZ)');
ylabel('amplitude');
grid on

% Adjusting the figure size to fit the subplots
set(gcf, 'Position', [100, 100, 800, 800])

sgtitle('ITPC Plots');