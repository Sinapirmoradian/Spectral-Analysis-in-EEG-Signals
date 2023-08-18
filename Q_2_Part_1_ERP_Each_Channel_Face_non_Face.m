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
%% PLOT ERP FOR FACE DATA FOR EACH CHANNEL
for j = 1:126
    erp_f(j, :) = mean(face_erp(j, :, :), 3);
    epc_f(j, :, :) = face_erp(j, :, :); % Uncomment this line to calculate error bars
end

% Plotting ERPs for each channel
figure;
for j = 1:126
    subplot(9, 14, j);
    plot(t, erp_f(j, :), 'k', 'LineWidth', 1.5);
    ylim([-5 7]);
    xlim([-120 1000]);
    ylabel('Voltage');
    title(sprintf('ERP-Ch%d', j));
    ax = gca;
    ax.YLim = [-5 7];
    line([0 0], ax.YLim, 'Color', 'k');
end
%% PLOT ERP FOR NON_FACE DATA FOR EACH CHANNEL
for j = 1:126
    erp_nf(j, :) = mean(n_face_erp(j, :, :), 3);
    epc_nf(j, :, :) = n_face_erp(j, :, :); % Uncomment this line to calculate error bars
end

% Plotting ERPs for each channel
figure;
for j = 1:126
    subplot(9, 14, j);
    plot(t, erp_nf(j, :), 'k', 'LineWidth', 1.5);
    ylim([-5 7]);
    xlim([-120 1000]);
    ylabel('Voltage');
    title(sprintf('ERP-Ch%d', j));
    ax = gca;
    ax.YLim = [-5 7];
    line([0 0], ax.YLim, 'Color', 'k');
end


