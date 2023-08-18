f = load('face.mat');
n_face_erp = f.face_data;

nf = load('nonface.mat');
n_nonface_erp = nf.nonface_data;

t = linspace(-100, 1000, 282);
alpha = 0.05; % Significance level

% Initialize arrays for p-values
p_values = zeros(1, 126);

% Perform Wilcoxon test for each channel
for j = 1:126
    [~, p_values(j)] = signrank(squeeze(n_face_erp(j, :, :)), squeeze(n_nonface_erp(j, :, :)));
end

% Plot p-values
figure
bar(p_values)
ax = gca;
ax.YLim = [0 1];
ax.XLim = [1 126];
xlabel('Channel');
ylabel('p-value');
title('Wilcoxon Test: Face vs Non-Face Conditions');

% Add significance asterisks
for j = 1:126
    if p_values(j) < alpha
        text(j, p_values(j), '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end
