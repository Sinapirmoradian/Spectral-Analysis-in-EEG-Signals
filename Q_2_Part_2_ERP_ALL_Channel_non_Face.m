nf = load('nonface.mat');
n_face_erp = nf.nonface_data;
epc_nf=[];
erp_nf=[];
err_nf=[];
t=linspace(-100,1000,282);
for j=1:126
    erp_nf(j,:)=mean(n_face_erp(j,:,:),3);
    hold on
    plot(t,erp_nf(j,:))
    ax.YLim = [-5 7];
    ax.XLim=[-100 1000];
    line([0 0], ax.YLim, 'Color', 'k');
    line(ax.XLim, [0 0], 'Color', 'k');
end
xlabel('duration time per trial(ms)-latency');
ylabel('voltage');
title(sprintf('ERP of nonface condition of all channels'))
m_nf_er= mean(erp_nf,1);
plot(t,m_nf_er,'k','LineWidth',1.5)
ylim([-5 7]);
xlim([-120 1000]);
xlabel('duration time per trial(ms)-latency');
ylabel('voltage');
ax.YLim = [-5 7];
line([0 0], ax.YLim, 'Color', 'k');
legend('Mean non-face ERP')
 %% non-face condition and confidence interval of allchannels
m_nf_er= mean(erp_nf,1);
sem_nf_amp=1.96*std(erp_nf)/sqrt(126);
figure;
t=linspace(-100,1000,282);
errorbar(t,m_nf_er,sem_nf_amp,'m')
hold on
plot(t,m_nf_er,'k','LineWidth',1.5)
ylim([-7 6]);
xlim([-120 1000]);
xlabel('duration time per trial(ms)-latency');
ylabel('voltage');
title(sprintf('ERP of non-face condition of all chanels'))
ax.YLim = [-7 6];
ax.XLim=[-100 1000];
line([0 0], ax.YLim, 'Color', 'k');
line(ax.XLim, [0 0], 'Color', 'k');
legend('confidence_interval-%95')
%% PLOT MEAN OF SIGNAL
plot(t,m_nf_er,'k','LineWidth',1.5)
ylim([-0.5 0.5]);
title(sprintf('ERP of non-face condition of all chanels'))
