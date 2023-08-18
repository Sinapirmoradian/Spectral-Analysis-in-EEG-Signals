f = load('face.mat');
face_erp = f.face_data;
epc_f=[];
erp_f=[];
err_f=[];
t=linspace(-100,1000,282);
for j=1:126
    erp_f(j,:)=mean(face_erp(j,:,:),3);
    hold on
    plot(t,erp_f(j,:))
end
line([0 0], ylim, 'Color', 'k');
line(xlim, [0 0], 'Color', 'k');
xlabel('duration time per trial(ms)-latency');
ylabel('voltage');
title(sprintf('ERP of face condition of all channels'))
m_f_er= mean(erp_f,1);
plot(t,m_f_er,'k','LineWidth',1.5)
ylim([-5 7]);
xlabel('duration time per trial(ms)-latency');
ylabel('voltage');
line([0 0], ylim, 'Color', 'k');
line(xlim, [0 0], 'Color', 'k');
legend('Mean face ERP')
 %% face condition and confidence interval of allchannels
m_f_er= mean(erp_f,1);
sem_f_amp=1.96*std(erp_f)/sqrt(126);
figure
t=linspace(-100,1000,282);
errorbar(t,m_f_er,sem_f_amp,'m')
hold on
plot(t,m_f_er,'k','LineWidth',1.5)
ylim([-7 6]);
xlabel('duration time per trial(ms)-latency');
ylabel('voltage');
title(sprintf('ERP of face condition of all chanels'))
ylim([-1 1]);
xlim=[-100 1000];
line([0 0], ylim, 'Color', 'k');
line(xlim, [0 0], 'Color', 'k');
legend('confidence_interval-%95')