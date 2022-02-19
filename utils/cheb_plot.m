clc
clear
close all


figure(1);clf;
% subplot(1,2,1);
hold on;
% grid on;
% box on;
% f = chebfun('sqrt(abs(x-3))',[0,6],'splitting','on');
f = chebfun('x*cos(x)',[-3,3],'splitting','on');
[p1,err1] = minimax(f,4);
[p3,err3] = minimax(f,80);

a2 = plot(p1,'m--');hold on;
set(a2,'linewidth',3,'linestyle',':');
set(a2,'color',[247 77 77]/255);

a4 = plot(p3,'k-');hold on;
set(a4,'linewidth',3,'linestyle','-');
set(a4,'color',[35 145 213]/255);hold on;

a1 = plot(f,'k');hold on;
set(a1,'linewidth',3,'linestyle','-.');


% set(a1,'linewidth',2.4,'linestyle','-');
% set(a2,'linestyle','--');
% set(a4,'linewidth',1.0,'linestyle','-');
legend([a1,a2,a4],{'$y$','$4^{th}$','$80^{th}$'}, 'Interpreter','latex','location','best','Fontsize',24);
n=0:1:2;
set(gca,'ytick',n)
% set(gca,'ytick',[-2,0,2]);
set(gca,'xtick',[3,6]);
set(gca,'FontSize',16,'Fontname','Times');
set(gca,'Box','on');
ax = gca;
ax.BoxStyle = 'full';
ax.LineWidth = 1.7;
xlabel('$x$','Interpreter','latex','Fontsize',22);
ylabel('$\xi$','Interpreter','latex','Fontsize',22);

%%
% subplot(1,2,2); 
figure(2);clf
b2 = plot(f-p3,'k'); hold on
set(b2,'linewidth',3,'linestyle','-');
set(b2,'color',[35 145 213]/255);hold on;

a5 = plot([0 6],err3*[1 1],'--k');
a6 = plot([0 6],-err3*[1 1],'--k');
set(a5,'linewidth',1,'linestyle','-.');
% set(a5,'color',[35 145 213]/255);hold on;
set(a6,'linewidth',1,'linestyle','-.');
% set(a6,'color',[35 145 213]/255);hold on;

b3 = plot(f-p1,'m'); hold on
set(b3,'linewidth',3,'linestyle',':');
set(b3,'color',[247 77 77]/255);hold on;

a7 = plot([0 6],err1*[1 1],'--k');
a8 = plot([0 6],-err1*[1 1],'--k');
set(a7,'linewidth',1,'linestyle','-.');
% set(a7,'color',[247 77 77]/255);hold on;
set(a8,'linewidth',1,'linestyle','-.');
% set(a8,'color',[247 77 77]/255);hold on;
ylim(1.3*err1*[-1,1.2])

% set(b3,'linestyle','--');
% grid on;
legend([b3,b2],{'$4^{th}$','$80^{th}$'}, 'Interpreter','latex','location','best','Fontsize',24,'Orientation','horizon');


set(gca,'ytick',[-err1,-err3,0,err3,err1])
set(gca,'xtick',[0,3,6]);
set(gca,'Box','on');
ax = gca;
ax.BoxStyle = 'full';
ax.LineWidth = 1.7;
set(gca,'FontSize',16,'Fontname','Times');
xlabel('$x$','Interpreter','latex','Fontsize',22);
ylabel('$\xi$','Interpreter','latex','Fontsize',22);