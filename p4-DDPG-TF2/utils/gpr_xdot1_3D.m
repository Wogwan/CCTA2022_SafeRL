function [mean1,hyp1,delta,rmse] = gpr_xdot1_3D(x,y,xtest,ytest,it,noise,poly_deg)

tic
syms x1 x2 x3;
dxdt2 = 0;
X = [1];
%% Mean function
m1 = {@meanPoly,poly_deg}; 
hyp_m1 = zeros([3*poly_deg 1]);
m2 = {@meanConst}; hyp_m2 = 0;
meanfunc = {'meanSum',{m2,m1}}; hyp.mean = [hyp_m2; hyp_m1]; 
%%
% deg = 3;
% m1 = {@meanConst}; hyp_m1 = 0;
% m2 = {@meanLinear}; hyp_m2 = [0;0];
% msu = {'meanSum',{m1,m2}};  hypsu = [hyp_m1; hyp_m2];    % 1+2*x1+3*x2
% mpo = {'meanPow',deg,msu};       hyppo = hypsu;          % third power
% meanfunc = mpo; hyp.mean = hyppo; 
%% Cov function
% sf = 0.4; ell = 0.1; 
% cov1 = {@covSEiso}; hyp_cov1 = log([ell;sf/2]); 
% cov2 = {@covConst}; hyp_cov2 = 0.1; 
% covfunc = {'covSum',{cov1,cov2}}; hyp.cov = [hyp_cov1; hyp_cov2]; 
%%
sf = 0.1; ell = 0.1; 
cov1 = {@covSEiso}; hyp_cov1 = log([ell;sf/2]); 
covfunc = cov1; hyp.cov = hyp_cov1; 
%% Lik function
likfunc = @likGauss; sn = noise; hyp.lik = log(sn);
%% Inf function
inf_func = {'infGaussLik'};
%%
hyp1 = minimize(hyp, @gp, -it, inf_func, meanfunc, covfunc, likfunc, x, y);
%%
[ymu, ys2, fmu, fs2] = gp(hyp1, inf_func, meanfunc, covfunc, likfunc, x, y, xtest);    
mean1 = hyp1.mean;   
toc
%%
pvar x1 x2 x3
X1 = p2s(monomials(x1,(1:poly_deg)));
X2 = p2s(monomials(x2,(1:poly_deg)));
X3 = p2s(monomials(x3,(1:poly_deg)));
for i = 1:length(X1)
    X = [X;X1(i);X2(i);X3(i)];
end
for i = 1:length(X)
    dxdt2 = dxdt2 + X(i)*mean1(i);
end

delta = 0;
%%
% sigma = 2*real(sqrt(mean(ys2)));
% c = 0;
% x_error = abs(y - ymu);    
% for n = 1:length(y)
%     % check if the event occurs
%     if x_error(n,1) >= sigma
%         % find the number of occurrences
%         c = c + 1;
%     end
% end    
% delta = c/length(y)
% delta = 0;
%%
n = size(xtest,1);
figure(801);clf;
Output = [ymu+2*sqrt(ys2);flipdim(ymu-2*sqrt(ys2),1)];
error = abs(ytest-ymu);
rmse = sqrt(sum(error.^2)/n) 
subplot(211);hold on;
fill([(1:n)'; flipdim((1:n)',1)], Output, [0 7 0]/8, 'EdgeColor', [0 7 0]/8);
plot((1:n)',ymu,'k.','LineWidth',1);
plot((1:n)',ytest, 'r+', 'LineWidth',1);
xlabel('x'); ylabel('y');
legend('2$\sigma$','GP Predict','Real Data', 'Interpreter','latex','Orientation','horizon');
xlim([0 n]); ylim([min(ytest) max(ytest)])
time1 = toc
txt = ['GPML Time: ' num2str(time1) 's RMSE:' num2str(rmse)];
text(n/5,max(ytest)/3,txt)
set(gca,'xtick',[0,round(n/3),round(n/3*2),n]);
set(gca,'ytick',[min(ytest),max(ytest)]);
set(gca,'FontSize',18,'Fontname','Times');
set(gca,'Box','on');
ax = gca;
ax.BoxStyle = 'full';
ax.LineWidth = 1.2;
xlabel('$x_1$','Interpreter','latex','Fontsize',22,'Fontname','Times');
ylabel('$d_{\xi}(x)$','Interpreter','latex','Fontsize',22,'Fontname','Times','Rotation',0);

%%
subplot(212);hold on;
syms x1 x2 x3
dXtr_3 = [];
for num = 1:length(xtest(:,1))
    num_0 = vpa(subs(dxdt2,[x1,x2,x3],[xtest(num,1),xtest(num,2),xtest(num,3)]));
    dXtr_3 = [dXtr_3; num_0];
end
dXtr_3 = reshape(double(dXtr_3),1,[])';
Output_2 = [dXtr_3+2*sqrt(ys2);flipdim(dXtr_3-2*sqrt(ys2),1)];
error_2 = abs(ytest-dXtr_3);
rmse_poly = sqrt(sum(error_2.^2)/n) 
fill([(1:n)'; flipdim((1:n)',1)], Output_2, [0 7 0]/8, 'EdgeColor', [0 7 0]/8);
hold on
plot((1:n)',dXtr_3,'k.','LineWidth',1);
plot((1:n)',ytest, 'r+', 'LineWidth',1);
xlabel('x'); ylabel('y');
legend('2$\sigma$','GP Predict','Real Data','Interpreter','latex','Orientation','horizon');
xlim([0 n]); ylim([min(ytest) max(ytest)])
txt = ['Polynomial mean GPML Time: ' num2str(time1) 's RMSE:' num2str(rmse_poly)];
text(n/5,max(ytest)/3,txt)
set(gca,'xtick',[0,round(n/3),round(n/3*2),n]);
set(gca,'ytick',[min(ytest),max(ytest)]);
set(gca,'FontSize',18,'Fontname','Times');
set(gca,'Box','on');
ax = gca;
ax.BoxStyle = 'full';
ax.LineWidth = 1.2;
xlabel('$x$','Interpreter','latex','Fontsize',22,'Fontname','Times');
ylabel('$d_{\xi}(x)$','Interpreter','latex','Fontsize',22,'Fontname','Times','Rotation',0);
refreshdata; drawnow; 
end