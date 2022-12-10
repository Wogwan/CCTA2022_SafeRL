function [mean1,hyp1,delta,rmse] = gpr_xdot2(x,y,xtest,ytest,it,noise,poly_deg)

%%SOS_GP_CHECK_XDOT2 generates the xdot2 value [MEAN] to construct the polynomial function
% In:
%     Xtr_0      double   400  x  2   Input training X
%     dXtr_0     double   400  x  2   Input training Y
% Out:
%     mean1      double   6  x  1    Learned kernel with polynomial m
% Copyright (c) by Huang Hejun (CUHK) under BSD License
% Last modified: Huang Hejun 2021-05
tic
syms x1 x2;
dxdt2 = 0;
X = [1];
%% Mean function
m1 = {@meanPoly,poly_deg};
hyp_m1 = zeros([2*poly_deg 1]);
m2 = {@meanConst}; 
hyp_m2 = 0.11;
% hyp_m2 = 0;
meanfunc = {'meanSum',{m2,m1}}; hyp.mean = [hyp_m2; hyp_m1];
%%
% sf = 0.1; ell = 0.2;
sf = 0.2; ell = 0.4;
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
pvar x1 x2
X1 = p2s(monomials(x1,(1:poly_deg)));
X2 = p2s(monomials(x2,(1:poly_deg)));
for i = 1:length(X1)
    X = [X; X1(i); X2(i)];
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
m = size(x,1);
n = size(xtest,1);
%% Gathered
figure(801);clf;hold on;
% set(gca,'Position',[300 300 1000 500]);
% Output = [ymu+2*sqrt(ys2);flipdim(ymu-2*sqrt(ys2),1)];
error = abs(ytest-ymu);
rmse = sqrt(sum(error.^2)/n)
% subplot(211);
% fill([(1:n)'; flipdim((1:n)',1)], Output, [0 7 0]/8, 'EdgeColor', [0 7 0]/8);
le1 = plot((1:n)',ymu,'k*','LineWidth',0.8);
le2 = plot((1:n)',ytest, 'r+', 'LineWidth',0.7);
% legend('2$\sigma$','GP Predict','Real Data', 'Interpreter','latex','Orientation','horizon');
xlim([0 length(x)]); ylim([0.098 0.152])
time1 = toc
% txt = ['GPML Time: ' num2str(time1) 's RMSE:' num2str(rmse)];
% text(40,0.1,txt)

%%
hold on;
syms x1 x2
dXtr_3 = [];
for num = 1:length(xtest(:,1))
    num_0 = vpa(subs(dxdt2,[x1,x2],[xtest(num,1),xtest(num,2)]));
    dXtr_3 = [dXtr_3; num_0];
end
dXtr_3 = reshape(double(dXtr_3),1,[])';
error_2 = abs(ytest-dXtr_3);
rmse_poly = sqrt(sum(error_2.^2)/n)
le3 = plot((1:n)',dXtr_3,'bo','LineWidth',0.7);

xlim([0 length(x)]); ylim([0.098 0.152])
set(gca,'xtick',[0,100,200,300]);
set(gca,'ytick',[0.10,0.11,0.12,0.13,0.14,0.15]);
set(gca,'FontSize',24,'Fontname','Times');
set(gca,'Box','on');
ax = gca;
ax.BoxStyle = 'full';
ax.LineWidth = 1.2;
xlabel('$x_1$','Interpreter','latex','Fontsize',24,'Fontname','Times');
ylabel('$d_{\xi}$','Interpreter','latex','Fontsize',24,'Fontname','Times');
legend([le2;le1;le3],{'Data','SE Kernel','Polynomial'},'Interpreter','latex','Orientation','horizon');

end