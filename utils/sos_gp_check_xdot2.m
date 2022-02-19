function [mean1,hyp1,delta] = gp_xdot2(Xtr_0,dXtr_0,it,noise)
%%SOS_GP_CHECK_XDOT2 generates the xdot2 value [MEAN] to construct the polynomial function
% In:
%     Xtr_0      double   400  x  2   Input training X
%     dXtr_0     double   400  x  2   Input training Y
% Out:
%     mean1      double   6  x  1    Learned kernel with polynomial m
% Copyright (c) by Huang Hejun (CUHK) under BSD License 
% Last modified: Huang Hejun 2021-05
    syms x1 x2;
    % Kernel design
%     m1 = {@meanPoly,3}; hyp_m1 = [0;0;0;0;0;0];
    m1 = {@meanPoly,2}; hyp_m1 = [0;0;0;0];
    m2 = {@meanConst}; hyp_m2 = 0;
    meanfunc = {'meanSum',{m1,m2}}; hyp.mean = [hyp_m1; hyp_m2]; 
%     meanfunc = {@meanConst}; hyp.mean = 0.1;
    % cp1 = {@covPoly,1}; cp2 = {@covPoly,2}; cp3 = {@covPoly,3}; sf = 0.1; c = 0;
    % hycp1 = log([c,sf]); hycp2 = log([c,sf]); hycp3 = log([c,sf]);
    % covfunc = {'covSum', {cp1,cp2,cp3}}; hyp.cov = [hycp1; hycp2; hycp3];
    %%
    sf = 0.1; ell = 0.2; 
    n = 5; D = 3; x = randn(n,D); xs = randn(3,D);  % create a data set
    al = 2; L = rand(D,1);
    covfunc = {@covSEiso}; hyp.cov = log([ell;sf/2]); 
%     covfunc = {'covRQard'}; hyp.cov = log([L;sf]); 
%     covfunc = {'covMaterniso',3}; hyp.cov = log([ell;sf]); % Matern class d=3
    
    %%
    a = 1; sn = 0; nu = 4;
    likfunc = @likGauss; sn = 0.01; hyp.lik = log(sn);

    
    %%
%     inf_list = {'infGaussLik','infLaplace','infEP','infVB','infKL'};% inference algs
    inf_func = {'infGaussLik'};
    
    %%
    x = Xtr_0;
    y2 = dXtr_0;
    hyp1 = minimize(hyp, @gp, -it, inf_func, meanfunc, covfunc, likfunc, x, y2);
    
   %%
    [ymu ys2 fmu fs2] = gp(hyp1, inf_func, meanfunc, covfunc, likfunc, x, y2, x);    
    mean1 = hyp1.mean;    
%     dxdt2 = mean1(1)*x1+mean1(2)*x2+mean1(3)*x1^2+mean1(4)*x2^2+mean1(5)*x1^3+mean1(6)*x2^3+mean1(7);   % GP-posterior result for polynomial dynamic systems
    dxdt2 = mean1(1)*x1+mean1(2)*x2+mean1(3)*x1^2+mean1(4)*x2^2+mean1(5);   % GP-posterior result for polynomial dynamic systems

%     [nlZ dnlZ          ] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, y2);
%     mean2 = dnlZ.mean;    
%     dxdt3 = mean2(1)*x1+mean2(2)*x2+mean2(3)*x1^2+mean2(4)*x2^2+mean2(5);   % GP-posterior result for polynomial dynamic systems
  
    %%
    %     delta = 0.02;
    %     lh = zeros(size(x,2)+2,size(y2,2));
    %     for i=1:size(x,2)
    %         lh(i) = log(ell);
    %     end
    %     lh(end,1) = log(noise);
    %     lh(end-1,1) = log(sf);
    %     Kxx = gprt(lh, x); 
    %     it = 1;
    %     beta_t = 2*func^2+300*max(max(gammaT))*(log(it/delta)^3);
    %     eta = sqrt(beta_t)*noise;

%     delta_0 = 0.02;
%     
%     X = y2;
%     nsq=sum(X.^2,2);
%     K=bsxfun(@minus,nsq,(2*X)*X.');
%     K=bsxfun(@plus,nsq.',K);
%     K=exp(-K);
%     I = eye(size(K));
%     gammaT = 0.5*log(I+K*noise^(-2));
    
%     d = size(y2,2);
%     k = size(y2,1);
%     norm_f = 2*norm(y2)^2;
%     error = abs(max(y2 - ymu));
%     sigma = 2*real(sqrt(mean(ys2)));
%     delta = real(compute_delta(d,k,norm_f,error,sigma,gammaT))
%     beta_sqrt_sigma = 0;

    sigma = 2*real(sqrt(mean(ys2)));
    c = 0;
    x_error = abs(y2 - ymu);    
    for n = 1:length(y2)
        % check if the event occurs
        if x_error(n,1) >= sigma
            % find the number of occurrences
            c = c + 1;
        end
    end    
    delta = c/length(y2)

    
    %%
%     figure;clf
    plot3(x(:,1),x(:,2),y2(:,1),'ro'); hold on;
    
    dXtr_3 = [];
    for num = 1:length(Xtr_0(:,1))
        num_0 = vpa(subs(dxdt2,[x1,x2],[Xtr_0(num,1),Xtr_0(num,2)]));
        dXtr_3 = [dXtr_3; num_0];
    end
    dXtr_3 = reshape(double(dXtr_3),1,[])';
    plot3(Xtr_0(:,1),Xtr_0(:,2),dXtr_3(:,1),'b*'); hold on;
    
    plot3(Xtr_0(:,1),Xtr_0(:,2),ymu(:,1),'g^'); hold on;
%     plot3(Xtr_0(:,1),Xtr_0(:,2),ymu(:,1)+fs2(:,1),'co'); hold on;
    
    %%
%     figure(67);clf;    
%     plot3(x(:,1),x(:,2),y2(:,1),'ro'); hold on;        
%     plot3(Xtr_0(:,1),Xtr_0(:,2),ymu(:,1)+fs2(:,1),'co'); hold on;
%     plot3(Xtr_0(:,1),Xtr_0(:,2),ymu(:,1)-fs2(:,1),'mo'); hold on;
%     plot3(x(:,1),x(:,2),ymu(:,1),'ko'); hold on;
    
%     %% 
%     figure(68);clf;
%     plot3(x(:,1),x(:,2),y2(:,1),'ro'); hold on;        
%     plot3(Xtr_0(:,1),Xtr_0(:,2),ymu(:,1),'co'); hold on;
%     plot3(Xtr_0(:,1),Xtr_0(:,2),dXtr_3(:,1),'b*'); hold on;
%     
%     dXtr_4 = [];
%     for num = 1:length(Xtr_0(:,1))
%         num_0 = vpa(subs(dxdt3,[x1,x2],[Xtr_0(num,1),Xtr_0(num,2)]));
%         dXtr_4 = [dXtr_4; num_0];
%     end
%     dXtr_4 = reshape(double(dXtr_4),1,[])';
%     plot3(Xtr_0(:,1),Xtr_0(:,2),dXtr_4(:,1),'ko'); hold on;
    
end