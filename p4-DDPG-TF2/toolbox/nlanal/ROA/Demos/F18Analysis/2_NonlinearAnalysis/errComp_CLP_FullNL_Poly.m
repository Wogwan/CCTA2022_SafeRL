%
% Evaluate F18 Full Nonlinear and Cubic polynomial Closed Loop Model at 
% different points. The sampling of the points are done in a gridded
% region of X'*N*X <= beta, where beta is from Monte Carlo results. The
% relative error norm is plotted as the points are sampled from increasingly
% larger ellipsoid region less than beta. 
% 
% Abhijit 02/16/10

clear all

%==========================================================================
% Set up Model 
baseline = 0;  
revised = 1; 
%==========================================================================
% ------------- Cases for Trimming Nonlinear Model 
 
Level       = 0; 
CoordTurn   = 1;
UnCoordTurn = 0;
opts.poly   = 0; 

% ---------- Generate Trim Values / Open-loop linear model
GenF18LinModel; 

%----------- Generate Linear Closed-Loop Models / Load SS for Controller
GenF18LinCLPModel; 


%==========================================================================
% ------------ Set Model and Controllers

if baseline 
    fB = load('2_Mar14_PolyF18BaselineModel_Phi_35');
    f = fB.xcldotB;  
    xx = fB.x; 

    % ------------ Set Controller 
    opts.K = CssB; 
    opts.Acl = fB.CLB.A; 
    
else
    fR = load('2_Mar14_PolyF18RevisedModel_Phi_35');
    f = fR.xcldotR;  
    xx = fR.x; 

    % ------------ Set Controller 
    opts.K = CssR;  
    opts.Acl = fR.CLR.A; 
end

%==========================================================================
% --- These output matrix can be formed to speed up 
Jy = jacobian(y,u);
Hy = y-Jy*u;


%--------------------- TEST out if Closed-loop realization is right or not

if 0 
    x0 = zeros(7,1); 
    % ------ Evaluate  output matrix and inverse of M 
    Jeval = double(subs(Jy,xx,x0 ));
    Heval = double(subs(Hy,xx,x0 )); 
    Meval  = double(eye(7) + Jeval*opts.K.D); 
    Minv  = inv(Meval); 

    % Build up Options 

    opts.J = Jeval; 
    opts.H = Heval; 
    opts.Minv = Minv; 
    opts.xtrim = xtrim;
    opts.utrim = utrim; 

    [xcldot2,f_full,f1,g1]  = f18Clpfull([],x0,opts)

end

%==========================================================================
% ---------------------- Sample and Evaluate at X'*N*X <= beta


%----------------------  Create Shape Function N
d2r = pi/180;
Dmax = diag([10*d2r 25*d2r 35*d2r 30*d2r 15*d2r 25*d2r 20*d2r]);  
N = inv(Dmax^2);   %  Scale by inverse of max state values


%----------------------  Initialize 

b       = 6;                            % MC beta Vlaue
Ndiv    = 600;                           % Number of divisions for creating a beta grid  
bgrid   = linspace(0.001,b ,Ndiv);      % beta grid 
Nsample = 50;                           % at each beta sample points
norm_relerVec = zeros(Ndiv, Nsample); 

for i1 = 1:Ndiv 
    
    for j1 = 1:Nsample

    % Find random initial condition and scale so that x0'*N*x0 = bgrid(i1)
%     x0 = randn(7,1);
%     x0(end) = 0;   % Special for the A/C example--set controller state to zero
%     scl = sqrt( (x0'*N*x0)/bgrid(i1) );
%     x0 = x0/scl;
    
% Find Random IC within the range and scale it back to 
    Xdata.range = [-10 10; -5 45; -35 35; -30 30; -15 15; 10 60; -20 20]*d2r; 
    x0 = Xdata.range(:,1) + rand(7,1).*(Xdata.range(:,2) - Xdata.range(:,1));
    x0(end) = 0;
    scl = sqrt( (x0'*N*x0)/bgrid(i1) );
    x0 = x0/scl;
    % ------ Evaluate  output matrix and inverse of M 
    Jeval = double(subs(Jy,xx,x0));
    Heval = double(subs(Hy,xx,x0)); 
    Meval  = double(eye(7) + Jeval*opts.K.D); 
    Minv  = inv(Meval); 

    % Build up Options 

    opts.J = Jeval; 
    opts.H = Heval; 
    opts.Minv = Minv; 
    opts.xtrim = xtrim;
    opts.utrim = utrim; 

    [xcldot2,f_full,f1,g1]  = f18Clpfull([],x0,opts);
    f_poly  = double(subs(f,xx,x0)); 
    
    % -------------------------------------------------------------------
    % --------------- Compute error statistics

    erVec = f_poly - f_full;
    norm_erVec = norm(erVec,2); %sqrt(sum((erVec.*erVec),1)); 
    norm_full =  norm(f_full,2); %sqrt(sum((f_full.*f_full),1));
    norm_relerVec(i1,j1) = norm_erVec/norm_full;
    
%      plot(bgrid(i1),norm_relerVec(i1,j1),'x'); 
%      hold on 
%     
    end
end

% box('on');grid('on'); hold('all');
% xlabel('||x^TNx||','FontSize',18)
% ylabel('relative error','FontSize',18)


return

figure5 =figure;
axes('Parent',figure5,'FontSize',18);
for i1 = 1:Ndiv 
    for j1 = 1:Nsample
       plot(bgrid(i1),norm_relerVec(i1,j1),'x'); 
       hold on 
    end
end



