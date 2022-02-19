%------------------------------------------------------------------
% Chesi example 3 with Global Lyapunov function
%   Compute a Lyapunov function for Chesi example 3 which proves x=0 is 
%   a globally stable equilibrium point.
%
% References: 
%  1) G. Chesi, A. Garulli, A. Tesi and A. Vicino. LMI-based
%     computation of optimal quadratic Lyapunov functions for
%     odd polynomial systems. Int. J. Robust Nonlinear Control
%     15:35–49, 2005.
%  2) U. Topcu, A. Packard, and P. Seiler, Local stability analysis 
%     using simulations and sum-of-squares programming," To appear in 
%     Automatica. 
%  3) U. Topcu.  Quantitative Local Analysis of Nonlinear Systems.
%     Ph.D. Thesis, Univ. of California, Berkeley, Fall 2008.
%     (See Section 3.3.2 for results on this example)
%  4) W. Tan and A. Packard. Stability region analysis using polynomial
%     and composite polynomial Lyapunov functions and sum of squares
%     programming. IEEE Transactions on Automatic Control. 
%     53:2:565-571, 2008.
%  
%------------------------------------------------------------------
clear all; close all;

% Create vector field
pvar x1 x2;
x = [x1;x2];
x1dot = -1.3*x1-1.4*x2-x1^5-1.8*x1^4*x2+1.7*x1^3*x2^2+3.2*x1^2*x2^3 ...
    - 0.4*x1*x2^4+0.9*x2^5;
x2dot = 2*x1-0.8*x2-4*x1^5+3.5*x1^4*x2-2.8*x1^3*x2^2-2.2*x1^2*x2^3 ...
        -0.1*x1*x2^4-0.8*x2^5;
f = [x1dot; x2dot];

% Use sosopt to find a Lyapunov function which proves
% x=0 is a Globally Exponentially Stable eq. pt.


% Define decision variable for quadratic Lyapunov function
zV = monomials(x,2:4);
V = polydecvar('c',zV,'vec');

% Constraint 1 : V(x) - L1 in SOS
L1 = 1e-6*(x1^2+x2^2);
sosconstr(1) = V-L1 >= 0;

% Constraint 2: -Vdot - L2 in SOS
% r is the exponential rate of convergence
L2 = 1e-6*(x1^2+x2^2);
Vdot = jacobian(V,x)*f;
sosconstr(2) =-Vdot-L2 >=0;

% Solve with feasibility problem
[info,dopt,sossol] = sosopt(sosconstr,x);

if info.feas
    disp('Global Stability Analysis: Feasible');
    Vs = subs(V,dopt);
    Vs = cleanpoly(Vs,1e-10)
else
    disp('Global Stability Analysis: Not Feasible');
    return
end

% Verify that Vs-L1 = z'*Q*z and Q>=0
z1 = sossol(1).z;
Q1 = sossol(1).Q;
e1 = (Vs-L1) - (z1'*Q1*z1);
fprintf('\nMax magnitude coefficient of s(1)-w1''*Q1*w1 is:')
disp(full(max(abs(e1.coefficient))))
fprintf('Minimum eigenvalue of Q1 is:')
disp(min(eig(Q1)));

% Verify that -(Vdot+L2)=w'*Q*w and Q>=0 
z2 = sossol(2).z;
Q2 = sossol(2).Q;
e2 = (-jacobian(Vs,x)*f) - (z2'*Q2*z2);
fprintf('\nMax magnitude coefficient of s(2)-w2''*Q2*w2 is:')
disp(full(max(abs(e2.coefficient))))
fprintf('Minimum eigenvalue of Q2 is:')
disp(min(eig(Q2)));

% Overlay trajectories of system and additional contours of V
domain=[-2 2 -2 2];
simdomain=[-2 2 -2 2];
Npts = 10;
x1ic = linspace(simdomain(1),simdomain(2),Npts);
x2ic = linspace(simdomain(3),simdomain(4),Npts);
x0 = [x1ic x1ic repmat(x1ic(1),1,Npts) repmat(x1ic(end),1,Npts);
    repmat(x2ic(1),1,Npts) repmat(x2ic(end),1,Npts) x2ic x2ic];
xtraj = psim(f,x,x0,20);

figure(1);
hold on;
[C,ph(1)]=pcontour(Vs,[16 8 4 2 1 0.5 0.1 0.01],simdomain,'k',[1e3 1e3]); 
for i1=1:size(x0,2)
    plot(xtraj{i1,2}(:,1),xtraj{i1,2}(:,2),'g');hold on;
end
hold off;
axis(domain);

%----------------
return;

% Plot Results
tmp = -2:0.5:2;
lt = length(tmp);
x0 = [tmp(:) repmat(-2,[lt 1]); tmp(:) repmat(+2,[lt 1]); ...
    repmat(-2,[lt 1]) tmp(:); repmat(+2,[lt 1]) tmp(:)];
tfinal = 10;
[xtraj,xconv]=psim(f,x,x0',tfinal);

for i1=1:length(xtraj)
    plot(xtraj{i1,2}(:,1),xtraj{i1,2}(:,2),'b'); hold on;
end

[X1,X2] = meshgrid(-2:.1:2, -2:.1:2);
Vgrid = double(subs(Vs,{x1 x2},{X1 X2}));
contour(X1,X2,Vgrid,[1 5 20 40 60])
xlabel('x1');
ylabel('x2');
hold off;
axis([min(tmp) max(tmp) min(tmp) max(tmp)]);
%----------------------


return
% Create shape function and vector of Lyapunov function monomials
% The shape function is from W. Tan's Thesis (Sec 3.1.4.2)
Vdeg = 2;
p = x'*x;
zV = monomials(x,2:Vdeg);    
z1maxd = ceil((Vdeg-p.maxdeg)/2);
z1 = monomials(x, 0:z1maxd ); 
z2 = monomials(x, 1:2 );
L2 = 1e-6*(x'*x);

% Construct Lyap function from linearization 
Vlin=linstab(f,x);

% Run V-s iteration
Nsteps = 20; 
ph = [];
opts = [];
opts.gmin = 0; 
opts.gmax = 50;
opts.L2 = L2;
for i1=1:Nsteps;    
    if i1==1
        V = Vlin;        
    else
        % roavstep
        [V,c] = roavstep(f,p,x,zV,b,g,s1,s2,opts);
        
        % Scale V by previous gamma.  This roughly normalizes V
        V = V/g;
    end
    
    % get gamma
    [gbnds,s2]=pcontain(jacobian(V,x)*f+L2,V,z2,opts);
    g = gbnds(1);

    % get beta
    [bbnds,s1]=pcontain(V-g,p,z1,opts);
    b = bbnds(1);
    
    fprintf('i1 = %d  \t beta = %4.3f\n',i1,b);

    % Plot V/p contours for ROA estimate 
    if  ~isempty(ph)
        delete(ph)
    end
    domain = [-2 2 -2 2];
    [C,ph(1)]=pcontour(V,g,domain,'r'); hold on;
    [C,ph(2)]=pcontour(p,b,domain,'b'); hold off;
    title(['Iteration = ' int2str(i1) ' beta = ' num2str(b)]);
    axis(domain)
    drawnow;    
end

fprintf('Volume of {x : V(x)<=g} = %4.3f*pi\n',pvolume(V,g)/pi);

% Overlay trajectories of system and additional contours of V
simdomain=[-2 2 -2 2];
Npts = 10;
x1ic = linspace(simdomain(1),simdomain(2),Npts);
x2ic = linspace(simdomain(3),simdomain(4),Npts);
x0 = [x1ic x1ic repmat(x1ic(1),1,Npts) repmat(x1ic(end),1,Npts);
    repmat(x2ic(1),1,Npts) repmat(x2ic(end),1,Npts) x2ic x2ic];
xtraj = psim(f,x,x0,20);

figure(1);
ax = axis;
hold on;
[C,ph(1)]=pcontour(V,g*[0.3 0.1 0.05 0.01],domain,'k',[1e3 1e3]); 
for i1=1:size(x0,2)
    plot(xtraj{i1,2}(:,1),xtraj{i1,2}(:,2),'g');hold on;
end
hold off;
axis(ax);