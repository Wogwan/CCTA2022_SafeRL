%------------------------------------------------------------------
% Hauser & Lai Example 2 with Bootstrap
%   Compute an estimate of the region of attraction for Hauser & 
%   Lai Example 2 using the bootstrap method.  The V-S iteration 
%   is run with a quadratic Lyapunov function and the shape function 
%   p(x)=x'*x.  The iteration is initialized with the Lyapunov function 
%   from linearized analysis.  Then the iteration is run again with the 
%   shape function set equal to the quadratic Lyapunov function obtained 
%   during the first iteration.
%
% References: 
%  1) J. Hauser and M.C. Lai. Estimating quadratic stability domains by 
%     nonsmooth optimization. Proceedings of the ACC. Chicago, IL. 
%     pp. 571–576, 1992.
%  2) U. Topcu, A. Packard, and P. Seiler, Local stability analysis 
%     using simulations and sum-of-squares programming," To appear in 
%     Automatica. 
%  3) U. Topcu.  Quantitative Local Analysis of Nonlinear Systems.
%     Ph.D. Thesis, Univ. of California, Berkeley, Fall 2008.
%     (See Section 3.3.2 for results on this example)
%
% 11/18/2010 Updated with roaest code.  
%------------------------------------------------------------------
clear all; close all;

% Create vector field
pvar x1 x2 x3;
x = [x1;x2;x3];
x1dot = x2+2*x2*x3;
x2dot = x3;
x3dot = -0.5*x1-2*x2-x3;
f = [x1dot; x2dot;x3dot];

% Create options for estimating region of attraction
Vdeg = 2;
p = x'*x;
z1maxd = ceil((Vdeg-p.maxdeg)/2);
z1 = monomials(x, 0:z1maxd ); 
z2 = monomials(x, 1:2 );
L2 = 1e-6*(x'*x);
Nsteps = 20; 
ropt = roaoptions(f,x,'Vdeg',Vdeg,'p',p,'z1',z1,'z2',z2,'L2',L2,'Nstep',Nsteps);

% Run the iteration code
info = roaest(f,x,ropt);

%Extract Results
V = info.iteration(Nsteps).V; 
b = info.iteration(Nsteps).beta; 
g = info.iteration(Nsteps).gamma; 

fprintf('Volume of {x : V(x)<=g} = %4.3f*(4/3*pi)\n',pvolume(V,g)/(4/3*pi));
fprintf('\n--------------------------------------------------\n');
fprintf('Quartic Lyapunov function \n');
fprintf('shape function & initial V initialized to final quadratic V  \n ');
fprintf('\n');

% Use Lyapunov function from previous iteration to "bootstrap" shape func.
p = V;
Vdeg = 4;
zV = monomials(x,2:Vdeg);    
z1maxd = ceil((Vdeg-p.maxdeg)/2);
z1 = monomials(x, 0:z1maxd ); 
Nsteps = 20; 
gamstepopts = gsosoptions('maxobj',50);
ropt = roaoptions(f,x,'Vdeg',Vdeg,'p',p,'z1',z1,'z2',z2,'L2',L2,'Nstep',...
                  Nsteps,'gammastepopts',gamstepopts,'Vin',V);

% Run the iteration code
info = roaest(f,x,ropt);

% extract results
V = info.iteration(Nsteps).V; 
b = info.iteration(Nsteps).beta; 
g = info.iteration(Nsteps).gamma;              

fprintf('Volume of {x : V(x)<=g} = %4.3f*(4/3*pi)\n',pvolume(V,g)/(4/3*pi));

% Plot 3d V/p contours for ROA estimate
figure(1)
domain = [-10 10 -10 10 -10 10];
npts = 50*[1 1 1];
ph1 = patch(pcontour3(V,g,domain,npts));
set(ph1, 'FaceColor', 'red', 'FaceAlpha', 0.25, 'EdgeColor', 'red' );
ph2 = patch(pcontour3(p,b,domain,npts));
xlabel('x1');ylabel('x2');zlabel('x3');
set(ph2, 'FaceColor', 'blue', 'EdgeColor', 'none' );
view(90,0)  % Set view for x2/x3 slice
ax = axis;

% Plot 2d slice of V/p contours for ROA estimate
figure(2);
V23 = subs(V,x1,0);
p23 = subs(p,x1,0);
[C,ph(1)]=pcontour(V23,g,domain(3:end),'r'); hold on;
[C,ph(2)]=pcontour(p23,b,domain(3:end),'b'); hold off;
axis(ax(3:end));

