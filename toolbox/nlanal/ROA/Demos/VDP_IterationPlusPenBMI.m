%------------------------------------------------------------------
% VDP with VS Iteration Initialized with Vlin
%   Compute an estimate of the region of attraction for the van 
%   der Pol oscillator using the V-S iteration.  The iteration is 
%   initialized with the Lyapunov function from linearized analysis.
%
% 11/19/2010 Updated but need to check PENBMI... 
%------------------------------------------------------------------
clear all; close all;

% Create vector field
pvar x1 x2;
x = [x1;x2];
x1dot = -x2;
x2dot = x1+(x1^2-1)*x2;
f = [x1dot; x2dot];

% Plot VDP limit cycle from backward simulation
X0 = [1/2; 1/2];
Tf = 100;
[xtraj,isconv]=psim(-f,x,X0,Tf);
t = xtraj{1};
xtraj = xtraj{2};
idx = min(find(t>=Tf*0.8));

figure(1)
hold off;
plot(xtraj(idx:end,1),xtraj(idx:end,2));
xlabel('x1');
ylabel('x2');
hold on;


% Create options for estimating region of attraction
Vdeg = 6;
p = x'*x;
zV = monomials(x,2:Vdeg);    
z1maxd = ceil((Vdeg-p.maxdeg)/2);
z1 = monomials(x, 0:z1maxd ); 
z2 = monomials(x, 1:2 );
L2 = 1e-6*(x'*x);
Nsteps = 5; 
gamstepopts = gsosoptions('maxobj',50);
ropt = roaoptions(f,x,'Vdeg',Vdeg,'p',p,'z1',z1,'z2',z2,'L2',L2,'Nstep',...
                  Nsteps,'gammastepopts',gamstepopts);

% Run the iteration code
info = roaest(f,x,ropt);

% plot sublevel sets of Lyapunov function and shape function
domain = [-3 3 -3 3];

V = info.iteration(Nsteps).V; 
b = info.iteration(Nsteps).beta; 
g = info.iteration(Nsteps).gamma; 
[C,ph(1)]=pcontour(V,g,domain,'r'); hold on;
[C,ph(2)]=pcontour(p,b,domain,'b'); hold off;
title(['Iteration = ' int2str(Nstep) ' beta = ' num2str(b)]);

% Call PenBMI
ic.beta = b;
ic.gamma = g;
ic.s1 = s1;
ic.s2 = s2;
ic.V = V;
[bopt,gopt,s1opt,s2opt,Vopt] = roabmi(f,x,p,zV,opts,ic);

set(ph,'LineStyle',':');
[C,ph(1)]=pcontour(Vopt,gopt,domain,'r');
[C,ph(2)]=pcontour(p,bopt,domain,'b');
set(ph,'LineWidth',2)
title(['PENBMI:  beta = ' num2str(bopt)]);
