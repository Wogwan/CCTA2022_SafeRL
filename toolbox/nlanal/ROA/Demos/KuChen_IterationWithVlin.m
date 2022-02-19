%------------------------------------------------------------------
% Ku/Chen example with VS Iteration Initialized with Vlin
%   Compute an estimate of the region of attraction for Ku/Chen
%   example using the V-S iteration.  The iteration is initialized
%   with the Lyapunov function from linearized analysis.
%
% References: 
%  1) H. Ku and C. Chen. Stability study of a third-order servomechanism
%     with multiplicative feedback control. AIEE Transactions, Part 1,
%     77:131-136, 1958.
%
%  2) W. Tan.  Nonlinear Control Analysis and Synthesis using
%     Sum-of-Squares Programming. Ph.D. Thesis, Univ. of California,
%     Berkeley, Spring 2006.
%     (See Section 3.1.4.3 for results on Ku/Chen example)
%
% 11/18/2010 Updated with roaest code.  
%------------------------------------------------------------------
clear all; close all;

% Create vector field
pvar x1 x2 x3;
x = [x1;x2;x3];
x1dot = -x2;
x2dot = -x3;
x3dot = -0.915*x1+(1-0.915*x1^2)*x2-x3;
f = [x1dot; x2dot; x3dot];


% Create options for estimating region of attraction
Vdeg = 4;
p = x'*[12.5 -8.1 3.0; -8.1 20.8 -8.5; 3.0 -8.5 13.4]*x;
zV = monomials(x,2:Vdeg);    
z1maxd = ceil((Vdeg-p.maxdeg)/2);
z1 = monomials(x, 0:z1maxd ); 
z2 = monomials(x, 1:2 );
L2 = 1e-6*(x'*x);
NstepMinTrace = 20;
roadisplay = 'on';
ropt = roaoptions(f,x,'zV',zV,'p',p,'z1',z1,'z2',z2,'L2',L2,...
    'NstepMinTrace',NstepMinTrace,'display',roadisplay);


% Run the iteration code
[b,V,g,s1,s2,iter] = roaest(f,x,ropt);

%Extract Results

% Plot 3d V/p contours for ROA estimate
figure(1)
domain = [-2 2 -2 2 -2 2];
npts = 50*[1 1 1];
ph1 = patch(pcontour3(V,g,domain,npts));
set(ph1, 'FaceColor', 'red', 'FaceAlpha', 0.25, 'EdgeColor', 'red' );
ph2 = patch(pcontour3(p,b,domain,npts));
xlabel('x1');ylabel('x2');zlabel('x3');
set(ph2, 'FaceColor', 'blue', 'EdgeColor', 'none' );
view(90,0)  % Set view for x2/x3 slice

% Plot 2d slice of V/p contours for ROA estimate
figure(2);
V23 = subs(V,x1,0);
p23 = subs(p,x1,0);
[C,ph(1)]=pcontour(V23,g,domain(3:end),'r'); hold on;
[C,ph(2)]=pcontour(p23,b,domain(3:end),'b'); hold off;
axis([-1 1 -1 1])
