%------------------------------------------------------------------
% van der Pol example
% This example demonstrates that the effect of multiplier order 
% for the gamma step on the V-S iteration.  We compute the
% quadratic Lyapunov function from the linearization of the VDP
% oscillator.  Then we compute a gamma step using second and
% fourth order multipliers. Both multipliers, s2 and s4, give the 
% same value of gamma.  We then hold the multiplier fixed at either 
% s2 or s4 and compute feasible region of the Lyapunov coefficients.
% in the V-step. The fourth order multiplier gives a much larger
% feasible region for the Lyapunov function in the V-step and hence
% leads to better performance in the V-S iteration.
%
% 11/19/2010 Checked and Work with new codes
%------------------------------------------------------------------
clear all; close all;

% Create vector field
pvar x1 x2;
x = [x1;x2];
x1dot = -x2;
x2dot = x1+(x1^2-1)*x2;
f = [x1dot; x2dot];

% Create shape function and vector of Lyapunov function monomials
Vdeg = 2;
p = x'*x;
zV = monomials(x,2:Vdeg);    
L2 = 1e-6*(x'*x);

% Construct Lyap function from linearization 
fprintf('\nConstructing Lyap. function from linearization');
Vlin=linstab(f,x);
V = Vlin;
    
% get gamma : 2nd order multiplier
fprintf('\nGamma step with second order multiplier');
z2 = monomials(x, 1 );
[g2,s2]=pcontain(jacobian(V,x)*f+L2,V,z2);
g2 = g2(1);

% get gamma : 4th order multiplier
fprintf('\nGamma step with fourth order multiplier');
z4 = monomials(x, 1:2 );
[g4,s4]=pcontain(jacobian(V,x)*f+L2,V,z4);
g4 = g4(1);

fprintf('\nBoth multipliers should give the same value of gamma:');
fprintf('\ng2=%4.3f \t g4 = %4.3f\n',g2,g4);

% Now fix multiplier s2 and determine space of feasible V's 
% [We hold the third coeff fixed and search over the remaing 2 coefs]
fprintf('\nComputing feasible region for V while holding s2 fixed\n');
Npts1a = 17;
Npts2a = 19;
c1data = linspace(1.498,1.504,Npts1a);
c2data = linspace(-1.012,-0.998,Npts2a);
c3a = 1;
feas2 = zeros(Npts2a,Npts1a);
for i1 = 1:Npts1a
    for i2 = 1:Npts2a;
        V = c1data(i1)*x1^2+c2data(i2)*x1*x2+c3a*x2^2;        
        Vdot = jacobian(V,x)*f;
        
        sosconstr =  -( (Vdot+L2)+(g2-V)*s2 );
        [info,dopt,sossol] = sosopt(sosconstr,x);        
        if ~isempty(sossol)
            feas2(i2,i1) = 1;
        end        
    end
end
hold off;


% Now fix multiplier s4 and determine space of feasible V's 
fprintf('Computing feasible region for V while holding s4 fixed\n');
Npts1c = 25;
Npts2c = 27;
c1datb = linspace(1.45,1.6,Npts1c);
c2datb = linspace(-1.15,-0.95,Npts2c);
c3c = 1;
feas4 = zeros(Npts2c,Npts1c);
for i1 = 1:Npts1c
    for i2 = 1:Npts2c;
        V = c1datb(i1)*x1^2+c2datb(i2)*x1*x2+c3c*x2^2;        
        Vdot = jacobian(V,x)*f;
        
        sosconstr =  -( (Vdot+L2)+(g4-V)*s4 );
        [info,dopt,sossol] = sosopt(sosconstr,x);        
        if ~isempty(sossol)
            feas4(i2,i1) = 1;
        end        
    end
end

% Plot contours of feasible regions
fprintf('Plotting feasible sets of Lyap. function coefs\n');
figure(1);
contour(c1data,c2data,feas2,0.1*[1 1],'b');hold on;
contour(c1datb,c2datb,feas4,0.1*[1 1],'r');hold off;
axis([c1datb(1) c1datb(end) c2datb(1) c2datb(end)]);
xlabel('c1');
ylabel('c2');
legend('Feas for s2','Feas for s4');
title('Feasible Region for Quadratic Lyapunov Fnc Coefs (with c3=1)')
