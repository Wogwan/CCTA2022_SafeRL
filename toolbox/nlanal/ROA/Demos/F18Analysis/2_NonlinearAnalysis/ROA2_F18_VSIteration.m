%==========================================================================
%   F18_ROA_Vs Iteration
%   Compute an estimate of the region of attraction for  
%   F/A-18 using the V-S Iteration.
%==========================================================================

clc; 

%==========================================================================
% Generate Model 
% State Ordering 
% x     = [beta  ;      alpha;   p;   q;  r; phi; xc];
%==========================================================================



%==========================================================================
% Construct Lyap function from linearization 
Vlin=linstab(f,x);
Vdeg = Vlin.maxdeg;


%==========================================================================
% Create multiplier and function to force dV/dt<0
z2 =  monomials(x, 1); 
L2 =  1e-6*(x'*x);
L1 =  1e-6*(x'*x);
z1 =  monomials(x,0:1);
Vdeg = 2;
zV = monomials(x,2:Vdeg); 


%==========================================================================
% Create Shape Function
d2r = pi/180;
Dmax = diag([10*d2r 25*d2r 35*d2r 30*d2r 15*d2r 25*d2r 20*d2r]);  
N = inv(Dmax^2);   %  Scale by inverse of max state values
%N = N/max(N(:));   %  Normalize by the largest entry
p = x'*N*x;



%==========================================================================
% Set options for V-s iteration
Nsteps = 40; 
ph = [];
opts = [];
opts.gmin = 0; 
opts.gmax = 50;
opts.L2 = L2;
opts.bistol = 1e-6; % For Baseline Model
opts.checkfeas = 'fast'; 
opts.fid = 0;

%==========================================================================
% Initialize Storage Variable
info.iteration = struct; 
info.iteration.V = cell(Nsteps,1); 
info.iteration.beta = cell(Nsteps,1); 
info.iteration.gamma = cell(Nsteps,1);
info.iteration.s1 = cell(Nsteps,1);
info.iteration.s2 = cell(Nsteps,1);
info.iteration.time = cell(Nsteps,1);

for i1=1:Nsteps;    
    
    tic; 
    %======================================================================
    % Find V step: 
    % Hold s1, s2, b, g fixed and solve the Gamma Step and Beta Step 
    % for Appropiate Lyapunov function with an additional constraint 
    % V-L1 in SOS
    %======================================================================
     
    if i1==1
        V = Vlin;        
    else
        [V,c] = roavstep(f,p,x,zV,b,g,s1,s2,opts);
    end
    
    
    
    %======================================================================
    % Gamma Step: Solve the following problem 
    % {x:V(x) <= gamma} is contained in {x:grad(V)*f <0}
    % max gamma subject to 
    %                 -[grad(V)*f + (gamma - V)*s2] in SOS, s2 in SOS
    %======================================================================

    [gbnds,s2]=pcontain(jacobian(V,x)*f+L2,V,z2,opts);
    g = gbnds(1);

    
    
    %======================================================================
    % Beta Step: Solve the following problem 
    % {x: p(x)) <= beta} is contained in {x: V(x) <= gamma}
    % max beta subject to 
    %                 -[(V - gamma) + (beta - p)*s1] in SOS, s1 in SOS
    %======================================================================

    [bbnds,s1]=pcontain(V-g,p,z1,opts);
    b = bbnds(1);
    
    fprintf('i1 = %d  \t beta = %4.6f\n',i1,b);
    
    info.iteration(i1,1).V = V; 
    info.iteration(i1,1).beta = b; 
    info.iteration(i1,1).gamma = g;
    info.iteration(i1,1).s1 = s1;
    info.iteration(i1,1).s2 = s2;
    info.iteration(i1,1).time = toc;
    
    save(V0,'info','xtrim','utrim','N','Dmax');
   
    opts.bistol = 1e-4; 
end



