load 2_Mar14_PolyF18BaseLineModel_Phi_35
f = cleanpoly(xcldotB,10^-06);  
d2r = pi/180;  r2d = 1/ d2r; 

load 1_Mar14_ROAVS_F18Baseline_2_35; 

V0 = strcat('1_Mar14_ROAVS_F18Baseline_3_',num2str(xtrim(7)*r2d));

%==========================================================================
% Create multiplier and function to force dV/dt<0
z2 =  monomials(x, 1); 
L2 =  1e-6*(x'*x);
L1 =  1e-6*(x'*x);
z1 =  monomials(x,0:1);
Vdeg = 4;
zV = monomials(x,2:Vdeg); 

%==========================================================================
% Create Shape Function
d2r = pi/180;  r2d = 1/ d2r; 
Dmax = diag([10*d2r 25*d2r 35*d2r 30*d2r 15*d2r 25*d2r 20*d2r]);  
N = inv(Dmax^2);   %  Scale by inverse of max state values
%N = N/max(N(:));   %  Normalize by the largest entry
p = x'*N*x;


opts.gmin = 0; 
opts.gmax = 50;
opts.L2 = L2;
opts.bistol = 1e-4; % For Baseline Model
opts.checkfeas = 'fast'; 
opts.fid = 0;

s1 = info.iteration(end).s1; 
s2 = info.iteration(end).s2; 
b = info.iteration(end).beta; 
g = info.iteration(end).gamma; 


for i1=66:80;    
    
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
    
    save(V0,'info','xtrim','utrim','N');
   
end

