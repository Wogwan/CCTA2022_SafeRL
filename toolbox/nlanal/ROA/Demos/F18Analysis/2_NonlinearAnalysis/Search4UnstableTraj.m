%  Search4UnstableTraj.m
%  This file attempts to solve via Monte Carlo search:
%      min  x0'*N*x0
%      subject to:  dx/dt = f(x) diverges from the initial condtion x0
%
%  The file randomly picks initial conditions on the ellipse x0'*N*x0 = g.
%  If a divergent trajectory is found then g is decreased by gfac and the
%  search continues.  A trajectory is considered "convergent" if it 
%  satisfies one of two conditions:
%     1) norm( x(t) ) < nsmall before Tfinal
%     2) abs( x_i(t) ) < xmall(i) for each state i before Tfinal 
%  Any trajectory not satisfying one of these two conditions is considered
%  "divergent".  The final solution should be checked to verify that the
%  trajectory actually does diverge from the initial condition.


%==========================================================================
% Create Shape Function
d2r = pi/180;
Dmax = diag([10*d2r 25*d2r 35*d2r 30*d2r 15*d2r 25*d2r 20*d2r]);  
N = inv(Dmax^2);   %  Scale by inverse of max state values
%N = N/max(N(:));   %  Normalize by the largest entry
p = x'*N*x;

%------------------------------------------------------------------
% Initialization of Simulation Parameter
%------------------------------------------------------------------

Ntry = 500000;         % Number of Monte carlo simulation runs 
Tfinal = 70;           % Simulation stopping time

gtry = 10.00;            % Size of the Initial Ellipsoid
gfac = 0.995;           % Decreasing Factor for the Ellipsoid

ns = length(x); 
x0min = zeros(ns,1);  

ind =1;                 % Index Initialization for Saving Data

%------------------------------------------------------------------
% Run Monte Carlo Simulation 
%------------------------------------------------------------------
tic
for i1=1:Ntry;

    % Find random initial condition and scale so that x0'*N*x0 = gtry
    x0 = randn(ns,1);
    x0(end) = 0;   % Special for the A/C example--set controller state to zero
    scl = sqrt( (x0'*N*x0)/gtry );
    x0 = x0/scl;
    
    % Simulate and check for convergence
    event_params.nbig = 1e6; event_params.nsmall = 1e-6;
    event_params.xbig = 300; event_params.xsmall = 1e-4;
    event_params.funchandle = [];
    [xtraj,xconv]= psim(f,x,x0,Tfinal,event_params);
    
     n0 = norm(x0);
     nf = norm(xtraj{2}(end,:));
     div = ( nf > 5*n0 );
     %div = (xconv==0);  
     % plot(xtraj{1},xtraj{2}); %title(['beta = ' num2str(gmin)]); 
    if div
        
        % Found divergent trajectory
        [gmin,minidx]=min(sum(xtraj{2}*N.*xtraj{2},2)); 
        fprintf('Div. traj. found at step=%d \t',i1);
        fprintf('gtry = %4.6f\t gmin= %4.6f\n',gtry,gmin);        
        plot(xtraj{1},xtraj{2}); title(['beta = ' num2str(gmin)]); 
        drawnow;
        
        % Save the Data once found divergent trajectories
        x0min = xtraj{2}(minidx,:)';
        xun(:,ind) = x0min; %x0;
        betaun(ind) = gmin; %gtry;
        divstep(ind) = i1;
        save(V0,'xun', 'betaun' ,'i1'); 
        ind = ind +1; 
        
        % Decrease the ellipsoid shape since divergent trajectories found 
        %gtry = gfac*gtry;
        gtry = gfac*gmin;
    end
    
    % Update Display every 100 iterations
    if rem(i1,1000)== 0
        fprintf('Current Iteration =%d \t gtry = %4.6f\n',i1,gtry)
    end
     
end
toc

%-----------------------------------------------------------------------
% Simulate smallest IC and verify divergence
%-----------------------------------------------------------------------
[xtraj,xconv] = psim(f,x,x0min,Tfinal,event_params);
figure(2)
plot(xtraj{1},xtraj{2}); title(['beta = ' num2str(gmin)]); 
betaub = double(subs(p,x,x0min))
title(['beta = ' num2str(betaub)]);

% Results
% Ntry = 10000
% beta = 3.131402
% t =  1721.4 sec

