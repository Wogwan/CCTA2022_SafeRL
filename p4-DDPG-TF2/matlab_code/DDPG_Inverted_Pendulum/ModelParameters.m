clear all, close all, clc

% Cart parameters
m = 1; % pendulum mass
M = 5; %cart mass
L = 2; %pendulum length
g = -9.81; % Gravity
d = 1; % dampping term opposing the force input 
% slows the pendulum down, the higher this is the faster the pendulum slows down 

s = 1; % pendulum up (s=1)


A = [0 1          0                0;               %x
     0 -d/M       -m*g/M           0;            %x_dot
     0 0          0                1;            % theta
     0 -s*d/(M*L) -s*(m+M)*g/(M*L) 0];          % theta_dot
                                                     %state space representation of the dependent forces
 
 
B = [0; 1/M; 0; s*1/(M*L)];                  %state space representation of the linear forces

eig(A) %% give the eigen values of A, where the poles are. 

C = eye(4); % sets an identity matrix. 
sys = ss(A,B,C,0*B); % converts a dynamic system, to space state representation

%%

tspan = 0:.001:10;
if(s==-1)
    y0 = [0; 0; 0; 1.5];
    %initial is used to plot the output of the system of equations.
    %yL is the output response, t is the time vector for the simulation,
    %and xL is the state trajectories. Uses sys to make this. 
    
    [yL,t,xL] = initial(sys,y0,tspan); 
    % solves the ode for the given syste, time and initial conditions.
    % returns the time frame, and the output vector. 
    [t,yNL] = ode45(@(t,y)cartpend(y,m,M,L,g,d,0),tspan,y0);
elseif(s==1)
    y0 = [0; 0; pi+.0001; 0];
    [yL,t,xL] = initial(sys,y0-[0; 0; pi; 0],tspan);
    [t,yNL] = ode45(@(t,y)cartpend(y,m,M,L,g,d,0),tspan,y0);
else    
end

figure (10);  
%plot(t,yL);
plot(t,yL+ones(10001,1)*[0; 0; pi; 0]'); % as far as I can see they are the same.
xlabel("x"); ylabel("y"); 
title("Plot for the Output Functions for the Pendulum Equations System");
hold on; 
%%hold off; 
pause;


figure (11);
for k=1:50:length(t)
    drawcartpend(yNL(k,:),m,M,L); 
    %hold on; % use if you want to see where the pendulum is at every
    %frame.   
end

%plotting the output data
%finding the graphical functions of the control parameters 
%Need to learn how to break these up into their respective components
figure (3); plot(yNL); legend("Cart Position", "Cart Acceleration", "Pendulum Angle", "Pendulum Rotational Speed"); 
xlabel("t sample 'frames'"); ylabel("Amplitude"); title ("y Parameter Impulse Response");  