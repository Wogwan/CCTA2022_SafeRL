%% Actor Critic DDPG code for PendulumModel

% If not done, run ModelParameters first. 
ModelParameters

% sets the random number generator to a seed, for reproducibility of the
% training
rng(0); 

% opens the relevant model. 
model = 'PendulumModel';
open_system(model);

%% Linking the Agent Block to MATLAB file
%linking the RL agent block to the code
agent = 'PendulumModel/RL Agent';

% Signal bus for the observation of the RL agent block 
% obsBus = Simulink.Bus(); % Creating the relevant signals within the Simulink Virtual Bus 
% obs(1) = Simulink.BusElement;
% obs(1).Name = 'x';
% obs(2) = Simulink.BusElement;
% obs(2).Name = 'velocity';
% obs(3) = Simulink.BusElement;
% obs(3).Name = 'sintheta';
% obs(4) = Simulink.BusElement;
% obs(4).Name = 'costheta';
% obs(5) = Simulink.BusElement;
% obs(5).Name = 'thetadot';
% obsBus.Elements = obs; 

obsInfo = rlNumericSpec([5 1]);
obsInfo.Name = 'observations';


% Action Signal bus, for the action of the RL agent into the control effort
% input of the pendulum model. 
% actBus = Simulink.Bus();% Linking Tau bus signal to the simulink 
% act(1) = Simulink.BusElement;
% act(1).Name = 'tau';
% act(1).Min = -10;
% act(1).Max = 10;
% actBus.Elements = act;
actInfo = rlNumericSpec((1),...
    'LowerLimit',(-15),...
    'UpperLimit',(15));
actInfo.Name = 'torque';

% Creating the objects within Simulink
% obsInfo = bus2RLSpec('obsBus','Model',model); 
%actInfo = bus2RLSpec('actBus','Model',model);

% Creates the end environment - TOP LAYER for the RL agent. 
env = rlSimulinkEnv(model,agent,obsInfo,actInfo);

% Reset for the RL agent after epoch termination - effectively model
% initial conditions 
%env.ResetFcn = @(in)setVariable(in,'y0',[0,0, pi+0.001,0] ,'Workspace', model);

% Simulation time and sample rate
Ts = 0.05;
Tf = 10;

% getting the number of observations for the size of the input.
obsInfo = rlNumericSpec([5, 1, 1]);
obsInfo.Name = 'observations';
%numObservations = obsInfo.Dimension(1);


%% Creating the DDPG Agent. 
% Here an actor neural network will be created. 
%% Create Array of Layers
% Actornetwork = layerGraph();
% 
% tempLayers = imageInputLayer([1,1,1],"Name","thetadot",'Normalization','none');
% Actornetwork = addLayers(Actornetwork,tempLayers);
% 
% tempLayers = imageInputLayer([1,1,1],"Name","velocity",'Normalization','none');
% Actornetwork = addLayers(Actornetwork,tempLayers);
% 
% tempLayers = imageInputLayer([1,1,1],"Name","sintheta",'Normalization','none');
% Actornetwork = addLayers(Actornetwork,tempLayers);
% 
% tempLayers = imageInputLayer([1,1,1],"Name","costheta",'Normalization','none');
% Actornetwork = addLayers(Actornetwork,tempLayers);
% 
% tempLayers = imageInputLayer([1,1,1],"Name","x",'Normalization','none');
% Actornetwork = addLayers(Actornetwork,tempLayers);
% 
% tempLayers = [
%     additionLayer(5,"Name","addition")
%     fullyConnectedLayer(10,"Name","ActorFC1")
%     reluLayer("Name","ActorRelu1")
%     fullyConnectedLayer(10,"Name","ActorFC2")
%     reluLayer("Name","ActorRelu2")
%     fullyConnectedLayer(1,"Name","ActorFC3")
%     tanhLayer("Name","ActorTanh1")
%     scalingLayer("Name","ActorScaling1","Scale",10)];
% Actornetwork = addLayers(Actornetwork,tempLayers);
% 
% % clean up helper variable
% clear tempLayers;
% %% Connect Layer Branches
% % Connect all the branches of the network to create the network graph.
% Actornetwork = connectLayers(Actornetwork,"thetadot","addition/in5");
% Actornetwork = connectLayers(Actornetwork,"costheta","addition/in4");
% Actornetwork = connectLayers(Actornetwork,"velocity","addition/in2");
% Actornetwork = connectLayers(Actornetwork,"sintheta","addition/in3");
% Actornetwork = connectLayers(Actornetwork,"x","addition/in1");
% 
%% Plot Layers

Actornetwork = [
    imageInputLayer([5,1,1],"Name","observations","Normalization","none")
    fullyConnectedLayer(200,"Name","ActorFC1")
    reluLayer("Name","ActorRelu1")
    fullyConnectedLayer(200,"Name","ActorFC2")
    reluLayer("Name","ActorRelu2")
    fullyConnectedLayer(1,"Name","ActorFC3")
    tanhLayer("Name","ActorTanh1")
    scalingLayer('Name','ActorScaling1','Scale', [15])];


actorOptions = rlRepresentationOptions('LearnRate',5e-4,'GradientThreshold',1);

actor = rlDeterministicActorRepresentation(Actornetwork,obsInfo,actInfo,...
    'Observation',{'observations'},'Action',{'ActorScaling1'},actorOptions);

%% Critic Network creation

Criticnetwork = layerGraph();
%% Add Layer Branches
% Add the branches of the network to the layer graph. Each branch is a linear 
% array of layers.

% tempLayers = imageInputLayer([1, 1, 1],"Name","thetadot",'Normalization','none');
% Criticnetwork = addLayers(Criticnetwork,tempLayers);
% 
% tempLayers = [
%     imageInputLayer([1 1 1],"Name","actioninput","Normalization","none")
%     fullyConnectedLayer(10,"Name","fcaction")];
% Criticnetwork = addLayers(Criticnetwork,tempLayers);
% 
% tempLayers = imageInputLayer([1, 1, 1],"Name","x",'Normalization','none');
% Criticnetwork = addLayers(Criticnetwork,tempLayers);
% 
% tempLayers = imageInputLayer([1, 1, 1],"Name","sintheta",'Normalization','none');
% Criticnetwork = addLayers(Criticnetwork,tempLayers);
% 
% tempLayers = imageInputLayer([1, 1, 1],"Name","costheta",'Normalization','none');
% Criticnetwork = addLayers(Criticnetwork,tempLayers);
% 
% tempLayers = imageInputLayer([1, 1, 1],"Name","velocity",'Normalization','none');
% Criticnetwork = addLayers(Criticnetwork,tempLayers);
% 
% tempLayers = [
%     additionLayer(5,"Name","addition")
%     fullyConnectedLayer(10,"Name","fcstate1")
%     reluLayer("Name","staterelu")
%     fullyConnectedLayer(10,"Name","fcstate2")];
% Criticnetwork = addLayers(Criticnetwork,tempLayers);
% 
% tempLayers = [
%     additionLayer(2,"Name","addition_2")
%     reluLayer("Name","relu")
%     fullyConnectedLayer(1,"Name","criticoutput")];
% Criticnetwork = addLayers(Criticnetwork,tempLayers);
% 
% clean up helper variable
% clear tempLayers;
 %% Connect Layer Branches
%  Connect all the branches of the network to create the network graph.
% 
% Criticnetwork = connectLayers(Criticnetwork,"thetadot","addition/in5");
% Criticnetwork = connectLayers(Criticnetwork,"x","addition/in1");
% Criticnetwork = connectLayers(Criticnetwork,"sintheta","addition/in3");
% Criticnetwork = connectLayers(Criticnetwork,"costheta","addition/in4");
% Criticnetwork = connectLayers(Criticnetwork,"fcaction","addition_2/in2");
% Criticnetwork = connectLayers(Criticnetwork,"velocity","addition/in2");
% Criticnetwork = connectLayers(Criticnetwork,"fcstate2","addition_2/in1");

%% Plot Layers

tempLayers = [
    imageInputLayer([1,1,1],"Name","actioninput","Normalization","none")
    fullyConnectedLayer(200,"Name","fcaction")];
Criticnetwork = addLayers(Criticnetwork,tempLayers);

tempLayers = [
    imageInputLayer([5,1,1],"Name","observations","Normalization","none")
    fullyConnectedLayer(200,"Name","fcstate1")
    reluLayer("Name","staterelu")
    fullyConnectedLayer(200,"Name","fcstate2")];
Criticnetwork = addLayers(Criticnetwork,tempLayers);

tempLayers = [
    additionLayer(2,"Name","addition_2")
    reluLayer("Name","relu")
    fullyConnectedLayer(1,"Name","criticoutput")];
Criticnetwork = addLayers(Criticnetwork,tempLayers);

% clean up helper variable
clear tempLayers;
%% Connect Layer Branches
% Connect all the branches of the network to create the network graph.

Criticnetwork = connectLayers(Criticnetwork,"fcaction","addition_2/in2");
Criticnetwork = connectLayers(Criticnetwork,"fcstate2","addition_2/in1");
plot(Criticnetwork);

%% Creating the Critic Network Options
criticOpts = rlRepresentationOptions('LearnRate',1e-03,'GradientThreshold',1);
%%
useGPU = 1;
if useGPU
   criticOpts.UseDevice = "gpu";
end
%%
critic = rlQValueRepresentation(Criticnetwork,obsInfo,actInfo, ...
                'Observation',{'observations'},'Action',{'actioninput'},criticOpts);
                     
%% Creation of DDPG options 

agentOpts = rlDDPGAgentOptions(... 
    'SampleTime',Ts,...
    'TargetSmoothFactor',1e-3,...
    'ExperienceBufferLength',1e6,...
    'MiniBatchSize',256);
agentOptions.NoiseOptions.Variance = 1;
agentOptions.NoiseOptions.VarianceDecayRate = 0.01;
%% Linking to the RL agent in Simulink 
%%
if useGPU
   agentOpts.UseDevice = "gpu";
end
agent = rlDDPGAgent(actor,critic,agentOpts);

%% Training the DDPG options

maxepisodes = 5000;
maxsteps = ceil(Tf/Ts);
trainOpts = rlTrainingOptions(...
    'MaxEpisodes',maxepisodes,...
    'MaxStepsPerEpisode',maxsteps,...
    'ScoreAveragingWindowLength',10,...
    'Verbose',false,...
    'Plots','training-progress',...
    'StopTrainingCriteria','AverageReward',...
    'StopTrainingValue',100);

%% Training 

    % Train the agent.
    trainingStats = train(agent,env,trainOpts);

%% Running Final Result 

simOptions = rlSimulationOptions('MaxSteps',200);
experience = sim(env,agent,simOptions);

%% plotting outputs
outputcosangle = experience.Observation.costheta;
outputsinangle = experience.Observation.sintheta;
outputx = experience.Observation.x;
outputvelocity = experience.Observation.velocity;
outputthetadot = experience.Observation.thetadot;

figure(20); plot(outputsinangle); grid on;

figure(21); plot(outputcosangle); grid on;

figure(22); plot(outputx); grid on;

figure(23); plot(outputvelocity); grid on;

figure(23); plot(outputthetadot); grid on; 
