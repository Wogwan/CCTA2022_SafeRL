%% Set parameters
% mdl = 'rlSimplePendulumModel';
mdl = 'rlDDPGmdl';
open_system(mdl)
useGPU = 1;

%%
env = rlPredefinedEnv('SimplePendulumModel-Continuous');
numObs = 3;
set_param('rlSimplePendulumModel/create observations','ThetaObservationHandling','sincos');
env.ResetFcn = @(in)setVariable(in,'theta0',pi,'Workspace',mdl);
Ts = 0.05;
Tf = 20;
rng(0);

statePath = [
    featureInputLayer(numObs,'Normalization','none','Name','observation')
    fullyConnectedLayer(400,'Name','CriticStateFC1')
    reluLayer('Name', 'CriticRelu1')
    fullyConnectedLayer(300,'Name','CriticStateFC2')];
actionPath = [
    featureInputLayer(1,'Normalization','none','Name','action')
    fullyConnectedLayer(300,'Name','CriticActionFC1','BiasLearnRateFactor',0)];
commonPath = [
    additionLayer(2,'Name','add')
    reluLayer('Name','CriticCommonRelu')
    fullyConnectedLayer(1,'Name','CriticOutput')];

criticNetwork = layerGraph();
criticNetwork = addLayers(criticNetwork,statePath);
criticNetwork = addLayers(criticNetwork,actionPath);
criticNetwork = addLayers(criticNetwork,commonPath);
% figure(1)
% plot(criticNetwork)
criticNetwork = connectLayers(criticNetwork,'CriticStateFC2','add/in1');
criticNetwork = connectLayers(criticNetwork,'CriticActionFC1','add/in2');
% figure(2)
% plot(criticNetwork)

criticOpts = rlRepresentationOptions('LearnRate',1e-03,'GradientThreshold',1);
if useGPU
    criticOpts.UseDevice = "gpu";
end

obsInfo = getObservationInfo(env);
actInfo = getActionInfo(env);
critic = rlQValueRepresentation(criticNetwork,obsInfo,actInfo,'Observation',{'observation'},'Action',{'action'},criticOpts);


actorNetwork = [
    featureInputLayer(numObs,'Normalization','none','Name','observation')
    fullyConnectedLayer(400,'Name','ActorFC1')
    reluLayer('Name','ActorRelu1')
    fullyConnectedLayer(300,'Name','ActorFC2')
    reluLayer('Name','ActorRelu2')
    fullyConnectedLayer(1,'Name','ActorFC3')
    tanhLayer('Name','ActorTanh')
    scalingLayer('Name','ActorScaling','Scale',max(actInfo.UpperLimit))];

actorNetwork_plot = [
    featureInputLayer(numObs,'Normalization','none','Name','observation')
    fullyConnectedLayer(400,'Name','ActorFC1')
    reluLayer('Name','ActorRelu1')
    fullyConnectedLayer(300,'Name','ActorFC2')
    reluLayer('Name','ActorRelu2')
    fullyConnectedLayer(1,'Name','ActorFC3')
    tanhLayer('Name','ActorTanh')
    scalingLayer('Name','ActorScaling','Scale',max(actInfo.UpperLimit))];
actorplotNetwork = layerGraph();
actorplotNetwork = addLayers(actorplotNetwork,actorNetwork_plot);

actorOpts = rlRepresentationOptions('LearnRate',1e-04,'GradientThreshold',1);

if useGPU
    actorOpts.UseDevice = "gpu";
end
actor = rlDeterministicActorRepresentation(actorNetwork,obsInfo,actInfo,'Observation',{'observation'},'Action',{'ActorScaling'},actorOpts);

agentOpts = rlDDPGAgentOptions(...
    'SampleTime',Ts,...
    'TargetSmoothFactor',1e-3,...
    'ExperienceBufferLength',1e6,...
    'DiscountFactor',0.99,...
    'MiniBatchSize',128);
agentOpts.NoiseOptions.StandardDeviation = 0.6;
agentOpts.NoiseOptions.StandardDeviationDecayRate = 1e-5;

agent = rlDDPGAgent(actor,critic,agentOpts);
maxepisodes = 5000;
maxsteps = ceil(Tf/Ts);
trainOpts = rlTrainingOptions(...
    'MaxEpisodes',maxepisodes,...
    'MaxStepsPerEpisode',maxsteps,...
    'ScoreAveragingWindowLength',5,...
    'Verbose',false,...
    'Plots','training-progress',...
    'StopTrainingCriteria','AverageReward',...
    'StopTrainingValue',-740,...
    'SaveAgentCriteria','EpisodeReward',...
    'SaveAgentValue',-740);

% figure(1);plot(actorplotNetwork);
% figure(2);plot(criticNetwork)

%% Do training
doTraining = false;
% doTraining = true;
if doTraining
    % Train the agent.
    trainingStats = train_sos(agent,env,trainOpts);
else
    % Load the pretrained agent for the example.
    load('SimulinkPendulumDDPG_test.mat','agent')
end

%% Run simulation
simOptions = rlSimulationOptions('MaxSteps',500);
experience = sim(env,agent,simOptions);