function TrainingStatistics = train_sos(this,env,varargin)
% TRAININGSTATISTICS = TRAIN(AGENT,ENVIRONMENT,OPTIONS)
%
% Train the AGENT against the ENVIRONMENT with training options OPTIONS
% created from rlTrainingOptions. 
%
% For multi-agents, AGENT is an array of reinforcement learning agents, and
% ENVIRONMENT is a Simulink environment.
%
% TRAININGSTATISTICS are returned when training is finished detailing the
% history of the training progress. TRAININGSTATISTICS is a 1xN structure 
% with the following fields, where N is the number of trained agents.
%
%   EpisodeIndex        : Episode number.
%   EpisodeReward       : Cumulative episode reward.
%   EpisodeSteps        : Steps taken during episode.
%   AverageReward       : Average cumulative reward of episodes in the
%                         averaging window (specified by
%                         ScoreAveragingWindowLength in the training
%                         options).
%   TotalAgentSteps     : Total number of steps taken.
%   EpisodeQ0           : Critic estimate of the long-term reward at the
%                         initial conditions of the environment (only for
%                         agents with a critic e.g. DQN and DDPG).
%
% See also: rlTrainingOptions, <a href="matlab: help rl.env.Abstract\sim">sim</a>

% Copyright 2018-2020 The MathWorks, Inc.

    parser = inputParser();
    addRequired(parser,'Environment',@(x)isa(x,'rl.env.AbstractEnv'));
    addOptional(parser,'TrainingOptions',rlTrainingOptions(),@(x)isa(x,'rl.option.rlTrainingOptions'));
    parse(parser,env,varargin{:});

    trainingOptions = parser.Results.TrainingOptions;
    
	%% check compatibility with deployment settings
    if isdeployed && strcmpi(trainingOptions.Plots,'training-progress')
        error(message('rl:general:TrainingPlotNotDeployable'))
    end
	
    %% perform validations
    validateStopTrainingFunction(trainingOptions);
    validateSaveAgentFunction(trainingOptions);

	% validate agents, environment and training options
    rl.util.validateAgentsWithEnv(this,env,trainingOptions);

    % validate training compatibility
    if numel(this) == 1
        trainingOptions = validateAgentTrainingCompatibility(this, trainingOptions);
        
        % validate parallel options
        if trainingOptions.UseParallel
            parOpts = trainingOptions.ParallelizationOptions;
            switch parOpts.DataToSendFromWorkers
                case "experiences"
                    if ~isa(this,'rl.agent.mixin.ExperienceParallelTrainable')
                        error(message('rl:general:errParallelSendExpNotSupport'));
                    end
                case "gradients"
                    if ~isa(this,'rl.agent.mixin.GradientParallelTrainable')
                        error(message('rl:general:errParallelSendGradNotSupport'));
                    end
                otherwise
                    % should never be hit
                    error('unsupported parallel option')
            end
        end
    end
    
    %% create the training manager
    trainMgr = rl.train.TrainingManager(env,this,trainingOptions);
    clnup    = onCleanup(@() delete(trainMgr));

    %% run the training
    TrainingStatistics = run(trainMgr);
end