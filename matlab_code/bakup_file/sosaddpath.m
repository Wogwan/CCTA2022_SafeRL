function sosaddpath

% Add multipoly

cm = computer;
k = 1;

if cm(1) == 'M' ||  cm(1)=='G'
    set(0,'DefaultFigureWindowStyle','docked')
    % Add chebfun-master
    addpath([pwd '/toolbox/chebfun-master']);

    % Add multipoly
    addpath([pwd '/toolbox/multipoly']);
    
    % Add nlanal
    addpath([pwd '/toolbox/nlanal']);
    
    % Add my version of SOSTools
    addpath([pwd '/toolbox/sosopt']);
    addpath([pwd '/toolbox/sosopt/Demos']);
    
    % Add polysys
    addpath([pwd '/toolbox/polysystems_1_0_3'])
    addpath([pwd '/toolbox/polysystems_1_0_3/demo'])
    
    % Add utils
    addpath([pwd '/utils'])
    
    % Add updates
    %     addpath([pwd '/updates'])

    
    % Add demo_3D
    %     addpath([pwd '/demo_3D'])
    %     addpath([pwd '/demo_3D/paper_example'])
    
    fprintf('Adding Path');
    
elseif cm(1) == 'P'
    set(0,'DefaultFigureWindowStyle','docked')
    %     set(0,'DefaultFigureWindowStyle','normal')
    % Add chebfun-master
    addpath([pwd '/toolbox/chebfun-master']);
    
    % Add gpml-matlab-master
%     addpath([pwd '/toolbox/gpml']);
%     path_gpml = '/toolbox/gpml/startup.m';
%     path_gp = append(pwd,path_gpml);
%     run path_gp
    run('/home/wogwan/Documents/safeRL/ASCC2022_SafeRL-dev/DDPG-TF2/toolbox/gpml/startup.m')
    
    % Add multipoly
    addpath([pwd '/toolbox/multipoly']);
    
    % Add nlanal
    addpath([pwd '/toolbox/nlanal']);
    
    % Add my version of SOSTools
    addpath([pwd '/toolbox/sosopt']);
    addpath([pwd '/toolbox/sosopt/Demos']);
    
    % Add polysys
    addpath([pwd '/toolbox/polysystems_1_0_3'])
    addpath([pwd '/toolbox/polysystems_1_0_3/demo'])
    
    % Add utils
    addpath([pwd '/utils'])
    
    % Add updates
    %     addpath([pwd '/updates'])

    fprintf('Adding Path');
    % Add demo_3D
    %     addpath([pwd '/demo_3D'])
    %     addpath([pwd '/demo_3D/paper_example'])
%     savepath
end