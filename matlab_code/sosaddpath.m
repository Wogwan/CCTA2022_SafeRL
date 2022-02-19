function sosaddpath

cm = computer;

if cm(1) == 'M' ||  cm(1)=='G'
    set(0,'DefaultFigureWindowStyle','docked');
    % Add chebfun-master
    addpath '../toolbox/chebfun-master';
    
    % Add multipoly
    addpath '../toolbox/multipoly';
    
    % Add nlanal
    addpath '../toolbox/nlanal';
    
    % Add my version of SOSTools
   addpath '../toolbox/sosopt';
    addpath '../toolbox/sosopt/Demos';
    
    % Add polysys
    addpath '../toolbox/polysystems_1_0_3';
    
    % Add demo_3D
    %     addpath([pwd '/demo_3D'])
    %     addpath([pwd '/demo_3D/paper_example'])
    
    fprintf('Adding Path\n');
    
elseif cm(1) == 'P'
    set(0,'DefaultFigureWindowStyle','docked');

    % Add chebfun-master
    addpath '../toolbox/chebfun-master';
    
    % Add gpml-matlab-master
    addpath '../toolbox/gpml';
    run '../toolbox/gpml/startup.m';
    %     run('/home/wogwan/Documents/safeRL/ASCC2022_SafeRL-dev/DDPG-TF2/toolbox/gpml/startup.m')
    
    % Add multipoly
    addpath '../toolbox/multipoly';
    
    % Add nlanal
    addpath '../toolbox/nlanal';
    
    % Add my version of SOSTools
    addpath '../toolbox/sosopt';
    addpath '../toolbox/sosopt/Demos';
    
    % Add polysys
    addpath '../toolbox/polysystems_1_0_3';
    addpath '../toolbox/polysystems_1_0_3/demo';
    
    fprintf('Adding Path\n');

end