function display(sys)
%DISPLAY Display the properties of a POLYSYS model in the Command Window.

% 7.20.2007: TJW - Initial coding.
% 7.22.2007: TJW - Changed class name from twpds to polysys.


%% Collect necessary information.

% Get format spacing from command window
isLooseDisplay = not(isequal(get(0,'FormatSpacing'),'compact'));

% See what kind of polysys object this is.
isContinuous = isequal(sys.sampleTime,0);

if isLooseDisplay
    disp(' ')
end


%% Display a simple message for empty systems
if isempty(sys)
    disp('Empty polysys object.')
    if isLooseDisplay
        disp(' ')
    end
else

    
%% Display the type of system.
    if sys.isDynamic
        if isContinuous
            disp('Continuous-time polynomial dynamic system.')
            if isLooseDisplay
                disp(' ')
            end 
        else
            disp('Discrete-time polynomial dynamic system.')
            if isLooseDisplay
                disp(' ')
            end 
            if sys.sampleTime == -1
                disp('Inherited sampling time.')
            else
                disp(['Sampling time: ',num2str(sys.sampleTime)])
            end
        end
    else
        disp('Static polynomial map.')
        if isLooseDisplay
            disp(' ')
        end         
    end
    
   
    
%% Display each of the variables.
    if not(isempty(sys.states))
        stateString = 'States: ';
        stateNames = sys.states.varname;
        for i = 1:length(stateNames)
            stateString = [stateString,stateNames{i}];
            if i<length(stateNames)
                stateString = [stateString,','];
            end
        end
        disp(stateString);
    end
    if not(isempty(sys.inputs))
        inputString = 'Inputs: ';
        inputNames = sys.inputs.varname;
        for i = 1:length(inputNames)
            inputString = [inputString,inputNames{i}];
            if i<length(inputNames)
                inputString = [inputString,','];
            end
        end
        disp(inputString)
    end

%% Display the maps.

    if sys.isDynamic
        if isContinuous
            disp('State transition map is x''=f(x,u) where')
        else
            disp('State transition map is x(k+1)=f(x(k),u(k)) where')
        end
        for i = 1:length(sys.stMap)
            stTextCell = char(sys.stMap(i));
            disp(['  f',num2str(i,1),' = ',stTextCell{1}])
        end
    end
    if isContinuous
        disp('Output response map is y=g(x,u) where')
    else
        disp('Output response map is y(k)=g(x(k),u(k)) where')
    end
    for j = 1:length(sys.orMap)
        orTextCell = char(sys.orMap(j));
        disp(['  g',num2str(j,1),' = ',orTextCell{1}])
    end
    
    if isLooseDisplay
        disp(' ')
    end

end % end if isempty(sys)
