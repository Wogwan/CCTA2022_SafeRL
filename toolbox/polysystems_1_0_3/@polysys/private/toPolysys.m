function sys = toPolysys(obj)

if isempty(obj)
    sys = polysys();
    return;
end

%% Input is an LTI system.
if isa(obj,'lti')
    
    
    [A,B,C,D,Ts] = ssdata(obj);

    for i = 1:size(A,1)
        states(i,1) = polynomial(1,1,{['x',num2str(i,1)]},[1 1]);
    end
    for j = 1:size(B,2)
        inputs(j,1) = polynomial(1,1,{['u',num2str(j,1)]},[1 1]);
    end

    stMap = A*states + B*inputs;
    orMap = C*states + D*inputs;
    name = get(obj,'Name');
    if Ts == 0
        sampleTime = [];
    else
        sampleTime = Ts;
    end

    sys = polysys(stMap,orMap,states,inputs,sampleTime);    
    sys = set(sys,'Name',name);

    

%% Input is an array of numbers.
elseif isnumeric(obj)
    
    % Else, create a static map.
    if ndims(obj) > 2
        error('POLYSYS:ndims',...
            'Can only convert two dimensional arrays to polysys objects.')
    else
        obj = double(obj);
        inputs = polynomial();
        for i = 1:size(obj,2);
            inputs(i,1) = polynomial(1,1,{['u',num2str(i,1)]},[1 1]);
        end

        orMap = polynomial();
        for j = 1:size(obj,1)
            orMap(j,1) = obj(j,:)*inputs;
        end
        sys = polysys([],orMap,[],inputs);
    end


%% Input is a polynomial object.
elseif isa(obj,'polynomial')
    
    if size(obj,2) == 1

        varNames = obj.varname;
        if isempty(varNames)
            sys = toPolysys(double(obj));
        else
            nInputs = length(varNames);
            inputCoeff  = eye(nInputs);
            inputDegMat = eye(nInputs);
            inputMatDim = [nInputs, 1];
            inputs = polynomial(inputCoeff,inputDegMat,varNames,inputMatDim);
            sys = polysys([],obj,[],inputs,'');
        end
    else
        error('POLYSYS:orMapWide',...
            'Static maps must be in column vector form.')
    end    
  
    
%% Cannot convert input.
else
    error('POLYSYS:inputClass',...
        'Cannot convert class ''%s'' to a polysys object.', class(obj));
end
