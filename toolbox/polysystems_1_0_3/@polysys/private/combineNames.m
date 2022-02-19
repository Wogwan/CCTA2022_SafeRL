function name = combineNames(sys1,sys2)
%Private utility function for the POLYSYS class.

% 08.22.2007: TJW - Intial coding.

if isempty(sys1.name)
    name = sys2.name;
elseif isempty(sys2.name)
    name = sys1.name;
else
    name = [sys1.name,',',sys2.name];
end