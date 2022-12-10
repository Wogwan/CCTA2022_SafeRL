function T = pevalIndex2(poly,vars)


% Get the cell arrays of the variable names.
polyNames = poly.varname;
varNames  = vars.varname;

% Determine where vars are used in poly
[logicalIndex,locationIndex] = ismember(varNames,polyNames);

% Change logical index into a subscript index.
nVars = length(vars);
varList = (1:nVars)';
index = varList(logicalIndex);

% Create the transformation matrix
n = length(index);
try
  T = full(sparse( locationIndex(logicalIndex),index,ones(n,1),length(polyNames),nVars));
catch
  keyboard
end