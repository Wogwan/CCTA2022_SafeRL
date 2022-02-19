function [newPoly,keepIndex] = robustSubs(poly,old,new)
%Private utility function for the POLYSYS class.
%
%ROBUSTSUBS  Symbolic substitution for POLYNOMIAL objects.
%
%   NEWP = robustSubs(P,OLDVARS,NEWVALUE)  finds which (if any) elements of
%   OLDVARS are being used by P and then substitutes values from NEWVALUE.
%   OLDVARS is a column vector of class POLYNOMIAL where each row is a
%   monomial with coefficient 1.  NEWVALUE is the same size as OLDVARS and
%   may be a vector of POLYNOMIAL objects or numerical values.
%
%   Note: this utility avoids some errors in the polynomial/subs method,
%   and allows us to "try" to substitute values for OLDVARS that are not
%   used in P.  For example, a state-transition map may not use all 
%   inputs to a system.

% 7.21.2007: TJW - Initial coding.
% 7.22.2007: TJW - Help documentation and notes.


%% Verify that this is a one-to-one substitution.

if length(old) ~= length(new)
    error('Old variables and new variables must have the same length.')
end


%% See what kind of substitution we're dealing with.

if isnumeric(new)
    hasCoeff1 = false(1);
    areMonomials = true(1);
else
    newCoeff = full(new.coefficient);
    hasCoeff1 = all( sum(newCoeff,1) == ones(1,size(newCoeff,2)) )...
                && all( sum(newCoeff,2) == ones(size(newCoeff,1),1) );
    newDeg = full(new.degmat);
    areMonomials = all( sum(newDeg,1) == ones(1,size(newDeg,2)) )...
                   && all( sum(newDeg,2) == ones(size(newDeg,1),1) );
end

%% Do the substitution.
           
% If new is a simple vector of monomial variables...
if hasCoeff1 && areMonomials
    
    % Do a quick and simple substitution using varname field.
    polyVars = poly.varname;
    newPolyVars = polyVars;
    keepIndex = [];
    for i = 1:length(old)
        thisVarname = get(old(i),'varname');
        isUsed = strcmp(thisVarname{1},polyVars);
        if any(isUsed)
            keepIndex = [keepIndex;i];
            thisUsedVarname = get(new(i),'varname');
            newPolyVars{isUsed} = thisUsedVarname{1};
        end
    end
    newPoly = polynomial(poly.coefficient,poly.degmat,newPolyVars,size(poly));

% Else, user is doing a 'fancier' substitution...
else
    
    polyVars = poly.varname;
    keepIndex = [];
    for i = 1:length(old)
        thisOldVarname = get(old(i),'varname');
        if ismember(thisOldVarname{1},polyVars);
            keepIndex = [keepIndex;i];
        end
    end
    if not(isempty(keepIndex))
        oldCell = mat2cell(old(keepIndex),ones(length(old(keepIndex)),1),1)';
        newCell = mat2cell(new(keepIndex),ones(length(new(keepIndex)),1),1)';
        newPoly = subs(poly,oldCell,newCell);

        % Sometimes subs will transponse the result (why?).
        if size(newPoly') == size(poly)
            newPoly = newPoly';
        end
    else
        newPoly = poly;
    end
        
end
