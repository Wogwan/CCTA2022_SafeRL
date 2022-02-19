function reply = isMonoVector(vec)
%reply = isMonoVector(vec)

% 7.2X.2007: TJW - Initial coding.
% 8.07.2007: TJW - Improved performace by doing unique() on varname, rather 
%                  than on polynomial vector itself.
%                - Instead of testing row- and col-sum, code now tests for
%                  all ones and zeros, and unitary.

if isempty(vec)
    reply = true(1);
elseif isa(vec,'polynomial')

    n = length(vec);
    haveUnitCoeffs = local_isPermOfIdentity(vec.coefficient,n);
    areMonomials   = local_isPermOfIdentity(vec.degmat,n);
    areUnique = length(unique(vec.varname)) == length(vec);
    
    reply = isvector(vec) && haveUnitCoeffs && areMonomials && areUnique;
    
else
    reply = false(1);
end



% A unitary matrix of all ones and zeros is a permutation of identity.
function reply = local_isPermOfIdentity(M,n)
I = eye(n);
isOnesAndZeros = isequal( M, double(logical(M)) );
isUnitary = isequal(I,M'*M) && isequal(M*M',I);
reply = isOnesAndZeros && isUnitary;
