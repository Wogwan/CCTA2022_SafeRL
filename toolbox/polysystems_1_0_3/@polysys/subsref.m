function value = subsref(sys,index)
%SUBSREF  Subscript reference into a POLYSYS model.

switch index(1).type
    
    case '.'
        value = sys.(index(1).subs);
    case '()'
        switch length(index(1).subs)
            case 1
                error('POLYSYS:subsref:singleIndex', ...
                    'Single index referencing is not allowed.')
            case 2
                
                outIndex = index(1).subs{1};
                inIndex = index(1).subs{2};
                value = subsystem(sys,outIndex,inIndex);
               
            otherwise
                error('POLYSYS:subsref:multipleIndex', ...
                    'Arrays of polysys objects are not supported.')
        end
                
    case '{}'
        error('POLYSYS:subsref:cellRef', ...
            'Subscript cell array referencing is not allowed.')
    otherwise
        error('POLYSYS:subsref:unknownType', ...
            'Invalid subscript referencing.')
end

if length(index) > 1
    value = subsref(value,index(2:end));
end