function output_list = repblank( varargin )
%REPBLANK Replaces blanks with underscores. This function receives any 
% string as input, and returns a modified string obtained by replacing all 
% white spaces " " by underscores "_". Leading and trailing white space 
% is removed rather than replaced with underscores. Multiple consecutive
% spacing is replaced by a single underscore.

% Check the number and type of input arguments
if nargin ~= 1
    
    error('repblank:IncorrectNumberInputs','Incorrect number of input arguments.')
    
else
    
    str = varargin{:};
    
    if iscell(str)
        
        nstrings = length(str);
        output_list = cell(nstrings,1);
        for i=1:nstrings
            
            output_list{i} = str_repblank(str{i});
            
        end
        
    elseif ischar(str)
        
        output_list = str_repblank(str);
        
    else
        
        error('repblank:InvalidDataType','Invalid data type.')
        
    end
    
end

end


%-------------------------------------------%
function output_str = str_repblank(input_str)

idx = isspace(input_str);

if all(idx)
    error('repblank:AllBlankString', ...
        'Input string cannot be all blanks.')
end

while idx(1)
    idx(1)=[];
    input_str(1)=[];
end

while idx(end)
    idx(end)=[];
    input_str(end)=[];
end

idx1 = find(idx);
idx2 = find(diff(idx1)==1)+1;
iDel = idx1(idx2);
input_str(iDel) = [];
idx(iDel) = [];

output_str = input_str;
output_str(idx) = '_';

end



