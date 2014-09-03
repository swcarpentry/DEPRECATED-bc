function d = strdist( varargin )
    % STRDIST: Computes the Levenshtein distance between two strings.
    % Inputs can be two strings, one string and a cell array of strings,
    % or two cell arrays containing strings.

if nargin == 1
    
    s = varargin{1};
    
    if ischar(s)
        
        d = 0;
        
    elseif iscell(s)
        
        n = length(s);
        b = nchoosek(1:n,2);
        d = zeros(1,size(b,1));
        for i=1:size(b,1)
            d(i) = levdist(s{b(i,1)},s{b(i,2)});
        end
        
    else
        
        error('strdist:InvalidDataType','Invalid data type.')
        
    end
    
elseif nargin == 2
    
    s1 = varargin{1};
    s2 = varargin{2};
        
    if ischar(s1) && ischar(s2)
        
        d = levdist(s1, s2);
        
    elseif ischar(s1) && iscell(s2)
        
        nstrings = length(s2);
        d = zeros(1,nstrings);
        for i=1:nstrings
            is2 = s2{i};
            d(i) = levdist(s1,is2);
        end
        
    elseif iscell(s1) && ischar(s2)
        
        n = length(s1);
        d = zeros(1,n);
        for i=1:n
            is1 = s1{i};
            d(i) = levdist(is1,s2);
        end
        
    elseif iscell(s1) && iscell(s2)
        
        n1 = length(s1);
        n2 = length(s2);
        d = zeros(n1,n2);
        for i=1:n1
            is1 = s1{i};
            for j=1:n2
                is2 = s2{j};
                d(i,j) = levdist(is1,is2);
            end
        end
    else
        
        error('strdist:InvalidDataType','Invalid data type.')
        
    end
    
else
    
    error('strdist:InvalidDataType','Incorrect number of input arguments.')
    
end


%-------------------------------%
function y = levdist(s1, s2)
% Calculates the Levenshtein distance

mtx = (1:length(s1)+1)'*(1:length(s2)+1)-1; 

for n=1:length(s1) 
    for k=1:length(s2) 
        temp = mtx(n,k); 
        if not(s1(n) == s2(k)) 
            temp = min([temp,mtx(n,k+1),mtx(n+1,k)])+1;
        end;
        mtx(n+1,k+1) = temp; 
    end;
end;

y = mtx(end,end); 