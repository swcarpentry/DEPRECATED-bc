function [ W , T ] = getwords( varargin )
    % GETWORDS: Performs a frequency count on a given text file. The first
    % output W is a column vector of words contained in the input text
    % file. The second output T is a table with two columns: column one
    % contains the unique words from the text file, and column two contains
    % the corresponding frequency counts. The table is sorted in descending
    % order of word frequencies.
    % Examples: add the Literature folder to the path, and then enter:
    % >> [W1, T1] = getwords('sherlock_holmes.txt');
    % >> [W2, T2] = getwords('don_quixote.txt');

if nargin ~= 1
    
    error('getwords:IncorrectNumberInputs','Incorrect number of input arguments.')
    
elseif ischar(varargin{:})
    
    
    try
        
        fileID = fopen(varargin{:});
        
        if fileID == -1
            error('getwords:noFile','no such file')
        end
        
        W = textscan(fileID, '%s');
        
        if fileID ~=-1
            fclose(fileID);
        end
        
        W = W{1};
        W = regexprep(W, '[^A-Za-z0-9]', '');
        L = cellfun(@isempty, W);
        W(L) = [];
        W = lower(W);
        
        vocabulary = unique(W);
        n = length(vocabulary);
        freq = zeros(n,1);
        for i=1:n
            freq(i) = nnz(ismember(W,vocabulary{i}));
        end
        [freq, idx] = sort(freq, 'descend');
        vocabulary = vocabulary(idx);
        T = table(vocabulary,freq,'VariableNames',{'Vocabulary' 'Freq'});
        
    catch Ex
        rethrow(Ex);
    end
    
    
else
    
    error('getwords:InvalidDataType','Invalid data type.')
    
end


