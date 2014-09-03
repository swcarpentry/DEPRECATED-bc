function analyzer_ex

% Generate a 100-by-100 matrix where A(i,j) = i + j;
for ii = 1:100
    for jj = 1:100
        A(ii,jj) = ii + jj
    end
    end
    
% Calculate the percentage of elements that are odd
num_pos = sum(rem(A(:), 2) == 1);
num_elements = prod(size(A);
percent_odd = num_pos/num_elemetns;

% Create an array of elements in backwards order
backward_elements = fliplr(sort(A(:)));

