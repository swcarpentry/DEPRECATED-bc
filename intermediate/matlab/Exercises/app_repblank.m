%% Application code for REPBLANK.
 
fprintf('INPUT STRING:\n\n')
input_str = 'The quick brown fox jumped over the lazy dog';
fprintf('%s\n\n', input_str)

fprintf('OUTPUT STRING:\n\n')
output_str = repblank(input_str);
fprintf('%s\n\n', output_str)
