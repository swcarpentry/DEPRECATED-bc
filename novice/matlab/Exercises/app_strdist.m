%% Application code for STRDIST.

%% TEST CHAR s
fprintf('TEST CHAR s:\n')
s = 'juan';
d = strdist(s);
disp(d)

%% TEST CELL ARRAY s
fprintf('TEST CELL ARRAY s:\n')
s = {'juan', '?uan'};
d = strdist(s);
disp(d)

%% TEST CHAR s1 AND CELL ARRAY s2
fprintf('TEST CHAR s1 AND CELL ARRAY s2:\n')
s1 = 'juan';
s2 = {'juan', 'john'};
d = strdist(s1,s2);
disp(d)

%% TEST CELL ARRAY s1 AND CELL ARRAY s2 
fprintf('TEST CELL ARRAY s1 AND CELL ARRAY s2:\n')
s1 = {'jos', 'juan', 'greg', 'ken'};
s2 = {'jos', 'juan', 'mike', 'aleksandra'};
d = strdist(s1,s2);
disp(d)