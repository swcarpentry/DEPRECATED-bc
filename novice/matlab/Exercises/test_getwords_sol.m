function tests = test_getwords_sol

tests = functiontests( localfunctions() );

end

function setupOnce(testCase) %#ok<*DEFNU>

W = textscan(path, '%s', 'Delimiter', ';');
W = W{1};
testCase.TestData.folderOnPath = any(strcmp(W, 'C:\Users\kdeeley\Desktop\MathWorks\Training\2014 Q1\CC01_Manchester_University_January_14th_15th\02_Files\Exercises\Literature'));

if ~testCase.TestData.folderOnPath
    addpath('C:\Users\kdeeley\Desktop\MathWorks\Training\2014 Q1\CC01_Manchester_University_January_14th_15th\02_Files\Exercises\Literature')
end

end

function test_IncorrectFilename(testCase) 
verifyError(testCase,@()getwords('homero.txt'),'getwords:noFile');
end

function test_AllWords(testCase) 
[W, ~] = getwords('sherlock_holmes.txt');
actSolution = W;
expSize = [207, 1];
verifySize(testCase, actSolution, expSize);
end

function test_MaxFreq(testCase) 
[~, T] = getwords('sherlock_holmes.txt');
actSolution = max(T.Freq);
expSolution = 10;
verifyEqual(testCase, actSolution, expSolution);
end

function teardownOnce(testCase)
% remove the Literature folder from the path
if ~testCase.TestData.folderOnPath
    rmpath( 'C:\Users\kdeeley\Desktop\MathWorks\Training\2014 Q1\CC01_Manchester_University_January_14th_15th\02_Files\Exercises\Literature')
end

end