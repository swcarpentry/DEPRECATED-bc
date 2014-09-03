function tests = test_strdist_sol
    
    tests = functiontests( localfunctions() );
    
end

function setupOnce(testCase) %#ok<*DEFNU>
    
    % Add the test data folder to the path.
    addpath('STRDIST_Test_Data')
    
    % Load data from the words.mat file.
    testCase.TestData = load('words.mat');
end

function test_OneWord(testCase)
    w = testCase.TestData.words{1};
    actSolution = strdist(w);
    expSolution = 0;
    verifyEqual(testCase, actSolution, expSolution);
end

function test_ThreeWords(testCase)
    w = testCase.TestData.words;
    actSolution = strdist(w(1:3));
    expSolution = [7,5,7];
    verifyEqual(testCase, actSolution, expSolution);
end

function test_DimensionsOutput(testCase)
    w = testCase.TestData.words;
    d = strdist(w(1:4));
    actSolution = size(d);
    expSolution = [1 6];
    verifyEqual(testCase, actSolution, expSolution);
end

function teardownOnce(~)
    % remove the folder from the path
    rmpath('STRDIST_Test_Data');
end