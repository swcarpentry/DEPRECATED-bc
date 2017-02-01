function tests = test_F03_fitQuadModel_003
    % TEST_FITQUADMODEL Main test function for the algorithm implemented in
    % fitQuadModel. To run the tests defined here, use RUNTESTS.
    % See also F02_fitQuadModel, runtests.
    
    % Test array constructed from local functions in this file.
    tests = functiontests( localfunctions() );
    
end % test_fitQuadModel

function setupOnce(testCase) %#ok<*DEFNU>
    
    addpath('../Test_Data')
    testCase.TestData = load('fitQuadModel_TestData.mat');
    testCase.TestData.currentRNG = rng;
    
end % setupOnce

function test_nargin(testCase)
    
    
end % test_nargin

function test_invalidInputs(testCase)    
    
    
end % test_invalidInputs

function test_validInputs(testCase)
   
    
end % test_validInputs

function test_basicFit(testCase)
    
    
end % test_basicFit

function test_allNaNs(testCase)
    
        
end % test_allNaNs

function test_deficientRank(testCase)
       
    
end % test_deficientRank

function teardownOnce(testCase)
    
    % Remove the test directory from the path.
    rmpath('../Test_Data')
    
    % Restore the random number generator settings.
    s = testCase.TestData.currentRNG;
    rng(s)
    
end % teardownOnce