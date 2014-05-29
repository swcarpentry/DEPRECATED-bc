function tests = test_F03_fitQuadModel_004
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
    
    % Ensure that fitQuadModel throws an exception when called with 0, 1 or
    % 4 or more input arguments.
    verifyError(testCase, @() F02_fitQuadModel(),...
        'MATLAB:narginchk:notEnoughInputs')
    verifyError(testCase, @() F02_fitQuadModel(0),...
        'MATLAB:narginchk:notEnoughInputs')
    verifyError(testCase, @() F02_fitQuadModel(1, 2, 3, 4), ...
        'MATLAB:TooManyInputs')
    
end % test_nargin

function test_invalidInputs(testCase)
    
    % The first input argument, X, should be a real, 2D, nonempty double
    % matrix with no infinite values.
    % The second input argument, y, should be a real, nonempty column vector
    % of type double with no infinite values.
    % The third input argument, showplot, should be a logical scalar value.
    % X should have at least 3 rows and at most 2 columns.
    % The number of rows of X and y should coincide. 
    
    rng('default') % For reproducibility of tests.
    verifyError(testCase, @() F02_fitQuadModel(rand(10,0), rand(10,1)), ...
        'MATLAB:F02_fitQuadModel:expectedNonempty')
    verifyError(testCase, @() F02_fitQuadModel(rand(3,2), rand(3,1), true(1,2)), ...
        'MATLAB:F02_fitQuadModel:expectedScalar')
    verifyError(testCase, @() F02_fitQuadModel(rand(2,2), rand(2,1)), ...
        'F02_fitQuadModel:XTooSmall')
    verifyError(testCase, @() F02_fitQuadModel(rand(3,4), rand(3,1)), ...
        'F02_fitQuadModel:TooManyXCols')
    verifyError(testCase, @() F02_fitQuadModel(rand(3,2), rand(2,1)), ...
        'F02_fitQuadModel:DimMismatch')
    verifyError(testCase, @() F02_fitQuadModel(rand(10, 1), Inf(10, 1)), ...
        'F02_fitQuadModel:expectedNoninfinite')
    verifyError(testCase, @() F02_fitQuadModel(Inf(10, 1), rand(10, 1)), ...
        'F02_fitQuadModel:expectedNoninfinite')
    
    
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