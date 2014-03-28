function tests = test_F03_fitQuadModel
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
    
    rng('default') % For reproducibility of test results.
    verifyInstanceOf(testCase, F02_fitQuadModel(rand(30,1), rand(30,1)), 'double')
    verifyInstanceOf(testCase, F02_fitQuadModel(rand(50,2), rand(50,1)), 'double')
    
    c = F02_fitQuadModel(rand(100,1), rand(100,1));
    verifySize(testCase, c, [3, 1])
    c = F02_fitQuadModel(rand(200,2), rand(200,1));
    verifySize(testCase, c, [6, 1])
    
end % test_validInputs

function test_basicFit(testCase)
    
    % Test the algorithm on data for which we know the exact solution.
    x = testCase.TestData.X1;
    y = testCase.TestData.y1;
    cExpected = testCase.TestData.c1;
    cActual = F02_fitQuadModel(x, y);
    verifyEqual(testCase, cActual, cExpected, 'AbsTol', 1e-10)
    
    x = testCase.TestData.X2;
    y = testCase.TestData.y2;
    cExpected = testCase.TestData.c2;
    cActual = F02_fitQuadModel(x, y);
    verifyEqual(testCase, cActual, cExpected, 'AbsTol', 1e-10)
    
    % Test the function against another solution method (LINSOLVE).
    rng('default')
    x1 = rand(500, 1); x2 = rand(500, 1);
    y = rand(500, 1);
    A = [ones(size(x1)), x1, x2, x1.^2, x2.^2, x1.*x2];
    c_linsolve = linsolve(A, y);
    cActual = F02_fitQuadModel([x1, x2], y);
    verifyEqual(testCase, cActual, c_linsolve, 'AbsTol', 1e-10)
    
    % Test the function against another solution method (POLYFIT).
    rng('default')
    x = rand(500, 1);
    y = 15 + 4*x - 7*x.^2 + randn(size(x));
    c_polyfit = polyfit(x, y, 2);
    c_polyfit = flipud( c_polyfit.' );
    cActual = F02_fitQuadModel(x, y);
    verifyEqual(testCase, cActual, c_polyfit, 'AbsTol', 1e-10)
    
    
end % test_basicFit

function test_allNaNs(testCase)
    
    rng('default') % Reproducibility.
    x = rand(10, 1);
    y = NaN(10, 1);
    verifyError(testCase, @() F02_fitQuadModel(x, y), 'F02_fitQuadModel:AllNaNs')
    
end % test_allNaNs

function test_deficientRank(testCase)
    
    % Test in the case when X has two columns.
    x = (1:10).'; X = [x, 2*x];
    y = (11:20).';
    verifyError(testCase, @() F02_fitQuadModel(X, y), ...
        'F02_fitQuadModel:fitModel:RankDeficientDesignMatrix')
    
    % Test in the case when X has one column.
    x = zeros(10, 1);
    y = (1:10).';
    verifyError(testCase, @() F02_fitQuadModel(x, y), ...
        'F02_fitQuadModel:fitModel:RankDeficientDesignMatrix')
    
end % test_deficientRank

function teardownOnce(testCase)
    
    % Remove the test directory from the path.
    rmpath('../Test_Data')
    
    % Restore the random number generator settings.
    s = testCase.TestData.currentRNG;
    rng(s)
    
end % teardownOnce