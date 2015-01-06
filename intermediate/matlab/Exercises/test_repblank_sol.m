function tests = test_repblank_sol

tests = functiontests( localfunctions() );

end

function test_OneBlank(testCase) %#ok<*DEFNU>
actSolution = repblank('this is');
expSolution = 'this_is';
verifyEqual(testCase, actSolution, expSolution);
end

function test_ManyBlank(testCase)
actSolution = repblank('this  is');
expSolution = 'this_is';
verifyEqual(testCase,actSolution,expSolution);
end

function test_FirstBlank(testCase)
actSolution = repblank(' this');
expSolution = '_this';
verifyNotEqual(testCase,actSolution,expSolution);
end

function test_LastBlank(testCase)
actSolution = repblank('this is  ');
expSolution = 'this_is';
verifyEqual(testCase,actSolution,expSolution);
end

function test_AllBlanksError(testCase)
    
    verifyError(testCase, @() repblank('   '), 'repblank:AllBlankString')
    
end



