function predictorRankings = F04_findBestPredictors_V2()
    % FINDBESTPREDICTORS: Function ranking pairs of predictor variables
    % in the medical data based on their ability to correctly classify an
    % individual's BMI category. We want to find the best two predictors
    % (excluding height, weight, BMI and IDNum) for classifying an
    % individual's BMI category. The output of this function is a table
    % containing all pairs of predictor variables and their corresponding
    % predictive accuracy (%).
    
    % Step 1: Load and preprocess the data.
    if exist('S02_MedData.mat', 'file') ~= 2
        error('findBestPredictors:MissingDataFile', ...
            'Unable to locate S02_MedData.mat')
    end
    X = load('S02_MedData.mat');
    if ~isfield(X, 'MedData')
        error('findBestPredictors:MissingVariable', ...
            'Unable to find variable MedData in MAT-file.')
    end
    paredData = preprocess(X.MedData);
    
    % Step 2: Partition the data into training and test sets.
    [trainingIdx, testIdx, yTrain, yTest, Xtab] = partition(paredData);
    
    % Step 3: Assemble a list of distinct pairs of variable names.
    [varPairs, tabNames] = getpairs(Xtab);
    
    % Step 4: Measure the performance of each pair of predictor variables,
    % and store the results in a table.
    rng('default') % Reproducible results.
    [predictorRankings, bestPredictors, predictedBMIs] = ...
        measurePerformance(Xtab, varPairs, tabNames, ...
        trainingIdx, testIdx, yTrain, yTest);
    
    % Step 5: Visualise the results.
    visualiseResults(paredData, bestPredictors, ...
        predictedBMIs, predictorRankings, yTest, testIdx)
    
end % findBestPredictors

function paredData = preprocess(MedData)
    % PREPROCESS: Accepts the full medical data table as input, and returns
    % a preprocessed table as output. Certain predictor variables are
    % removed, missing data is standardised and categorical variables are
    % converted to their double equivalents (apart from the response BMI
    % classification variable).
    varsToRemove = {'IDNum', 'Height', 'Weight',...
        'BMI', 'ReportHeight', 'ReportWeight'};
    paredNames = setdiff(MedData.Properties.VariableNames, varsToRemove);
    paredData = MedData(:, paredNames);
    paredData = standardizeMissing(paredData, '<undefined>');
    tf = varfun(@iscell, paredData, 'OutputFormat', 'uniform') & ...
        ~strcmp(paredData.Properties.VariableNames, 'BMIClass');
    f = @(V) double(categorical(V));
    tmp = varfun(f, paredData(:, tf));
    idx = find(tf);
    nms = paredData.Properties.VariableNames(idx);
    for k = 1:numel(idx)
        paredData.(nms{k}) = tmp{:, k};
    end
    
    BMICats = {'Anorexic', 'Underweight', 'Optimal', ...
        'Overweight', 'Obese', 'Morbidly Obese'};
    paredData.BMIClass = categorical(paredData.BMIClass, BMICats, ...
        'Ordinal', true);
    
end % preprocess

function [trainingIdx, testIdx, yTrain, yTest, Xtab] = ...
        partition(paredData)
    % PARTITION: Function taking the reduced table data as input and
    % returning the following outputs:
    % trainingIdx - indices of training observations
    % testIdx - indices of test observations
    % yTrain - response training data, as a cell array of strings
    % yTest - response test data, as a cell array of strings
    % Xtab - table of remaining predictor variables after removing the
    % response.
    
    nObs = size(paredData, 1);
    m = floor(0.85*nObs);
    trainingIdx = 1:m;
    testIdx = (m+1):nObs;
    
    allButBMIClass = setdiff(paredData.Properties.VariableNames, 'BMIClass');
    yTrain = cellstr(paredData.BMIClass( trainingIdx ));
    yTest = cellstr(paredData.BMIClass( testIdx ));
    
    Xtab = paredData(:, allButBMIClass);
    
end % partition

function [varPairs, tabNames] = getPairs(Xtab) %#ok<DEFNU>
    % GETPAIRS: Assembles a nx2 cell array of pairs of variable names from
    % the input table Xtab. Also returns tabNames, a cell array of variable
    % names from the input table Xtab.
    
    tabNames = Xtab.Properties.VariableNames;
    varPairs = cell( nchoosek(size(Xtab, 2), 2), 2 );
    k = 1;
    for k1 = 2:size(Xtab, 2)
        for k2 = 1:k1-1
            varPairs{k, 1} = tabNames{k1};
            varPairs{k, 2} = tabNames{k2};
            k = k + 1;
        end
    end
    
end % getPairs

function [predictorRankings, bestPredictors, predictedBMIs] = ...
        measurePerformance(Xtab, varPairs, tabNames, ...
        trainingIdx, testIdx, yTrain, yTest)
    % MEASUREPERFORMANCE: Function measuring the performance of each pair
    % of predictor variables in classifying an individual's BMI class.
    % Input arguments are the preprocessed table with BMI class removed,
    % Xtab, the cell array of all pairs of variable names varPairs, and the
    % cell array of variable names tabNames. We also pass in the training
    % and test indices computed previously, as well as the response data.
    % Output arguments are as follows:
    % predictorRankings - table containing each pair of predictor variables
    % and the predictive accuracy of each pair (%).
    % bestPredictors - cell array containing the names of the two best
    % predictor variables.
    % predictedBMIs - cell array of strings containing the predicted BMI
    % classes from the two best predictor variables.
    
    X = table2array(Xtab);
    
    k = 1;
    PredictiveAccuracy = NaN(nchoosek(size(X, 2), 2), 1);
    for k1 = 1:size(X, 2)
        for k2 = 1:(k1-1)
            Xiter = X(trainingIdx, [k1, k2]);
            F = TreeBagger(5, Xiter, yTrain);
            predictedBMIs = predict(F, X(testIdx, [k1, k2]));
            matches = strcmp(predictedBMIs, yTest);
            predAcc = 100*nnz(matches)/numel(matches);
            fprintf('Accuracy for predictors %s and %s (%%): %.2f\n', ...
                tabNames{k1}, tabNames{k2}, ...
                predAcc)
            PredictiveAccuracy(k) = predAcc;
            k = k + 1;
        end
    end
    
    predictorRankings = table( varPairs(:, 1), varPairs(:, 2), PredictiveAccuracy );
    predictorRankings.Properties.VariableNames(1:2) = {'Predictor1', 'Predictor2'};
    predictorRankings = sortrows(predictorRankings, 'PredictiveAccuracy', 'descend');
    bestPredictors = predictorRankings{1, 1:2};
    bestPredictorsIdx = ismember(Xtab.Properties.VariableNames, bestPredictors);
    F = TreeBagger(500, X(trainingIdx, bestPredictorsIdx), yTrain );
    predictedBMIs = predict(F, X(testIdx, bestPredictorsIdx));
    
end % measurePerformance

function visualiseResults(paredData, bestPredictors, ...
        predictedBMIs, predictorRankings, yTest, testIdx)
    % VISUALISERESULTS: Show the results of the analysis.
    
    BMICats = categories(paredData.BMIClass);
    nCats = numel(BMICats);
    screenPos = get(0, 'ScreenSize');
    dx = screenPos(3); dy = screenPos(4);
    centrePos = [dx/4, dy/4, dx/2, dy/2];
    fh = figure('Renderer', 'Painters', 'Color', 'w', 'Position', centrePos);
    
    axesOpts = {'Parent', fh, ...
        'XTick', 1:nCats, 'XTickLabel', 1:nCats, ...
        'XLimMode', 'manual', 'XLim', [0, 7], 'NextPlot', 'add', ...
        'YLim', [-10, 320], 'YLimMode', 'manual', 'DrawMode', 'fast', 'Color', 0.9*ones(1, 3), ...
        'FontWeight', 'Bold', 'FontSize', 12, 'FontName', 'Century'};
    ah(1) = subplot(2, 2, 1, axesOpts{:} );
    ah(2) = subplot(2, 2, 3, axesOpts{:} );
    
    grid(ah(1), 'on')
    grid(ah(2), 'on')
    textOpts = {'FontSize', 12, 'FontWeight', 'Bold', 'FontName', 'Century'};
    th(1) = title(ah(1), 'Predicted', textOpts{:});
    th(2) = title(ah(2), 'Actual', textOpts{:});
    
    uthResults = uitable(fh, 'Data', NaN(numel(yTest), 2), ...
        'ColumnName', {'Predicted', 'Actual'}, ...
        'RowName', testIdx, ...
        'Units', 'Normalized', ...
        'Position', [0.55, 0.50, 0.3, 0.4]) ;
    
    uitable(fh, 'ColumnName', {'Category', 'BMI Class'}, ...
        'RowName', '', ...
        'Units', 'Normalized', ...
        'Position', [0.55, 0.15, 0.19, 0.22], ...
        'Data', [num2cell((1:6).'), BMICats], ...
        'ColumnWidth', {75, 100});
    
    BMIMap = containers.Map( BMICats, ...
        num2cell(1:nCats) );
    
    predAcc = predictorRankings{1, end};
    set(th(1), 'String', ...
        sprintf('Predictions from %s and %s\n Accuracy (%%): %.2f',...
        bestPredictors{:}, predAcc))
    
    binCountsActual = zeros(1, nCats);
    binCountsPredicted = zeros(1, nCats);
    
    actualPlotOpts = {'Marker', 'o', 'MarkerSize', 5, 'Color', 'b'};
    predictedPlotOpts = {'Marker', 'o', 'MarkerSize', 5, 'Color', 'g'};
    
    
    for k = 1:numel(yTest)
        
        actualClass = BMIMap( yTest{k} );
        plot(ah(2), actualClass, binCountsActual(actualClass), actualPlotOpts{:} )
        binCountsActual(actualClass) = binCountsActual(actualClass)+1;
        
        predictedClass = BMIMap( predictedBMIs{k} );
        if predictedClass ~= actualClass, predictedPlotOpts{end} = 'r'; end
        plot(ah(1), ...
            predictedClass, binCountsPredicted(predictedClass), predictedPlotOpts{:})
        binCountsPredicted(predictedClass) = binCountsPredicted(predictedClass)+1;
        predictedPlotOpts{end} = 'g';
        
        currentData = get(uthResults, 'Data');
        currentData(k, :) = [predictedClass, actualClass];
        set(uthResults, 'Data', currentData)
        
    end % for
    
end % visualiseResults