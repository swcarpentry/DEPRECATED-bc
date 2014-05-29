%% Model fitting.

%% Load data.
load('S02_MedData.mat')

%% Fit a quadratic model for pulse vs. age.

% Clean data.
x = MedData.Age;
y = MedData.BPDiff;
missingIdx = isnan(x) | isnan(y);
xClean = x(~missingIdx);
yClean = y(~missingIdx);

% Fit model.
designMat = [ones(size(xClean)), xClean, xClean.^2];
modelCoeffs = designMat\yClean;
model = designMat*modelCoeffs;

% Visualise results.
figure
scatter(x, y, 'kx')
hold on
plot(xClean, model, 'r*')
hold off
xlabel('Age')
ylabel('Pulse pressure')

%% Fit a quadratic model for weight vs. height and waist.

% Clean data.
x1 = MedData.Height;
x2 = MedData.Waist;
y = MedData.Weight;
missingIdx = isnan(x1) | isnan(x2) | isnan(y);
x1Clean = x1(~missingIdx); 
x2Clean = x2(~missingIdx);
yClean = y(~missingIdx);

% Fit model.
designMat = [ones(size(x1Clean)), x1Clean, x2Clean, x1Clean.^2, ...
    x2Clean.^2, x1Clean.*x2Clean];
modelCoeffs = designMat\yClean;
model = designMat*modelCoeffs;

% Visualise results.
figure
scatter3(x1, x2, y, 'kx')
hold on
x1Vec = linspace(min(x1Clean), max(x1Clean));
x2Vec = linspace(min(x2Clean), max(x2Clean));
[X1Grid, X2Grid] = meshgrid(x1Vec, x2Vec);
modelFun = @(c, x1, x2) ...
    c(1) + c(2)*x1 + c(3)*x2 + c(4)*x1.^2 + ...
    c(5)*x2.^2 + c(6)*x1.*x2;
YGrid = modelFun(modelCoeffs, X1Grid, X2Grid);
surf(X1Grid, X2Grid, YGrid, 'FaceColor', 'interp', 'EdgeAlpha', 0) 
hold off
xlabel('Height (cm)'), ylabel('Waist (cm)'), zlabel('Weight (kg)')