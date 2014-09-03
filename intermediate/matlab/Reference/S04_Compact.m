%% Vectorisation Examples.

%% Load data.
load('S02_MedData.mat')

%% Compute approximations to mean arterial pressure (MAP) and 
% body surface area (BSA).
%
% $$ MAP \approx P_{dias} + \frac{1}{3}\left(P_{sys}-P_{dias}\right)$$
%
% $$ BSA \approx 0.007184\times W^{0.425}\times H^{0.725} $$ 

% Extract variables.
P_dias = MedData.BPDias1;
P_sys = MedData.BPSyst1;
W = MedData.Weight;
H = MedData.Height;

% % Preallocate for the results.
% nObs = height(MedData);
% MAP = NaN(nObs, 1);
% BSA = NaN(nObs, 1);
% 
% % Compute the approximations.
% for k = 1:nObs
%     MAP(k) = P_dias(k) + (1/3)*(P_sys(k)-P_dias(k));
% end
% 
% for k = 1:nObs
%     BSA(k) = 0.007184*(W(k)^0.425)*(H(k)^0.725);
% end

MAP = P_dias + (1/3)*(P_sys-P_dias);
BSA = 0.007184*(W.^0.425).*(H.^0.725);

%% Compute the mean and std height of male individuals by Ethnicity.

% ethCats = unique(MedData.Ethnicity);
% males = strcmp(MedData.Sex, 'M');
% 
% heightStats = NaN(numel(ethCats), 2);
% 
% for k = 1:numel(ethCats)
%     L = males & strcmp(MedData.Ethnicity, ethCats{k});
%     heightStats(k, 1) = mean(MedData.Height(L));
%     heightStats(k, 2) = std(MedData.Height(L));    
% end
% 
% heightStatsTable = table(ethCats, heightStats);

requiredStatsFun = @(V) [mean(V), std(V)];
heightStats = varfun(requiredStatsFun, MedData, ...
    'InputVariables', 'Height', ...
    'GroupingVariables', {'Sex', 'Ethnicity'});

%% Operations on cell arrays.
load('S04_CellData.mat')
disp(LTW_cats)

% % Compute BMI statistics for individuals classified 
% % by their "LikeToWeigh" status.
requiredStats = {@numel, @median, @iqr, @skewness, @kurtosis};
% BMI_stats = NaN(numel(LTW_cats), numel(requiredStats));
% 
% for k1 = 1:numel(LTW_cats)
%     for k2 = 1:numel(requiredStats)
%         BMI_stats(k1, k2) = requiredStats{k2}( C{k1} );
%     end    
% end
% 
% BMI_stats = [[ {'Category/Statistic'}, requiredStats ] ; ...
%     LTW_cats.', num2cell(BMI_stats)];

requiredStatsFun = @(V) [numel(V), median(V), iqr(V), skewness(V), kurtosis(V)];
BMI_stats = cellfun(requiredStatsFun, C, 'UniformOutput', false);

BMI_stats = [ [ {'Category/Statistic'}, requiredStats]; ...
    LTW_cats.', num2cell( cat(1, BMI_stats{:}) )];
disp(BMI_stats)

%% Operations on structures.
load('S04_StructData.mat')
% 
% % Same calculations as above.
requiredStats = {@numel, @median, @iqr, @skewness, @kurtosis};
% S_fields = fieldnames(S);
% BMI_stats = NaN(numel(S_fields), numel(requiredStats));
% 
% for k1 = 1:numel(S_fields)
%     for k2 = 1:numel(requiredStats)
%         BMI_stats(k1, k2) = requiredStats{k2}( S.(S_fields{k1}) );
%     end
% end
% 
% % Format the results into a structure.
% for k = 1:numel(S_fields)
%     T.(S_fields{k}) = BMI_stats(k, :);
% end
% T.statNames = requiredStats;
% disp(T)

requiredStatsFun = @(V) [numel(V), median(V), iqr(V), skewness(V), kurtosis(V)];
T = structfun(requiredStatsFun, S, 'UniformOutput', false);
T.statNames = requiredStats;
disp(T)