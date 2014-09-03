%% Load the data.
load('Ex04_House.mat')

%% Extract the region names.
regions = pricesToWages.Properties.RowNames;

%% Identify Manchester in the data.
L = strcmp(regions, 'Manchester');
ManchesterData = pricesToWages{L, :};
Yrs = 1997:2012;
figure
plot(Yrs, ManchesterData, 'LineWidth', 2, 'Marker', '*')
title(sprintf('Data for %s', regions{L}), 'FontWeight', 'Bold')
xlabel('Year')
ylabel('House Price to Median Earnings Ratio')
grid

%% Test the function.
plotPrices_sol(find(L))