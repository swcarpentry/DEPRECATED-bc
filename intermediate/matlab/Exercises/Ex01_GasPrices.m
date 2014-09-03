%% Import the data into a table.
G = readtable('gasprices.csv', 'Delimiter', ',', 'Headerlines', 4);

%% Extract the data for Japan, and compute mean and std.
JP = G.Japan;
JP_mean = mean(JP);
JP_std = std(JP);

%% Extract the data for Europe, and compute the mean European price in each year.
Europe = [G.France, G.Germany, G.Italy, G.UK];
annualEuroMeans = mean(Europe, 2);

%% Compute the return series for Europe.
Euro_Returns = log( Europe(2:end, :) ./ Europe(1:end-1, :) );

%% Compute and visualise correlation.
C = corrcoef(Euro_Returns);
figure
imagesc(C, [-1, 1])
set(gca, 'XTick', 1:4, 'XTickLabel', {'France', 'Germany', 'Italy', 'UK'}, ...
    'YTick', 1:4, 'YTickLabel', {'France', 'Germany', 'Italy', 'UK'})
colorbar


