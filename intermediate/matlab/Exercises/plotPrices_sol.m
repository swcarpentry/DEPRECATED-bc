function plotPrices_sol(k)
    
    S = load('Ex04_House.mat');
    Yrs = 1997:2012;
    figure
    plot(Yrs, S.pricesToWages{k, :}, 'LineWidth', 2, 'Marker', '*')
    xlabel('Year')
    ylabel('House Price to Median Earnings Ratio')
    title(sprintf('Data for %s', ...
        S.pricesToWages.Properties.RowNames{k}), 'FontWeight', 'Bold')
    grid
    
end