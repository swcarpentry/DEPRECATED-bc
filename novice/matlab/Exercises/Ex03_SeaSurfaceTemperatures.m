%% Load the data.
load('Ex03_SST.mat')

%% Visualise the data for the first time observation.
figure
scatter(lon, lat, 5, sst(:, 1))
xlabel('Longitude')
ylabel('Latitude')

figure
scatter3(lon, lat, sst(:, 1), 'k.')
xlabel('Longitude')
ylabel('Latitude')
zlabel('Surface Sea Temperature')

%% Fit a surface to the first set of data.
x = linspace(min(lon), max(lon));
y = linspace(min(lat), max(lat));
[X, Y] = meshgrid(x, y);
Z = griddata(lon, lat, sst(:, 1), X, Y);
figure
surf(X, Y, Z, 'FaceColor', 'interp', 'EdgeAlpha', 0)
hold on
scatter3(lon, lat, sst(:, 1), 'k.')
hold off
xlabel('Longitude')
ylabel('Latitude')
zlabel('Surface Sea Temperature')

%% Call the function for each set of temperature data.
figure
for k = 1:size(sst, 2)
    subplot(4, 6, k)
    [X, Y, Z] = SST_grid_sol(lon, lat, sst(:, k));
    surf(X, Y, Z, 'FaceColor', 'interp', 'EdgeAlpha', 0)
    hold on
    scatter3(lon, lat, sst(:, k), 'k.')
    hold off
    title(sprintf('Data series %d', k))
end