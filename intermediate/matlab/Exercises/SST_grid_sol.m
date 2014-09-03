function [lonGrid, latGrid, tempGrid] = SST_grid_sol(lon, lat, seaTemp)
    
    x = linspace(min(lon), max(lon));
    y = linspace(min(lat), max(lat));
    [lonGrid, latGrid] = meshgrid(x, y);
    tempGrid = griddata(lon, lat, seaTemp, lonGrid, latGrid);
    
end