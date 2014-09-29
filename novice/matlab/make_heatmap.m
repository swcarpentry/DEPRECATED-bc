
function heatmap = make_heatmap(data, low_band, high_band)

    % Make a 3-colored heatmap from
    % a 2D array of data.

    [height, width] = size(data);

    heatmap = zeros(height, width);
    center = mean(data(:));

    for y = 1:height
        for x = 1:width

            if data(y, x) > high_band*center
                heatmap(y, x) = 1;
            elseif data(y, x) < low_band*mean(data(:))
                heatmap(y, x) = -1;
            else
                heatmap(y, x) = 0;
            end
        end
    end

end
