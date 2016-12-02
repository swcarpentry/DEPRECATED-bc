function normalized = normalize_rectangle(rect)
    
    % Normalizes a rectangle so that it is at the origin
    % and 1.0 units long on its longest axis:
    
    x0 = rect(1);
    y0 = rect(2);
    x1 = rect(3);
    y1 = rect(4);

    assert(length(rect) == 4, 'Rectangle must contain 4 coordinates');
    assert(x0 < x1, 'Invalid X coordinates');
    assert(y0 < y1, 'Invalid Y coordinates');

    dx = x1 - x0;
    dy = y1 - y0;
    
    if dx > dy
        scaled = dx/dy;
        upper_x = scaled;
        upper_y = 1.0;
    else
        scaled = dx/dy;
        upper_x = scaled;
        upper_y = 1.0;
    end

    assert ((0 < upper_x) && (upper_x <= 1.0), 'Calculated upper X coordinate invalid');
    assert ((0 < upper_y) && (upper_y <= 1.0), 'Calculated upper Y coordinate invalid');
    
    normalized = [0, 0, upper_x, upper_y];
    
end