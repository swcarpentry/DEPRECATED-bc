def border(grid, color):
    assert grid.width > 1, 'Must have at least two columns to draw border.'
    assert grid.height > 1, 'Must have at least two rows to draw border.'

    grid[0,  :] = color
    grid[-1, :] = color
    grid[:,  0] = color
    grid[:, -1] = color
