#!/usr/bin/env python
"""
Conway's game of life, final version.

Rules:
 * Any cell with fewer than two neighbors dies (underpopulation)
 * Any cell with more than three neighbors dies (overpopulation)
 * Any empty spot with three neighbors becomes a live cell (reproduction)

Check out the wikipedia article for more information, including interesting
starting populations.
"""


def conway(population, generations=100):
    """Runs Conway's game of life on an initial population."""
    population = set(population)
    for i in range(generations):
        population = evolve(population)
    return list(population) 


def evolve(population):
    """Evolves the population by one generation."""
    # Get a unique set of discrete cells that need to be checked
    active_cells = population | set([neighbor for p in population
                                    for neighbor in neighbors(p)])
    # For each cell in the set, test if it lives or dies
    new_population = set()
    for cell in active_cells:
        count = sum([neighbor in population for neighbor in neighbors(cell)])
        if count == 3 or (count == 2 and cell in population):
            new_population.add(cell)
    # Return the new surviving population
    return new_population


def neighbors(cell):
    """Returns the neighbors of a given cell."""
    x, y = cell
    return [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1), (x + 1, y + 1),
            (x + 1, y - 1), (x - 1, y + 1), (x - 1, y - 1)]

glider = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 2)]
print conway(glider)
