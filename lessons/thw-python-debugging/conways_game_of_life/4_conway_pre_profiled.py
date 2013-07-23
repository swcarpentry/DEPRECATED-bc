#!/usr/bin/env python
"""
Conway's game of life example, part four.

This code is functional, but it doesn't make use of sets in a scenario where 
they're the perfect data type.
"""


def conway(population, generations=100):
    """Runs Conway's game of life on an initial population."""
    for i in range(generations):
        population = evolve(population)
    return population


def evolve(population):
    """Evolves the population by one generation."""
    # Get a unique set of discrete cells that need to be checked
    active_cells = population[:]
    for cell in population:
        for neighbor in neighbors(cell):
            if neighbor not in active_cells:
                active_cells.append(neighbor)
    # For each cell in the set, test if it lives or dies
    new_population = []
    for cell in active_cells:
        count = sum(neighbor in population for neighbor in neighbors(cell))
        if count == 3 or (count == 2 and cell in population):
            if cell not in new_population:
                new_population.append(cell)
    # Return the new surviving population
    return new_population


def neighbors(cell):
    """Returns the neighbors of a given cell."""
    x, y = cell
    return [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1), (x + 1, y + 1),
            (x + 1, y - 1), (x - 1, y + 1), (x - 1, y - 1)]

glider = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 2)]
print conway(glider)

