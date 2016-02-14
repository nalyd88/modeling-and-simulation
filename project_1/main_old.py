#!/usr/bin/env python
# This file is the main script for the simulations.
# The program simulates the flow of people leaving a stadium and heading
# to their desired locations. The simulation is based on Cellular Automata.

#import matplotlib
#matplotlib.use('macosx')
from enum import IntEnum
from pylab import *
import time


class CellType(IntEnum):
    """Define cell types with an enum."""
    Sidewalk = 1,
    Street = 2,
    Destination = 3,
    Start = 4,
    Prohibited = 0


class Cell(object):
    """This class is used to represent a walking_grid cell."""

    def __init__(self, cell_type, cell_size, i, j):
        self.occupied = False
        self.occupant = None
        self.type = cell_type
        self.traffic = 0
        self.safe = False
        self.size = cell_size
        self.coordinate = (i, j)
        self.id = 0


class Grid(object):
    """This class is used to represent the simulation walking_grid."""

    def __init__(self):
        self.size = (0, 0)
        self.cells = None

    def build(self, map_file_path=None):
        if map_file_path is None:
            return

        with open(map_file_path, 'r') as fileID:

            # Read the header and get the number of nodes and edges.
            rows, cols, size = [int(x) for x in fileID.readline().split()]
            self.size = (rows, cols)

            # Create a walking_grid of prohibited cells
            self.cells = [[Cell(CellType.Prohibited, size, i, j)
                           for j in range(cols)] for i in range(rows)]

            for cellID, line in enumerate(fileID):
                sr, sc, er, ec, t = [int(x) for x in line.split()]
                for i in range(sr, er + 1):
                    for j in range(sc, ec + 1):
                        if t > int(self.cells[i][j].type):
                            if t == int(CellType.Sidewalk):
                                self.cells[i][j].type = CellType.Sidewalk
                            elif t == int(CellType.Street):
                                self.cells[i][j].type = CellType.Street
                            elif t == int(CellType.Destination):
                                self.cells[i][j].type = CellType.Destination
                            elif t == int(CellType.Start):
                                self.cells[i][j].type = CellType.Start
                            else:
                                self.cells[i][j].type = CellType.Prohibited
                        self.cells[i][j].id = cellID
                        self.cells[i][j].size = size


class Pedestrian(object):
    """This class represents a pedestrian traveling in the walking_grid."""

    def __init__(self, speed, destination):
        self.speed = speed
        self.destination = destination


class PedestrianExitSimulation(object):

    def __init__(self, map_file_path, num_peds):
        self.map_file_path = map_file_path
        self.grid = Grid()
        self.pedestrians = None
        self.num_peds = num_peds
        self.map_plot = None
        self.fig = None

    def initialize(self):
        """Initialize system state"""
        self.pedestrians = [Pedestrian(1, 1) for i in range(self.num_peds)]
        self.grid.build(self.map_file_path)

    def observe(self):
        """Visualize system state"""

        # Print state to the console
        for n in range(0, self.grid.size[0]):
            for m in range(0, self.grid.size[1]):
                if self.grid.cells[n][m].occupied:
                    print()

        plot_grid = zeros(self.grid.size, dtype=int)
        for n in range(0, self.grid.size[0]):
            for m in range(0, self.grid.size[1]):
                plot_grid[n, m] = int(self.grid.cells[n][m].type)

        if self.map_plot is None:
            self.fig = figure()
            self.map_plot = imshow(plot_grid, interpolation='nearest', origin='lower')

        else:
            self.map_plot.set_data(plot_grid)
        #ion()
        show()#block=False)
        #draw()
        print('Observe')
        print(plot_grid)


    def update(self):
        """Update system states for one discrete time step"""
        pass

    def run(self, num_iter):

        print("Begin Simulation")

        loop_update = 1.  # Loop update time (seconds)

        self.initialize()
        self.observe()

        for t in xrange(1, num_iter):
            start_time = time.time()
            self.update()
            self.observe()
            elapsed_time = time.time() - start_time
            if elapsed_time < loop_update:
                time.sleep(loop_update - elapsed_time)

        print("End Simulation")
        show()


if __name__ == "__main__":

    # Path to the map file
    map_path = "map2.txt"

    # Number of pedestrians to simulate
    pedestrian_count = 10

    # Create a simulator object.
    sim = PedestrianExitSimulation(map_path, pedestrian_count)

    # Number of iterations
    niter = 10

    # Run the simulation
    sim.run(niter)
