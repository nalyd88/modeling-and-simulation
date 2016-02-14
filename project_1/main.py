#!/usr/bin/env python
# This program simulates the flow of people leaving a starting point and heading
# to their desired locations. The simulation is based on Cellular Automata (CA).
# Dylan Crocker
# CSE 6730 Project 1
# Due: Feb. 26, 2016


from pylab import *
import collections
import random


# Globals -------------------------------------------------------------------- #
prohibited_id = 0   # Prohibited (non-walkway)
walkway_id = 1      # Designated walkway
street_id = 2       # Street (prohibited for walking)
crossing_id = 3     # Transient walkway (e.g., street crossing)
exit_id = 4         # Cells representing exit points (destinations)
entry_id = 5        # Entry cells (where pedestrians enter simulation_grid)
exit_offset = 100   # Offset used for integer ids of exit points (destinations)
entry_offset = 200  # Offset used for integer ids of entry points
# ---------------------------------------------------------------------------- #


class Cell(object):
    """This class is used to represent a walking_grid cell.

    Types are denoted by integer values:
    0 - Prohibited (non-walkway)
    1 - Designated walkway
    2 - Street (prohibited for walking)
    3 - Transient walkway (e.g., street crossing)
    """

    def __init__(self, coordinates, cell_type=0, is_open=0, group_id=999):
        """

        :param coordinates: The x,y coordinates in the walking_grid
        :param cell_type: Integer value indicating the cell type
        :param is_open: For use when the cell represents a transient path
        :param group_id: Used to indicate membership to a group of cells
        """
        self.coordinates = coordinates
        self.type = cell_type
        self.open = is_open
        self.group_id = group_id


class Grid(object):
    """This class is used to represent the simulation walking_grid."""

    def __init__(self):
        """Object initialization"""
        self.size = (0, 0)
        self.cells = []
        self.exit_group_ids = set()
        self.exit_coordinates = []
        self.walking_grid = None
        self.map = None
        self.floor_fields = []

    def build_from_map_file(self, map_file_path):
        """This function builds the simulation gid from the given map file.

        :param map_file_path: Path to the map file
        """
        global prohibited_id, walkway_id, street_id, crossing_id,\
            exit_offset, entry_offset, exit_id, entry_id

        with open(map_file_path, 'r') as fid:

            # Skip comment lines
            map_line = fid.readline()
            while map_line[0] == '#':
                map_line = fid.readline()

            # Read the header and get the number of nodes and edges.
            r, c, a = [int(v) for v in map_line.split()]

            # Create a walking_grid of prohibited cells
            self.cells = [[Cell((i, j)) for j in range(c)] for i in range(r)]
            self.size = (r, c)

            # initialize the walking_grid
            self.walking_grid = zeros((r, c), dtype=int)
            self.map = zeros((r, c), dtype=int)

            # Read each line of the map file
            for map_line in fid:

                # Skip comment lines
                if map_line[0] == '#':
                    continue

                # Parse the values from the line of text
                sr, sc, er, ec, t = [int(v) for v in map_line.split()]

                # All the walking grid values default to prohibited. Update the
                # cells to the types specified by the map file.
                for i in range(min((sr, er)), max((sr, er)) + 1):
                    for j in range(min((sc, ec)), max((sc, ec)) + 1):

                        # Branch to account for the different cell types
                        # Note: Higher order types override
                        if t == walkway_id:  # Walkway
                            self.walking_grid[i, j] = walkway_id
                            # Check for street crossing
                            if self.cells[i][j].type == street_id:
                                self.cells[i][j].type = crossing_id
                                self.map[i, j] = crossing_id
                            else:
                                self.cells[i][j].type = walkway_id
                                self.map[i, j] = walkway_id

                        elif t == street_id:  # Street
                            if self.walking_grid[i, j] == walkway_id:
                                self.walking_grid[i, j] = walkway_id
                                self.cells[i][j].type = crossing_id
                                self.map[i, j] = crossing_id
                            else:
                                self.map[i, j] = street_id
                                self.cells[i][j].type = street_id

                        elif t == crossing_id:  # Street crossing
                            self.walking_grid[i, j] = walkway_id
                            self.cells[i][j].type = crossing_id
                            self.map[i, j] = crossing_id

                        # Destination cell (exit simulation area)
                        elif exit_offset <= t < entry_offset:
                            self.walking_grid[i, j] = walkway_id
                            self.map[i, j] = exit_id
                            self.cells[i][j].type = exit_id
                            group_id = t % exit_offset
                            self.exit_group_ids.add(group_id)
                            self.cells[i][j].group_id = group_id
                            self.exit_coordinates.append((i, j))
                            # Log the start and destination locations
                            # dest_ij.append((i, j))

                        # Simulation entry point
                        elif t >= entry_offset:
                            self.walking_grid[i, j] = walkway_id
                            self.map[i, j] = entry_id
                            self.cells[i][j].type = entry_id
                            group_id = t % entry_offset
                            self.cells[i][j].group_id = group_id
                            # Log the start and destination locations
                            # start_ij.append((i, j))

    def calculate_floor_fields(self):
        """Calculate the static floor fields for each destination (exit)"""

        # Create floor fields for each destination
        for xid in self.exit_group_ids:

            # Get a list of the cell locations relating to this exit group
            exit_loc = []
            for loc in self.exit_coordinates:
                if self.cells[loc[0]][loc[1]].group_id == xid:
                    exit_loc.append(loc)

            # Initialize the floor field matrix
            floor_field = ones(shape(self.walking_grid))*inf

            # Keep a list of cells already mapped into the floor field
            mapped = set()

            # Find the (allowed) cells that neighbor the exit cells.
            neighbor_list = []
            for d in exit_loc:
                floor_field[d[0], d[1]] = 1
                mapped.add(d)
                for n in self._find_neighbor_cells(d):
                    if n not in exit_loc:
                        neighbor_list.append(n)
                        floor_field[n[0], n[1]] = 2
                        mapped.add(n)

            # create a queue for the neighbors
            neighbor_queue = collections.deque(neighbor_list)

            level = 3  # current level of neighbors to find.
            level_count = len(neighbor_queue)
            count = 0
            while len(neighbor_queue) > 0:

                # Calculate which level to mark the neighboring cells. When
                # |count| is equal to |level_count|, all the next neighbors
                # will be at the next distance level.
                if count == level_count:
                    level_count = len(neighbor_queue)
                    level += 1
                    count = 0

                cell_coordinates = neighbor_queue.pop()
                temp = self._find_neighbor_cells(cell_coordinates)
                for n in temp:
                    if n not in mapped:
                        neighbor_queue.appendleft(n)
                        mapped.add(n)
                        floor_field[n[0], n[1]] = level
                count += 1  # Increment the count

            self.floor_fields.append(floor_field)

    def _find_neighbor_cells(self, coordinate):
        """Helper method that finds (allowed) neighboring cells."""
        global prohibited_id
        x, y = coordinate[0], coordinate[1]
        max_x = self.size[0] - 1
        max_y = self.size[1] - 1
        min_x, min_y = 0, 0
        # print(x, y, max_x, max_y, min_x, min_y)
        neighbors = []
        if x + 1 <= max_x:
            if self.walking_grid[x + 1, y] != prohibited_id:
                neighbors.append((x + 1, y))
        if x - 1 >= min_x:
            if self.walking_grid[x - 1, y] != prohibited_id:
                neighbors.append((x - 1, y))
        if y + 1 <= max_y:
            if self.walking_grid[x, y + 1] != prohibited_id:
                neighbors.append((x, y + 1))
        if y - 1 >= min_y:
            if self.walking_grid[x, y - 1] != prohibited_id:
                neighbors.append((x, y - 1))
        return neighbors

    def show_map(self):
        """Show a visual representation of the map."""
        print(self.map)
        imshow(self.map, interpolation='nearest', origin='lower')
        show()

    def show_walking_map(self):
        """Show a visual representation of the walking map."""
        print(self.walking_grid)
        imshow(self.walking_grid, interpolation='nearest', origin='lower')
        show()

    def show_floor(self):
        """Show a visual representation of the floor field(s)."""
        for field in self.floor_fields:
            print(field)
            imshow(field, interpolation='nearest', origin='lower')
            show()


class Pedestrian(object):
    """This class represents a pedestrian traveling in the walking_grid."""

    def __init__(self, destination, speed=1):
        self.speed = speed
        self.destination = destination


class PedestrianExitSimulation(object):
    """This object encapsulates the entire simulation program."""

    def __init__(self, map_file_path, pedestrian_count, rng_seed):
        random.seed(rng_seed)  # Seed the RNG so results can be reproduced
        self.map_file_path = map_file_path
        self.pedestrian_count = pedestrian_count
        self._pedestrian_queue = None
        self._grid = Grid()

    def _queue_pedestrians(self):
        destinations = list(self._grid.exit_group_ids)
        self._pedestrian_queue = collections.deque()
        for n in range(0, self.pedestrian_count - 1):
            destination = destinations[random.randint(0, len(destinations) - 1)]
            self._pedestrian_queue.appendleft(Pedestrian(destination))

    def initialize(self):
        """Initialize system state"""
        self._grid.build_from_map_file(self.map_file_path)
        self._grid.calculate_floor_fields()
        self._queue_pedestrians()
        print(self._pedestrian_queue)

    def observe(self):
        pass
        #self._grid.show_map()
        #self._grid.show_floor()

    def run(self):
        self.initialize()
        self.observe()


if __name__ == "__main__":

    # Path to the map file
    map_path = "map_test_02.txt"

    # Number of pedestrians to simulate
    n_pedestrians = 10

    # Seed to RNG
    seed = 1

    # Create a simulator object.
    sim = PedestrianExitSimulation(map_path, n_pedestrians, seed)

    # Number of iterations
    # niter = 10

    # Run the simulation
    sim.run()
