#!/usr/bin/env python
# This program simulates the flow of people leaving a starting point and heading
# to their desired locations. The simulation is based on Cellular Automata (CA).
# Dylan A. Crocker
# CSE 6730 Project 1
# Due: March 4, 2016


# Imports -------------------------------------------------------------------- #


from pylab import *
import collections


# Constants ------------------------------------------------------------------ #


PROHIBITED_ID = 0   # Prohibited (non-walkway)
WALKWAY_ID = 1      # Designated walkway
STREET_ID = 2       # Street (prohibited for walking)
CROSSING_ID = 3     # Transient walkway (e.g., street crossing)
EXIT_ID = 4         # Cells representing exit points (destinations)
ENTRY_ID = 5        # Entry cells (where pedestrians enter simulation_grid)
EXIT_OFFSET = 100   # Offset used for integer ids of exit points (destinations)
ENTRY_OFFSET = 200  # Offset used for integer ids of entry points


# Functions ------------------------------------------------------------------ #


def min_rand(idum=0):
    """Minimal random number generator of Park and Miller. Returns a uniform
    random deviate between 0.0 and 1.0. Set or reset idum to any integer value
    (except the unlikely value MASK) to initialize the sequence; idum must not
    be altered between calls for successive deviates in a sequence.
    Ref. Numerical Recipes in C 2nd ed."""

    # Define constants
    ia = 16807
    im = 2147483647
    am = (1.0/im)
    iq = 127773
    ir = 2836
    mask = 123459876

    # Only allow the generator to be seeded once
    if "seed" not in min_rand.__dict__:
        # XORing with MASK allows use of zero and other simple bit patterns for
        # seed.
        min_rand.seed = idum ^ mask

    # Compute idum=(IA*idum) % IM without over-flows by Schrage's method.
    k = min_rand.seed/iq
    min_rand.seed = ia*(min_rand.seed - k*iq) - ir*k
    if min_rand.seed  < 0:
        min_rand.seed += im
    ans = am*min_rand.seed  # Convert to a floating result.
    min_rand.seed ^= mask   # Unmask before return.

    return ans


def min_rand_int(min_int=0, max_int=10):

    # Ensure inputs are integers
    min_int = int(min_int)
    max_int = int(max_int)

    # Since int() always rounds down add one to the delta
    delta = max_int - min_int + 1

    # Scale the random number generated
    random_number = min_rand()*delta

    # In the rare case where the random number exactly equals delta, decrement
    # enough so that the result will round to the next integer down.
    if random_number == delta:
        random_number -= 0.1

    # Scale back to the desired range before returning
    return int(random_number + min_int)


def gaussian(mu=0, sigma=1):
    """Generate a random number from a Gaussian distribution.

    Adapted from Numerical Recipes in C 2nd ed.
    """

    if "iset" not in gaussian.__dict__:
        gaussian.iset = 0
    if "gset" not in gaussian.__dict__:
        gaussian.gset = 0

    if gaussian.iset == 0:
        v1, v2, rsq = 0.0, 0.0, 0.0
        while rsq >= 1.0 or rsq == 0.0:
            v1 = 2.0*min_rand() - 1
            v2 = 2.0*min_rand() - 1
            rsq = v1*v1 + v2*v2
        fac = np.sqrt(-2.0*np.log(rsq)/rsq)
        gaussian.iset = 1
        gaussian.gset = v1*fac
        xi = v2*fac
    else:
        gaussian.iset = 0
        xi = gaussian.gset

    # Scale from standard normal distribution to that desired.
    # See page 242 in Statistics for Engineers and Scientists 2nd ed. by Navidi
    return mu + xi*sigma


def gaussian_distribution(mu=0., sigma=1., n=100):
    """Create a list of random numbers from a Gaussian distribution."""
    return [gaussian(mu, sigma) for i in range(n)]


def generate_pedestrian_speeds(n_ped):
    """Generate a normally distributed list of pedestrian walking speeds.

    This function utilizes data obtained from "Evaluation of Pedestrian Walking
    Speeds in Airport Terminals" by Seth B. Young

    :param n_ped: The number of pedestrians for which to generate speeds.
    :return: List of normally distributed walking speeds.
    """
    speed_mu = 1.34        # Average pedestrian walking speed (m/s)
    speed_sigma = 0.87     # Std. Dev. of walking speeds (m/s)
    speed_minimum = 0.153  # Minimum walking speed (m/s)
    speeds = gaussian_distribution(speed_mu, speed_sigma, n_ped)
    for idx, speed in enumerate(speeds):
        if speed < speed_minimum:
            speeds[idx] = speed_minimum
    return speeds


# Classes -------------------------------------------------------------------- #


class Pedestrian(object):
    """This class represents a pedestrian traveling in the walking_grid."""

    def __init__(self, destination, speed=1):
        self.speed = speed
        self.destination = destination
        self.coordinates = [0, 0]
        self.moves = []
        self.to_move = 0
        self.uid = -1

    def normalize_moves(self):
        sum_probs = sum([m[0] for m in self.moves])
        for i, m in enumerate(self.moves):
            self.moves[i] = (m[0]/sum_probs, m[1], m[2], m[3])

    def move(self, location):
        self.coordinates = list(location)
        self.moves = []


class Cell(object):
    """This class is used to represent a walking_grid cell.

    Types are denoted by integer values:
    0 - Prohibited (non-walkway)
    1 - Designated walkway
    2 - Street (prohibited for walking)
    3 - Transient walkway (e.g., street crossing)
    """

    def __init__(self, coordinates, cell_type=0, is_open=0, group_id=999):
        """Initialize cell object

        :param coordinates: The x,y coordinates in the walking_grid
        :param cell_type: Integer value indicating the cell type
        :param is_open: For use when the cell represents a transient path
        :param group_id: Used to indicate membership to a group of cells
        """
        self.coordinates = coordinates
        self.type = cell_type
        self.open = is_open
        self.group_id = group_id
        self.occupant = None


class Grid(object):
    """This class is used to represent the simulation walking_grid."""

    def __init__(self):
        """Object initialization"""
        self.size = (0, 0)
        self.cells = []
        self.exit_group_ids = set()
        self.exit_coordinates = []
        self.entry_coordinates = []
        self.walking_grid = None
        self.map = None
        self.floor_fields = []
        self.pedestrians = []
        self.cell_len = 0.4  # Length of a cell (m)

    def setup(self, map_file_path):
        self._build_from_map_file(map_file_path)
        self._calculate_floor_fields()

    def find_reachable_cells(self, loc, steps=1):
        """Find reachable cells given the walking distance (number of cells).

        This function utilizes an algorithm similar to the breadth first search
        in graph theory to find all the cells reachable by a pedestrian willing
        to travel a certain amount of steps. The algorithm avoids obstacles
        (prohibited cells, occupied cells, etc.) by design.

        :param loc: Current location or starting point (i,j)
        :param steps: The maximum number of steps (cells) allowed to transverse
        :return: List of cells that are reachable
        """

        # Create a hash table of cells to store only unique moves. Add the
        # starting location to the since it will be seen as a neighbor when
        # evaluating its neighbors.
        start_cell = self.cells[loc[0]][loc[1]]
        cells = set([start_cell])

        # create a queue for the neighbors
        starts = collections.deque([loc])

        # Loop through neighbors for each step distance away.
        for n in range(steps):

            # A start is a neighbor that has not yet had its neighbors explored
            for m in range(len(starts)):

                # Find neighbor cells
                neighbor_locations = self._find_neighbor_cells(starts.pop())
                for coor in neighbor_locations:
                    cell = self.cells[coor[0]][coor[1]]
                    if cell not in cells:  # Only add unique neighbor cells
                        cells.add(cell)
                starts.extendleft(neighbor_locations)

        # Remove the original location?
        #cells.remove(start_cell)  # What if the pedestrian wants to stand still?

        return list(cells)  # Return a list - not a hash table

    def _build_from_map_file(self, map_file_path):
        """This function builds the simulation gid from the given map file.

        Note: The map file is specified in (row, column) format. All the data
        arrays are likewise indexed (row, column).

        :param map_file_path: Path to the map file
        """
        global PROHIBITED_ID, WALKWAY_ID, STREET_ID, CROSSING_ID,\
            EXIT_OFFSET, ENTRY_OFFSET, EXIT_ID, ENTRY_ID

        with open(map_file_path, 'r') as fid:

            # Skip comment lines
            map_line = fid.readline()
            while map_line[0] == '#':
                map_line = fid.readline()

            # Read the header and get the number of nodes and edges.
            r, c, a = [int(v) for v in map_line.split()]

            # Create a walking_grid of prohibited cells
            self.cells = [[Cell((j, i)) for i in range(c)] for j in range(r)]
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
                for row in range(min((sr, er)), max((sr, er)) + 1):
                    for col in range(min((sc, ec)), max((sc, ec)) + 1):

                        # Branch to account for the different cell types
                        # Note: Higher order types override
                        if t == WALKWAY_ID:  # Walkway
                            self.walking_grid[row, col] = WALKWAY_ID
                            # Check for street crossing
                            if self.cells[row][col].type == STREET_ID:
                                self.cells[row][col].type = CROSSING_ID
                                self.map[row, col] = CROSSING_ID
                            else:
                                self.cells[row][col].type = WALKWAY_ID
                                self.map[row, col] = WALKWAY_ID

                        elif t == STREET_ID:  # Street
                            if self.walking_grid[row, col] == WALKWAY_ID:
                                self.walking_grid[row, col] = WALKWAY_ID
                                self.cells[row][col].type = CROSSING_ID
                                self.map[row, col] = CROSSING_ID
                            else:
                                self.map[row, col] = STREET_ID
                                self.cells[row][col].type = STREET_ID

                        elif t == CROSSING_ID:  # Street crossing
                            self.walking_grid[row, col] = WALKWAY_ID
                            self.cells[row][col].type = CROSSING_ID
                            self.map[row, col] = CROSSING_ID

                        # Destination cell (exit simulation area)
                        elif EXIT_OFFSET <= t < ENTRY_OFFSET:
                            self.walking_grid[row, col] = WALKWAY_ID
                            self.map[row, col] = EXIT_ID
                            self.cells[row][col].type = EXIT_ID
                            group_id = t % EXIT_OFFSET
                            self.exit_group_ids.add(group_id)
                            self.cells[row][col].group_id = group_id
                            # Log the start and destination locations
                            self.exit_coordinates.append((row, col))

                        # Simulation entry point
                        elif t >= ENTRY_OFFSET:
                            self.walking_grid[row, col] = WALKWAY_ID
                            self.map[row, col] = ENTRY_ID
                            self.cells[row][col].type = ENTRY_ID
                            group_id = t % ENTRY_OFFSET
                            self.cells[row][col].group_id = group_id
                            # Log the start and destination locations
                            self.entry_coordinates.append((row, col))

    def _calculate_floor_fields(self):
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
        global PROHIBITED_ID
        i, j = coordinate[0], coordinate[1]
        max_i = self.size[0] - 1
        max_j = self.size[1] - 1
        min_i, min_j = 0, 0
        neighbors = []
        if i + 1 <= max_i:
            if self.walking_grid[i + 1, j] != PROHIBITED_ID:
                neighbors.append((i + 1, j))
        if i - 1 >= min_i:
            if self.walking_grid[i - 1, j] != PROHIBITED_ID:
                neighbors.append((i - 1, j))
        if j + 1 <= max_j:
            if self.walking_grid[i, j + 1] != PROHIBITED_ID:
                neighbors.append((i, j + 1))
        if j - 1 >= min_j:
            if self.walking_grid[i, j - 1] != PROHIBITED_ID:
                neighbors.append((i, j - 1))
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


class PedestrianExitSimulation(object):
    """This object encapsulates the entire simulation program."""

    def __init__(self, map_file_path, pedestrian_count, rng_seed=0):
        min_rand(rng_seed)  # Seed the RNG so results can be reproduced
        self.map_file_path = map_file_path
        self.pedestrian_count = pedestrian_count
        self._pedestrian_queue = None
        self._grid = Grid()

    def _queue_pedestrians(self):
        """Build a queue of pedestrian entities."""

        # Get a list of possible destinations
        destinations = list(self._grid.exit_group_ids)

        # Create a distribution of walking speeds.
        n = self.pedestrian_count
        speeds = generate_pedestrian_speeds(n)

        # Create the queue
        self._pedestrian_queue = collections.deque()
        for pidx in range(0, n):

            # Randomly choose a destination
            destination = destinations[min_rand_int(0, len(destinations) - 1)]

            # Enqueue a new pedestrian object
            self._pedestrian_queue.appendleft(
                Pedestrian(destination, speed=speeds[pidx])
            )

    def _calculate_moves(self, time_elapsed):
        """Calculate the pedestrian movements within the grid"""

        # Loop over each pedestrian and calculate the move probabilities
        for p in self._grid.pedestrians:

            # clear the list of moves
            p.moves = []

            # Calculate the distance the pedestrian can move
            move_distance = p.speed*time_elapsed + p.to_move

            # convert distance to grid steps
            steps = int(np.floor(move_distance/self._grid.cell_len))

            # Save any truncated distance
            p.to_move = move_distance % self._grid.cell_len

            # Find possible cells given the distance the pedestrian can travel
            cells = self._grid.find_reachable_cells(p.coordinates, steps=steps)

            # Calculate movement probability for each reachable cell
            dest_id = p.destination
            prow, pcol = p.coordinates
            for cell in cells:
                i, j = cell.coordinates
                prev_potential = self._grid.floor_fields[dest_id][prow, pcol]
                new_potential = self._grid.floor_fields[dest_id][i, j]
                prob = np.exp(prev_potential - new_potential)

                # Add a random number in the second element to break ties in
                # movement probabilities during sorting. Prevents bias towards
                # largest row value.
                p.moves.append((prob, min_rand(), i, j))

            # Normalize the move probabilities
            p.normalize_moves()

            # Sort by priority!
            p.moves.sort(reverse=True)  # This is important!

        # Initialize the moves and probabilities grids
        nrows, ncols = self._grid.size
        moves = [[None for j in range(ncols)] for i in range(nrows)]
        probabilities = zeros(self._grid.size, dtype=float)

        # Implement an iterative process to determine moves
        pedestrians_to_move = collections.deque(self._grid.pedestrians)
        while len(pedestrians_to_move) > 0:

            pedestrian = pedestrians_to_move.pop()  # Pedestrian to move

            # Loop over all the possible moves and determine which will not
            # result in conflicts.
            for potential_move in pedestrian.moves:
                prob, rn, i, j = potential_move

                # If no pedestrian has yet claimed the spot...
                # Note: This will always happen for the current cell (if it gets
                # that low in the move probabilities.
                if moves[i][j] is None:
                    probabilities[i, j] = prob
                    moves[i][j] = pedestrian
                    break

                # If someone else has claimed the spot but with a lower
                # probability, place them back in the queue and take the spot.
                if probabilities[i, j] < prob:
                    pedestrians_to_move.appendleft(moves[i][j])
                    probabilities[i, j] = prob
                    moves[i][j] = pedestrian
                    break

        return moves

    def _move_pedestrian(self, pedestrian, new_location):

        # Set the old cell to empty
        i, j = pedestrian.coordinates
        self._grid.cells[i][j].occupant = None

        # Move the pedestrian
        pedestrian.move(new_location)

        # Set the new cell to occupied
        i, j = pedestrian.coordinates
        self._grid.cells[i][j].occupant = pedestrian

    def initialize(self):
        """Initialize system state"""
        self._grid.setup(self.map_file_path)
        self._queue_pedestrians()
        # self._grid.show_map()
        # self._grid.show_floor()

    def observe(self, visual=False):
        """Observe"""

        # Get the current pedestrian locations
        xs, ys = [], []
        for ped in self._grid.pedestrians:
            xs.append(ped.coordinates[1])  # column
            ys.append(ped.coordinates[0])  # row

        # Write the state to file

        if visual:
            # Plot the map with the pedestrians shown as dots over top
            imshow(self._grid.map, interpolation='nearest', origin='lower')
            scatter(x=xs, y=ys, c='r', s=100)
            show()

    def update(self, time_step=1):
        """Update the simulation state."""

        # Check for pedestrians who have arrived at the destination and remove
        # them from the list. Use a set for faster membership tests.
        exits = set(self._grid.exit_coordinates)
        self._grid.pedestrians[:] = [p for p in self._grid.pedestrians
                                     if tuple(p.coordinates) not in exits]

        # Determine the pedestrian movements
        moves = self._calculate_moves(time_step)

        # Execute the moves
        nrows, ncols = self._grid.size
        for i in range(nrows):
            for j in range(ncols):
                pedestrian = moves[i][j]
                if pedestrian is None:
                    continue
                self._move_pedestrian(pedestrian, (i, j))

        # Add pedestrians to the grid
        if len(self._pedestrian_queue) > 0:

            # For every entry coordinate make a random choice whether or not to
            # add a pedestrian to that cell.
            for entry_coordinate in self._grid.entry_coordinates:
                i, j = entry_coordinate
                if min_rand() > 0.5 and len(self._pedestrian_queue) > 0:
                    entry_cell = self._grid.cells[i][j]
                    if entry_cell.occupant is None:
                        new_pedestrian = self._pedestrian_queue.pop()
                        entry_cell.occupant = new_pedestrian
                        new_pedestrian.coordinates = [i, j]
                        self._grid.pedestrians.append(new_pedestrian)

    def run(self, max_time_steps=100):
        """Execute the simulation loop."""
        self.initialize()
        self.observe()
        for t in range(1, max_time_steps):  # Time steps are once a second
            self.update()
            if len(self._grid.pedestrians) == 0:
                break
            self.observe(visual=True)


# ---------------------------------------------------------------------------- #


def main():

    # Path to the map file
    map_path = "map_test_02.txt"

    # Number of pedestrians to simulate
    n_pedestrians = 100

    # Seed the RNG
    seed_rng = 34643

    # Create a simulator object.
    sim = PedestrianExitSimulation(map_path, n_pedestrians, seed_rng)

    # Number of iterations
    niter = 200

    # Run the simulation
    sim.run(niter)

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)
