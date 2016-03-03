import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


PROHIBITED_ID = 0   # Prohibited (non-walkway)
WALKWAY_ID = 1      # Designated walkway
STREET_ID = 2       # Street (prohibited for walking)
CROSSING_ID = 3     # Transient walkway (e.g., street crossing)
EXIT_ID = 4         # Cells representing exit points (destinations)
ENTRY_ID = 5        # Entry cells (where pedestrians enter simulation_grid)
EXIT_OFFSET = 100   # Offset used for integer ids of exit points (destinations)
ENTRY_OFFSET = 200  # Offset used for integer ids of entry points


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


def build_from_map_file(map_file_path):
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
        r, c, a = [float(v) for v in map_line.split()]
        r = int(r)
        c = int(c)

        # initialize the walking_grid
        walking_grid = np.zeros((r, c), dtype=int)
        map_grid = np.zeros((r, c), dtype=int)

        # Read each line of the map file
        for map_line in fid:

            # Skip comment lines
            if map_line[0] == '#':
                continue

            values = [int(v) for v in map_line.split() if v.isdigit()]
            if len(values) != 5:
                continue

            # Parse the values from the line of text
            sr, sc, er, ec, t = values

            # All the walking grid values default to prohibited. Update the
            # cells to the types specified by the map file.
            for row in range(min((sr, er)), max((sr, er)) + 1):
                for col in range(min((sc, ec)), max((sc, ec)) + 1):

                    # Branch to account for the different cell types
                    # Note: Higher order types override
                    if t == WALKWAY_ID:  # Walkway
                        walking_grid[row, col] = WALKWAY_ID
                        # Check for street crossing
                        if map_grid[row, col] == STREET_ID:
                            map_grid[row, col] = CROSSING_ID
                        else:
                            map_grid[row, col] = WALKWAY_ID

                    elif t == STREET_ID:  # Street
                        if walking_grid[row, col] == WALKWAY_ID:
                            walking_grid[row, col] = WALKWAY_ID
                            map_grid[row, col] = CROSSING_ID
                        else:
                            map_grid[row, col] = STREET_ID

                    elif t == CROSSING_ID:  # Street crossing
                        walking_grid[row, col] = WALKWAY_ID
                        map_grid[row, col] = CROSSING_ID

                    # Destination cell (exit simulation area)
                    elif EXIT_OFFSET <= t < ENTRY_OFFSET:
                        walking_grid[row, col] = WALKWAY_ID
                        map_grid[row, col] = EXIT_ID
                        group_id = t % EXIT_OFFSET

                    # Simulation entry point
                    elif t >= ENTRY_OFFSET:
                        walking_grid[row, col] = WALKWAY_ID
                        map_grid[row, col] = ENTRY_ID
                        group_id = t % ENTRY_OFFSET

                    # Unused grid cell ID
                    else:
                        if map_grid[row, col] == PROHIBITED_ID:
                            map_grid[row, col] = t

    return (map_grid, walking_grid)


def show_map(mapgrid):
    """Show a visual representation of the map."""
    print(mapgrid)
    map_colors = np.array(
        [[0, 0, 0], [255, 255, 255], [0, 100, 255], [255, 150, 0],
         [255, 0, 0], [227, 0, 252], [0, 255, 50]], dtype=float
    )
    fig, a = plt.subplots()
    a.xaxis.set_visible(False)
    a.yaxis.set_visible(False)
    plt.imshow(
        mapgrid, interpolation='nearest', origin='lower', aspect='equal',
        cmap=mpl.colors.ListedColormap(map_colors/255)
    )
    plt.show()


if __name__ == "__main__":
    import os
    print(os.getcwd())
    show_map(build_from_map_file("GT_Boddy_Dodd_South_Map.txt")[0])
