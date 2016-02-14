from pylab import *
import collections

def find_neighbor_cells(cells, coordinate):
    x, y = coordinate[0], coordinate[1]
    max_x = shape(cells)[0] - 1
    max_y = shape(cells)[1] - 1
    min_x, min_y = 0, 0
    neighbors = []
    if x + 1 <= max_x:
        if cells[x + 1, y] > 0:
            neighbors.append((x + 1, y))
    if x - 1 >= min_x:
        if cells[x - 1, y] > 0:
            neighbors.append((x - 1, y))
    if y + 1 <= max_y:
        if cells[x, y + 1] > 0:
            neighbors.append((x, y + 1))
    if y - 1 >= min_y:
        if cells[x, y - 1] > 0:
            neighbors.append((x, y - 1))
    # print(coordinate, neighbors)
    return neighbors

# Build walking_grid map
if __name__ == "__main__":

    start_ij = []
    dest_ij = []

    plot_grid = []
    with open("map.txt", 'r') as fileID:

        # Read the header and get the number of nodes and edges.
        rows, cols, size = [int(value) for value in fileID.readline().split()]

        plot_grid = zeros((rows, cols), dtype=int)

        for cellID, line in enumerate(fileID):
            sr, sc, er, ec, t = [int(value) for value in line.split()]
            #print(line)
            for i in range(min((sr, er)), max((sr, er)) + 1):
                for j in range(min((sc, ec)), max((sc, ec)) + 1):
                    if plot_grid[i, j] <= t:
                        plot_grid[i, j] = t
                        if t == 3:
                            dest_ij.append((i, j))
                        if t == 4:
                            start_ij.append((i, j))
            #print(plot_grid)

        print(start_ij)
        print(dest_ij)

    imshow(plot_grid, interpolation='nearest', origin='lower')
    show()

    # Initialize the floor field
    floor_field = ones(shape(plot_grid))*inf

    # Keep a list of cells already mapped into the floor field
    mapped = set()

    # Find the (allowed) cells that neighbor the destination cells.
    neighbor_list = []
    for d in dest_ij:
        floor_field[d[0], d[1]] = 1
        mapped.add(d)
        for n in find_neighbor_cells(plot_grid, d):
            if n not in dest_ij:
                neighbor_list.append(n)
                floor_field[n[0], n[1]] = 2
                mapped.add(n)

    print(neighbor_list)
    print(mapped)

    # create a queue (it currently holds max level 1)
    neighbor_queue = collections.deque(neighbor_list)
    print(neighbor_queue)

    #pause(20)

    #

    level = 3
    level_count = len(neighbor_queue)

    count = 0
    while len(neighbor_queue) > 0:
        #print(count)

        if count == level_count:
            level_count = len(neighbor_queue)
            level += 1
            count = 0

        cell_coordinates = neighbor_queue.pop()
        temp = find_neighbor_cells(plot_grid, cell_coordinates)
        for n in temp:
            if n not in mapped:
                neighbor_queue.appendleft(n)
                mapped.add(n)
                floor_field[n[0], n[1]] = level

        #print(len(temp), len(neighbor_queue))
        count += 1
        #pause(1)

    print(floor_field)
    imshow(floor_field, interpolation='nearest', origin='lower')
    show()

