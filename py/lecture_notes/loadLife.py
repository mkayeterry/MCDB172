import matplotlib.pyplot as plt
import numpy as np
from IPython.display import clear_output
from time import sleep

def loadlife(filename):
    
    def read_grid(j, lines):
        grid = []
        
        while not lines[j].startswith('#P'):
            
            line = list(lines[j])[:-1]
            
            tx_line = [0 if item == '.' else item for item in line]
            tx_line = [1 if item == '*' else item for item in tx_line]
            
            grid.append(tx_line)
            
            j += 1
            if j == len(lines): break
                
        return grid
    
    file = open(filename, 'r')
    
    lines = file.readlines()
    
    point_map = {}
    points = []
    
    for i, line in enumerate(lines):
        if line.startswith('#P'):
            point = (int(line.split()[1]), int(line.split()[2]))
            points.append(point)
            grid = read_grid(i + 1, lines)
            point_map[point] = grid
            
    for point, grid in point_map.items():
        # find max length if grid
        max_length = max([len(row) for row in grid])

        #Add leftovers
        for row in grid:
            row.extend([0]*(max_length - len(row)))
        
    most_neg_x = min([point[0] for point in points])
    most_pos_y = max([point[1] for point in points])
    
    most_pos_x = max([point[0] + len(point_map[point][0]) - 1 for point in points])
    most_neg_y = min([point[1] - len(point_map[point]) + 1 for point in points])
    
    x_len = most_pos_x - most_neg_x + 1
    y_len = most_pos_y - most_neg_y + 1
    
    universe = [[0 for _ in range(x_len)] for _ in range(y_len)]
    
    normalized_numpized_point_map = {}
    
    for point, grid in point_map.items():
        #normalize point
        x, y = point
        
        x_n = x - most_neg_x
        y_n = (-y) + most_pos_y
        
        point_n = x_n, y_n
        normalized_numpized_point_map[point_n] = np.asarray(grid)
        
    universe = np.asarray(universe)
    universe = np.zeros((y_len,x_len))
    
    for point, grid in normalized_numpized_point_map.items():
        
        rowu = point[1]
        rowl = int(point[1] + grid.shape[0])
        coll = point[0]
        colr = int(point[0] + grid.shape[1])
        
    
#         print('col dif ', colr-coll)
#         print('cols in universe ', universe.shape[1])
#         print('row dif ', rowl-rowu)
#         print('rows in universe ', universe.shape[0])
#         print('grid shape', grid.shape)
#         print('rowu ', rowu, 'row l ', rowl)
        
        universe[rowu:rowl, coll:colr] = grid
    
    return universe
