import numpy as np


A = np.array([[10,5,0,0], [5,10,-4,0],[0,-4,8,-1],[0,0,-1,5]])


B = np.array([6,25,-11,-11])


solution = np.linalg.solve(A, B)

print("Solution of part 3 Q. 11:", solution)