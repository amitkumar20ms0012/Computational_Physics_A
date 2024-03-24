import numpy as np


A = np.array([[10,-1,0], [-1,10,-2],[0,-2,10]])


B = np.array([9,7,6])


solution = np.linalg.solve(A, B)

print("Solution of part 2 Q. 11:", solution)