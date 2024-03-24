import numpy as np


A = np.array([[3, -1, 1], [3, 6, 2],[3, 3, 7]])


B = np.array([1, 0, 4])


solution = np.linalg.solve(A, B)

print("Solution of part 1 Q. 11:", solution)