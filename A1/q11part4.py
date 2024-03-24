import numpy as np

A = np.array([[4, 1, 1, 0, 1], [-1, -3, 1, 1, 0], [2, 1, 5, -1, -1], [-1, -1, -1, 4, 0], [0, 2, -1, 1, 4]])

B = np.array([6,6,6,6,6])

solution = np.linalg.solve(A, B)

print("Solution of part 4 Q. 11:", solution)