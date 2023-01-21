import numpy as np

# Define the x and y values of the data points
x_values = np.array([1, 2, 3, 4])
y_values = np.array([2, 4, 8, 16])

# Create the Vandermonde matrix
n = len(x_values)
V = np.array([x_values**(n-1-i) for i in range(n)]).T

print(V)

# Solve for the polynomial coefficients using the Vandermonde matrix
coefficients = np.linalg.solve(V, y_values)

# Print the polynomial equation
print("Polynomial equation:")
for i, c in enumerate(coefficients):
    print(f"{c}x^{i}", end=" + ")

"""A = np.array([
    [0, 1, 1, 0, 0, 0, 1, 0, 0, 1],
    [1, 0, 0, 1, 1, 1, 0, 0, 0, 1],
    [1, 0, 0, 1, 0, 0, 1, 0, 1, 0],
    [0, 1, 1, 0, 1, 0, 1, 0, 1, 0],
    [0, 1, 0, 1, 0, 1, 0, 1, 0, 0],
    [0, 1, 0, 0, 1, 0, 0, 0, 0, 0],
    [1, 0, 1, 1, 0, 0, 0, 1, 1, 0],
    [0, 0, 0, 0, 1, 0, 1, 0, 0, 0],
    [0, 0, 1, 1, 0, 0, 1, 0, 0, 0],
    [1, 1, 0, 0, 0, 0, 0, 0, 0, 0]
])
B = np.dot(np.dot(A, A), np.dot(A, A))

print(B)
print(B[4,8])
print(B[8,4])"""
