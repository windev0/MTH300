import numpy as np


# Définir les coefficients du système d'équations linéaires sous forme de matrice augmentée
A = np.array([[-1, 1, 2],[1, -1, 3],[0, 1, 1]])
B = np.array([2,3,2])

# Concaténer la matrice A et le vecteur B en une seule matrice augmentée
AB = np.column_stack((A, B))

# Trouver la taille de la matrice augmentée
n, m = AB.shape

# Appliquer la méthode de Gauss-Jordan
for i in range(n):
    # Diviser la i-ème ligne par l'élément diagonal
    div = AB[i, i]
    AB[i, :] /= div
    
    # Soustraire la i-ème ligne multipliée par l'élément non diagonal des autres lignes
    for j in range(n):
        if i != j:
            mult = AB[j, i]
            AB[j, :] -= mult * AB[i, :]
            
# Extraire la solution x à partir de la dernière colonne de la matrice augmentée
x = AB[:, -1]

print("La solution du système d'équations linéaires est:", x)
print(np.linalg.solve(A, B))


# Print the polynomial equation
"""print("Polynomial equation:")
for i, c in enumerate(coefficients):
    print(f"{c}x^{i}", end=" + ")"""

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
