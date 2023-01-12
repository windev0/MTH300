import numpy as np
import matplotlib.pyplot as plt
# from math import sin

def lagrange(x, y):
    xplt = np.linspace(x[0], x[-1])
    yplt = np.array([], float)
    
    for xp in xplt:
        yp = 0
        for xi, yi in zip(x, y):
            yp += yi * (np.prod((xp - x[x != xi]) / (xi - x[x != xi])))
        yplt = np.append(yplt, yp)
    
    return xplt, yplt
        
    # # courbe
    # plt.plot(x, y, 'go', xplt, yplt, 'r--')
    # plt.title('Lagrange')
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.show()

# ===============================================================================================
def interpolation_newton(x, y, x_new):
    
    # xplt = np.linspace(x[0], x[-1])
    # yplt = np.array([], float)
    
    # Initialiser la matrice de différences divisées
    f = [[0 for _ in range(len(x))] for _ in range(len(x))]
    # Remplir la première colonne de la matrice avec les valeurs de y
    for i in range(len(x)):
        f[i][0] = y[i]
    # Calculer les différences divisées
    for i in range(1, len(x)):
        for j in range(i, len(x)):
            f[j][i] = (f[j][i-1] - f[j-1][i-1]) / (x[j] - x[j-i])
    # Calculer le polynôme d'interpolation de Newton
    
    p = 0
    for i in range(len(x)):
        term = f[i][i]
        for j in range(i):
            term *= (x_new - x[j])
        p += term
    return p
    
    # # courbe
    # plt.plot(x, y, 'ro', x_new, p, 'b--')
    # plt.title('Newton')
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.show()
# ===============================================================================================
def computeNewtonCoefs(x, y):
    n = len(x) # nb de points
    c = np.zeros(x, 'float')
    
    z = np.zeros(len(x), 'float')
    c[0] = y[0]
    for j in range(1, n):
        z[0] = 1.0
        for k in range(j):
            z[k+1] = x[j]-x[k]
        pm1 = 0.0
        w = np.zeros(j, 'float')
        for k in range(j):
            w[k] = np.prod(z[:k+1])
            pm1 += (c[k]*w[k])
        c[j] = (y[j]-pm1)/np.prod(z[:j+1])
    return c
# ===============================================================================================

# ===============================================================================================

    
x = np.array([0, 20, 40, 60, 80, 100])
y = np.array([26.0, 48.0, 61.6, 71.2, 74.8, 75.2], float)
xplt, yplt = lagrange(x, y)
x_new = np.linspace(x[0], x[-1], 50)
p = interpolation_newton(x,y,x_new)

# courbe
plt.plot(x, y, 'go', xplt, yplt, 'r--', label= "Lagrange")
plt.plot(x, y, 'yo', x_new, p, 'b.', label= 'Newton')
plt.legend()
plt.title('INTERPOLATION')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
# lagrange(x,y)
# interpolation_newton(x,y, x_new)
