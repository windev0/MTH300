import numpy as np
import matplotlib.pyplot as plt
from math import sin

def lagrange(x, y):
    xplt = np.linspace(x[0], x[-1])
    yplt = np.array([], float)
    
    for xp in xplt:
        yp = 0
        for xi, yi in zip(x, y):
            yp += yi * (np.prod((xp - x[x != xi]) / (xi - x[x != xi])))
        yplt = np.append(yplt, yp)
        
    # courbe
    plt.plot(x, y, 'ro', xplt, yplt, 'b-')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()


x = np.array([0, 20, 40, 60, 80, 100])
y = np.array([26.0, 48.0, 61.6, 71.2, 74.8, 75.2], float)
lagrange(x,y)
