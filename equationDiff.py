import numpy as np
import matplotlib.pyplot as plt

def euler(x0, y0, pas, nb_etapes):
    # Les tableaux pour stocker les résultats
    x = np.linspace(x0, x0 + nb_etapes*pas, nb_etapes+1)
    y = np.zeros(nb_etapes+1)
    y[0] = y0

    # Boucle pour itérer la méthode d'Euler
    for i in range(nb_etapes):
        y[i+1] = y[i] + pas * y[i]
        
    return x, y # les couples de points (x,y)
# ================================================================================================
def runge_kunta4(x0, y0, pas, nb_etapes):
    # Les tableaux pour stocker les résultats
    x = np.linspace(x0, x0 + nb_etapes*pas, nb_etapes+1)
    y = np.zeros(nb_etapes+1)
    y[0] = y0

    # Boucle pour itérer la méthode de runge kunta
    for i in range(nb_etapes):
        k1 = pas * y[i]
        k2 = pas * (y[i] + k1/2)
        k3 = pas * (y[i] + k2/2)
        k4 = pas * (y[i] + k3)
        y[i+1] = y[i] + (k1 + 2*k2 + 2*k3 + k4) / 6
        
    return x, y # les couples de points (x,y)
# =========================================================================================================
def euler_modified(f, y0, t0, tf, h):  
    t = np.arange(t0, tf+h, h)
    y = np.zeros(len(t))
    y[0] = y0
    for i in range(1, len(t)):
        y_pred = y[i-1] + h*f(t[i-1], y[i-1])
        y[i] = y[i-1] + (h/2)*(f(t[i-1], y[i-1]) + f(t[i], y_pred))
    return t, y

# définir la fonction f(t, y) de l'équation différentielle
def f(x, y):
    return -y + x + 1
# =========================================================================================================

# Affichage des résultats
x1, y1  = euler(0, 1, 0.1, 10)
x2, y2  = runge_kunta4(0, 1, 0.1, 10)
x3, y3 = euler_modified(f, 1, 0, 2, 0.1)

# courbe des fonctions
plt.plot(x1, y1, 'b--', label='Euler')
plt.plot(x2, y2, 'r--', label='Runge kunta 4')
plt.plot(x3, y3, 'y--', label='Euler modifé')
plt.legend()
plt.title('EQUATION DIFFERENTIELLE')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
