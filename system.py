import numpy as np
# from numpy.linalg import norm, inv, inf
from equaNonLineaire import effacer_console

def matrice_diagonale(A, b):
    try:  
        n, m, p = A.shape[0], A.shape[1], len(b)     # recueil du nb de lignes et de colonnes
        assert n == m == p                           # verifier si la matrice est carree
        for i in range(n):                           # on s'assure de A est diagonale
            for j in range(n):
                if (i == j) and A[i, j] == 0 :       # si la diagonale contient 0
                    raise PermissionError
                if (i != j) and A[i, j] != 0:        # si les autres éléments sont non nuls
                    raise PermissionError
                               
        sol = np.zeros(shape=(n, 1))                 # matrice de n ligne, 1 colonnes initialisé à 0
        for i in range(n):
            sol[i] = b[i] / A[i, i]                  # les éléments de la matrice solution
        print("\n\t La solution est S = \n")
        # print(sol)                  
        return sol                                   # renvoie de la solution                       
    
    except AssertionError:
        print("\n\tERREUR: la matrice n'est pas carree ou la matrice b n'est pas un vecteur") 
    except TypeError or ValueError or NameError:
        print("\n\tUne ou plusieurs donnée(s) entrée(s) est (sont) non valide(s)")   
    except PermissionError:
        print("\n\tERREUR: Matrice non diagonale")
# ==================================================================================================
def matrice_triangulaire_sup(A, b):
    try:
        n, m, p = A.shape[0], A.shape[1], len(b)   # recueil du nb de lignes et de colonnes
        assert n == m == p                         # verifier si la matrice est carree
       
        for i in range(n):                         # verifier si la matrice est triangulaire sup
            for j in range(n):
                if (i > j) and A[i, j] != 0:             # si les autres éléments en dessous de la diagonale sont non nuls
                    raise PermissionError
                
        sol = np.zeros(shape=(n, 1))    # initialisation de la matrice vecteur solution
        sol[n-1] = b[n-1] / A[n-1, n-1]     # la solution z 
        # la remontée
        for i in range(n-2, -1, -1):
            sol[i] = b[i]
            for j in range(i+1, n):
                sol[i] -= A[i, j]*sol[j]
            sol[i] /= A[i, i]
            
        # print("\n\t La solution est S = \n")
        # print(sol) 
        return  sol                            # renvoie de la solution
    except AssertionError:
        print("\n\tERREUR: la matrice n'est pas carree ou la matrice b n'est pas un vecteur")
    except TypeError or ValueError or NameError:
        print("\n\tUne ou plusieurs donnée(s) entrée(s) est (sont) non valide(s)")
    except PermissionError:
        print("\n\tERREUR: Matrice non triangulaire supérieure !")
# ==================================================================================================
def matrice_triangulaire_inf(A, b):
    try:
        n, m, p = A.shape[0], A.shape[1], len(b)   # recueil du nb de lignes et de colonnes
        assert n == m == p                 # verifier si la matrice est carree
        
        for i in range(n):                         # verifier si la matrice est triangulaire sup
            for j in range(n):
                if (i < j) and A[i, j] != 0:             # si les autres éléments en dessous de la diagonale sont non nuls
                    raise PermissionError
                
        sol = np.zeros(shape=(n, 1))    # initialisation de la matrice vecteur solution
        sol[0] = b[0] / A[0, 0]         # l'element x[0,0] de la solution x
        
        # la decsente
        for i in range(1, n):
            sol[i] = b[i]
            for j in range(i):
                sol[i] -= A[i, j]*sol[j]
            sol[i] /= A[i, i]
            
        return sol
    except AssertionError:
        print("\n\tERREUR: la matrice n'est pas carree ou pas triangulaire inférieure ou la matrice b n'est pas un vecteur")  
    except TypeError or ValueError or NameError:
        print("\n\tUne ou plusieurs donnée(s) entrée(s) est (sont) non valide(s)")  
    except PermissionError:
        print("\n\tERREUR: Matrice non triangulaire inférieure !")
# ==============================================================================================
def Gauss_elimination(A, B):
    try:
        n, m = A.shape[0], A.shape[1]   # nb de lignes et de colonnes
        assert n == m                   # verifier si la matrice est carrée
        Aw = np.zeros([n, n+1])         # matrice de travail, matrice augmenté
        Aw[:, 0:n] = A                  # la matrice de travail, de 0 à m-1 colonnes contient A
        Aw[:, n] = B                    # la derniere colonne est le vecteur B 
        sol = np.zeros(n)               # la solution est un vecteur de dimension le nb de ligne de A que la matrice A
        
        # echelonnement
        # boucle en i(pivots) de 1 a n-1 <==> 0, 1, 2, ..., n-2 (ici n-1 a cause de arrage)
        for i in np.arange(0, n-1):    
            # boucle en k (lignes) de i+1 a n <==> 0, 1, 2, ..., n-1 (ici n a cause de arrage) 
            for k in np.arange(i+1, n):
                # formule : Aw[k, :] =  Aw[k, :] - Aw[k, i] / Aw[i, i] * Aw[i, :]
                Aw[k, :] =  Aw[k, :] - Aw[k, i] / Aw[i, i] * Aw[i, :]
        
        # remontee
        # sol[n-1] = Aw[n-1, n] / Aw[n-1, n-1]    # initialisation de la solution
        # for i in np.arange(n-2, -1, -1):        # boucle de n-2 a 0 <==> n-1 a 1 en python
        #     sol[i] = (Aw[i, n] - np.sum(Aw[i,i+1:n]*sol[i+1:n])) / Aw[i, i]
            
        # print("\nX =", sol)
        # print(np.linalg.solve(A, B))
        
        sol = matrice_triangulaire_sup(Aw[:, 0:n], Aw[:, n]) # appel à la procedure pour la remontee
        return sol                              # renvoie de la solution
            
    except AssertionError:
        print("\nLa matrice n'est pas carrée !")
    except ValueError or TypeError or NameError:
        print("\n\nERREUR: Une erreur s'est produite, veuillez reprendre avec d'autres données")
    except ZeroDivisionError:
        print("\n\tERREUR: Division par zéro !")
# ==============================================================================================
def Gauss_jordanV2(A, B):
    try:
        # verifier si la matrice A est carrée
        m, n = np.shape(A)[0], np.shape(A)[1]   # respectivemnt nb de lignes et de colonnes
        assert m == n                           # matrice carree
        
        for i in range(n):
            if A[i, i] == 0:    # le pivot est nul
                # permutation des lignes et colonnes de A et de B
                # on peut aussi juste remplacer la matrice par sa transposée car il s'agit d'une matrice carree
                for k in range(n):  # on parcourt toutes les lignes
                    while A[i, i] == 0: # pour tout pivot nul
                        A[i, k] = A[i+1, k] # nous échangeons les eléments de la ligne courante avec ceux de la ligne suivante
                        B[i] = B[i+1]       # nous échangeons les eléments de la ligne courante avec ceux de la ligne suivante
                   
            # divison par le pivot
            B[i] = B[i] / A[i, i]
            A[i, :] = A[i, :] / A[i, i]
            
            for j in range(n):  # pour chaque colonne de la ligne i
                if j != i:  # seulement les éléments hors de la diagonale
                    B[j] = B[j] - B[i] * A[j, i]    # nvlle valeur de bj
                    A[j, i:n] = A[j, i:n] - (A[i, i:n] * A[j, i])
            
    except AssertionError: 
        print("\nLa matrice n'est pas carrée !")
        
    return B, A      # B est le vecteur solution X, et A la nouvelle matrice dont les pivot sont à 1     
# ==============================================================================================
def Gauss_jordan(A, B):
    A = np.array(A, int)
    B = np.array(B, int)
    n = len(B)
    
    # boucle principale
    for k in range(n):
        # partial pivoting
        if np.fabs(A[k, k]) < 1.0e-12:
            for i in range(k+1, n):
                if np.fabs(A[i, k]) > np.fabs(A[k, k]):
                    for j in range(k, n):
                        A[k, j], A[i, j] = A[i, j], A[k, j]
                    B[k], B[i] = B[i], B[k]
                    break
        # Division of the pivot row
        pivot = A[k, k]
        for j in range(k, n):
            A[k, j] /= pivot
        B[k] /= pivot
        # Elimination loop
        for i in range(n):
            if i == k or A[i, k] == 0: continue
            factor = A[i, k]
            for j in range(k, n):
                A[i, j] -= factor * A[k, j]
            B[i] -= factor * B[k]
            
    return B, A
# ==============================================================================================

# ==============================================================================================
def decomposition_LU(A, m = 1):
    try:
         # verifier si la matrice A est carrée
        m, n = np.shape(A)[0], np.shape(A)[1]   # respectivemnt nb de lignes et de colonnes
        assert m == n
        
        L = np.zeros(shape = (n, n))  # initilaisation de L
        U = np.zeros(shape = (n, n))  # initilaisation de U
        
        if m == 1:          # pour Doolite
            for i in range(n):
                # calcul des élément de la matrice L
                for j in range(i-1):
                    L[i, j] = A[i, j]
                    for k in range(j-1):
                        L[i, j] = L[i, j] - L[i, k] * U[k, j]
                    L[i, j] = L[i, j] - U[j, j]
                    
                # calcul des élément de la matrice U
                for j in range(i, n):
                    U[i, j] = A[i, j]
                    for k in range(i-1):
                        U[i, j] = U[i, j] - L[i, k] * U[k, j]
                
            for i in range(n):
                # rendre la diagonle de L à 1
                L[i, i] = 1
        elif m == 2:        # pour crout
            pass
        else:               # pour cholesky
            pass
            
    except AssertionError:
        print("\nLa matrice A n'est pas carrée")
    except ZeroDivisionError:
        print("\n\tERREUR: Division par zéro")
    except Exception:
        print("\n\tERREUR: Une erreur s'est produite, réessayez avec d'autres données !")
    return L, U     
# ----------------------------------------------------------------------------------------
def calcul_Y_pour_LU(L, B):
    # n = np.size(B)  # la taille du vecteur B
    # Y = np.zeros(shape = (n, 1))        # initialisation de y
    # for i in range(n):
    #     Y[i] = B[i]
    #     for k in range(i-1):
    #         Y[i] = Y[i] - L[i, k] * Y[k]
    
    Y = matrice_triangulaire_inf(L, B) # appel à procedure pour la remontée
    return Y            # renvoie de la valeur de y calculé
# ----------------------------------------------------------------------------------------
def vecteur_solution(U, Y): # permet de déterminer le x à partir du y
    # n = np.size(Y)          # recueil de la taille de y
    # X = np.zeros(shape = (n, 1))    # initialisation du vecteur x à n
    # i = n-1
    # X[i] = Y[i]     
    # for k in range(i+1, n):
    #     X[i] = X[i] - U[i, k] * X[k]
        
    # X[i] = X[i] / U[i, i]
    X = matrice_triangulaire_sup(U, Y) # appel à procedure pour la descente
    return X
# ----------------------------------------------------------------------------------------
def methode_LU_Doolite(A, B):
    L, U = decomposition_LU(A, 1)
    Y = calcul_Y_pour_LU(L, B)
    X = vecteur_solution(U, Y)
    return L, U, Y, X
# ==============================================================================================
def methode_LU_Crout(A, B):
    L, U = decomposition_LU(A, 2)
    Y = calcul_Y_pour_LU(L, B)
    X = vecteur_solution(U, Y)
    return L, U, Y, X
# ==============================================================================================
def methode_LU_Cholesky(A, B):
    L, U = decomposition_LU(A)
    Y = calcul_Y_pour_LU(L, B)
    X = vecteur_solution(U, Y)
    return L, U, Y, X
# ==============================================================================================


A = np.array([
    [4, -1, 1],
    [4, -8, 1],
    [-2, 1, 5]
])
# C = np.diag([1, 5, 9]) # pour la deckaration d'une matrice diagonale
# b = np.array([[7],[-21],[15]])
b = np.array([7,-21,15])

L, U, Y, X = methode_LU_Doolite(A, b)
print(L)
print(U)
print(Y)
print("\n", X)



# x, D = Gauss_jordanV2(A, b)
# x2 = Gauss_elimination(A,b)
# print(D)
# print(x)
# print(x2)
# matrice_diagonale(A, b)
# matrice_diagonale(C, b)
# matrice_triangulaire_sup(A, b)

# TEST
"""A = np.array([
    [1, 3, 1, 1],
    [2, -1, 0, 2],
    [5, 4, 3, -3],
    [0, 1, 4, 2]
])

B = np.array([8, 6, -2, 0])
print(B)"""

"""A = np.array([
    [1, 3, 1, 1],
    [2, -1, 0, 2],
    [5, 0, 3, -3],
    [0, 1, 4, 2]
])

B = np.array([8, 6, 2, 0])"""


# Gauss_elimination(A, B) # appel de la fonction pivot
# X, A = Gauss_jordan(A, B)
# X, A = Gauss_jordanV2(A, B)
# x, n = jaccobi()


# print("cc\n", A[0:3, 1:3])
# for i in range(4):
    # print(A[i][i])
# print(np.shape(A)[0], np.shape(A)[1])
# A, B = B, A
# print(np.shape(A)[0])
# print('---------')

# print('================================')
# print("L = \n\n", L)
# print('================================')
# print("U = \n\n", U)
# print('================================')
# print("Y = \n\n", Y)
# print('================================')
# print("X = \n\n", X)

# effacer_console() # on rend au propre la console

# Test des différentes méthodes
# print('\n')
# print('\n========================== E Q U A T I O N   D U  T Y P E   A * X  =  B ==========================')

# if 1 == 1:
    
#     continuer = 'o' # variable servant de condition de continuation
#     while continuer == 'o':
        
#         # affichage du menu
#         print("\n\t M E N U ")
        
#         print("\n\t 1- ELIMINATION DE GAUSS ")
#         print("\n\t 2- GAUSS JORDAN ")
#         print("\n\t 3- LU - CROUT ")
#         print("\n\t 4- LU - CROUT ")
#         try:
#             choix = int(input("\nFaire un choix entre {1} {2} {3} {4}: "))
#             if choix in [1, 2, 3, 4]:
#                 # appel aux procedures des methodes
#                 if choix == 1:
#                     print("\n\t================ M E T H O D E ================")
#                     # appel fonction
#                 elif choix == 2:
#                     print("\n\t================ M E T H O D E ================")
#                     # appel fonction
#                 elif choix == 2:
#                     print("\n\t================ M E T H O D E ================")
#                     # appel fonction
#                 elif choix == 3:
#                     print("\n\t================ M E T H O D E ================")
#                     # appel fonction
#                 elif choix == 4:
#                     print("\n\t================ M E T H O D E ================")
#                     # appel fonction
#                 else:
#                     pass
#             else:
#                 raise PermissionError
    
#             continuer = str.lower(input("\nVoulez vous continuer (o/n)? ")) # forcer la valeur du choix a etre en minuscule
#             if continuer != 'n' and continuer != 'o':
#                 print("\n\tErreur: Saisie invalide ")
            
#         except PermissionError:
#             print("\nChoix invalide !")
            
#         effacer_console() # on rend au propre la console

# else:
#     print("\nAucune solution pour cette équation !")

# print('\n')
