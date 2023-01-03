import numpy as np
from math import sqrt
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
def calcul_Y_pour_LU(L, B):
    Y = matrice_triangulaire_inf(L, B) # appel à procedure pour la remontée
    return Y            # renvoie de la valeur de y calculé
# ----------------------------------------------------------------------------------------
def vecteur_solution(U, Y): # permet de déterminer le x à partir du y
    X = matrice_triangulaire_sup(U, Y) # appel à procedure pour la descente
    return X
# ==============================================================================================
def doolittle(A, B):
    try:
         # verifier si la matrice A est carrée
        m, n, p = np.shape(A)[0], np.shape(A)[1], len(B)   # respectivemnt nb de lignes et de colonnes
        assert m == n == p
        
        L = np.zeros(shape = (n, n))  # initilaisation de L
        U = np.zeros(shape = (n, n))  # initilaisation de U
    
        for j in range(n):
            U[0][j] = A[0][j]
            L[j][0] = A[j][0] / U[0][0]
        for i in range(1, n):
            for j in range(i, n):
                s1 = sum(U[k][j] * L[i][k] for k in range(i))
                U[i][j] = A[i][j] - s1
            for j in range(i + 1, n):
                s2 = sum(U[k][i] * L[j][k] for k in range(i))
                L[j][i] = (A[j][i] - s2) / U[i][i]
        y = [0.0] * n
        for i in range(n):
            s3 = sum(L[i][k] * y[k] for k in range(i))
            y[i] = B[i] - s3
        x = [0.0] * n
        for i in range(n - 1, -1, -1):
            s4 = sum(U[i][k] * x[k] for k in range(i + 1, n))
            x[i] = (y[i] - s4) / U[i][i]
            
    except AssertionError:
        print("\nLa matrice A n'est pas carrée ou b n'est pas une matrice colonne")
    except ZeroDivisionError:
        print("\n\tERREUR: Division par zéro")
    except Exception:
        print("\n\tERREUR: Une erreur s'est produite, réessayez avec d'autres données !")
    return L, U, y, x
# ==========================================================================================
def crout(A,b):
    try:
         # verifier si la matrice A est carrée
        m, n, p = np.shape(A)[0], np.shape(A)[1], len(b)  # respectivemnt nb de lignes et de colonnes
        assert m == n == p
        
        L = np.zeros(shape = (n, n))  # initilaisation de L
        U = np.zeros(shape = (n, n))  # initilaisation de U
   
        for j in range(n):
            U[j][j] = 1.0
            for i in range(j, n):
                s1 = sum(U[k][j] * L[i][k] for k in range(j))
                L[i][j] = A[i][j] - s1
            for i in range(j + 1, n):
                s2 = sum(U[k][j] * L[j][k] for k in range(j))
                U[j][i] = (A[j][i] - s2) / L[j][j]
        Y = calcul_Y_pour_LU(L, b)
        X = vecteur_solution(U, Y)
        return L, U, Y, X
                
    except AssertionError:
        print("\nLa matrice A n'est pas carrée ou b n'est pas une matrice colonne")
    except ZeroDivisionError:
        print("\n\tERREUR: Division par zéro")
    except Exception:
        print("\n\tERREUR: Une erreur s'est produite, réessayez avec d'autres données !")
# ======================================================================================
def cholesky(A, b):
    try:
        n = len(A)
        L = np.zeros(shape = (n, n)) 
        # vérifier si la matrice est symétrique
        for i in range(n):
            for j in range(n):
                if A[i,j] != np.transpose(A)[i,j]:
                    raise AssertionError
        # vérification de la définition positive
        sum_items = 0
        for i in range(n):
            for j in range(n):
                if i != j:
                    sum_items = abs(A[i,j]) # somme des termes hors de la diagonale
            if A[i, i] <= sum_items:
                raise PermissionError
        # calcul
        for i in range(n):
            for j in range(i + 1):
                s = sum(L[i][k] * L[j][k] for k in range(j))
                if (i == j):
                    L[i][j] = sqrt(A[i][i] - s)
                else:
                    L[i][j] = (1.0 / L[j][j] * (A[i][j] - s))
        Y = calcul_Y_pour_LU(L, b)
        X = vecteur_solution(np.transpose(L), Y)
        return L, np.transpose(L), Y, X
    except AssertionError:
        print("\n\tERREUR: Matrice non symétrique")
    except PermissionError:
        print("\n\tERREUR: Matrice non définie positive")
# ===========================================================================================
A = np.array([
    [4, 1, 1],
    [1, 5, 2],
    [1, 2, 6]
])
b = np.array([7,-21,15])
L, transposedL, Y, X2 = crout(A, b)
print(L)
print(transposedL)
print(Y)
print(X2)


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
