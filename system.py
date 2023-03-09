import numpy as np
import sys
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
        sys.exit('Fin du programme') 
    except TypeError or ValueError or NameError:
        print("\n\tUne ou plusieurs donnée(s) entrée(s) est (sont) non valide(s)")   
        sys.exit('Fin du programme') 
    except PermissionError:
        print("\n\tERREUR: Matrice non diagonale")
        sys.exit('Fin du programme') 
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
        sys.exit('Fin du programme') 
    except TypeError or ValueError or NameError:
        print("\n\tUne ou plusieurs donnée(s) entrée(s) est (sont) non valide(s)")
        sys.exit('Fin du programme')
    except PermissionError:
        print("\n\tERREUR: Matrice non triangulaire supérieure !")
        sys.exit('Fin du programme')
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
        sys.exit('Fin du programme')  
    except TypeError or ValueError or NameError:
        print("\n\tUne ou plusieurs donnée(s) entrée(s) est (sont) non valide(s)")
        sys.exit('Fin du programme')    
    except PermissionError:
        print("\n\tERREUR: Matrice non triangulaire inférieure !")
        sys.exit('Fin du programme') 
# ==============================================================================================
def Gauss_elimination(A, B):
    try:
        n, m = A.shape[0], A.shape[1]   # nb de lignes et de colonnes
        assert n == m                   # verifier si la matrice est carrée
        Aw = np.zeros([n, n+1])         # matrice de travail, matrice augmenté
        Aw[:, 0:n] = A                  # la matrice de travail, de 0 à m-1 colonnes contient A
        Aw[:, n] = B                    # la derniere colonne est le vecteur B 
        sol = np.zeros(n)               # la solution est un vecteur de dimension le nb de ligne de A que la matrice A
        
        # if Aw == None:
        #     sys.exit()
            
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
        sys.exit('Fin du programme')  
    except ValueError or TypeError or NameError:
        print("\n\nERREUR: Une erreur s'est produite, veuillez reprendre avec d'autres données")
        sys.exit('Fin du programme')  
    except ZeroDivisionError:
        print("\n\tERREUR: Division par zéro !")
        sys.exit('Fin du programme')  
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
        sys.exit('Fin du programme')  
        
    return B, A      # B est le vecteur solution X, et A la nouvelle matrice dont les pivot sont à 1     
# ==============================================================================================
def gauss_jordan(A, b):
    n = len(b)
    Ab = np.hstack([A, b.reshape(-1, 1)]) # matrice augmentee
    # trianguler la matrice augmentee
    for i in range(n):
        # recherche du pivot maximal dans la colonne i
        pivot = i
        for j in range(i+1, n):
            if abs(Ab[j, i]) > abs(Ab[pivot, i]):
                pivot = j
        # échanger les lignes i et pivot
        Ab[[i, pivot]] = Ab[[pivot, i]]
        # normaliser la ligne i
        Ab[i] = Ab[i] / Ab[i, i]
        # soustraire la ligne i des autres lignes
        for j in range(n):
            if j != i:
                Ab[j] = Ab[j] - Ab[j, i] * Ab[i]
    # séparer la matrice étendue en deux matrices
    x = Ab[:, n]
    A = Ab[:, :n]
    return A, x
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
        sys.exit('Fin du programme')  
    except ZeroDivisionError:
        print("\n\tERREUR: Division par zéro")
        sys.exit('Fin du programme')  
    except Exception:
        print("\n\tERREUR: Une erreur s'est produite, réessayez avec d'autres données !")
        sys.exit('Fin du programme')  
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
        sys.exit('Fin du programme')  
    except ZeroDivisionError:
        print("\n\tERREUR: Division par zéro")
        sys.exit('Fin du programme')  
    except Exception:
        print("\n\tERREUR: Une erreur s'est produite, réessayez avec d'autres données !")
        sys.exit('Fin du programme') 
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
        sys.exit('Fin du programme')
    except PermissionError:
        print("\n\tERREUR: Matrice non définie positive")
        sys.exit('Fin du programme')
# ===========================================================================================
def jacobi(A, b, nmax = 100, tol = 1e-5):
    try:
        # verifier si la matrice A est carrée
        m, n, p = np.shape(A)[0], np.shape(A)[1], len(b)   # respectivemnt nb de lignes et de colonnes
        assert m == n == p
        
        D = np.diag(np.diag(A))      # D est la diagonale de A
        E = -np.tril(A,-1)
        F = - np.triu(A, 1)
        
        x = np.zeros((n), 'float')
        cpt = 0
        # critere d'arret : ||AX-b|| ~~ 0
        while cpt < nmax and np.max(abs(np.matmul(A,x)-b)) > tol:
            x = np.dot(np.linalg.inv(D), np.dot(E+F, x) + b )
            cpt += 1
        if cpt < nmax:
            print("\n\tConvergence atteinte")
        else:
            print("\n\tNombre maximal d'itérations atteinte")
        return x
            
    except AssertionError:
        print("\nLa matrice A n'est pas carrée ou b n'est pas une matrice colonne")
        sys.exit('Fin du programme')
    except ZeroDivisionError:
        print("\n\tERREUR: Division par zéro")
        sys.exit('Fin du programme')
    except Exception:
        print("\n\tERREUR: Une erreur s'est produite, réessayez avec d'autres données !")
        sys.exit('Fin du programme')
# ===========================================================================================
def gauss_seidel(A, b, x = [0, 0, 0, 0], max_iter=100, tol=1e-5):
    try:
        # verifier si la matrice A est carrée
        m, n, p = np.shape(A)[0], np.shape(A)[1], len(b)   # respectivemnt nb de lignes et de colonnes
        assert m == n == p
        cpt = 0
        while cpt < max_iter:
            x_new = x.copy()
            for j in range(n):
                s1 = sum(A[j][k] * x_new[k] for k in range(j))
                s2 = sum(A[j][k] * x[k] for k in range(j+1, n))
                x_new[j] = (b[j] - s1 - s2) / A[j][j]
            if all(abs(x_new[i] - x[i]) < tol for i in range(n)):
                return x_new
            x = x_new
        return x
    except AssertionError:
        print("\nLa matrice A n'est pas carrée ou b n'est pas une matrice colonne")
        sys.exit('Fin du programme')
    except ZeroDivisionError:
        print("\n\tERREUR: Division par zéro")
        sys.exit('Fin du programme')
    except Exception:
        print("\n\tERREUR: Une erreur s'est produite, réessayez avec d'autres données !")
        sys.exit('Fin du programme')
# ===========================================================================================
def thomas(A, d):
    try:
        # print(np.diag(A, k=-1))
        a = np.diag(A, k=-1)     # diagonale  inférieure
        b = np.diag(A, k=0)      # diagonale principale
        c = np.diag(A, k=1)      # diagonale supérieure
        m, p = matrice_tri_diag(b, c, a)
        assert m == True
        n = len(d)
        # étape de substitution avant
        for i in range(1, n):
            b[i] = b[i] - (a[i-1] * c[i-1]) / b[i-1]
            d[i] = d[i] - (d[i-1] * c[i-1]) / b[i-1]
        # résolution du système
        x = np.zeros(n)
        x[-1] = d[-1] / b[-1]
        for i in range(n-2, -1, -1):
            x[i] = (d[i] - a[i] * x[i+1]) / b[i]
        return x
    except AssertionError:
        print("\n\tErreur:Matrice non tri-diagonale")
        sys.exit('Fin du programme')
    except IndexError: 
        print("\n\tErreur: Demande d'accès à un élément de la matrice qui n'existe pas")
        sys.exit('Fin du programme')
# =================================================================================================================
def matrice_tri_diag(main_diag, upper_diag, lower_diag):
    
    # créer la matrice tri-diagonale
    tri_diag = np.diagflat(main_diag) + np.diag(upper_diag, k=1) + np.diag(lower_diag, k=-1)

    return True, tri_diag # retourner la matrice tri-diagonale
# ===================================================================================================================
def decompoThomas(A,b):
    try:
        # verifier si la matrice A est carrée
        m, n, p = np.shape(A)[0], np.shape(A)[1], len(b)   # respectivemnt nb de lignes et de colonnes
        if m != n and n != p: raise PermissionError
        
        a = np.diag(A, k=-1)     # diagonale  inférieure
        b = np.diag(A, k=0)      # diagonale principale
        c = np.diag(A, k=1)      # diagonale supérieure
        m, h = matrice_tri_diag(b, c, a)
        assert m == True
        
        # n = len(b)
        B, C = np.zeros(shape=(n,n)), np.zeros(shape=(n,n))
        
        C[0,0] = A[0,0]
        
        for i in range(n-1):
            C[i,i+1] = A[i,i+1]
        for i in range(n):
            B[i,i] = 1
            if(i >= 1):
                B[i,i-1] = A[i,i-1]/C[i-1,i-1]
                C[i,i] = A[i,i] - B[i,i-1]*A[i-1,i]
        return B, C
    except AssertionError:
        print("\n\tErreur:Matrice non tri-diagonale")
        sys.exit('Fin du programme')
    except PermissionError:
        print("\n\tErreur:Matrice non tri-diagonale")
        sys.exit('Fin du programme')
    except IndexError: 
        print("\n\tErreur: Demande d'accès à un élément de la matrice qui n'existe pas")
        sys.exit('Fin du programme')
# ===================================================================================================================
def resolveThomas(A,b):
    n = np.shape(A)[0]
    B, C = decompoThomas(A,b)
    D, X = np.zeros(shape=(n,1)), np.zeros(shape=(n,1))
    
    D[0,0] = b[0]
    for i in range(1,n):
        D[i,0] = b[i] - B[i,i-1]*D[i-1,0]
    X[n-1,0] = D[n-1,0]/C[n-1,n-1]
    
    for i in range(n-2, -1, -1):
        X[i,0] = (D[i,0] - C[i,i+1]*X[i+1,0])/C[i,i]
        
    return X
# ===================================================================================================================


if __name__ == "__main__":
    effacer_console() # on rend au propre la console

    # Test des différentes méthodes
    print('\n')
    print('\n========================== E Q U A T I O N   D U  T Y P E   A * X  =  B ==========================')

    # SAISIE DES MATRICES
    # A = np.array([[-1, 1, 2],[1, -1, 3],[0, 1, 1]])
    # b = np.array([2,3,2])
    
    A = np.array([[7, 1, 2, 3],[1, 7, 2, 3],[1, 2, 7, 3],[1, 2, 3, 7]])
    b = np.array([5,-7,3,-5])
    
    # A = np.array([[1, 4, -3],[2, -3, -1],[3, 2, 1]])
    # b = np.array([1, 2, 3])
    
    if int(np.linalg.det(A)) != 0: # cas ou la matrice est inversible

        continuer = 'o' # variable servant de condition de continuation
        valide = True
        while continuer == 'o':
            
            # A = np.array([[-1, 1, 2],[1, -1, 3],[0, 1, 1]])
            # b = np.array([2,3,2])
            
            A = np.array([[7, 1, 2, 3],[1, 7, 2, 3],[1, 2, 7, 3],[1, 2, 3, 7]])
            b = np.array([5,-7,3,-5])
            
            # A = np.array([[1, 4, -3],[2, -3, -1],[3, 2, 1]])
            # b = np.array([1, 2, 3])
            
            # affichage du menu
            print("\n\t M E N U ")

            print("\n\t 1- ELIMINATION DE GAUSS ")
            print("\n\t 2- GAUSS JORDAN ")
            print("\n\t 3- LU - DOOLITE ")
            print("\n\t 4- LU - CROUT ")
            print("\n\t 5- LU - CHOLESKY ")
            print("\n\t 6- JACOBI ")
            print("\n\t 7- GAUSS - SEIDEL ")
            print("\n\t 8- THOMAS ")
            print("\n\t 9- QUITTER ")
            
            try:
                choix = int(input("\nFaire un choix entre {1} {2} {3} {4} {5} {6} {7} {8} {9}: "))
                if choix in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
                    # appel aux procedures des methodes
                    if choix == 1:
                        print("\n\t================ M E T H O D E   D E   G A U S S ================")
                        print("\nX = \n", Gauss_elimination(A, b))
                        # print("\nX = \n", np.linalg.solve(A, b))
                    elif choix == 2:
                        print("\n\t================ M E T H O D E   D E   G A U S S - J O R D A N ================")
                        x, newA = Gauss_jordanV2(A,b)
                        print("\nA = \n", newA)
                        print("\nX = ", x)
                        # print("\nX = ", np.linalg.solve(A, b))
                    elif choix == 3:
                        print("\n\t================ M E T H O D E   L U   ( D O O L I T E ) ================")
                        L, U, y, x = doolittle(A, b)
                        print("\nL = \n", L)
                        print("\nU = \n", U)
                        print("\nY = ", y)
                        print("\nX = ", x)
                        # print("\nX = \n", np.linalg.solve(A, b))
                    elif choix == 4:
                        print("\n\t================ M E T H O D E   L U   ( C R O U T ) ================")
                        L, U, y, x = crout(A, b)
                        print("\nL = \n", L)
                        print("\nU = \n", U)
                        print("\nY = \n", y)
                        print("\nX = \n", x)
                        # print("\nX = \n", np.linalg.solve(A, b))
                    elif choix == 5:
                        print("\n\t================ M E T H O D E   L U   ( C H O L E S K Y ) ================")
                        L, transposedL, y, x = cholesky(A, b)
                        print("\nL = \n", L)
                        print("\ntranspL = \n", transposedL)
                        print("\nY = \n", y)
                        print("\nX = \n", x)
                        print("\nX = \n", np.linalg.solve(A, b))
                    elif choix == 6:
                        print("\n\t================ M E T H O D E    D E   J A C O B I  ================")
                        print("\nX = \n", jacobi(A, b))
                        # print("\nX = \n", np.linalg.solve(A, b))
                    elif choix == 7:
                        print("\n\t================ M E T H O D E    D E   G A U S S - S E I D E L  ================")
                        print("\nX = \n", gauss_seidel(A, b))
                        # print("\nX = \n", np.linalg.solve(A, b))
                    elif choix == 8:
                        print("\n\t================ M E T H O D E    D E   T H O M A S  ================")
                        print("\nX = \n", resolveThomas(A,b))
                        # print("\nX = \n", thomas(A,b))
                        # print("\nX = \n", np.linalg.solve(A, b))
                    elif choix == 9:
                        continuer = 'n'
                    else:
                        pass
                else:
                    raise PermissionError
                
                if continuer == 'o':
                    continuer = str.lower(input("\nVoulez vous continuer (o/n)? ")) # forcer la valeur du choix a etre en minuscule
                    valide = True
                    if continuer != 'n' and continuer != 'o':
                        valide = False
                        print("\n\tErreur: Saisie invalide ")
                
                    while valide == False:        
                        continuer = str.lower(input("\nVoulez vous continuer (o/n)? ")) # forcer la valeur du choix a etre en minuscule
                        if continuer != 'n' and continuer != 'o':
                            valide = False
                            print("\n\tErreur: Saisie invalide ")
                        else:
                            valide = True
                    
            except PermissionError:
                print("\nChoix invalide !")
            except ValueError:
                print("\n\tERREUR: Saisie invalide !")
                    
                # effacer_console() # on rend au propre la console
        print('\n')
    else:
        print("\n\tERREUR: Matrice non inversible !")