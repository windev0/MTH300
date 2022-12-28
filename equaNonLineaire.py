from math import ceil, exp, fabs, log, log10, sqrt
import os

def nbre_solutions(x1, x2, pas):    # cette fonction retourne le nb de solution apres balayage
    a = x1 + pas                    # initialisation de a, et servira de borne sup temporaire
    cpt = 0                         # initialisation de cpt, pour compter les solutions
    while a <= x2:                  # tant que nous ne sommes pas arrivé à la borne supérieure
        x1 += pas                   # x1 est augmente d'un pas vers la gauche
        a += pas                    # a egalement est augmente d'un pas vers la gauche
        if f(x1)*f(a) < 0:          # si la fonction change de signe sur [x1, a], alors...
            cpt += 1                # il y a existence d'une solution, on incrémente cpt
    return cpt                      # on retourne le nombre de solutions trouvées
#---------------------------------------------------------------------------------------------
def effacer_console():
    cmd = 'clear'
    if os.name in ('nt', 'dos'):
        cmd = 'cls' # cas de windows, la commande sera cls et non clear
    os.system(cmd) # effacer la console
#---------------------------------------------------------------------------------------------

def saisie_bornes():
    try:
        x1 = int(input("\n\tVeuillez saisir la borne inférieure X1 = "))
        x2 = int(input("\n\tVeuillez saisir la borne supérieure X2 = "))
        assert x1 < x2
        if x1 < x2:
            return x1, x2
    except AssertionError:
        print('\nErreur: la borne inférieure doit etre strictement inférieure à la borne supérieure')
    except ValueError or TypeError:
        print("\nErreur: Veuillez bien choisir un entier !")
#---------------------------------------------------------------------------------------------
def saisie_X0():
    # cette fonction permet de saisir la valeur initiale
    try:
        x0 = float(input("\n\tVeuillez saisir une valeur initiale: "))
        return x0
    except ValueError or TypeError:
        print("\n\tErreur: Saisie invalide, choisir une autre valeur entière!")
#---------------------------------------------------------------------------------------------
def calcul_N(x1, x2, e):
    # cette fonction permet de calculer le nombre maximal d'iteration avec dichotomie
    try:
        N = ceil((log(x2-x1) - log(2*e))/log(2))
        print('\n\tLe nombre maximal d\'itération est N =',N)
        return N
    except ValueError or TypeError:
        print("\n\tErreur: données invalide !")
#---------------------------------------------------------------------------------------------
def saisie_tolerance():
    try:
        e = float(input("\n\tQuelle est la tolérence e = "))
        return e
    except ValueError or TypeError:
        print("\n\tErreur: Saisie invalide, essayer une autre valeur !")
#---------------------------------------------------------------------------------------------
def convergence(x1, x2, e, nb, xm = 1):
    try:
        if nb == 1: # 1 s'il s'agit de la convergence de la methode dichotomie
            return (fabs(x2-x1) / 2*fabs(xm)) < e
        else:   # autre entier pour la convergence des methodes de points fixes, NEWTON
            return (fabs(x2-x1) / fabs(xm)) < e
    except ZeroDivisionError:
        print("\n\tErreur: Division par zéro, essayer une autre valeur car", xm, ' = 0')
    except TypeError or ValueError:
        print("\n\tErreur: donnée(s) invalide(s) !")
#---------------------------------------------------------------------------------------------
f = lambda x : exp(-x) - x
    # return x**3 + x**2 - 3*x - 3
    # return pow(x,2) - 2
    # return sqrt(2*x + 3)
    # return pow(x,2) - 2*x - 3
    # return pow(x,2) - 3*x + 2
#---------------------------------------------------------------------------------------------
g = lambda x : exp(-x) + x
#---------------------------------------------------------------------------------------------

devf = lambda x : -exp(-x) - 1
    # reperesente la fonction derivee 
#---------------------------------------------------------------------------------------------
# permet de tronquer un nombre nb a un m chiffres apres la virgule
def tronquer(nb, m):
    nb = str(nb)
    partie_ent, partie_decim = nb.split('.')
    nb = '.'.join([partie_ent, partie_decim[:m]])
    return float(nb) 

#=============================================================================================

def secante():
    try:
        x0 = saisie_X0()
        x1 = saisie_X0()
        assert x0 < x1
        
        if f(x0) == 0:
            print("\tla solution est: X* = {}".format(x0))
        elif f(x1) == 0:
            print("\tla solution est: X* = {}".format(x1))
        else:
            e = saisie_tolerance()
            N = int(input("\n\tVeuillez saisir le nombre total d'itération N = "))
            x2 = x1 - ((f(x1)*(x1-x0))/(f(x1)-f(x0)))  # initialisation de x2
            i = 1
            while not(convergence(x1, x2, e, 2, x2)) and i < N:
                # mise a jour des valeurs 
                x0 = x1
                x1 = x2
                x2 = x1 - ((f(x1)*(x1-x0))/(f(x1)-f(x0)))
                i += 1
            # Saisie des résultats
            print('\n-----------------------------------------------')
            if convergence(x1, x2, e, 2, x2):
                print("\tconvergence atteinte en ", i, 'itérations')
            if i == N and (not convergence(x1, x2, e, 2, x2)):
                print("\tConvergence non atteinte en {} itérations".format(i))
            print("\tla solution est: X* = {}".format(x2))
            print('-----------------------------------------------')
    except TypeError or ValueError or NameError:
        print("\n\tUne ou plusieurs donnée(s) entrée(s) est (sont) non valide(s)")
    except AssertionError:
        print("\n\tErreur: La premiere valeur initiale doit etre inférieure à la seconde!")
    except ZeroDivisionError:
        print("\n\tErreur: Division par zéro, car f(", x1, ') - f(', x0, ') = 0; essayez avec d\'autres valeurs')

#---------------------------------------------------------------------------------------------

def newton():
    try:
        e = saisie_tolerance()
        N = int(input("\n\tVeuillez saisir le nombre total d'itération N = "))
        x0 = saisie_X0()
        x1 = x0 - (f(x0) / devf(x0)) # initialisation de x1
        
        if f(x0) == 0:
                print("\tla solution est: X* = {}".format(x0))
        elif f(x1) == 0:
            print("\tla solution est: X* = {}".format(x1))
        else:
            i = 1
            while not(convergence(x0, x1, e, 2, x1)) and i < N:
                # mise a jour des valeurs 
                x0 = x1
                x1 = x0 - (f(x0) / devf(x0))
                i += 1
            # Saisie des résultats
            print('\n-----------------------------------------------')
            if convergence(x0, x1, e, 2, x1):
                print("\tconvergence atteinte en ", i, 'itérations')
            if i == N and (not convergence(x0, x1, e, 2, x1)):
                print("\tConvergence non atteinte en {} itérations".format(i))
            print("\tla solution est: X* = {}".format(x1))
            print('-----------------------------------------------')
            
    except TypeError or ValueError and NameError:
        print("\n\tErreur: Une ou plusieurs donnée(s) entrée(s) est (sont) non valide(s)")
    except ZeroDivisionError:
        print("\n\tErreur: Division par zéro, car f'(", x0, ') = 0', "essayer avec d'autres valeurs !")

# -----------------------------------------------------------------------------------------------------------------------------------------

def point_fixe():
    try:
        x0 = saisie_X0()
        x1 = g(x0)  # initialisation de x1
        e = saisie_tolerance()
        N = int(input("\n\tVeuillez saisir le nombre total d'itération N = "))

        if f(x0) == 0:
            print("\tla solution est: X* = {}".format(x0))
        elif f(x1) == 0:
            print("\tla solution est: X* = {}".format(x1))
        else:
            i = 1
            while not(convergence(x0, x1, e, 2, x1)) and i < N:
                # mise a jour des valeurs 
                x0 = x1
                x1 = f(x0)
                i += 1
            # Saisie des résultats
            print('\n-----------------------------------------------')
            if convergence(x0, x1, e, 2, x1):
                print("\tconvergence atteinte en", i, 'itérations')
            if i == N-1 and (not convergence(x0, x1, e, 2, x1)):
                print("\tConvergence non atteinte en {} itérations".format(i))
            print("\tla solution est: X* = {}".format(x1))
            print('-----------------------------------------------')
        
    except TypeError or ValueError:
        print("\n\tUne ou plusieurs donnée(s) entrée(s) est (sont) non valide(s)")    
    except ZeroDivisionError:
        print("\n\tErreur: Division par zéro, essayer avec d'autres valeurs !")
    except OverflowError:
        print("\n\tErreur: Dépassement de capacité !")
#---------------------------------------------------------------------------------------------
        
def dichotomie():
    try:
        x1, x2 = saisie_bornes()    # saisie des bornes de l'intervale
        if f(x1) == 0:
            print('\n-------------------------------------------------------------------------------------------')
            print("\tLa solution est x =", x1)
            print('-------------------------------------------------------------------------------------------')
        elif f(x2) == 0:
            print('\n-------------------------------------------------------------------------------------------')
            print("\tLa solution est x =", x2)
            print('-------------------------------------------------------------------------------------------')
        else:
            if f(x1)*f(x2) < 0:         # la fonction change de signe sur l'intervalle donné
                print("\n\tLa fonction est monotone sur l'intervalle [{},{}]".format(x1, x2))
                print("\ton a {} solution(s)".format(nbre_solutions(x1, x2, 0.01))) # affichage du nombre de solutions trouvées
                xm = (x1 + x2 )/2           # initialisation du milieu de l'intervalle
                if f(xm) == 0:              # xm est la racine cherchee
                    print("\n\t X1 = {}\t\t\t X2 = {}\t\t\t Xm = {}\n f(X1) = {}\t\t f(X2) = {}\t\t f(Xm) = {}".format(x1, x2, xm, f(x1), f(x2), f(xm)))
                else:
                    e = saisie_tolerance()   # saisie du critère d'arret et le nb max d'itération
                    N = calcul_N(x1, x2, e)  # calcul du nb max d'iteration a l'aide de la tolerance
                    i = 1
                    while (not (convergence(x1, x2, e, 1, xm))) and i < N and f(x1)*f(x2) < 0:
                        if f(x1)*f(xm) < 0:
                            x2 = xm             # changement de la borne supérieure
                        if f(xm)*f(x2) < 0:
                            x1 = xm             # changement de la borne inférieure
                        xm = (x1 + x2 )/2       # mise à jour du milieu de l'intervalle
                        i += 1                  # incrémentation du nb d'itérations
                    print('\n-------------------------------------------------------------------------------------------')
                    if i == N and (not convergence(x1, x2, e, 1, xm)):    # nb d'itération atteinte
                        print("\tConvergence non atteinte en {} itérations".format(i))
                    if convergence(x1, x2, e, 1, xm):           # si la convergence est atteinte, alors...
                        print("\tConvergence atteinte en", i, 'itérations')
                    # Saisie des résultats
                    print("\n\t X1 = {}\t\t\t X2 = {}\t\t\t Xm = {}\n\t f(X1) = {}\t\t f(X2) = {}\t\t f(Xm) = {}".format(\
                    tronquer(x1, 6), tronquer(x2, 6), tronquer(xm, 6), tronquer(f(x1), 6), tronquer(f(x2), 6), tronquer(f(xm), 6)))
                print('-------------------------------------------------------------------------------------------')
            else:
                print("\n\tla fonction ne change pas de signe sur l'intervalle [{},{}]".format(x1, x2))

    except TypeError or ValueError or NameError:
        print("\n\tErreur: une ou plusieurs donnée(s) entrée(s) est (sont) non valide(s) !")
    except ZeroDivisionError:
        print("\n\tErreur: Division par zéro, essayer avec d'autres valeurs !")
# -----------------------------------------------------------------------------------------------------------------------------------------
def main():
    effacer_console() # on rend au propre la console

    # Test des différentes méthodes
    print('\n')
    print('\n========================== E Q U A T I O N   D U  T Y P E   F ( X ) = 0 ==========================')

    a = nbre_solutions(0, 8, 0.01) # détermination du nombre de solutions
    if a > 0:
        print("\n Il existe {} solution (s)".format(a)) # affichage du nombre de solutions après balayage
        
        continuer = 'o' # variable servant de condition de continuation
        while continuer == 'o':
            
            # affichage du menu
            print("\n\t M E N U ")
            
            print("\n\t 1- DICHOTOMIE ")
            print("\n\t 2- SECANTE ")
            print("\n\t 3- NEWTON ")
            print("\n\t 4- POINTS FIXES ")
            try:
                choix = int(input("\nFaire un choix entre {1} {2} {3} {4}: "))
                if choix in [1, 2, 3, 4]:
                    # appel aux procedures des methodes
                    if choix == 1:
                        print("\n\t================ D I C H O T O M I E ================")
                        dichotomie()
                    elif choix == 2:
                        print("\n\t================ L A  S E C A N T E  ================")
                        secante()
                    elif choix == 3:
                        print("\n\t================  N  E  W  T  O  N  ================")
                        newton()
                    elif choix == 4:
                        print("\n\t================ P O I N T  F I X E ================")
                        point_fixe()
                    else:
                        pass
                else:
                    raise PermissionError
        
                continuer = str.lower(input("\nVoulez vous continuer (o/n)? ")) # forcer la valeur du choix a etre en minuscule
                if continuer != 'n' and continuer != 'o':
                    print("\n\tErreur: Saisie invalide ")
                
            except PermissionError:
                print("\nChoix invalide !")
                
            effacer_console() # on rend au propre la console

    else:
        print("\nAucune solution pour cette équation !")

    print('\n')

# main()
# dichotomie()
# point_fixe()
# newton()
# secante()
"""print(tronquer(3.99999999999998, 5))
liste = ['', 1, 0, 1, '\n', 4, 5, 1, "\n", 8, 1, 2]
print(*liste)
mot_de_passe = getpass()
print(mot_de_passe)"""
