'''
David Bernier

Mai 2021

Ce fichier contient la classe Complexe et une fonction servant au calcul d'homologie du complexe.  On retrouve les
méthodes pour l'initialisation, pour le calcul des rang d'homologie et pour le calcule du champs gradient avec les
valeurs sur les sommets.

'''
from Cellule import Cellule
import numpy as np
from scipy.sparse import lil_matrix
from bisect import bisect_left
from Champ_Gradient import Champ_Gradient




class Complexe :
    """
    Classe modélisant un complexe simplicial ou cubique.

    ...

     Attributes
    ----------
    id_max : list
        la liste des nombres de cellules dans le complexe.

    IJ : tuple
        Dimention de l'image (utiliser pour l'affichage).

    Dim : int
        La dimension de la plus grande des cellules.

    Liste_Cellules : list
        Une liste (selon la dimention) de liste des indexe des (dim)-cellules.

    Methods
    -------

    Méthode pour l'initialisation :
    -----
    Init_Rectangle(I, J)
        Initialise un complexe en forme de grille uniforme 2D I par J.

    Init_Tore_Exemple()
        Initialise un tore comme la figure 11.2 (page 381) de Computational Homology.

    Add_Face_Tore(Id_Face, Liste_Id_Arete, Liste_orientation)
        Ajoute une 2-cellule dans le complexe.

    Add_Arete_Tore(Id_Arete, Id_moins, Id_plus)
        Ajoute une 1-cellule dans le complexe.

    Méthode sur le calcul de l'homologie :
    -----
    Get_Frontiere()
        Calcule les applications de fontiere du complexe en forme matriciel.

    Calcul_CCR()
        Procède à la réduction du complexe de chaîne.

    Méthode du calcul du champs gradient :
    -----
    Calcul_Champs_Gradient(Valeur_0)
         Calcul le champ gradient avec l'algo lowerstar de Robins peut importe la dim du complexe.

    Etoile(cube, dim)
        Calcule l'étoile  d'une cellule.

    Face_Non_Appariees(index, dim, etoile)
        Donne les faces non appariées de la index ème (dim)-cellule et dans etoile.

    Pair_Cellule(dim, index_0, index_1)
        Change les cellules pour indiquer leur pair.

    Pair_Cellule_seul(dim, index_0)
        Change la cellule pour indiquer qu'elle est critique.


    Autres méthode de la classe :
    -----
    Get_Cellules(dim,index)
        Retourne la index ème (dim)-cellule du complexe.

    Get_IJ()
        Retourne les valeurs des dimension de l'image soit I et J.

    Get_Nombre_cellule()
        Retourne une liste avec le nombre de chaque (dim)-cellules.

    Get_dim()
        Retourne la dimension de la plus grande des cellules.

    """

    def __init__(self):
        """
        Initialise la classe Complexe avec les valeurs nulles

        """
        self.id_max = []
        self.Dim = 0
        self.Liste_Cellules = []
        self.IJ = (0,0)

    '''
    Méthode pour l'initialisation :
    '''

    def Init_Rectangle(self, I, J) :
        """
        Initialise un complexe en forme de grille uniforme 2D I par J.

        Parameters
        ----------
        I : int
            Dimension de la grille horizontale.

        J : int
             Dimension de la grille verticale.

        """
        self.IJ = (I, J)
        self.Dim = 2

        # Dim 0
        Compteur = 0
        self.Liste_Cellules.append([])          # Listes des 0-celulles
        for i in range(I):
            for j in range(J):
                self.Liste_Cellules[0].append(Cellule(Compteur, 0, [], []))
                Compteur = Compteur + 1
        self.id_max.append(Compteur)


        # Dim 1
        self.Liste_Cellules.append([])          # Listes des 1-celulles
        Compteur = 0
        # Ligne du haut
        for j in range(J - 1):
            # Ajout de l'arete
            self.Liste_Cellules[1].append(Cellule(Compteur, 1, [j, j+1], [-1, 1]))
            self.Liste_Cellules[0][j].Add_Coface(Compteur) # Ajustement des cofaces
            self.Liste_Cellules[0][j+1].Add_Coface(Compteur)
            Compteur = Compteur + 1

        # Autres lignes
        for i in np.arange(1, I):
            for j in np.arange(0, J - 1):
                # Ajout horizontale
                self.Liste_Cellules[1].append(Cellule(Compteur, 1, [j + 1 + i * J, j + i * J], [1, -1]))
                self.Liste_Cellules[0][j + 1 + i * J].Add_Coface(Compteur)
                self.Liste_Cellules[0][j + i * J].Add_Coface(Compteur)
                Compteur = Compteur + 1

                # Ajout verticale
                self.Liste_Cellules[1].append(Cellule(Compteur, 1, [j + (i - 1) * J, j + i * J], [1, -1]))
                self.Liste_Cellules[0][j + (i - 1) * J].Add_Coface(Compteur)
                self.Liste_Cellules[0][j + i * J].Add_Coface(Compteur)
                Compteur = Compteur + 1


            # Ajout verticale last
            self.Liste_Cellules[1].append(Cellule(Compteur, 1, [i * J - 1, (i + 1) * J - 1], [1, -1]))
            self.Liste_Cellules[0][i * J - 1].Add_Coface(Compteur)
            self.Liste_Cellules[0][(i + 1) * J - 1].Add_Coface(Compteur)
            Compteur = Compteur + 1
        self.id_max.append(Compteur)


        # Dim 2
        self.Liste_Cellules.append([])                   # Listes des 2-celulles
        Compteur = 0

        for j in range(J - 2):
            # Ordre des faces : bas, droite, gauche, haut
            self.Liste_Cellules[2].append(Cellule(Compteur, 2, [(J - 1) + 2 * j, (J + 2) + 2 * j,
                                                                J + 2 * j, j],[1, 1, -1, -1]))
            self.Liste_Cellules[1][(J - 1) + 2 * j].Add_Coface(Compteur)
            self.Liste_Cellules[1][(J + 2) + 2 * j].Add_Coface(Compteur)
            self.Liste_Cellules[1][J + 2 * j].Add_Coface(Compteur)
            self.Liste_Cellules[1][j].Add_Coface(Compteur)
            Compteur = Compteur + 1

        self.Liste_Cellules[2].append(Cellule(Compteur, 2, [(J - 1) + 2 * (J - 2),(J + 2) + 2 * (J - 2) - 1,
                                                            J + 2 * (J - 2), J - 2], [1, 1, -1, -1]))
        self.Liste_Cellules[1][(J - 1) + 2 * (J - 2)].Add_Coface(Compteur)
        self.Liste_Cellules[1][(J + 2) + 2 * (J - 2) - 1].Add_Coface(Compteur)
        self.Liste_Cellules[1][J + 2 * (J - 2)].Add_Coface(Compteur)
        self.Liste_Cellules[1][(J - 2)].Add_Coface(Compteur)
        Compteur = Compteur + 1


        for i in np.arange(1, I - 1):
            for j in range(J - 2):
                self.Liste_Cellules[2].append(Cellule(Compteur, 2, [(J - 1) + (2 * J - 1) * i + 2 * j,
                                                                    (J + 2) + (2 * J - 1) * i + 2 * j,
                                                                    J + (2 * J - 1) * i + 2 * j,
                                                                    (J - 1) + (2 * J - 1) * (i - 1) + 2 * j],
                                                      [1, 1, -1, -1]))
                self.Liste_Cellules[1][(J - 1) + (2 * J - 1) * i + 2 * j].Add_Coface(Compteur)
                self.Liste_Cellules[1][(J + 2) + (2 * J - 1) * i + 2 * j ].Add_Coface(Compteur)
                self.Liste_Cellules[1][J + (2 * J - 1) * i + 2 * j].Add_Coface(Compteur)
                self.Liste_Cellules[1][(J - 1) + (2 * J - 1) * (i - 1) + 2 * j].Add_Coface(Compteur)
                Compteur = Compteur + 1

            self.Liste_Cellules[2].append(Cellule(Compteur, 2, [(J - 1) + (2 * J - 1) * i + 2 * (J - 2),
                                                                (J + 2) + (2 * J - 1) * i + 2 * (J - 2) - 1,
                                                                J + (2 * J - 1) * i + 2 * (J - 2),
                                                                (J - 1) + (2 * J - 1) * (i - 1) + 2 * (J - 2)],
                                                  [1, 1, -1, -1]))
            self.Liste_Cellules[1][(J - 1) + (2 * J - 1) * i + 2 * (J - 2)].Add_Coface(Compteur)
            self.Liste_Cellules[1][(J + 2) + (2 * J - 1) * i + 2 * (J - 2) - 1].Add_Coface(Compteur)
            self.Liste_Cellules[1][J + (2 * J - 1) * i + 2 * (J - 2)].Add_Coface(Compteur)
            self.Liste_Cellules[1][(J - 1) + (2 * J - 1) * (i - 1) + 2 * (J - 2)].Add_Coface(Compteur)
            Compteur = Compteur + 1

        self.id_max.append(Compteur)

    def Init_Tore_Exemple(self):
        """
        Initialise un tore comme la figure 11.2 (page 381) de Computational Homology.

        Orientation issues de mes notes de cours page 22.

        """

        self.Dim = 2

        # Dim 0
        Compteur = 0
        self.Liste_Cellules.append([])  # Listes des 0-celulles
        for i in range(9):
                self.Liste_Cellules[0].append(Cellule(Compteur, 0, [], []))
                Compteur = Compteur + 1
        self.id_max.append(Compteur)


        # Dim 1
        self.Liste_Cellules.append([])  # Listes des 1-celulles
        self.Add_Arete_Tore(0,0,1)
        self.Add_Arete_Tore(1,1,2)
        self.Add_Arete_Tore(2,2,0)
        self.Add_Arete_Tore(3,3,0)
        self.Add_Arete_Tore(4,6,3)
        self.Add_Arete_Tore(5,0,6)
        self.Add_Arete_Tore(6,3,4)
        self.Add_Arete_Tore(7,6,7)
        self.Add_Arete_Tore(8,4,5)
        self.Add_Arete_Tore(9,7,8)
        self.Add_Arete_Tore(10,5,3)
        self.Add_Arete_Tore(11,8,6)
        self.Add_Arete_Tore(12,4,1)
        self.Add_Arete_Tore(13,5,2)
        self.Add_Arete_Tore(14,7,4)
        self.Add_Arete_Tore(15,8,5)
        self.Add_Arete_Tore(16,1,7)
        self.Add_Arete_Tore(17,2,8)
        self.Add_Arete_Tore(18,3,1)
        self.Add_Arete_Tore(19,6,4)
        self.Add_Arete_Tore(20,0,7)
        self.Add_Arete_Tore(21,1,8)
        self.Add_Arete_Tore(22,7,5)
        self.Add_Arete_Tore(23,4,2)
        self.Add_Arete_Tore(24,5,0)
        self.Add_Arete_Tore(25,8,3)
        self.Add_Arete_Tore(26,2,6)
        self.id_max.append(27)

        # Dim 2
        self.Liste_Cellules.append([])  # Listes des 2-celulles
        self.Add_Face_Tore(0, [18,3,0], [-1,1,1])
        self.Add_Face_Tore(1, [6,18,12], [-1,1,-1])
        self.Add_Face_Tore(2, [12,1,23], [1,1,-1])
        self.Add_Face_Tore(3, [8,23,13], [-1,1,-1])
        self.Add_Face_Tore(4, [13,2,24], [1,1,-1])
        self.Add_Face_Tore(5, [10,24,3], [-1,1,-1])
        self.Add_Face_Tore(6, [19,4,6], [-1,1,1])
        self.Add_Face_Tore(7, [7,19,14], [-1,1,-1])
        self.Add_Face_Tore(8, [22,14,8], [-1,1,1])
        self.Add_Face_Tore(9, [22,15,9], [1,-1,-1])
        self.Add_Face_Tore(10, [15,10,25], [1,1,-1])
        self.Add_Face_Tore(11, [11,25,4], [-1,1,-1])
        self.Add_Face_Tore(12, [5,7,20], [1,1,-1])
        self.Add_Face_Tore(13, [0,20,16], [-1,1,-1])
        self.Add_Face_Tore(14, [16,9,21], [1,1,-1])
        self.Add_Face_Tore(15, [1,21,17], [-1,1,-1])
        self.Add_Face_Tore(16, [26,17,11], [-1,1,1])
        self.Add_Face_Tore(17, [2,26,5], [-1,1,-1])
        self.id_max.append(18)

    def Add_Face_Tore(self, Id_Face, Liste_Id_Arete, Liste_orientation):
        """
        Ajoute une 2-cellule dans le complexe.

        Parameters
        ----------
        Id_Face : int
            Index de la face que l'on ajoute.

        Liste_Id_Arete : list
            Liste des index (int) des 1-cellules composant la face.

        Liste_orientation : list
             Liste des orientations (int) respectives .

        """

        self.Liste_Cellules[2].append(Cellule(Id_Face, 2, Liste_Id_Arete, Liste_orientation))
        for i in range(len(Liste_Id_Arete)):
            self.Liste_Cellules[1][Liste_Id_Arete[i]].Add_Coface(Id_Face)

    def Add_Arete_Tore(self, Id_Arete, Id_moins, Id_plus):
        """
        Ajoute une 1-cellule dans le complexe.

        Parameters
        ----------
        Id_Arete : int
            Index de l'arete que l'on ajoute.

        Id_moins : int
            Index de la 0-cellule initiale (orientation).

        Id_plus : int
            Index de la 0-cellule finale (orientation).

        """

        self.Liste_Cellules[1].append(Cellule(Id_Arete, 1, [Id_moins, Id_plus], [-1, 1]))
        self.Liste_Cellules[0][Id_moins].Add_Coface(Id_Arete)
        self.Liste_Cellules[0][Id_plus].Add_Coface(Id_Arete)


    '''
    Méthode sur le calcul de l'homologie :
    '''

    def Get_Frontiere(self):
        """
        Calcule les applications de fontiere du complexe en forme matriciel.

        Returns
        -------
        Frontiere : list
            liste des matrices (lil_matrix) représentant les applications de frontère.

        """

        Frontiere = []

        if len(self.Liste_Cellules) != 0:  # Dim 0
            Bd_0 = lil_matrix(np.zeros((len(self.Liste_Cellules[0]), 1)))
            Frontiere.append(Bd_0)

        for h in range(self.Dim):  # Pour chaque Dim > 0
            temp = lil_matrix(np.zeros((len(self.Liste_Cellules[h + 1]), len(self.Liste_Cellules[h]))), dtype=np.int8)
            for i in range(len(self.Liste_Cellules[h + 1])):  # Pour chaque cellule
                faces = self.Liste_Cellules[h + 1][i].Get_Face()
                Orientation = self.Liste_Cellules[h + 1][i].Get_Orientation()
                for j in range(len(faces)):  # Pour chaque face d'une cellules
                    temp[i, faces[j]] = Orientation[j]
            Frontiere.append(temp)
        return Frontiere

    def Calcul_CCR(self):
        """
        Procède à la réduction du complexe de chaîne.

        Affiche le nombre de générateur restant.

        Returns
        -------
        BD : list
            liste des matrices (lil_matrix) représentant les applications de frontère.

        IND : list
            Liste de l'index (int) des faces restantes.

        """

        IND = [np.arange(len(self.Liste_Cellules[i])) for i in range(len(self.Liste_Cellules))]  # Liste de id
        BD =self.Get_Frontiere()
        for i in range(len(IND) - 1, 0, -1):
            found = False
            while found == False:
                found = True
                Stop = False
                for b in range(IND[i].shape[0]):
                    for a in range(IND[i - 1].shape[0]):
                        if Stop == False:
                            if abs(BD[i][b, a]) == 1:
                                IND, BD = Reduc_Pair(IND, BD, i, a, b)
                                Stop = True
                                found = False
        print('nb gen 0 = ', len(IND[0]), ', nb gen 1 = ', len(IND[1]), ' et nb gen 2 = ', len(IND[2]))
        return BD, IND


    '''
    Méthode du calcul du champs gradient :
    '''
    def Calcul_Champs_Gradient(self,Valeur_0):
        """
        Calcul le champ gradient avec l'algo lowerstar de Robins peut importe la dim du complexe.

        On trouve Pair_Vecteur, Cellule_critique et Valeur critique tel que decrit dans la classe Champs_Gradient

        Parameters
        ----------
        Valeur_0 : list de float
            La liste des valeurs sur les sommets initiales ligne par ligne.

        Returns
        -------

        Champ_gradient : class Champ_Gradient
            Objet de type Champ_ Gradient contenant les information necesaires au complexe de Morse.

        """

        # Initialiser les listes
        Cellule_critique = []   # liste des cubes critiques
        Valeur_critique = []    # liste des valeurs critiques
        Pair_Vecteur = []       # liste des fleches 0-1 et 1-2
        for d in range(self.Get_dim()) :
            Cellule_critique.append([])
            Valeur_critique.append([])
        for d in range(self.Get_dim()-1):
            Pair_Vecteur.append([])

        # File d'attente pour la filtration lowerstar
        PQ_0 = []
        PQ_1 = []

        # Liste des listes [selon la dim] des indexs des cellules de l'etoile inferieure du i eme sommet
        Ls_inf = []
        for i in range(len(Valeur_0)):
            Ls_inf.append([])

        # Ordre croissant des valeurs sur les somments
        sorted_index_0 = np.argsort(Valeur_0)


        # Calcul de g
        # Fontion_g est liste [selon la dim] de vecteurs contenant toutes les valeurs d'une dim donnee
        Fontion_g = []
        Fontion_g.append(Valeur_0)

        # Initialisation des vecteurs
        for i in np.arange(1, self.Get_dim()):         # Pour chaque dimention > 0
            temp = np.empty(self.Get_Nombre_cellule()[i])
            temp[:] = -10000000                                # Valeur Initial
            Fontion_g.append(temp)

        # En ordre decroissant sur les valeur_0
        for i in np.arange(len(Valeur_0) - 1, -1, -1):

            # On considere le 0-cube en position x
            x = int(sorted_index_0[i])

            # On Calcul l'etoile de la cellule
            Ls = self.Etoile(x, 0)

            # Intialisation de l'etoile inferieure de la cellule
            for j in range(len(Ls)):
                Ls_inf[x].append([])
            Ls_inf[x][0].append(x)

            # On associe la valeur a chaque cellule de l'etoie du sorted_index_0[i] eme sommet si la cellule en n'a pas
            for d in np.arange(1, len(Ls)):                  # Pour chaque dim,
                for j in range(len(Ls[d])):                   # Pour chaque d-cellules de l'etoile
                    if Fontion_g[d][Ls[d][j]] == -10000000:    # Si pas de valeur associee
                        Fontion_g[d][Ls[d][j]] = Valeur_0[x]    # On associe la valeur du sommet no x
                        Ls_inf[x][d].append(Ls[d][j])           # Ajout dans l'etoile inferieure

        # FIN Calcul de g

        # En ondre croisant sur les valeurs de somment
        for y in range(len(sorted_index_0)):
            x = int(sorted_index_0[y])
            # On considere le 0-cube en position x

            # Si x est un minimum donc si etoile inferieure vide
            if len(Ls_inf[x][1]) == 0:  # Si x est minimum local,
                Cellule_critique[0].append(x)  # On ajoute x au cellule critique
                Valeur_critique[0].append(Valeur_0[x])
                self.Pair_Cellule_seul(0,x)

            else:  # Sinon
                # Ajouter une fleche
                liste_g = [] # valeur des sommet autres que x
                for j in range(len(Ls_inf[x][1])):
                    for c in self.Get_Cellules(1,Ls_inf[x][1][j]).Get_Face():
                        if c != x :
                            liste_g.append(Fontion_g[0][c])
                index_delta = np.argmin(np.array(liste_g))    # On trouve la plus petite dans les 1-cellules

                delta = Ls_inf[x][1][index_delta]
                Pair_Vecteur[0].append([x, delta])  # fleche de x vers delta
                self.Pair_Cellule(0, x, delta)

                # Ajouter les autres 1-cube a PQ_0
                for j in range(len(Ls_inf[x][1])):
                    if j != index_delta:
                        PQ_0.append([1,Ls_inf[x][1][j]])

                # Ajouter les cubes  a PQ_1
                for d in np.arange(2,len(Ls_inf[x])):
                    for j in range(len(Ls_inf[x][d])):
                        list_face_non_app = self.FaceNonAppariees(Ls_inf[x][d][j], d, Ls_inf[x][d-1])
                        if len(list_face_non_app) == 1:
                            PQ_1.append([d, Ls_inf[x][d][j]])


                # On vide les files d'attente
                while (len(PQ_1) != 0 or len(PQ_0) != 0):
                    while len(PQ_1) != 0:

                        alpha = PQ_1.pop(0)
                        pair_alpha = self.FaceNonAppariees(alpha[1], alpha[0], Ls_inf[x][alpha[0]-1])

                        if len(pair_alpha) == 0:
                            PQ_0.append(alpha)
                        else:
                            # On ajoute la fleche
                            Pair_Vecteur[alpha[0]-1].append([pair_alpha[0], alpha[1]])
                            self.Pair_Cellule(alpha[0]-1, pair_alpha[0], alpha[1])
                            PQ_0.remove([alpha[0]-1, pair_alpha[0]])

                            # ajouter a PQ_1 les new 1 faces non appariees
                            for d in np.arange(2,len(Ls_inf[x])):
                                for j in range(len(Ls_inf[x][d])):

                                    # A modifier
                                    list_face_non_app = self.FaceNonAppariees(Ls_inf[x][d][j], d, Ls_inf[x][d - 1])
                                    if len(list_face_non_app) == 1:
                                        # On verifie que ce n'est pas dans la liste
                                        Test = False
                                        for w in range(len(PQ_1)):
                                            if PQ_1[w][1] == Ls_inf[x][d][j]:
                                                Test = True
                                        if Test == False :
                                            PQ_1.append([d, Ls_inf[x][d][j]])

                    if len(PQ_0) != 0:
                        gamma = PQ_0.pop(0)
                        Cellule_critique[gamma[0]].append(gamma[1])
                        Valeur_critique[gamma[0]].append(Fontion_g[gamma[0]][gamma[1]])
                        self.Pair_Cellule_seul(gamma[0], gamma[1])


                        # ajouter a PQ_1 les new 1 faces non appariees
                        for d in np.arange(2, len(Ls_inf[x])):
                            for j in range(len(Ls_inf[x][d])):
                                list_face_non_app = self.FaceNonAppariees(Ls_inf[x][d][j], d, Ls_inf[x][d - 1])
                                if len(list_face_non_app) == 1:
                                    # On verifie que ce n'est pas dans la liste
                                    Test = False
                                    for w in range(len(PQ_1)):
                                        if PQ_1[w][1] == Ls_inf[x][d][j]:
                                            Test = True
                                    if Test == False:
                                        PQ_1.append([d, Ls_inf[x][d][j]])

        Champ_gradient = Champ_Gradient(self, Pair_Vecteur, Cellule_critique, Valeur_critique)

        return Champ_gradient

    def Etoile(self, cube, dim):
        """
        Calcule l'étoile  d'une cellule.

        Parameters
        ----------
        cube : int
            l'index de la cellule ou l'on veut calculer l'étoile.

        dim : int
            la dimension de la cellule.

        Returns
        -------
        LS : list
            Liste (selon la dimension) de liste des index de chaques (dim)-cellule de l'étoile.

        """

        LS = []
        for i in range(dim):
            LS.append([])
        LS.append([cube])
        for i in np.arange(dim + 1, len(self.Liste_Cellules)):  # inclus dim+1
            temp = []
            for j in range(len(self.Liste_Cellules[i])):  # j eme i-cube
                stop = False
                face = self.Liste_Cellules[i][j].Get_Face()
                for x in face:  # x est l'indexe des (i-1)-cube qui compose le jeme cube de dimention i
                    for y in LS[i - 1]:
                        if x == y:
                            test = True
                            for h in temp:
                                if j == h:
                                    test = False
                            if test:
                                temp.append(j)
            LS.append(temp)

        return LS

    def FaceNonAppariees(self, index, dim, etoile):
        """
        Donne les faces non appariées de la index ème (dim)-cellule.

        Parameters
        ----------
        index : int
            l'index du cube ou l'on veux calculer le nombre de face non appariée.

        dim : int (>=1)
            La dimension de cube que l'on veut calculer le nombre de face non appariée.

        etoile : list
            La liste des indexs dans laquelle on limite notre recherche.

        Returns
        -------
        Liste_Face_Non_Appariees_Dans_Etoile : list
            La liste des faces non appariées de la cellule dans la restriction etoile.

        """

        if dim == 0:
            return 0, []

        Liste_Face = self.Get_Cellules(dim, index).Get_Face()
        Liste_Face_Non_Appariees = []
        for face in Liste_Face:
            if not self.Get_Cellules(dim - 1, face).Est_Appariee() :
                Liste_Face_Non_Appariees.append(face)


        Liste_Face_Non_Appariees_Dans_Etoile = []
        for face in Liste_Face_Non_Appariees:
            for et in etoile :
                if face == et:
                    Liste_Face_Non_Appariees_Dans_Etoile.append(face)


        return Liste_Face_Non_Appariees_Dans_Etoile

    def Pair_Cellule(self, dim, index_0, index_1):
        """
        Change les cellules pour indiquer leur pair.

        La index_0 eme cellule de dimension dim ce voit associer avec la index_0 eme cellule de dimemsion dim + 1

        Parameters
        ----------
        dim : int
            La dimention de la plus petite cellule (base de la fleche).

        index_0 : int
            L'index de la cellule de plus petite dimension.

        index_1 : int
            L'index de la cellule de dimension supérieure.

        """

        self.Liste_Cellules[dim][index_0].Set_Pair([dim + 1, index_1])
        self.Liste_Cellules[dim + 1][index_1].Set_Pair([dim, index_0])

    def Pair_Cellule_seul(self, dim, index_0):
        """
        Change la cellule pour indiquer qu'elle est critique.

        la index_0 eme cellule de dimention dim ce voit associer la valeur [index_0]

        Parameters
        ----------
        dim : int
            Dimention de la cellule.

        index_0 : int
            Index de la cellule

        """

        self.Liste_Cellules[dim][index_0].Set_Pair([index_0])


    '''
    Autres méthode de la classe :
    '''
    def Get_Cellules(self,dim,index):
        val = -1
        if dim >= len(self.Liste_Cellules):
            print('La dim est trop grande')

        elif index >= len(self.Liste_Cellules[dim]):
            print('L index est trop grand')

        else:
            val = self.Liste_Cellules[dim][index]

        return val

    def Get_IJ(self):
        return self.IJ

    def Get_Nombre_cellule(self):
        a = [len(self.Liste_Cellules[i]) for i in range(len(self.Liste_Cellules))]
        return a

    def Get_dim(self):
        return len(self.Liste_Cellules)




def Reduc_Pair(E,Bd,i,a,b):
    """
    Procède à la réduction d'une pair dans un complexe de chaîne.

    UtilisÉ dans Calcul_CCR).

    Parameters
    ----------
    E : list
        Liste (selon la dim) de liste des indexs (int) des faces restantes.

    Bd : list
        liste des matrices (lil_matrix) représentant les applications de frontère restantes.

    i : int
        Dimention de la pair.

    a : int
        Première élémment de la pair.

    b : int
        Second élément de la pair.

    Returns
    -------
    BD : list
            liste des matrices (lil_matrix) représentant les applications de frontère.

    IND : list
        Liste de l'index (int) des faces restantes.

    """

    for e in range(np.shape(Bd[i])[0]) :
        if (e != b and Bd[i][e,a] != 0) :
            for j in range(np.shape(Bd[i])[1]) :
                if (Bd[i][b,j] != 0 and j!=a) :
                    Bd[i][e,j] = Bd[i][e,j] - Bd[i][e,a]*Bd[i][b,a]*Bd[i][b,j]

    # Retrait des id
    E[i] = np.delete(E[i],b)
    E[i - 1] = np.delete(E[i - 1], a)

    # Retrait de la ligne b
    Bd[i].rows = np.delete(Bd[i].rows, b)
    Bd[i].data = np.delete(Bd[i].data, b)
    Bd[i]._shape = (Bd[i]._shape[0] - 1, Bd[i]._shape[1])

    # Retrait de la colonne a
    rows = Bd[i].rows
    data = Bd[i].data
    for j in np.arange(Bd[i].shape[0]):
        pos = bisect_left(rows[j], a)
        if pos == len(rows[j]):
            continue
        elif rows[j][pos] == a:
            rows[j].pop(pos)
            data[j].pop(pos)
            if pos == len(rows[j]):
                continue
        for pos2 in np.arange(pos, len(rows[j])):
            rows[j][pos2] -= 1
    Bd[i].rows =  rows
    Bd[i].data = data
    Bd[i]._shape = (Bd[i]._shape[0], Bd[i]._shape[1] - 1)

    # Retrait de la ligne a
    Bd[i - 1].rows = np.delete(Bd[i - 1].rows, a)
    Bd[i - 1].data = np.delete(Bd[i - 1].data, a)
    Bd[i-1]._shape = (Bd[i-1]._shape[0] - 1, Bd[i-1]._shape[1])
    return(E,Bd)
