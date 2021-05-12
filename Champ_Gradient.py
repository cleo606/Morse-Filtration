'''
David Bernier

Mai 2021

Ce fichier contient la classe champs grandient. Cette classe est utiliser lors du calcul de champs gradient dans la
classe Complexe.

'''
from MS_Complexe import MS_Complexe
import numpy as np
import matplotlib.pyplot as plt




class Champ_Gradient :
    """
    Stucture pour du complexe de Morse-Smale.

    ...

    Attributes
    ----------
    complexe : class Complexe
        Le complexe sur lequel on calcule le champs grandiant.

    Pair_Vecteur : list
        La liste (selon la dimension de la plus petite cellule) qui contient les vecteurs du champs gradient.

    Cellule_critique : list
        la liste (selon la dimention) qui condient l'index des cellules critiques.

    Valeur_critique : Liste
        La liste des listes de valeur des d-cellules critiques.

    Methods
    -------
    Affichage()
        Affiche le champs gradient avec la librairie pyplot.

    Calcul_Complexe_De_Morse()
        Calcul du complexe de Morse-Smale en trouvant les arcs de Morse.

    Calcul_Variete(dim, index, ascendant, M)
        Calcul les v-chemins ascendants ou descendants de la index ème (dim)-cellule.

    Get_Cell_Crit()
        Retoune Cellule_critique.

    Get_Pair_Cell()
        Retourne Pair_Vecteur.

    Get_IJ()
        Retourne les valeurs des dimension de l'image soit I et J.

    Get_Cellules(dim,index)
        Retourne la index ème (dim)-cellule du complexe.

    Get_Valeur_Critique()
        Retourne Valeur_critique.

    """

    def __init__(self, complexe, Pair_Vecteur, Cellule_critique, Valeur_critique): #
        """
        Initialise un objet de type champ gradient avec un complexe.

        Parameters
        ----------
        complexe : class Complexe
            Le complexe sur lequel on calcule le champs grandiant.

        Pair_Vecteur : list
            La liste (selon la dimension de la plus petite cellule) qui contient les vecteurs du champs gradient.

        Cellule_critique : list
            la liste (selon la dimention) qui condient l'index des cellules critiques.

        Valeur_critique : Liste
            La liste des listes de valeur des d-cellules critiques.

        """
        self.complexe = complexe
        self.Pair_Vecteur = Pair_Vecteur
        self.Cellule_critique = Cellule_critique
        self.Valeur_critique = Valeur_critique

    def Affichage(self, image=None):
        """
        Affiche le champs gradient avec la librairie pyplot.

        Parameters
        ----------
        image : Bool
            Est-ce que l'on affiche li'mage en background ?

        """

        ax = plt.axes()
        I,J = self.complexe.Get_IJ()

        # Pour les 0-cellules critiques

        list_0 = [[],[]]
        for i in self.Cellule_critique[0]:
            list_0[0].append(i % I)           # modulo
            list_0[1].append(i // I)          # divison entiere

        # Ajoutes les points
        ax.scatter(list_0[0],list_0[1], s=150, label='Point critique', color = 'm')

        # Pour les 1-cellules critiques

        for i in self.Cellule_critique[1]:
            sommet = self.complexe.Get_Cellules(1,i).Get_Face()
            p_0 = [sommet[0] % I, sommet[0] // I]
            p_1 = [sommet[1] % I, sommet[1] // I]
            # Affichage de l'arette P_0 -> p_1
            ax.arrow(p_0[0], p_0[1], p_1[0]-p_0[0], p_1[1]-p_0[1], head_width=0.01, head_length=0.01, color = 'g',
                     width = 0.1)

        # Pour les 2-cellules critiques

        list_2 = [[],[]]
        for i in self.Cellule_critique[2]:
            list_2[0].append(i % (I-1) + 0.5)  # modulo
            list_2[1].append((i // (I-1)) + 0.5)

        # Ajoutes des face en fome de caree bleu
        ax.scatter(list_2[0], list_2[1], s=60, label='Face critique', marker='s', color='b')


        # Pour les 0-1 pairs

        for i in self.Pair_Vecteur[0] :
            # index de la base du vecteur
            id_0 = i[0]
            for j in self.complexe.Get_Cellules(1,i[1]).Get_Face():
                if j != id_0 :
                    # index au bout du vecteur
                    id_1 = j


            # Calculc des points
            p_0 = [id_0 % I, id_0 // I]
            p_1 = [id_1 % I, id_1 // I]

            # Ajout du vecteur p_0 -> p_1
            if  p_1[0] - p_0[0] == 0 :
                # Vertical
                if p_1[1] - p_0[1] > 0:
                    ax.arrow(p_0[0], p_0[1], 0, 0.7, head_width=0.3, head_length=0.4, color = 'r', alpha=0.5)

                else :
                    ax.arrow(p_0[0], p_0[1], 0, -0.7, head_width=0.3, head_length=0.4, color='r', alpha=0.5)


            else:
                # Horizontal
                if p_1[0] - p_0[0] > 0:
                    ax.arrow(p_0[0], p_0[1], 0.7, 0, head_width=0.4, head_length=0.4, color = 'r', alpha=0.5)

                else:
                    ax.arrow(p_0[0], p_0[1], -0.7, 0, head_width=0.4, head_length=0.4, color='r', alpha=0.5)


        # Pour les  1-2 pair

        for i in self.Pair_Vecteur[1]:
            # Centre du cube
            p_1 = [i[1] % (I-1) + 0.5, (i[1] // (I-1)) + 0.5]
            # Centre de l'arette
            id = self.complexe.Get_Cellules(1,i[0]).Get_Face()
            # Calcule de la moyenne des deux sommet de l'arc
            p_0 = [(id[0] % I + id[1] % I)/2, ((id[0] // I) + (id[1] // I))/2]

            # Ajout du vecteur p_0 -> p_1
            if  p_1[0] - p_0[0] == 0 :
                # Vertical
                if p_1[1] - p_0[1] > 0:
                    ax.arrow(p_0[0], p_0[1], 0, 0.42, head_width=0.3, head_length=0.4, color = 'y', alpha=0.5)

                else :
                    ax.arrow(p_0[0], p_0[1], 0, -0.42, head_width=0.3, head_length=0.4, color='y', alpha=0.5)


            else:
                # Horizontal
                if p_1[0] - p_0[0] > 0:
                    ax.arrow(p_0[0], p_0[1], 0.42, 0, head_width=0.4, head_length=0.4, color = 'y', alpha=0.5)

                else:
                    ax.arrow(p_0[0], p_0[1], -0.42, 0, head_width=0.4, head_length=0.4, color='y', alpha=0.5)




        # Ajout de la grille
        ax.set_xticks(np.arange(0, I))
        ax.set_yticks(np.arange(0, J))
        ax.grid(True)

        # Ajout de l'image
        if type(image) != type(None) :
            ax.imshow(image, cmap='gray')
        plt.title('Champ gradient')
        plt.show()

    def Calcul_Complexe_De_Morse(self):
        """
        Calcul du complexe de Morse-Smale en trouvant les arcs de Morse.

        On doit avoir calculer le champs gradient.

        Returns
        -------
        MS : class MS_Complexe
            le complexe de Morse-Smale résultant

        """

        # Initialise la variable MS_Complexe
        MS = MS_Complexe(self)

        for d in np.arange(len(self.Cellule_critique) - 1, 0, -1):  # Pour chaque dim en ordre decroissant de dim à 1

            # Note : On procede pour le couple (d, d-1)

            # Initialisation des listes
            Liste_ascendant = []
            Liste_descendant = []
            Liste_Coface = [0]*self.complexe.Get_Nombre_cellule()[d-1]

            # Calcul des d-varietes descendantes
            for i in range(len(self.Cellule_critique[d])):  # Pour chaque celule critique de dim d
                # Ajoute des listes de chemins descendants dans Liste_descendant
                Liste_descendant.append(self.Calcul_Variete(d, self.Cellule_critique[d][i], False, [])[0])


            # Calcul des (d-1)-varietes ascendantes
            for j in range(len(self.Cellule_critique[d-1])):  # Pour chaque celule critique de dim d-1

                # Ajoute des listes de chemins ascendants dans Liste_ascendant
                Liste_ascendant.append(self.Calcul_Variete(d-1, self.Cellule_critique[d-1][j], True, [])[0])

                for k in range(len(Liste_ascendant[j])):  # Pour chaque (d-1)-chemin ascendant
                    # Comparaison avec les chemins des varietes descendantes
                    for i in range(len(self.Cellule_critique[d])):  # Pour chaque celule critique de dim d
                        for l in range(len(Liste_descendant[i])):  # Pour chaque d-chemin descendant
                            # Ici : Liste_descendant[j][l] = (d)-chemin  et Liste_ascendant[j][k] = (d-1)-chemin

                            # On test chaque cellule du (d-1)-chemin
                            z = 0
                            test = True
                            stop = False
                            for a in range(len(Liste_descendant[i][l])):
                                face = self.Get_Cellules(d, Liste_descendant[i][l][a]).Get_Face()
                                for id in face:
                                    if id == self.Cellule_critique[d - 1][j]:
                                        # Le chemin ascendant passe sur la d-cellule critique
                                        stop = True
                                        index = a
                            while(test and z < len( Liste_ascendant[j][k])and not stop):
                                # Si nouvelle cellule
                                if Liste_Coface[Liste_ascendant[j][k][z]] == 0:
                                    # Calcul des ces cofaces
                                    Liste_Coface[Liste_ascendant[j][k][z]]= self.complexe.\
                                        Get_Cellules(d-1, Liste_ascendant[j][k][z]).Get_Coface()


                                ''' Ancien
                                for a in range(len(Liste_Coface[Liste_ascendant[j][k][z]])):
                                    if Liste_Coface[Liste_ascendant[j][k][z]][a] == self.Cellule_critique[d][i]:
                                        # Le chemin ascendant passe sur la d-cellule critique
                                        stop = True
                                        index = a'''




                                # Recherche de l'element Liste_ascendant[j][k][z] dans  Liste_descendant[j][l]
                                # On regarde si un element de Liste_Coface[Liste_ascendant[j][k][z]] est dans
                                #   Liste_descendant[i][l][1:]
                                test = len(list(set(Liste_Coface[Liste_ascendant[j][k][z]]) -
                                                set(Liste_descendant[i][l][1:]))) !=\
                                       len(Liste_Coface[Liste_ascendant[j][k][z]])

                                z = z + 1

                            '''if test:
                                # Ajout d'un arc de Morse
                                arc = Liste_descendant[i][l]
                                MS.Add_Arc(d - 1, j, i, arc)'''

                            if stop:
                                # Ajout d'un arc de Morse
                                arc = Liste_descendant[i][l][:index+1]
                                MS.Add_Arc(d-1, j, i, arc)
                            elif test:
                                # Ajout d'un arc de Morse
                                arc = Liste_descendant[i][l]
                                MS.Add_Arc(d-1, j, i, arc)




        return MS

    def Calcul_Variete(self, dim, index, ascendant, M):
        """
        Calcul les v-chemins ascendants ou descendants de la index ème (dim)-cellule.

        Implementation inspirée de GetManifold de Sousbie P.17

        Parameters
        ----------
        dim : int
            La dimention de la cellule.

        index : int
            L'indice de la cellule à traiter.
        ascendant : bool

            True pour ascendant et False pour descendant.

        M : list
            La liste des (dim)-cellules déjà visitées.

        Returns
        -------
        arc_curr : list
            La liste des v-chemins. Un v-chemin est une liste de d'index des (dim)-cellules composant le v-chemin.

        M : list
            La liste des (dim)-cellules déjà visitées.

        """

        if ascendant :
            A_curr = self.complexe.Get_Cellules(dim, index).Get_Coface()
            CoDim = dim + 1
        else :
            A_curr = self.complexe.Get_Cellules(dim, index).Get_Face()
            CoDim = dim - 1
        A_temp = [] # liste des dim-cellules en attentes
        arc_curr = [] # Liste de chemin, un chemin est une liste de dim-cellule

        Compteur = 0
        for c in range(len(A_curr)) :
            # On regarde l'admisibilite
            p = self.complexe.Get_Cellules(CoDim, A_curr[c]).Get_Pair()
            # Si la pair de c est de dimention dim
            if dim == p[0]:
                # on regarde si il est deja la
                Pas_la = True
                for i in M :
                    if i==A_curr[c] :
                        Pas_la = False

                # Si on n'a pas deja couvert et c n'est pas uen cellule critique
                if Pas_la and not self.complexe.Get_Cellules(CoDim, A_curr[c]).Est_Critique():

                    # on note la cellule c comme couverte
                    M.append(A_curr[c])
                    # on ajoute un chemin
                    arc_curr.append([])
                    # On iterera sur la cellule d'indice p recursivement
                    A_temp.append([Compteur,p[1]])
                    Compteur = Compteur + 1
                # si c est critique on stop

        M_base = []
        for i in M:
            M_base.append(i)

        for i in range(len(A_temp)):
            # Calcul recursif du chemin
            M = M_base
            # M=M_base
            new_arcs, M2 = self.Calcul_Variete(dim, A_temp[i][1], ascendant, M)

            # Si le chemin continue
            if len(new_arcs) > 0 :
                # Mettre a jours les arcs

                # Cas len(new_arcs) > 1
                for j in np.arange(1,len(new_arcs)):      # pour chaque nouvelles arcs sauf la premiere,
                    # On ajoute une nouvelle arc a la fin
                    arc_curr.append([])
                    # Copie des cellules existantes du chemin
                    for k in arc_curr[A_temp[i][0]]:
                        arc_curr[len(arc_curr)-1].append(k)
                    # Ajout des nouvelles cellules
                    for k in new_arcs[j] :
                        arc_curr[len(arc_curr)-1].append(k)

                # Cas j = 0
                # Ajout des nouvelles cellules
                for k in new_arcs[0]:
                    arc_curr[A_temp[i][0]].append(k)

        # On ajoute la cellule d'entree au debut de chaque arcs
        for i in range(len(arc_curr)):
            arc_curr[i] = [index] + arc_curr[i]

        if len(arc_curr) == 0 :
            arc_curr = [[index]]

        return arc_curr, M

    def Epsilon_Simplification(self,seuil = 0):
        """
        Procède à la simplification des cellules critiques qui sont collées ensemble.
        Le parametre de seuil sert à mettre une limite selon les valeurs des cellules.

        Parameters
        ----------
        seuil : float
            Le seuil selon laquel on garde la apir critique. Le valeur critique varie entre 0 et 1 dans le cas d'image
            en teinte de gris.--

        """

        for d in np.arange(1,len(self.Valeur_critique)):
            index_to_pop = []
            for i in range(len(self.Valeur_critique[d])):
                face = self.complexe.Get_Cellules(d,self.Cellule_critique[d][i]).Get_Face()
                for j in range(len(self.Valeur_critique[d-1])):
                    for cell in face :
                        if cell == self.Cellule_critique[d-1][j]:
                            if seuil <= self.Valeur_critique[d][i] - self.Valeur_critique[d-1][j] :
                                # epsilon simplicifation entre la i eme d-cell et la j eme (d-1)-cell

                                # ajout de la fleche dans le complexe
                                self.complexe.Pair_Cellule(d-1, self.Cellule_critique[d-1][j], self.Cellule_critique[d][i])

                                # ajout de la fleche dans les paires gradientes
                                self.Pair_Vecteur[d-1].append([self.Cellule_critique[d-1][j], self.Cellule_critique[d][i]])

                                # retrait des cellules dans les liste critique
                                self.Cellule_critique[d-1].pop(j)
                                self.Valeur_critique[d - 1].pop(j)
                                index_to_pop.append(i)
                                break


                    else:
                        continue
                    break

                continue

            # Retrait des d-cellules
            index_to_pop.sort(reverse=True)
            for i in index_to_pop:
                self.Valeur_critique[d].pop(i)
                self.Cellule_critique[d].pop(i)

    def Get_Cell_Crit(self):
        return self.Cellule_critique

    def Get_Pair_Cell(self):
        return self.Pair_Vecteur

    def Get_IJ(self):
        return self.complexe.Get_IJ()

    def Get_Cellules(self,dim,index):
        return self.complexe.Get_Cellules(dim,index)

    def Get_Valeur_Critique(self):
        return self.Valeur_critique

