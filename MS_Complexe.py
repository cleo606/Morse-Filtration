'''
David Bernier

Mai 2021

Ce fichier contion la classe MS_complexe. Cette classe est utilisé lors du calcul du complexe de Morse dans la
classe Champ_Gradient.

'''
from Cellule import Cellule
import numpy as np
import operator
import matplotlib.pyplot as plt



class MS_Complexe:
    """Complexe de Morse-Smale.

    ...

    Attributes
    ----------
    liste_cellule : List
        La liste des listes de d-cellules critiques.

    valeur_cellule : List
        La liste des listes de valeur des d-cellules critiques.

    arcs : list
        La liste v-chemins descendant. arc[dim][no cell][no chemin] = liste de (dim)-cellule.
        Même ordre que les faces le la cellule (dimno cell)

    IJ : tuple
        Dimention de l'image (utiliser dans Affichage).

    Methods
    -------
    Add_Arc(dim, index, pair)
        Ajoute une arc dans le MS_complexe.

    Affichage()
        Affichage du complexe de Morse avec la librairie pyplot.

    Resultat_MS_Complexe()
        Affiche le nombre arc de Morse a l'écran.

    Get_Arc_Multiplicity(dim, index, pair)
        Compte le nombre d'arc entre index et pair.

    Filtration()
        Procéde à la filtration du complexe de Morse.

    """

    def __init__(self, Champ_gradient):
        """Initialise la variable MS_complexe avec les valeurs sur les sommets et un ensemble d'arc nul.

        Parameters
        ----------
        Champ_gradient : class Champ_Gradient
            Champs gradient sur laquel est basée le complexe de Morse.

        """

        self.CG = Champ_gradient
        self.Liste_Cellule = []
        self.valeur_cellule =  Champ_gradient.Get_Valeur_Critique()
        self.arcs = []

        # Ajout des cellules critiques dans Liste_Cellule et initialiser
        Cellule_Crititque = Champ_gradient.Get_Cell_Crit()
        for d in range(len(Cellule_Crititque)):
            self.Liste_Cellule.append([])
            self.arcs.append([])
            for i in range(len(Cellule_Crititque[d])):
                self.Liste_Cellule[d].append(Cellule(Cellule_Crititque[d][i], d, [], []))
                self.arcs[d].append([])

    def Add_Arc(self, dim, index, pair, arc):
        """Ajoute une arc dans le MS_complexe

        L'arc est associer avec la dim-cellule en position index mais contient des (dim+1)-cellules.

        Parameters
        ----------
        dim : int
              La dimension de la plus petite cellule de l'arc.

        index : int
            l'index de la cellule de l'arc de plus petite dimension.

        pair : int
            L'index de la cellule de l'arc de plus grande dimension.

        arc : list
            Liste des index des (dim+1)-cellules représentant le v-chemin entre les deux cellules critiques

        """

        #On test si l'arette des deja la

        deja_la = False
        coface = self.Liste_Cellule[dim][index].Get_Coface()
        for i in range(len(self.arcs[dim][index])):
            if self.arcs[dim][index][i] == arc and  coface[i]==pair:
                deja_la = True


        if not deja_la:
            # arriver
            self.Liste_Cellule[dim][index].Add_Coface(pair)
            # depart
            self.Liste_Cellule[dim+1][pair].Add_Face(index)
            self.arcs[dim][index].append(arc)

    def Affichage(self, image=None):
        """Affichage du complexe de Morse avec la librairie pyplot.

        Parameters
        ----------
        image : narray
            Si la valeur est donner, on affiche l'image en background.

        """

        ax = plt.axes()
        I, J = self.CG.Get_IJ()

        # Pour les 0-cellules critiques
        list_0 = [[], []]
        for j in self.Liste_Cellule[0]:
            i = j.Get_id()
            list_0[0].append(i % I)  # modulo
            list_0[1].append(i // I)  # divison entiere

        # Ajoutes les points
        ax.scatter(list_0[0], list_0[1], s=150, label='Point critique', color='m')

        # Pour les 1-cellules critiques

        for j in self.Liste_Cellule[1]:
            i = j.Get_id()
            sommet = self.CG.Get_Cellules(1, i).Get_Face()
            p_0 = [sommet[0] % I, sommet[0] // I]
            p_1 = [sommet[1] % I, sommet[1] // I]
            # Affichage de l'arette P_0 -> p_1
            ax.arrow(p_0[0], p_0[1], p_1[0] - p_0[0], p_1[1] - p_0[1], head_width=0.01, head_length=0.01, color='g',
                     width=0.1)

        # Pour les 2-cellules critiques

        list_2 = [[], []]
        for j in self.Liste_Cellule[2]:
            i = j.Get_id()
            list_2[0].append(i % (I - 1) + 0.5)  # modulo
            list_2[1].append((i // (I - 1)) + 0.5)

        # Ajoutes des face en fome de caree bleu
        ax.scatter(list_2[0], list_2[1], s=60, label='Face critique', marker='s', color='b')




        # Ajout de la grille
        ax.set_xticks(np.arange(0, I))
        ax.set_yticks(np.arange(0, J))
        ax.grid(True)

        # Ajout de l'image
        if type(image) != type(None):
            ax.imshow(image, cmap='gray')



        # Ajout des arcs de Morse
        # v-chemin 0-1
        for liste_chemin in self.arcs[0] :
            for chemin in liste_chemin :
                for arete in chemin[1 :] :
                    # arete est l'index de la 1-cellule
                    sommet = self.CG.Get_Cellules(1, arete).Get_Face()
                    d, pair = self.CG.Get_Cellules(1, arete).Get_Pair()
                    if sommet[0] == pair :
                        p_0 = [sommet[0] % I, sommet[0] // I]
                        p_1 = [sommet[1] % I, sommet[1] // I]
                    else:
                        p_1 = [sommet[0] % I, sommet[0] // I]
                        p_0 = [sommet[1] % I, sommet[1] // I]

                    # Ajout du vecteur p_0 -> p_1
                    if p_1[0] - p_0[0] == 0:
                        # Vertical
                        if p_1[1] - p_0[1] > 0:
                            # Vers le haut
                            ax.arrow(p_0[0], p_0[1], 0, 0.7, head_width=0.3, head_length=0.4, color='r', alpha=0.5)

                        else:
                            # Vers le bas
                            ax.arrow(p_0[0], p_0[1], 0, -0.7, head_width=0.3, head_length=0.4, color='r', alpha=0.5)


                    else:
                        # Horizontal
                        if p_1[0] - p_0[0] > 0:
                            # Vers la droite
                            ax.arrow(p_0[0], p_0[1], 0.7, 0, head_width=0.4, head_length=0.4, color='r', alpha=0.5)

                        else:
                            # Vers la gauche
                            ax.arrow(p_0[0], p_0[1], -0.7, 0, head_width=0.4, head_length=0.4, color='r', alpha=0.5)


        # Pour les v-chemin 1-2

        for liste_chemin in self.arcs[1] :
            for chemin in liste_chemin :
                for face in chemin[1 :] :
                    d, arette_pair = self.CG.Get_Cellules(2, face).Get_Pair()
                    id = self.CG.Get_Cellules(1, arette_pair).Get_Face()
                    p_0 = [(id[0] % I + id[1] % I) / 2, ((id[0] // I) + (id[1] // I)) / 2]

                    #Centre du cube
                    p_1 = [face % (I-1) + 0.5, (face // (I-1)) + 0.5]

                    # Ajout du vecteur p_0 -> p_1
                    if p_1[0] - p_0[0] == 0:
                        # Vertical
                        if p_1[1] - p_0[1] > 0:
                            ax.arrow(p_0[0], p_0[1], 0, 0.42, head_width=0.3, head_length=0.4, color='y', alpha=0.5)

                        else:
                            ax.arrow(p_0[0], p_0[1], 0, -0.42, head_width=0.3, head_length=0.4, color='y', alpha=0.5)


                    else:
                        # Horizontal
                        if p_1[0] - p_0[0] > 0:
                            ax.arrow(p_0[0], p_0[1], 0.42, 0, head_width=0.4, head_length=0.4, color='y', alpha=0.5)

                        else:
                            ax.arrow(p_0[0], p_0[1], -0.42, 0, head_width=0.4, head_length=0.4, color='y', alpha=0.5)

        plt.title('Complexe de Morse-Smale')
        plt.show()

    def Resultat_MS_Complexe(self):
        """Affiche les arc de Morse a l'ecran.

        """

        print('Complexe de Morse :')
        for d in np.arange(0, len(self.Liste_Cellule) - 1):
            for i in range(len(self.Liste_Cellule[d])):
                for j in range(len(self.Liste_Cellule[d+1])):
                    mul = self.Get_Arc_Multiplicity(d,i,j)
                    if mul !=0 :
                        id_0 = self.Liste_Cellule[d][i].Get_id()
                        id_1 = self.Liste_Cellule[d+1][j].Get_id()
                        print('Arrete de multiplicite ', mul, 'entre', [d,id_0], 'et', [d+1,id_1])

    def Get_Arc_Multiplicity(self, dim, index, pair):
        """Compte le nombre d'arc entre index et pair.

        Parameters
        ----------
        dim : int
              La dimention de la plus petit cellule de l'arc.

        index : int
                l'index de la cellule de l'arc de plus petite dimension.

        pair : int
               L'index de la cellule de l'arc de plus grande dimension.

        Returns
        -------
        val : int
            Nombre d'arc entre index et pair

        """

        '''
        # Option d'affichage des chemins
        val = 0
        temp = self.Liste_Cellule[dim][index].Get_Coface()
        I,J = self.CG.Get_IJ()
        for i in range(len(temp)):
            if temp[i] == pair:
                val = val + 1
                print('chemin = ',self.arcs[dim][index][i])
                for j in self.arcs[dim][index][i]:
                    sommet = self.CG.Get_Cellules(1, j).Get_Face()
                    p_0 = [sommet[0] % I, sommet[0] // I]
                    p_1 = [sommet[1] % I, sommet[1] // I]
                    print('p_0 = ', p_0, ', p_1 = ', p_1)

        return val
        '''

        val = 0
        for i in self.Liste_Cellule[dim][index].Get_Coface():
            if i == pair :
                val = val + 1
        return val

    def Filtration(self):
        """
        Procède au calcul des pairs de cellules critiuqes selon l'article de Sousbie.
        Algo 2D.

        Returns
        -------
        Liste_pair : List
            Liste des pairs critiques. Chaque element de la liste est un tuple contenant la dimession de la plus petite
            cellule et un tuple des valeurs critiques de chaque cellule de la pair.

        """

        # STEP 1 : Ordonancement des cellules pour la filtration
        Liste = []
        for d in range(len(self.valeur_cellule)):
            for i in range(len(self.valeur_cellule[d])):
                Liste.append((d, i, self.valeur_cellule[d][i]))

        # On filtre sur la valeur critique puis sur la dim
        # Plusieur equalites entre les valeurs critique car defini avec le maximum sur les sommets.
        # On faut changer l'ordre pour consider un ordre lexicographique sur les sommet de chanque cellules
        # Il faut definir un operateur itemgetter.
        list_sorted = sorted(Liste, key=operator.itemgetter(2, 0))

        # Classer les cellule entre positive et negative.
        Liste = []
        Liste_negatif = []
        for cell in list_sorted:
            dim = cell[0]
            index = cell[1]
            if dim == 0 :
                # les 0-cellules sont toujours positives.
                self.Liste_Cellule[cell[0]][cell[1]].Set_Signe(1)
                Liste.append([cell])

            elif dim == 2 :
                # Les 2-cellule sont toujours negatives.
                self.Liste_Cellule[cell[0]][cell[1]].Set_Signe(-1)
                Liste_negatif.append(cell)

            else :
                # dim = 1
                arc_inf = self.Liste_Cellule[dim][index].Get_Face()
                # On cherche le nombre d'atache de cell dans Liste
                element_trouvee = []
                # pour chaque arc de dimention inferieur,
                for j in arc_inf :
                    # pour chaque element deja listée
                    # find j,
                    for l in range(len(Liste)):
                        for i in range(len(Liste[l])):
                            # si on trouve le bon élément
                            if j == Liste[l][i][1] and Liste[l][i][0] == dim - 1:
                                element_trouvee.append(l)

                if element_trouvee[0] == element_trouvee[1]:
                    # La cellule est positive.
                    Liste[element_trouvee[0]].append(cell)
                    self.Liste_Cellule[cell[0]][cell[1]].Set_Signe(1)

                else:
                    # La cellule est negative.
                    Liste[element_trouvee[0]] = Liste[element_trouvee[0]] + Liste[element_trouvee[1]]
                    Liste[element_trouvee[0]].append(cell)
                    Liste.pop(element_trouvee[1])
                    self.Liste_Cellule[cell[0]][cell[1]].Set_Signe(-1)
                    Liste_negatif.append(cell)

                if len(element_trouvee) != 2:
                    print('PROBLEME Calcul MS_complexe.Filtration')
                    # il faut changer la méthode de determiner signe selon le calcul de l'homoloogie
                    # on peut regarder les changement dans l'homologie mais cest plus long que juste faire des
                    # calcul de cycle dans un graphe comme maintenant. En 3D, on peut faire ce calule du signe des
                    # 1-cells critique comme présentement et calculer le signe des 2-cells avec la filtration de -f.
                    # C'est ce qui proposer dans soubies dans le cas 3D.


        # STEP 2 : Trouver la pair des cellules negatives (Algo 2 sousbie).

        # Initialisation des listes
        ppair=[[],[]]
        ppair_val = [[], []]
        cycles = []
        for d in range(len(self.Liste_Cellule)):
            cycles.append([])
            for i in range(len(self.Liste_Cellule[d])):
                cycles[d].append([])

        for cell in Liste_negatif:
            cur_set = []
            val_set = []
            a_nei = self.Liste_Cellule[cell[0]][cell[1]].Get_Face()
            for b in a_nei:
                # Si b est l'index d'une cellule positive :
                if self.Liste_Cellule[cell[0]-1][b].Get_Signe() == 1 :
                    cur_set.append(b)
                    val_set.append(self.valeur_cellule[cell[0]-1][b])

            while len(cur_set) != 0:
                id_cur = np.argmax(val_set) # index de la liste cur_set de l'élement maximum
                sigma_cur = cur_set[id_cur] # l'éléement maximum
                if len(cycles[cell[0]-1][sigma_cur])==0:
                    # On paire cell et sigma_cur
                    cycles[cell[0]-1][sigma_cur] = cur_set
                    cycles[cell[0]][cell[1]] = cur_set
                    ppair_val[cell[0]-1].append((val_set[id_cur], self.valeur_cellule[cell[0]][cell[1]]))
                    ppair[cell[0] - 1].append((sigma_cur, cell[1]))
                    break
                else:
                    cur_set.pop(id_cur)
                    val_set.pop(id_cur)
                    for b in cycles[cell[0]-1][sigma_cur]:
                        if b != sigma_cur:
                            cur_set.append(b)
                            val_set.append(self.valeur_cellule[cell[0] - 1][b])
            continue

        # STEP 3 : Traitement des cellules critiques restantes
        for d in range(len(self.Liste_Cellule)):
            for i in range(len(self.Liste_Cellule[d])):
                if len(cycles[d][i])==0:
                    ppair[d].append((i,'None'))
                    ppair_val[d].append((self.valeur_cellule[d][i],float('inf')))


        Liste_pair = []
        for d in range(len(ppair_val)):
            for pair in ppair_val[d]:
                Liste_pair.append((d,pair))

        return Liste_pair




