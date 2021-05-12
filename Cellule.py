'''
David Bernier

Mai 2021

Ce fichier contient la classe Cellule. Cette classe est utilisé lors de la constructuon d'un Complexe et sert à
modéliser la notion de simplexe ou de cube élémentaire.

'''


class Cellule :
    """
        Objet representant une celulle (simplicial ou cubique) et contient les informations pouvant etre utlie lors des
        calcul.

        ...

        Attributes
        ----------
        id : int
            l'index de la cellule dans le complexe.

        Dim : int
            la dimention de la cellule.

        Face : list
            la liste des ids des (dim-1)-faces de la cellule.

        Orientation : list
            La liste de l'orientation des faces de la cellule.

        Coface : list
            La liste des ids des cofaces de codimension 1 de la cellule.

        pair : list
            Liste comprenant la dimension et l'indexe de l'appariment fait dans le champs gradient.

        signe : int
            Signe de la cellule lors de la filtration de Morse (MS_complexe).

        Methods
        -------
        Add_Coface(id)
            Ajoute une coface à la cellules.

        Add_Face(id)
            Ajoute une face à la cellules.

        Est_Appariee()
            Vérifie si la cellule est appariée.

        Est_Critique()
            Vérifie si la cellule est critique.

        Get_Coface()
            Retourne la liste des cofaces.

        Get_Face()
            Retourne la liste des faces.

        Get_id()
            Retourne l'index de la cellules.

        Get_Orientation()
            Retourne la liste d'orientation.

        Get_Pair()
            Retourne une Liste comprenant la dimension et l'indexe de l'appariment fait dans le champs gradient.

        Get_Signe()
            Retourne le signe de la cellule.

        Set_Pair(pair)
            Fixe la valeur de la pair gradientes de la celulle.

        Set_Signe(signe)
            Fixe le signe de la cellule.
        """

    def __init__(self, id=-1, Dim=-1, Face=[], Orientation=[]):
        """
        Initialise une cellule simplicial ou cubique.

        Parameters
        ----------
        id : int
             l'index de la cellule dans le complexe.

        Dim : int
              la dimention de la cellule.

        Face : list
               la liste des ids des (dim-1)-faces de la cellule.

        Orientation : list
                      La liste de l'orientation des faces de la cellule.

        """

        self.id = id
        self.Dim = Dim
        self.Face = Face
        self.Orientation = Orientation
        self.Coface = []
        self.pair = []
        self.signe = 0

    def Add_Coface(self, id):
        """
        Ajoute une coface à la cellules.

        Parameters
        ----------
        id : int
            l'index de la cellule à ajouter dans les cofaces.

        """
        self.Coface.append(id)

    def Add_Face(self, id):
        """
        Ajoute une face à la cellules.

        Parameters
        ----------
        id : int
            l'index de la (dim-1)-cellule à ajouter dans les faces.

        """
        self.Face.append(id)

    def Est_Appariee(self):
        """
        Vérifie si la cellule est appariée.

        Returns
        -------
        val : bool
            vrai si la cellule est appariée.

        """

        if len(self.pair) == 0 :
            val = False
        else :
            val = True
        return val

    def Est_Critique(self):
        """
        Vérifie si la cellule est critique.

        Returns
        -------
         val : bool
            vrai si la cellule est critique.

        """

        if  len(self.pair) == 1:
            val = True
        else :
            val = False
        return val

    def Get_Coface(self):
        return self.Coface

    def Get_Face(self):
        return self.Face

    def Get_id(self):
        return self.id

    def Get_Orientation(self):
        return self.Orientation

    def Get_Signe(self):
        return self.signe

    def Get_Pair(self):
        return self.pair

    def Set_Pair(self, pair):
        """
        Fixe la valeur de la pair gradientes de la celulle.

        Parameters
        ----------
        pair : list
            Liste comprenant la dimension et l'indexe de l'appariment fait dans le champs gradient.

        """
        self.pair = pair

    def Set_Signe(self, signe):
        self.signe = signe

