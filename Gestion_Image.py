'''
David Bernier

Mai 2021

Ce fichier contient les méthodes que je me sert afin d'importer les images en teinte de gris et mettre un format
qui me permet de faire l'initialisation du complexe.

'''
import numpy as np
from skimage import io
from skimage.transform import resize
from skimage import exposure


def Import_Image_Fixe_Dim(dim, file = "spot_blanc.png", bruit = 0):
    """
    Import une image et change le format pour dim[0] par dim[1].

    Parameters
    ----------
    dim : list
        Liste des dimetion voulue de l'image.

    file : string
        Nom du fichier de l'image. Le fichier doit etre dans le dossier data. # Changer le Data_Path.

    bruit : float
        Valeur servant à ajouter du bruit issue d'une loi normal.
        Pour référance, les valeurs de l'image sont 'scale' en entre 0 et 1.

    Returns
    -------
    Valeur_0 : narray
        Liste des valeurs de l'image selon l'index du sommet corespondant.

    """

    # Repertoire des images
    Data_Path = "C:/Users/cleo6/Desktop/ecole/Maitrise/Cours/Cours topo computationelle/Projet/data/"

    # Differentes images
    # Data = "hibiscus.bmp"
    # Data = "sharingan.png"
    # Data = "spot_blanc.png"
    # Data = "image_o.png"
    # Data = "image_o2.png"
    # Data = "spot_blanc_3x.png"
    # Data = "spot_noir_3x.png"



    Data = file

    # Import de l'image
    File_Path = Data_Path + Data
    image = io.imread(File_Path)

    # On traite l'image en noir et blanc
    image_black = np.sum(image, axis=2)

    # On change le format de l'image
    image_black_scaled = resize(image_black, dim)

    # On normalise l'image entre 0 et 1.
    max = np.max(image_black_scaled)
    min = np.min(image_black_scaled)
    image_black_scaled_nomalize = (image_black_scaled - min)/(max-min)
    #image_black_scaled_nomalize = exposure.equalize_adapthist(image_black_scaled, clip_limit=0.03)

    # Ajout du bruit
    image_black_scaled_noisy = image_black_scaled_nomalize + bruit * np.random.normal(0,1,dim)

    return image_black_scaled_noisy


def Calcul_Valeur_Sommet(image_black_scaled):
    """

    Parameters
    ----------
    image_black_scaled : narray
        double vecteur numpy représantant l'image.

    Returns
    -------

    Valeur_0 : list
        La liste des valeurs de l'image en mettant les lignes de l'image boute à boute.

    """
    dim = np.shape(image_black_scaled)
    # Met les valeurs des lignes de pixels bout a bout.
    Valeur_0 = np.resize(image_black_scaled, (1, dim[0] * dim[1]))[0, :]

    return Valeur_0
