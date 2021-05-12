'''
David Bernier

Mai 2021

Ce fichier est le main. Il sert a faire les appeles de fonctions.

'''
import time

from Complexe import Complexe
from Gestion_Image import Calcul_Valeur_Sommet, Import_Image_Fixe_Dim
# pour l'affichage de la filtration
import gudhi as gd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
import gudhi.representations



# Debut du compteur de temps
tps1 = time.perf_counter()





'''
STEP 1 : Initialisation du complexe.
'''

### Grille uniforme avec image
I, J = 16,16
comp = Complexe()
comp.Init_Rectangle(I,J)

# on donne les valeurs de l'images
image_black_scaled = Import_Image_Fixe_Dim([I,J])
Valeur_0 = Calcul_Valeur_Sommet(image_black_scaled)


'''# Affichage de l'image
plt.imshow(image_black_scaled, cmap='gray')
plt.title('Image originale en gris')
plt.show()'''



'''### Tore
comp = Complexe()
comp.Init_Tore_Exemple()
Valeur_0=[0,8,5,7,2,4,1,6,3]
'''


'''# Test Methode Calcul_CCR : Sortie a l'ecran
comp.Calcul_CCR()'''




'''
STEP 2 : Calcul du Champs gradient.
'''

CG = comp.Calcul_Champs_Gradient(Valeur_0)


# Filtration des cellules critiques qui ce touchent
CG.Epsilon_Simplification()



# Affichage des resultats
'''
Liste_Cell_Crit = CG.Get_Cell_Crit()
for d in range(len(Liste_Cell_Crit)):
    print(len(Liste_Cell_Crit[d]), ' cellules critiques de dimension ', d , ' : ', Liste_Cell_Crit[d])
#print ('pair = ', CG.Get_Pair_Cell())
'''
# CG.Affichage()
# ou
# CG.Affichage(image_black_scaled)




'''
STEP 3 : Calcul du complexe de Morse-Smale.
'''


MS = CG.Calcul_Complexe_De_Morse()


# MS.Resultat_MS_Complexe()

# MS.Affichage(image_black_scaled)


dgm = MS.Filtration()
print(dgm)
gd.plot_persistence_diagram(dgm)





'''# Test Methode Calcul_Variete
arc, m = CG.Calcul_Variete(0,0,True,[])
print(arc)'''


tps2 = time.perf_counter()
print('Temps total : ', tps2 - tps1)