###############################################################################
#
#			PROJET 7 RESYS 2018 - M2 BIM
#
###############################################################################

#Auteurs
- Clotilde Garrido
- Adèle Weber

#Code
Codé en R version 3.4.3

#Dépendances :
- igraph
- optparse

#But:
Outil qui donne différentes mesures d'un réseau donné en entrée.

#Ligne de commande sur le terminal:
Rscript --vanilla projet.R -f <network-file> -e <format> [-t]

#Options:
- -f/--file : suivi du nom du fichier du réseau, obligatoire.
- -e/--extension : suivi du format du fichier du réseau, au choix parmis : "pajek", "ncol", "lgl", "edgelist", "graphml",
              "dimacs", "graphdb", "gml", "dl". Obligatoire.
- -t/--traceback : si l'option est présente, le programme fera un troisième output avec les plus courts chemins par paires de points

#Input:
- fichier du réseau, au format "pajek", "ncol", "lgl", "graphml",
              "dimacs", "graphdb", "gml", "dl", "edgelist"

#Output:
- fichier output.txt : différentes mesures du graphe
- image degree_distribution.png : graphique avec la distribution des degrés
- avec l'option -t, output_traceback.txt : fichier avec les plus courts chemins trouvés par une voire les 3 méthodes de calcul

#Exemple de commande :
Rscript --vanilla projet.R -f Graph.graphml -e graphml -t
