## Ginger/Soproner: Produits/Analyses poissons
# ** Code central pour lancer analyses densité/abondance/diversité
# Time-stamp: <2013-07-15 15:41:44 Laura>

########################################################
########################################################
if(!exists("poissons.tableau.brut")) source("GS_CodesPoissons_BioDens.r")

Run.poissons.all <- function() {

    # Tableau brut, tableau des données de comptage poissons reformatté
    # pour inclure les densités et autres infos
    poissons.brut <<- poissons.tableau.brut(save=TRUE)

    # Tableau densité + biodiversité par Campagne/Transect/Espèce
    BDtable <<- BioDens.sp.poissons() # tableau intermédiaire
    # requis par les fonctions poissons.ts1() et poissons.ts2()

    # Tableau synthèse 1:
    # moyenne et écart type pour: densité (d), biomasse (b),
    # richesse spécifique (rs), taille moyenne (tm)
    # calculé sur:
    # toutes espèces confondues (tot), commerciales (com), carnivore (car),
    # herbivore (her), piscivore (pis), planctonophage (pla),
    # sédentaire (sed), territoriale (ter), mobile (mo), très mobile (tmo)
    # cible pêche nouvelle-calédonie, ou non
    TS1 <- poissons.ts1(save=TRUE)

    # Tableau synthèse 2:
    # pour chaque espèce, moyenne de densité pour chaque année et totale,
    # pour chaque station et totale
    TS2 <- poissons.ts2(save=TRUE)

    # Graphiques -- non lancés pour le moment
    graph.poissons <- FALSE # TRUE pour produire les graphs poissons
    if(graph.poissons) poissons.p3(save=TRUE)
}
