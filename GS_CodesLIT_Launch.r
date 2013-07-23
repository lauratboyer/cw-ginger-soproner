## Ginger/Soproner: Produits/Analyses poissons
# ** Code central pour lancer analyses couvertures LIT moy/SE
# Time-stamp: <2013-07-22 16:33:37 Laura>

########################################################
########################################################

Run.LIT.all <- function() {

    if(!exists("LIT.tableau.brut")) source("GS_CodesLIT_Couvrt.r")

    # Tableaux seulement:
    ####################

    # Tableau brut formatté des données par transect:
    dmm <- LIT.tableau.brut(save=TRUE)

    # Données par transect, avec filtre annuel et semestriel
    dmm <- LIT.tableau.brut(save=TRUE, AS="A")
    dmm <- LIT.tableau.brut(save=TRUE, AS="S")

    # Moyenne + SE par type de coraux dans une categorie par geomorphologie
    Dmm <- LIT.resume(save=TRUE) # categorie par defaut: Coraux_Gen
    dmm <- LIT.resume(ff="Coraux_Acro",save=TRUE) # categorie acroporidae


    # Graphiques + tableaux:
    ########################

    LIT.bp1() #graphique des valeurs obtenues en lit.tb.1
    LIT.bp1(ff2="Coraux_Acro")

    # Statistiques calculees sur geomorpho x impact (tableau + graph)
    # couverture moyenne par geomorpho, comparaison impact/non impact
    LIT.ts1() # annuel
    LIT.ts1("S") # semestriel

    # Comparaison couverture moyenne par station, par geomorpho/impact/type de corail
    LIT.ts2() # annuel
    LIT.ts2("S") # semestriel
}


