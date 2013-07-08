## Ginger/Soproner: Produits/Analyses poissons
# ** Code central pour lancer analyses couvertures LIT moy/SE
# Time-stamp: <2013-07-04 11:27:11 Laura>

########################################################
########################################################

Run.LIT.all <- function() {

    # Tableaux seulement:
    ####################

    # Tableau brut formatté des données par transect:
    dmm <- LIT.tableau.brut(save=TRUE)
#    dmm <- TB.lit(save=TRUE)

    # Données par transect, avec filtre annuel et semestriel
    dmm <- TB.lit(save=TRUE, AS="A", filtre=TRUE)
    dmm <- TB.lit(save=TRUE, AS="S", filtre=TRUE)

    # Moyenne + SE par type de coraux dans une categorie par geomorphologie
    Dmm <- lit.tb.1(save=TRUE) # categorie par defaut: Coraux_Gen
    dmm <- lit.tb.1(ff="Coraux_Acro",save=TRUE) # categorie acroporidae


    # Graphiques + tableaux:
    ########################

    lit.BP.1() #graphique des valeurs obtenues en lit.tb.1
    lit.BP.1(ff2="Coraux_Acro")

    # Statistiques calculees sur geomorpho x impact (tableau + graph)
    # couverture moyenne par geomorpho, comparaison impact/non impact
    lit.TS.1() # annuel
    lit.TS.1("S") # semestriel

    # Comparaison couverture moyenne par station, par geomorpho/impact/type de corail
    lit.TS.2() # annuel
    lit.TS.2("S") # semestriel
}


