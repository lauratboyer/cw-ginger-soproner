## Ginger/Soproner: Produits/Analyses invertébrés
# ** Code central pour lancer analyses densité/abondance/diversité
# *et* spécifier groupes taxonomiques à analyser
# Time-stamp: <2013-07-25 09:18:17 Laura>

########################################################
########################################################
Run.INV.biodiv <- function() {

  source("GS_CodesInvertebres_Biodiv.r")
  all.bio <- inv.biodiv(save=TRUE)
  all.bio.T <- inv.biodiv(qunit="T",save=TRUE)
  all.bio.geom <- inv.biodiv.geom(AS="A",save=TRUE)
  all.bio.geom <- inv.biodiv.geom(AS="S",save=TRUE)
  invsr.A <- inv.sprich.tbl(AS="A",save=TRUE)
  invsr.S <- inv.sprich.tbl(AS="S",save=TRUE)

  # Richesse spécifique par transect et aggrégation taxonomique:
  dd <- sapply(c("Groupe","S_Groupe","Famille"), function(gt)
               sprich.by.aggrtaxo(AS="A", grtax=gt, save=TRUE))
  dd <- sapply(c("Groupe","S_Groupe","Famille"), function(gt)
               sprich.by.aggrtaxo(AS="S", grtax=gt, save=TRUE))
}

Run.INV.densite <- function(tabl.seulement = TRUE) {

  source("GS_CodesInvertebres_Densite.r")
  ###########################
  ######## Tableaux #########

  # Boucle par niveaux taxonomiques et filtres (incluant filtre absent)
  # Densité par transect:
  dd <- sapply(c("Genre","G_Sp"),
               function(tt) sapply(c("A","S","Absent"), function(ff)
                                   inv.dens.tbl(grtax=tt, smpl.unit="T", save=TRUE, AS=ff)))

  # Densité par transect, incluant les densité nulles, avec filtre A et S seulement:
  dd <- sapply(c("Genre","G_Sp"),
               function(tt) sapply(c("A","S","Absent"), function(ff)
                                   inv.dens.tbl.parT(grtax=tt, smpl.unit="T", save=TRUE, AS=ff, wZeroAll=TRUE)))

  # Densité par station:
  dd <- sapply(c("Groupe","S_Groupe","Famille","Genre","G_Sp"),
               function(tt) sapply(c("A","S","Absent"), function(ff)
                                   inv.dens.tbl(grtax=tt, save=TRUE, AS=ff)))

  # Densité moyenne par géomorphologie (et impact) et groupe/sous-groupe,
  # toutes espèces et pour les 10 espèces les plus abondantes en (1) 2006 et (2) 2006-2011
  # Tourner avec spttcampagnes FALSE et TRUE
  # spttcampagnes -> utilise seulement les espèces observées sur toutes les campagnes
  # lors du calcul des top10 abondance
  dd <- sapply(c(FALSE, TRUE), function(ii)
               sapply(c("A","S"), function(ff)
                      sapply(c(FALSE, TRUE), function(ss)
                             inv.dens.geom(AS=ff, aj.impact=ii, spttcampagnes=ss))))

  #############################
  ######## Graphiques #########

  if(!tabl.seulement) {
    # Graphs de séries temporelles: densité moyenne par géomorphologie/groupe
    inv.graph.TS()
    inv.graph.TS(wff="S_Groupe")

    # Graphs de séries temporelles: densité moyenne des 10 espèces les
    # plus abondantes par géomorphologie/groupe
    inv.graph.TS.top10()
    inv.graph.TS.top10(wff="S_Groupe")
    inv.graph.TS.top10(top10year=2006)
    inv.graph.TS.top10(wff="S_Groupe", top10year=2006)
}
}



