## Ginger/Soproner: Produits/Analyses invertébrés
# ** Code central pour lancer analyses densité/abondance/diversité
# *et* spécifier groupes taxonomiques à analyser
# Time-stamp: <2014-03-12 11:44:47 Laura>

########################################################
########################################################
Run.INV.biodiv <- function(wfiltre=c("A","S")) {

  source("GS_CodesInvertebres_Biodiv.r")
  all.bio <- inv.biodiv(save=TRUE)
  all.bio.T <- inv.biodiv(qunit="T",save=TRUE)
  dmm <- sapply(wfiltre, function(fltre) inv.biodiv.geom(AS=fltre,save=TRUE))
  dmm <- sapply(wfiltre, function(fltre) inv.sprich.tbl(AS=fltre,save=TRUE))

  # Richesse spécifique par transect et aggrégation taxonomique:
  dd <- sapply(wfiltre, function(fltre)
               sapply(c("Groupe","S_Groupe","Famille"), function(gt)
               sprich.by.aggrtaxo(AS=fltre, grtax=gt, save=TRUE)))
}

Run.INV.densite <- function(wfiltre=c("A","S","Absent"), tabl.seulement = TRUE) {

  source("GS_CodesInvertebres_Densite.r")
  ###########################
  ######## Tableaux #########

  # Boucle par niveaux taxonomiques et filtres (incluant filtre absent)
  # Densité par transect:
  dd <- sapply(c("Genre","G_Sp"),
               function(tt) sapply(wfiltre, function(ff)
                                   inv.dens.tbl(grtax=tt, smpl.unit="T", save=TRUE, AS=ff)))

  # Densité par transect, incluant les densité nulles, avec filtre A et S seulement:
  dd <- sapply(c("Genre","G_Sp"),
               function(tt) sapply(wfiltre, function(ff)
                                   inv.dens.tbl(grtax=tt, smpl.unit="T", save=TRUE, AS=ff, wZeroAll=TRUE)))

  # Densité par station:
  dd <- sapply(c("Groupe","S_Groupe","Famille","Genre","G_Sp"),
               function(tt) sapply(wfiltre, function(ff)
                                   inv.dens.tbl(grtax=tt, save=TRUE, AS=ff)))

  # Densité moyenne par géomorphologie (et impact) et groupe/sous-groupe,
  # toutes espèces et pour les 10 espèces les plus abondantes en (1) 2006 et (2) 2006-2011
  # Tourner avec spttcampagnes FALSE et TRUE
  # spttcampagnes -> utilise seulement les espèces observées sur toutes les campagnes
  # lors du calcul des top10 abondance
  dd <- sapply(c(FALSE, TRUE), function(ii)
               sapply(wfiltre, function(ff)
                      sapply(c(FALSE, TRUE), function(ss)
                             inv.dens.geom(AS=ff, aj.impact=ii, spttcampagnes=ss, save=TRUE))))

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
