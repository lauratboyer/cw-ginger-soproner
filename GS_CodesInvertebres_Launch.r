## Ginger/Soproner: Produits/Analyses invert�br�s
# ** Code central pour lancer analyses densit�/abondance/diversit�
# *et* sp�cifier groupes taxonomiques � analyser
# Time-stamp: <2013-08-06 17:24:08 Laura>

########################################################
########################################################
Run.INV.biodiv <- function(wfiltre=c("A","S")) {

  source("GS_CodesInvertebres_Biodiv.r")
  all.bio <- inv.biodiv(save=TRUE)
  all.bio.T <- inv.biodiv(qunit="T",save=TRUE)
  dmm <- sapply(wfiltre, function(fltre) inv.biodiv.geom(AS=fltre,save=TRUE))
  dmm <- sapply(wfiltre, function(fltre) inv.sprich.tbl(AS=fltre,save=TRUE))

  # Richesse sp�cifique par transect et aggr�gation taxonomique:
  dd <- sapply(wfiltre, function(fltre)
               sapply(c("Groupe","S_Groupe","Famille"), function(gt)
               sprich.by.aggrtaxo(AS=fltre, grtax=gt, save=TRUE)))
}

Run.INV.densite <- function(wfiltre=c("A","S","Absent"), tabl.seulement = TRUE) {

  source("GS_CodesInvertebres_Densite.r")
  ###########################
  ######## Tableaux #########

  # Boucle par niveaux taxonomiques et filtres (incluant filtre absent)
  # Densit� par transect:
  dd <- sapply(c("Genre","G_Sp"),
               function(tt) sapply(wfiltre, function(ff)
                                   inv.dens.tbl(grtax=tt, smpl.unit="T", save=TRUE, AS=ff)))

  # Densit� par transect, incluant les densit� nulles, avec filtre A et S seulement:
  dd <- sapply(c("Genre","G_Sp"),
               function(tt) sapply(wfiltre, function(ff)
                                   inv.dens.tbl(grtax=tt, smpl.unit="T", save=TRUE, AS=ff, wZeroAll=TRUE)))

  # Densit� par station:
  dd <- sapply(c("Groupe","S_Groupe","Famille","Genre","G_Sp"),
               function(tt) sapply(wfiltre, function(ff)
                                   inv.dens.tbl(grtax=tt, save=TRUE, AS=ff)))

  # Densit� moyenne par g�omorphologie (et impact) et groupe/sous-groupe,
  # toutes esp�ces et pour les 10 esp�ces les plus abondantes en (1) 2006 et (2) 2006-2011
  # Tourner avec spttcampagnes FALSE et TRUE
  # spttcampagnes -> utilise seulement les esp�ces observ�es sur toutes les campagnes
  # lors du calcul des top10 abondance
  dd <- sapply(c(FALSE, TRUE), function(ii)
               sapply(wfiltre, function(ff)
                      sapply(c(FALSE, TRUE), function(ss)
                             inv.dens.geom(AS=ff, aj.impact=ii, spttcampagnes=ss, save=TRUE))))

  #############################
  ######## Graphiques #########

  if(!tabl.seulement) {
    # Graphs de s�ries temporelles: densit� moyenne par g�omorphologie/groupe
    inv.graph.TS()
    inv.graph.TS(wff="S_Groupe")

    # Graphs de s�ries temporelles: densit� moyenne des 10 esp�ces les
    # plus abondantes par g�omorphologie/groupe
    inv.graph.TS.top10()
    inv.graph.TS.top10(wff="S_Groupe")
    inv.graph.TS.top10(top10year=2006)
    inv.graph.TS.top10(wff="S_Groupe", top10year=2006)
}
}



