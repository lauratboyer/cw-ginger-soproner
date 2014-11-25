## GS_new-code_draft.r
## Cleaning up invertebrates et al functions for new version of the Soproner code
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: November 25, 2014
## Time-stamp: <2014-11-25 08:40:36 Laura>

## A faire
## tableaux noms de fichiers?
## calcul de densité moy + SD par:
## facteur temporel (Campagne par défaut)
fact.temp <- "Campagne"
## facteur spatial (Geomorpho, Station, Transect)
fact.spat <- c("Geomorpho","St","T")
## facteur taxonomique (Groupe, S_Groupe, Famille, Genre, G_Sp)
fact.taxo <- c("Groupe","S_Groupe","Famille","Genre","G_Sp")
## facteur expl (N_Impact)

## General version

inv.dens.gnrl <- function(filt.camp="A", fspat="Geomorpho",
                          ftaxo="Groupe", smpl.unit="St", grtax="G_Sp",
                         save=FALSE, wZeroAll=FALSE) {

    departFunk() # message de depart
    on.exit(EM())

    ### 1. ########################################################
    ### Appliquer filtres campagnes et taxonomiques ###############
    # Filtre campagnes -- défini dans les arguments de la fonction
    wf <- paste("T",filt.camp,"inv",sep="_") # formatter nom du filtre
    ta.rawF <- filtreTable(dbio, wf)

    # Filtre espèces -- défini dans les options générales
    ta.rawF <- filtreTaxo(ta.rawF, action=taxoF.incl, taxtype=taxoF.utaxo, taxnom=taxoF.nom)


  # Densité par transect: Nombre d'individus/Aire du transect -> N/(50*Ltrans)

  ### 1. #################################
  ### Densité par transect ###############
  ff <- c("Campagne","St","T","Groupe","S_Groupe","Famille","Genre","G_Sp")
  ff <- ff[1:(which(ff==grtax))] # sélectionne la partie du vecteur allant jusqu'à grtax

    ff <- c(ftemp, fspat, ftaxo)
    print(ff)

  ## 3a. ########################################
  ### Densité de chaque grtax sur le transect ###
  dens.tb <- aggregate(list("D"=ta.rawF$D), as.list(ta.rawF[,ff]), sum)
  ta.raw <- merge(dens.tb, info.transect[,c("St","Geomorpho")])


  if(smpl.unit!="T") {

      ### 3b. ######################################################
      ### Densité moyenne(SD) by Campagne/St/Espèce ###############
      ff <- ff[!(ff=="T")] # ôter le facteur transect pour calculer sur stations
      print(ff)

      tb.1 <- aggregate(list("dens.moy"=ta.raw$D),
                        as.list(ta.raw[,ff]), mean)
      tb.1.b <- aggregate(list("dens.sd"=ta.raw$D),
                          as.list(ta.raw[,ff]), sd)
      tb.all <- merge(tb.1,tb.1.b,by=ff)

      tb.all[,c("dens.moy","dens.sd")] <- round(tb.all[,c("dens.moy","dens.sd")],5)

      } else {tb.all <- ta.raw # si smpl.unit == T, on continue avec le tableau transect seulement
          names(tb.all)[names(tb.all)=="D"] <- "dens" # renomme colonne densité
          }


  # Rajouter infos additionelles
  iti <- info.transect.INV[,c("Annee","Mois","Mission","Campagne",
                                "Geomorpho", "St","N_Impact")]
  tb.all <- merge(iti, tb.all)

  # ordonner colonnes:
  ord.ff <- c("Annee", "Mois", "Geomorpho", "N_Impact", "St")
  if(smpl.unit == "T") ord.ff <- c(ord.ff,"T") # ordonner par transect si applicable
  tb.all <- tb.all[do.call(order, tb.all[,ord.ff]),]

  ### 5. ###############################
  ### Sauvegarde des tableaux ############
  if(save) {
    ftag <- paste("_Filtre-",filt.camp,"_",sep="")
    taxotag <- taxotagFunk()
    write.csv(tb.all,file=paste(tabl.dir,"Inv_DensitePar_", smpl.unit, ftag, grtax, taxotag,
                            Sys.Date(),".csv",sep=""),row.names=FALSE)
  }

  return(tb.all)
}

## General version of inv.dens.geom()

## (changed my mind and working from inv.dens.tbl() now, above
inv.dens.by.fact <- function(AS="A", ftaxo="Groupe", fspat,
                        spttcampagnes=FALSE, save=FALSE) {

    departFunk() # message de depart
    on.exit(EM())

  ### 1. ###############################
  ### Définir fonction d'aggrégation ###

  aggr.funk <- function(ff) {
  # Requiert fonction *inv.dens.tbl*
  ### (i) #############################################################################
  ### Extraire tableau des densités moyennes pour le niveau taxonomique ###############

    # valeur filtre définie dans l'argument de la fonction inv.dens.geom
    # inv.dens.tbl retourne moy + SD de densité à l'échelle taxonomique 'grtax'
    # pour chaque site
    ta.raw <- inv.dens.tbl(grtax=ftaxo, AS=AS, save=FALSE)
    ta.raw <- ta.raw[,c(ff, "dens.moy")] # ote les colonnes superflues pour l'analyse

  ### (ii) #############################################################################
  ### Calculer densité moyenne sur la géomorphologie pour le niveau taxonomique ########

    # option de rajouter les zéros
    raj.0 <- FALSE
    if(raj.0) {
    # (1) rajouter abondances nulles pour géomorphologie/zone d'impact où
    # l'espèce (grtax) n'est pas observée sur toutes les stations
    # identification des espéces uniques observées sur le grspatial
    # par campagne / géomorphologie (impact) / St / grtax
    ind1 <- unique(ta.raw[,ff])
    ind2 <- nst.par.gs(AS=AS, impact=aj.impact)$nomsSt
    # rajoute les stations où l'espèce n'a pas été observée sur une géomorpho donnée
    ind3 <- merge(ind1, ind2, all=TRUE)

    ta.raw.0 <- merge(ta.raw, ind3, all=TRUE)
    ta.raw.0$dens.moy[is.na(ta.raw.0$dens.moy)] <-  0 # mettre les D à zéro si NA après le merge
  }else{ta.raw.0 <- ta.raw}
    idnow <<- idtax


    # (2) maintenant calculer moyenne/SD par groupement spatial
    tb.1 <- aggregate(list("dens.moy"=ta.raw.0$dens.moy),
                      as.list(ta.raw.0[,ff]), mean)
    tb.1.sd <- aggregate(list("dens.sd"=ta.raw.0$dens.moy),
                         as.list(ta.raw.0[,ff]), sd)
    tb.all <- merge(tb.1,tb.1.sd,by=ff)

  ### (iii) #####################################
  ### Formattage des tableaux de sortie ########
  ### (1) ordonner les valeurs par colonnes:
    tb.all <- tb.all[order(tb.all[,ff[1]],tb.all[,ff[2]],tb.all[,ff[3]]),]
    if(length(ff)>3) {
      tb.all <- tb.all[order(tb.all[,ff[1]],tb.all[,ff[2]],tb.all[,ff[3]],tb.all[,ff[4]]),]}

  ### (2) rajout colonnes info sur la campagne/géo
    tb.all <- merge(info.transect.INV.geo, tb.all, by=c("Campagne","Geomorpho"), sort=FALSE)
    tb.all <- tb.all[,c("Annee","Mois","Mission",ff, "dens.moy","dens.sd")]
    tb.all <- with(tb.all,tb.all[order(Annee, Mois),])

  ### (3) sauvegarde des fichiers
    if(save) {
      ftag <- paste("_Filtre-",AS,"_",sep="")
      top10tag <- ifelse(sum(grepl("G_Sp",ff))==1,paste("top10_",top10year,"_",sep=""),"")
      taxotag <- taxotagFunk()
      if(spttcampagnes) top10tag <- ifelse(sum(grepl("G_Sp",ff))==1,
                                           paste(top10tag,"SpObsTTcamp_", sep=""), "")

      write.csv(tb.all,file=paste(tabl.dir,"Inv_DensiteTS_",
                         gsub("N_","",paste(ff[-2],collapse="-")), ftag, top10tag,taxotag,
                            Sys.Date(),".csv",sep=""),row.names=FALSE) }
    return(tb.all)
  }

  ### 3. ###############################################
  ### Densité moyenne(SD) par facteurs X ###############
  ff.list <- list(c("Geomorpho","Campagne","Groupe"),
                  c("Geomorpho","Campagne","S_Groupe"))

  # Rajouter facteur impact si spécifié
  #if(aj.impact) ff.list <- lapply(ff.list, function(i) c(i[1:2],"N_Impact",i[3]))

  # Appliquer les facteurs à la fonction d'aggrégation
  tbl.list <- lapply(ff.list, aggr.funk)
  names(tbl.list) <- c("Groupe","S_Groupe") # deux tableaux produits: groupe/sous-gr

  ### 4. ########################################################################
  ### Densité moyenne(SD) pour les 10 espèces les plus abondantes ###############

  # Calcul de la densité moyenne de chaque espèce par morphologie
  # Rajouter G_Sp à la liste des facteurs et aggréger...
  ff.list.sp <- lapply(ff.list, function(i) c(i,"G_Sp"))


  ########################################################
  ########################################################
  return(list("allsp"=tbl.list))
  }
