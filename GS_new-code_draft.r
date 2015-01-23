## GS_new-code_draft.r
## Cleaning up invertebrates et al functions for new version of the Soproner code
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: November 25, 2014
## Time-stamp: <2015-01-13 10:23:29 Laura>

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
fact.expl <- c("N_Impact","Cote","Lieu")
ftaxo.defaut <- "Groupe"
fspat.defaut <- "Geomorpho"

## General version of inv.dens.tbl
INV.dens.gnrl <- function(filt.camp="X", fspat=fspat.defaut,
                          ftaxo=ftaxo.defaut, ftemp="Campagne",
                          fexpl, smpl.unit="St", grtax="G_Sp",
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

    ## Densité par transect: Nombre d'individus/Aire du transect -> N/(50*Ltrans)
    ## Calculée lors du formattage initial de dbio

    ### 1. #################################
    ### Densité par transect ###############

    # ****** remove this at some point
    #ff <- c("Campagne","St","T","Groupe","S_Groupe","Famille","Genre","G_Sp")
    #ff <- ff[1:(which(ff==grtax))] # sélectionne la partie du vecteur allant jusqu'à grtax

    ## Définir les facteurs d'aggrégation à partir des arguments donnés à la fonction
    ff <- c(ftemp, fspat, ftaxo)
    if(!missing(fexpl)) ff <- c(ff, fexpl)

    ## 3a. ########################################
    ### Densité de chaque grtax sur le transect ###
    ff.transect <- unique(c(ff, "St", "T"))
    dens.par.t <- aggregate(list("D"=ta.rawF$D), as.list(ta.rawF[,ff.transect]), sum)

    # on continue sur l'aggrégation spatiale n'est pas
    # au niveau du transect
    if(!("T" %in% fspat)) {

      ### Densité moyenne(SD) by Campagne/St/Espèce ###############
      tb.1 <- aggregate(list("dens.moy"=dens.par.t$D),
                        as.list(dens.par.t[,ff]), mean)
      tb.1.b <- aggregate(list("dens.sd"=dens.par.t$D),
                          as.list(dens.par.t[,ff]), sd)
      tb.all <- merge(tb.1,tb.1.b,by=ff)

      # ... et on arrondi à 5 valeurs décimales la moyenne + SD
      tb.all[,c("dens.moy","dens.sd")] <- round(tb.all[,c("dens.moy","dens.sd")],5)

      } else {
        tb.all <- dens.par.t # si smpl.unit == T, on continue avec le tableau transect seulement
        names(tb.all)[names(tb.all)=="D"] <- "dens" # renomme colonne densité, pas de moyenne
        }


    # Rajouter infos additionelles (ôté temporairement, A_2006 a deux mois)
    # iti <- unique(info.transect.INV[,c("Campagne","Annee","Mois","Mission")])
    # tb.all <- merge(iti, tb.all, all.y=TRUE)
    message("need to figure out which month to add for A_2006 given sampling in Nov 2006 and Feb 2007")

    # ordonner colonnes:
    ord.ff <- c("Campagne", fspat)
    tb.all <- tb.all[do.call(order, tb.all[,ord.ff]),]

    ### 5. ###############################
    ### Sauvegarde des tableaux ##########
    if(save) {
    ftag <- paste("_Filtre-",filt.camp,"_",sep="")
    taxotag <- taxotagFunk()
    expltag <- ifelse(missing(fexpl),"",paste("x",fexpl,"_",sep=""))

    write.csv(tb.all,file=paste(tabl.dir,"Inv_DensitePar_", smpl.unit, expltag,
                       ftag, ftaxo, taxotag,"_",
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
  return(list("allsp"=tbl.list))
  }

#################################################
## Tableau Indices de Biodiversité par station ##
#################################################

inv.biodiv <- function(AS="pas de filtre",qunit="St",wC="all",save=FALSE) {
    departFunk() # message de depart
    on.exit(EM())

    mf <- c("Campagne","St","T") # calcul des valeurs par transect
    # ôte les abondances = 0
    dbn <- dbio[dbio$N > 0,]

    # Filtre sur les stations
    wf <- paste("T",AS,"inv",sep="_") # colonne du filtre
    dbn <- filtreTable(dbn, wf)

  # Filtre sur les espèces au besoin
  dbn <- filtreTaxo(dbn, action=taxoF.incl, taxtype=taxoF.utaxo, taxnom=taxoF.nom)

  # Indice de Shannon #
  # nombre d'individus observés sur le transect
  N.St <- aggregate(dbn$N, as.list(dbn[,mf]),sum)
  names(N.St) <- c(mf,"N.Tot")
  # nombre d'individus *par espèce* observés sur le transect
  N.taxon <- aggregate(dbn$N, as.list(dbn[,c(mf,"G_Sp")]),sum)
  names(N.taxon) <- c(mf,"G_Sp","N.Tax")

  # calcul final de l'indice de Shannon
  shannon.df <- merge(N.taxon,N.St,by=mf)
  shannon.df$Nratio <- shannon.df$N.Tax/shannon.df$N.Tot
  shannon.df$logNratio <- log(shannon.df$N.Tax/shannon.df$N.Tot)
  shannon.df$NratioMult <- shannon.df$Nratio * shannon.df$logNratio
  shannon.final <- aggregate(list("H"=shannon.df$NratioMult), as.list(shannon.df[,mf]),sum)
  shannon.final$H <- -shannon.final$H # Final value Shannon Index by Campagne/Station

  # nombre d'espèces uniques observées sur le transect
  bio.St.i <- aggregate(dbn$G_Sp, as.list(dbn[,mf]), unique)
  bio.St.i$bd <- sapply(bio.St.i$x, length)
  bio.St <- bio.St.i[,c(mf,"bd")] # bio.St$bd = richesse spécifique sur le transect

  # Assemblage des tableaux:
  all.bio <- merge(bio.St,N.St,by=mf)
  all.bio <- merge(all.bio, shannon.final, by=mf)

  # Calcul de l'index *J*
  all.bio$J <- all.bio$H/log(all.bio$bd)

  # Calcul de l'index *d*
  all.bio$d <- (all.bio$bd-1)/log(all.bio$N.Tot)

  # lorsque diversité = 0, J+D -> infini
  val.inf <- sum(!is.finite(all.bio$J)) + sum(!is.finite(all.bio$d))
  all.bio[!is.finite(all.bio$J),"J"] <- NA
  all.bio[!is.finite(all.bio$d),"d"] <- NA
  if(val.inf != sum(is.na(all.bio[,c("H","J","d")]))) {
    print("attention certaines valeurs NA ne sont pas causees par N=1 ou S=1 sur un transect")}

  if(qunit=="St") {
  # Moyenne et écart type des indices par station/Campagne:
  all.bio.mean <- aggregate(list("Moy"=all.bio[,c("H","J","d")]),
                            as.list(all.bio[,c("Campagne","St")]), mean, na.rm=TRUE)
  all.bio.sd <- aggregate(list("ET"=all.bio[,c("H","J","d")]),
                          as.list(all.bio[,c("Campagne","St")]), sd, na.rm=TRUE)

  # Rajout des colonnes additionelles
  all.bioFN <- merge(all.bio.mean, all.bio.sd, by=c("Campagne","St"))
  } else { all.bioFN <- all.bio }

  all.bioFN <- merge(unique(info.transect[,c("St","Geomorpho")]),all.bioFN)
 # all.bioFN <- merge(info.transect.INV, all.bioFN, by=c("Campagne","St","Geomorpho"))
print("hop!")
  # Rajout de la richesse taxonomique:
  rsdf <- inv.RichSpecifique(qunit="T",fspat="St", AS="Absent")
  ff <- c("Campagne","Geomorpho","St","T")
  ff <- ff[1:which(ff==qunit)]

  all.bioFN <- merge(all.bioFN, rsdf, by=ff)

  # Ordonner le tableau par campagne
#  if(qunit == "T") {
#  all.bioFN <- with(all.bioFN, all.bioFN[order(Annee, Mois, Geomorpho, N_Impact, St, T),])
#  } else {
#    all.bioFN <- with(all.bioFN, all.bioFN[order(Annee, Mois, Geomorpho, N_Impact, St),])}

  if(save) {
      ftag <- ifelse(AS=="pas de filtre", "_Filtre-Absent_", paste("_Filtre-",AS,"_",sep=""))
      taxotag <- taxotagFunk()
      write.csv(all.bioFN,file=paste(tabl.dir,"Inv_IndexBiodivPar",
                          qunit,ftag,taxotag,
                        Sys.Date(),".csv",sep=""),row.names=FALSE) }

  return(all.bioFN)
}

##########################################################
##########################################################
##########################################################
inv.RichSpecifique <- function(AS="A", qunit="St", fspat=NA) {

    departFunk() # message de depart
    on.exit(EM())

    ### 1. #############################
    ### Appliquer filtres ##############
    # on commence par ôter les observations N=0
    ta.raw <- dbio[dbio$N > 0,]

    wf <- paste("T",AS,"inv",sep="_") # colonne du filtre
    ta.raw <- filtreTable(ta.raw, wf)

    # Filtre les espèces au besoin
    ta.raw <- filtreTaxo(ta.raw, action=taxoF.incl, taxtype=taxoF.utaxo, taxnom=taxoF.nom)

    # Définition des facteurs spatiaux d'aggrégation
    # Ordonnée de grand à petit
    ff <- c("Campagne","Geomorpho","St","T")
    # on rajoute facteur explicatif 'fspat' (ou qunit) entre Geomorpho et St au besoin
    if(!is.na(fspat) | (qunit %in% fact.expl)) {
      ff <- na.omit(unique(c("Campagne","Geomorpho",fspat,qunit,"St","T"))) }
    ff <- ff[1:which(ff==qunit)] # élimine les unités plus petites que qunit

    ### 2. ######################################################################
    ### Nombre d'espèces par Campagne/GéomorphologieOuSite/Groupe_Taxonomique ###
    e1 <- new.env() # créer environnment pour assembler les tableaux au fur et à mesure

    rsfunk <- function(grtax) {

    # Richesse spécifique totale par aggrégation spatiale:
    # Calcul du nombre ~total~ d'espèces par unité spatiale spécifiée sous 'qunit'
    count <- function(x) length(unique(x)) # cette fonction compte le nombre d'espèces uniques
    tb.1 <- aggregate(list(RS=ta.raw[,grtax]),
                      as.list(ta.raw[,ff]), count)

    # **Moyenne** de la richesse spécifique par Géomorphologie
    if(!is.na(fspat)) {
      ff <- ff[1:which(ff==fspat)] # on ôte les unités spatiales plus petites que fspat
      # moyenne
      tb.2 <- aggregate(list("Moy.RS"=tb.1$RS),
                        as.list(tb.1[,ff]), mean)
      # écart type
      tb.2.sd <- aggregate(list("ET.RS"=tb.1$RS),
                           as.list(tb.1[,ff]), sd)
      tb.all <- merge(tb.2, tb.2.sd, by=ff)
      # on nomme les colonnes
      names(tb.all) <- c(ff, paste(c("Moy.RS","ET.RS"),grtax,sep="."))

    }else{ # si on ne fait pas d'aggrégations statistiques, on garde RS seulement
      tb.all <- tb.1
      names(tb.all) <- c(ff, paste("RS",grtax,sep="."))
    }

    # ... et on rajoute au tableau principal de RS
    if(grtax != "G_Sp") e1$rs.all <- merge(e1$rs.all, tb.all, by=ff)

    return(tb.all)
    }

  e1$rs.all <- rsfunk("G_Sp")
  dmm <- lapply(c("Genre","Famille","S_Groupe","Groupe"), rsfunk)

  invisible(e1$rs.all)
}
