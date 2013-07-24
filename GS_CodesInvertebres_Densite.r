# Ginger/Soproner: Produits/Analyses invertebres
## ** Tableaux/graphiques de densite ** ##

# Time-stamp: <2013-07-24 12:28:12 Laura>

setwd(dossier.R)
fig.dir <- paste(dossier.R,"//Graphiques//",sep='')
tabl.dir <- paste(dossier.R,"//Tableaux//",sep='')

 # Extraction tableaux des bases de donnees
if(!exists("data.read")) source("GS_ExtractionDonnees.r")

print("TO DO: inv.dens.tbl should be the same as inv.dens.TS - not the case right now")
print("TO DO: make sure a single month is entered in DB")

#########################################################################
# Calcul de la densite moyenne des especes
# observees par transect/station/campagne

inv.dens.tbl.parT <- function(AS="A", grtax="G_Sp", wZeroAll=FALSE, save=FALSE) {

  ### 1. #############################
  ### Appliquer filtre ###############
  # Campagnes:
  wf <- paste("T",AS,"inv",sep="_") # colonne du filtre

  ta.raw <- filtreTable(dbio, wf)
  # Espèces:
  ta.raw <- filtreTaxo(ta.raw, action=taxoF.incl, taxtype=taxoF.utaxo, taxnom=taxoF.nom)

  # Densite par transect: Nombre d'individus/Aire du transect -> N/(50*Ltrans)
  # (confirme par valeurs "D" existantes)
  ### 2. #################################
  ### Densite par transect ###############
  ff <- c("Campagne","St","T",grtax) # facteurs d'aggrégation

  # Densité de chaque groupe taxonomique sur le transect
  ta.raw <- aggregate(list("D"=ta.raw$D), as.list(ta.raw[,ff]), sum)

  ### 3. Rajouter zéros partout si spécifié #######
  #################################################
  if(wZeroAll) {

  norm0 <- expand.grid(unique(ta.raw$St),
                       unique(ta.raw$T), unique(ta.raw[,grtax]))
  names(norm0) <- c("St","T",grtax)
  # Rattacher seulement les combinaisons de Campagne/St qui ont été échantillonnées
    filt.tag <- ifelse(AS %in% c("A","S"), AS, "noFiltre")
  norm0 <- merge(nst.par.gs(AS=filt.tag)$nomsSt[,c("Campagne","St")], norm0, all=TRUE)

  ta.raw2 <- merge(ta.raw,norm0,all=TRUE)

  ta.raw2$D[is.na(ta.raw2$D)] <- 0 # densités NA sont non-existantes, donc 0
  ta.raw <- ta.raw2
}

  # Formattage:
  # rajouter info taxonomique (info plus générale seulement, i.e. Espèce -> ... -> Groupe)
  ta.raw <- merge(index.invSp[,grep(grtax,names(index.invSp)):ncol(index.invSp)],ta.raw)

  # rajouter info sur les sites
  iti <- info.transect.INV[,c("Annee","Mois","Mission","Campagne",
                                "Geomorpho", "N_Impact", "St")]

  tb <- merge(iti, ta.raw, by=c("Campagne","St"))

  # ordonner colonnes:
  tb <- with(tb, tb[order(Annee, Mois, Geomorpho, N_Impact, St, T),])

  ### 4. ###############################
  ### Sauvegarde des tableaux ############
  if(save){
    ftag <- paste("_Filtre_",AS,"_",sep="")
    taxotag <- taxotagFunk()
    all0tag <- c("","AvecZeros_")[wZeroAll + 1]
    write.csv(tb,file=paste(tabl.dir,"Inv_DensiteParTransect_",grtax,ftag,all0tag,taxotag,
                   Sys.Date(),".csv",sep=""),row.names=FALSE)
  }
  return(tb)
}


#########################################################################
# Calcul de la densité moyenne des espèces observées par station/campagne

inv.dens.tbl <- function(AS="A", smpl.unit="St", grtax="G_Sp", save=FALSE, wZeroAll=FALSE) {

  # Densité par transect: Nombre d'individus/Aire du transect -> N/(50*Ltrans)

  ### 1. #################################
  ### Densité par transect ###############
  ff <- c("Campagne","St","T","Groupe","S_Groupe","Famille","Genre","G_Sp")
  ff <- ff[1:(which(ff==grtax))] # sélectionne la partie du vecteur allant jusqu'à grtax

  # Densité de chaque grtax sur le transect
  dens.tb <- aggregate(list("D"=dbio$D), as.list(dbio[,ff]), sum)
  ta.raw <- merge(dens.tb, info.transect[,c("St","Geomorpho")])

  ### 2. #############################
  ### Appliquer filtres ###############
  # Campagnes:
  wf <- paste("T",AS,"inv",sep="_") # colonne du filtre
  ta.raw <- filtreTable(ta.raw, wf)

  # Espèces:
  ta.raw <- filtreTaxo(ta.raw, action=taxoF.incl, taxtype=taxoF.utaxo, taxnom=taxoF.nom)

  if(smpl.unit!="T") {

      ### 3. ######################################################
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
  if(smpl.unit) ord.ff <- c(ord.ff,"T") # ordonner par transect si applicable
  tb.all <- tb.all[do.call(order, tb.all[,ord.ff]),]

  ### 4. ###############################
  ### Sauvegarde des tableaux ############
  if(save) {
    ftag <- paste("_Filtre-",AS,"_",sep="")
    taxotag <- taxotagFunk()
    write.csv(tb.all,file=paste(tabl.dir,"Inv_DensitePar_",smpl.unit,"_",ftag, grtax, taxotag,
                            Sys.Date(),".csv",sep=""),row.names=FALSE)
  }

  ### 5. ###############################
  ### Si spécifié, rajouter abondances nulles sur les sites non-occupés
  if(wZeroAll) {
      dlab <- ifelse(smpl.unit == "T", "dens", "dens.moy") # spécifier type de densité
      mfacts <- c("Campagne","St", grtax)

      if(smpl.unit == "T") {
          idtbl <- expand.grid(unique(tb.all$Campagne), unique(tb.all$St),
                               unique(tb.all[,grtax]), T=c("A","B","C"))
          mfacts <- c(mfacts, "T") # rajouter facteur Transect pour merge zéros
      }else{
          idtbl <- expand.grid(unique(tb.all$Campagne), unique(tb.all$St), unique(tb.all[,grtax]))}

      # nommer colonnes
      names(idtbl)[1:3] <- c("Campagne","St", grtax)

      # fusionner avec tableau-toutes-combinaisons (idtbl) pour rajouter les zéros
      tb.all2 <- merge(tb.all, idtbl, by=mfacts, all=TRUE)
      tb.all2[,dlab][is.na(tb.all2[,dlab])] <- 0 # mettre valeurs NA à zéro
      tb.all <- tb.all2 }

  return(tb.all)
}

#####################################################################################
## Calcul de la densité moyenne des espèces sur 4 types de description taxonomique ##
## (1) groupe (2) sous-groupe; ######################################################
## détails des 10 espèces les plus abondantes sur le (3) groupe; (4) sgroupe ########
#####################################################################################
## Option de rajouter l'effect Impact (aj.impact=TRUE) ##############################

inv.dens.geom <- function(AS="A", aj.impact=FALSE,
                        spttcampagnes=FALSE, save=FALSE) {

  ### 1. ###############################
  ### Définir fonction d'aggrégation ###

  aggr.funk <- function(ff, top10year=NA) {
  # Requiert fonction *inv.dens.tbl*
  ### (i) #############################################################################
  ### Extraire tableau des densités moyennes pour le niveau taxonomique ###############
    if(!is.na(top10year)) {idtax <- "G_Sp"
                          }else{ idtax <- ff[ff %in% c("Groupe","S_Groupe")] }

    # valeur filtre définie dans l'argument de la fonction inv.dens.geom
    ta.raw <- inv.dens.tbl(grtax=idtax, AS=AS, save=FALSE)

  ### (ii) #############################################################################
  ### Calculer densité moyenne sur la géomorphologie pour le niveau taxonomique ########

    # (1) rajouter abondances nulles pour géomorphologie/zone d'impact où
    # l'espèce (grtax) n'est pas observée sur toutes les stationsn
    # identification des espéces uniques observées sur le grspatial
    # par campagne / géomorphologie (impact) / St / grta
    ind1 <- unique(ta.raw[,ff])
    ind2 <- nst.par.gs(AS=AS, impact=aj.impact)$nomsSt
    ind3 <- merge(ind1, ind2, all=TRUE)
    ta.raw.0 <- merge(ta.raw, ind3, all=TRUE)
    ta.raw.0$dens.moy[is.na(ta.raw.0$dens.moy)] <-  0 # mettre les D à zéro si NA après le merge

    idnow <<- idtax
    # (2) maintenant calculer moyenne/SD par groupement spatial
    tb.1 <- aggregate(list("dens.moy"=ta.raw.0$dens.moy),
                      as.list(ta.raw.0[,ff]), mean)
    tb.1.sd <- aggregate(list("dens.sd"=ta.raw.0$dens.moy),
                         as.list(ta.raw.0[,ff]), sd)
    tb.all <- merge(tb.1,tb.1.sd,by=ff)

  ### (iii) ##################################################
  ### Extraire densité des espèces top 10 si spécifié ########
    if(sum(grepl("G_Sp",ff))==1) {
      grnow <- ff[ff %in% c("Groupe","S_Groupe")] # identification du groupement
      sptop10 <- top10sp(AS=AS, wyear=top10year, wff=grnow,
                         impact=aj.impact, spttcampagnes=spttcampagnes)
      tb.all <- merge(sptop10, tb.all, by=ff[ff != "Campagne"]) }

  ### (iv) #####################################
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
  if(aj.impact) ff.list <- lapply(ff.list, function(i) c(i[1:2],"N_Impact",i[3]))

  # Appliquer les facteurs à la fonction d'aggrégation
  tbl.list <- lapply(ff.list, aggr.funk)
  names(tbl.list) <- c("Groupe","S_Groupe") # deux tableaux produits: groupe/sous-gr

  ### 4. ########################################################################
  ### Densité moyenne(SD) pour les 10 espèces les plus abondantes ###############

  # Calcul de la densité moyenne de chaque espèce par morphologie
  # Rajouter G_Sp à la liste des facteurs et aggréger...
  ff.list.sp <- lapply(ff.list, function(i) c(i,"G_Sp"))
  # ... pour les 10 espèces les plus abondantes sur toute la série temporelle
  tbl.list.sp.1 <- lapply(ff.list.sp, function(x) aggr.funk(x, top10year="all"))
  names(tbl.list.sp.1) <- c("Groupe","S_Groupe")
  # ... et pour les 10 espèces les plus abondantes en 2006
  tbl.list.sp.2 <- lapply(ff.list.sp, function(x) aggr.funk(x, top10year="1ere"))
  names(tbl.list.sp.2) <- c("Groupe","S_Groupe")
  # ... et pour les 10 espèces les plus abondantes en 201X (dernière année)
  tbl.list.sp.3 <- lapply(ff.list.sp, function(x) aggr.funk(x, top10year="derniere"))
  names(tbl.list.sp.3) <- c("Groupe","S_Groupe")

  ########################################################
  ########################################################
  return(list("allsp"=tbl.list, "top10all"=tbl.list.sp.1,
              "top102006"=tbl.list.sp.2,"top10201X"=tbl.list.sp.3))
  }

top10sp <- function(AS="A", wyear="all", wff="Groupe", grtax="G_Sp",
                    impact=FALSE, spttcampagnes=FALSE) {

  # wyear set to "1ere" ou "derniere"

  # Cette function sélectionne les 10 espèces les plus abondantes ("wval")
  # à l'intérieur d'un facteur "wff", et par géomorphologie
  # dans le tableau "inv.Ab$AbdStEspece" (abondance totale espèce par station/campagne)
  wdf <- inv.dens.tbl(AS=AS, grtax=grtax)
  wdf <- wdf[wdf$dens.moy > 0,] # ote les espèces avec densité nulle sur une station

  # si spttcampagnes est TRUE, garder seulement les espèces observées sur toutes
  # les campagnes de l'unité spatiale spécifiée
  if(spttcampagnes) {
    ni <- top10sp.recc(AS=AS, grtax=grtax, impact=impact)
    fi <- "Geomorpho"
    if(impact) fi <- c(fi,"N_Impact")
    wdf <- merge(wdf, ni[,c(grtax, fi)], by=c(grtax,fi)) }

  # sélectionner la première ou dernière campagne si spécifié
  if(wyear != "all") {
    all.yr <- unique(wdf[,c("Annee","Mois")])
    all.yr <- all.yr[order(all.yr$Annee, all.yr$Mois),]
    if(wyear=="1ere") amfl <- all.yr[1,]
    if(wyear=="derniere") amfl <- all.yr[nrow(all.yr),]
    # garder seulement les observations sur la campagne spécifiée par "amfl"
    wdf <- wdf[wdf$Annee == amfl$Annee & wdf$Mois == amfl$Mois,]
  }

  # Calculer abondance moyenne par année/géomorpho/facteur/espèce + impact si spécifié
  ff <- c("Campagne","Geomorpho","N_Impact",wff, grtax)
  if(!impact) ff <- ff[ff != "N_Impact"]
  # abondance totale d'une espèce par campagne et géomorphologie (et impact si spécifié)
  tb <- aggregate(list("N"=wdf$dens.moy),as.list(wdf[,ff]),sum)

  # Abondance totale / nombre de sites appartenant au groupe
  ff.num <- c("Geomorpho","N_Impact")[1:ifelse(impact,2,1)]
  nsdf <- nst.par.gs(AS=AS, impact=impact)$numSt

  # ajouter colonne "N.st" - nombre de sites par Campagne & grspatial géo (+ impact si spécifié)
  tb.1 <- merge(tb, nsdf, by=c("Campagne",ff.num))
  tb.1$AbdMoy <- tb.1$N/tb.1$N.u

  # Si toutes les campagnes sont utilisées, recalculer la moyenne sur toute la série temp
  if(wyear == "all") {
    tb.2 <- aggregate(list("AbTot"=tb.1$AbdMoy), as.list(tb.1[,ff[ff!="Campagne"]]), sum)
    nsdf <- nst.par.gs(AS=AS, impact=impact, allyr=TRUE)$numSt
    tb <- merge(tb.2, nsdf, by=ff.num)
    tb$AbdMoy <- tb$AbTot/tb$N.u
  } else {tb <- tb.1}

  # Ordonner les facteurs par ordre croissant
  # et la valeur d'intérêt par ordre décroissant:
  # Les espèces sont ordonnées alphatiquement, donc quand les densités sont égales
  # le rang est assigné a > b > ... > z (pour assurer une consistence dans les choix)
  if(impact){ tb <- tb[order(tb$Geomorpho, tb$N_Impact, tb[,wff], -tb$AbdMoy, tb[,grtax]),]
           }else{ tb <- tb[order(tb$Geomorpho, tb[,wff], -tb$AbdMoy, tb[,grtax]),] }

  # l'ordre des colonnes est défini pour correspondre à Geomorpho/(N_Impact)/wff:
  ff2 <- ff[length(ff):1]; ff2 <- ff2[!(ff2 %in% c(grtax,"Campagne"))]
  #ff2 <- ff[!(ff %in% c(grtax,"Campagne"))]

  rdf <- aggregate(tb[,"AbdMoy"], by=as.list(tb[,ff2]), length)
  tb$rank <- unlist(sapply(rdf$x, seq))

  # Garde seulement les espèces dans les 10 plus abondantes:
  tb <- tb[tb$rank <= 10, c(ff2,grtax,"AbdMoy","rank")]

  return(tb)
}

top10sp.recc <- function(AS="A", grtax="G_Sp", filtre=FALSE, impact=FALSE) {

  # identify species that are observed in all campagnes, on the géomorpho/impact
  wdf <- inv.dens.tbl(AS=AS, grtax=grtax)
  wdf <- wdf[wdf$dens.moy > 0,] # ote les espèces avec densité nulle sur une station

  # filre
  if(filtre) wdf <- filtreTable(wdf, paste("T",AS,"inv",sep="_"))

  d1 <- merge(wdf, info.transect[,c("St","Geomorpho","N_Impact")])

  ff <- c("Geomorpho","N_Impact")
  if(!impact) ff <- ff[ff!="N_Impact"]

  d2 <- unique(d1[,c("Campagne",ff,grtax)])
  if(impact) {d2$SpatialGr <- paste(d2$Geomorpho, d2$N_Impact, sep="_")
            } else {d2$SpatialGr <- d2$Geomorpho}

  d3 <- tapply(d2[,grtax], as.list(d2[,c(grtax,"Campagne","SpatialGr")]), length)
  campf <- dimnames(d3)$Campagne
  geof <- dimnames(d3)$SpatialGr

  subf <- function(wg) {

    # Calcule le nombre de campagnes où chaque espèce a été observée
    s1 <- apply(d3[,,wg], 1, sum, na.rm=TRUE)
    # sélectionne les espèces observées sur toutes les camp
    s1 <- s1[s1 == length(campf)]

    if(length(s1)>0) {
    s2 <- data.frame(wg, names(s1))
    names(s2) <- c("SpatialGr",grtax)
    return(s2)}
}

  d4 <- do.call(rbind,lapply(geof, subf))
  if(impact) {
    d4col <- t(sapply(strsplit(d4$SpatialGr,"_"),function(ss) c(ss[1],ss[2])))
    d4$Geomorpho <- d4col[,1]; d4$N_Impact <- d4col[,2]
  } else { d4$Geomorpho <- d4$SpatialGr }

  d5 <- merge(d4[,names(d4) != "SpatialGr"], index.invSp, sort=FALSE)
  ni <- names(index.invSp)[which(names(index.invSp)==grtax):ncol(index.invSp)]
  d5 <- d5[,c(ff,ni)]
}

inv.graph.TS <- function(AS="A", wtype="allsp", wff="Groupe", top10year="", save=TRUE) {

  # option "wtype" peut-être "allsp" ou "top10" (pour les 10 espèces les plus abondantes)

  # Données:
  ts.data <- inv.dens.geom(AS=AS)

  # Série temporelle, groupe (ou sous-groupe) par géomorphologie
  # ts.data est une liste contenant plusieurs versions des tableaux de séries temporelles:
  # "allsp","top10all","top102006"
  # le tableau approprié est sélectionné à partir des arguments wtype et top10year
  dd <- ts.data[[ifelse(wtype=="allsp", "allsp", paste(wtype,top10year,sep=""))]][[wff]]
  dd$Campagne <- as.year(dd$Campagne)
  dd <- split(dd, dd$Geomorpho)

  # Fonction figure:
  graphfunk <- function(i) {

    rds <- round(dev.size(),2)
    if(rds[1] != 8 | rds[2] != 6.5) quartz(width=8, height=6.5)
    par(family="serif",mai=c(0.8,0.1,0.5,0.1), omi=c(0,0.75,0,0.4))
    layout(cbind(1,1,1,1,1,2,2))

    dfig <- dd[[i]]
    xl <- range(dfig$Campagne)
    yl <- c(0, max(dfig$dens.moy+dfig$dens.sd,na.rm=TRUE))
    lab.y <- c(expression(paste("Densite moyenne (individus/",m^2,")",sep="")))

    dfig <- split(dfig,dfig[,wff])
    leg.lab <- names(dfig)
    lab.y <- c(expression(paste("Densite moyenne (individus/",m^2,")",sep="")))

    plot(1,1,type="n",xlim=xl,ylim=yl,las=1,ylab="",xlab="Annee",
         cex.lab=1.5, cex.axis=1.25)
    mtext(paste(names(dd)[i],":",sep=""), adj=0, cex=1.2)
    mtext(lab.y, side=2, outer=TRUE, line=3.5)

    subfunk <- function(i) {
      dnow <- dfig[[i]]
      lines(dnow$Campagne, dnow$dens.moy, type="b", lty=i, col=i)
      sd.mat <- data.frame(dnow$Campagne,
                           dnow$dens.moy-dnow$dens.sd, dnow$dens.moy+dnow$dens.sd)
      draw.SE(sd.mat, couleur=i, typel=i)
    }

    lapply(1:length(dfig),subfunk)
    plot(1,1,type="n",ann=FALSE, axes=FALSE)
    legend("topleft",leg.lab, bty="n",lty=1:length(dfig), col=1:length(dfig),
          xpd=NA, ncol=1, cex=1.25)

    if(save) {
      ftag <- paste("_Filtre_",AS,"_",sep="")
      cat <- paste(ftag, wtype, top10year, wff, names(dd)[i])
      dev.copy2pdf(file=paste(fig.dir,"InvDensMoyTS_",make.names(cat),".pdf",sep=""))
    }
  }

  dmm <- sapply(1:length(dd), graphfunk)

}

inv.graph.TS.top10 <- function(AS="A", wff="Groupe", top10year="all", save=TRUE) {

  # option "wtype" peut-être "allsp" ou "top10" (pour les 10 espèces les plus abondantes)

  # Données:
 ts.data <- inv.dens.geom(AS=AS)

  # Série temporelle, groupe (ou sous-groupe) par géomorphologie
  # ts.data est une liste contenant plusieurs versions des tableaux de séries temporelles:
  # "allsp","top10all","top102006"
  # le tableau approprié est sélectionné à partir des arguments wtype et top10year
  dd <- ts.data[[paste("top10",top10year,sep="")]][[wff]]
  dd$Campagne <- as.year(dd$Campagne)
  jfact <- unique(dd[,wff])

  dd <- split(dd, dd$Geomorpho)

  # Fonction figure:
  graphfunk <- function(geo, gtype) {

    rds <- round(dev.size(),2)
    if(rds[1] != 8 | rds[2] != 6.5) quartz(width=8, height=6.5)
    layout(cbind(1,1,1,1,1,2,2))
    par(family="serif",mai=c(0.8,0.1,0.5,0.1), omi=c(0,0.75,0,0.4))

    dfig <- dd[[geo]][dd[[geo]][,wff]==gtype,]
    # faire seulement un graph s'il y a des données pour cette géomorpho/groupe
    if(nrow(dfig)>0){

      xl <- c(2006,2011)
    dfig$dens.sd[is.na(dfig$dens.sd)] <- 0
    yl <- c(0, max(dfig$dens.moy+dfig$dens.sd,na.rm=TRUE))

    dfig <- split(dfig,dfig[,"G_Sp"])
    leg.lab <- names(dfig)
    lab.y <- c(expression(paste("Densite moyenne (individus/",m^2,")",sep="")))

    plot(1,1,type="n",xlim=xl,ylim=yl,las=1,xlab="Annee",
         cex.lab=1.5, cex.axis=1.25)
      mtext(paste(names(dd)[geo],"/",gtype,":",sep=""), adj=0, cex=1.2)
      mtext(lab.y, side=2, outer=TRUE, line=3.5)

    subfunk <- function(i) {
      dnow <- dfig[[i]]
      lines(dnow$Campagne, dnow$dens.moy, type="b", lty=i, col=i)
      sd.mat <- data.frame(dnow$Campagne,
                           dnow$dens.moy-dnow$dens.sd, dnow$dens.moy+dnow$dens.sd)
      draw.SE(sd.mat, couleur=i, typel=i)
    }

      lapply(1:length(dfig),subfunk)
      plot(1,1,type="n",ann=FALSE, axes=FALSE)

    legend("topleft",leg.lab, bty="n",lty=1:length(dfig), col=1:length(dfig),
           xpd=NA, ncol=1, cex=1.25)

    if(save) {
      ftag <- paste("_Filtre_",AS,"_",sep="")
      cat <- paste(ftag, "top10",top10year, wff, names(dd)[geo], gtype)
      dev.copy2pdf(file=paste(fig.dir,"InvDensMoyTS_",make.names(cat),".pdf",sep=""))
    }
  }
  }

  dmm <- sapply(1:length(dd), function(x)
                sapply(jfact, function(y) graphfunk(geo=x,gtype=y)))
}





####################################################################################################
####################################################################################################
####################################################################################################

