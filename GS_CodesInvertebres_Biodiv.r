# Ginger/Soproner: Produits/Analyses invertébrés
# ** Indices de biodiversité ** 
# Time-stamp: <2013-01-21 14:45:51 Laura>

setwd(dossier.R)
fig.dir <- paste(dossier.R,"//Graphiques//",sep='')
tabl.dir <- paste(dossier.R,"//Tableaux//",sep='')

 # Extraction tableaux des bases de données
if(!exists("data.read")) source("GS_ExtractionDonnees.r")

#################################################
## Tableau Indices de Biodiversité par station ##
#################################################

inv.biodiv <- function(qunit="St",wC="all",save=FALSE) {
  print("Donnees non-filtrees par campagne")
  
  mf <- c("Campagne","St","T") # calcul des valeurs par transect
  # ôte les abondances = 0
  dbn <- dbio[dbio$N > 0,]

  # Filtre les espèces au besoin
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
  
  all.bioFN <- merge(info.transect[,c("St","Geomorpho")],all.bioFN)
  all.bioFN <- merge(info.transect.INV, all.bioFN, by=c("Campagne","St","Geomorpho"))

  # Rajout de la richesse taxonomique:
  rsdf <- inv.RichSpecifique(aggr=qunit,AS="Absent")
  ff <- c("Campagne","Geomorpho","St","T")
  ff <- ff[1:which(ff==qunit)]

  all.bioFN <- merge(all.bioFN, rsdf, by=ff)

  # Ordonner le tableau par campagne
  if(qunit == "T") {
  all.bioFN <- with(all.bioFN, all.bioFN[order(Annee, Mois, Geomorpho, N_Impact, St, T),])
  } else {
    all.bioFN <- with(all.bioFN, all.bioFN[order(Annee, Mois, Geomorpho, N_Impact, St),])}
  
  if(save) {
    taxotag <- taxotagFunk()
    write.csv(all.bioFN,file=paste(tabl.dir,"Inv_IndexBiodivPar",qunit,"_",taxotag,
                          Sys.Date(),".csv",sep=""),row.names=FALSE) }
  return(all.bioFN)
}

########################################################
## Tableau Indices de Biodiversité par Géomorphologie ##
########################################################

inv.biodiv.geom <- function(AS="A", save=FALSE) {
   # Campagnes "A"nnuelles ou "S"emestrielles

  all.bio <- inv.biodiv() # requiert tableau all.bio (~ 1 sec)

  ### 1. #############################
  ### Appliquer filtre ###############
  wf <- paste("T",AS,"inv",sep="_") # colonne du filtre
  dd.filt <- filtreTable(all.bio, wf)

  ### 2. #################################
  ### Définir fonction d'aggrégation #####
  
  aggr.funk <- function(ff) {
    print("for now ignoring J and d values that -> inf")
    dd.geo <- aggregate(dd.filt[,c("Moy.H","Moy.J","Moy.d")],
                      as.list(dd.filt[,ff]),mean,na.rm=TRUE)
    dd.geo.SE <- aggregate(dd.filt[,c("Moy.H","Moy.J","Moy.d")],
                           as.list(dd.filt[,ff]),sd,na.rm=TRUE)
    names(dd.geo.SE)[grep("Moy",names(dd.geo.SE))] <-
    sub("Moy","ET",names(dd.geo.SE)[grep("Moy",names(dd.geo.SE))])
    dd.geo.both <- merge(dd.geo, dd.geo.SE, by=ff)

    print("some geo/campagne surveyed in two different months")
    IT.tmp <- aggregate(list("Mois"=info.transect.INV$Mois),
                        as.list(info.transect.INV[,c("Annee","Mission","Campagne","Geomorpho")]),
                                last)

    dd.geo.both <- merge(IT.tmp,dd.geo.both,by=ff[1:2])
    dd.geo.both <- dd.geo.both[,c("Annee","Mission","Mois",ff,
                                  "Moy.H","ET.H","Moy.J","ET.J","Moy.d","ET.d")]

    # Rajouter info richesse spécifique
    if(length(ff)==2) {
      rsdf <- inv.RichSpecifique(aggr="geom", AS=AS)
      dd.geo.both <- merge(dd.geo.both, rsdf, by=c("Campagne","Geomorpho"))
      dd.geo.both <- dd.geo.both[order(dd.geo.both$Annee, dd.geo.both$Mois,
                                       dd.geo.both$Geomorpho),]
  } else {
      rsdf <- inv.RichSpecifique(aggr="geom", aj.impact=TRUE, AS=AS)
      dd.geo.both <- merge(dd.geo.both, rsdf, by=c("Campagne","Geomorpho","N_Impact"))
      dd.geo.both <- dd.geo.both[order(dd.geo.both$Annee, dd.geo.both$Mois,
                                       dd.geo.both$Geomorpho,
                                       dd.geo.both$N_Impact),]
      }

    ##### Sauvegarde ##### 

    if(save) {
      # format nom de fichier
      ftag <- paste("_Filtre_",AS,"_",sep="")
      gtag <- gsub("N_","",paste(ff[-1],collapse="-"))
      taxotag <- taxotagFunk()  
      write.csv(dd.geo.both,file=paste(tabl.dir,"Inv_IndexBiodiv_",gtag,ftag,taxotag,
                         Sys.Date(),".csv",sep=""),row.names=FALSE) }
  }

  ### 3. ########################################################
  ### Lancer la fonction par géomorphologie & géom+n.impact #####
  dd <- aggr.funk(c("Campagne","Geomorpho"))
  dd <- aggr.funk(c("Campagne","Geomorpho","N_Impact"))

}

###################################################
## Tableau: Richesse spécifique ###################
###################################################

inv.RichSpecifique <- function(AS="A", aj.impact=FALSE, aggr="St") {

  ### 1. #############################
  ### Appliquer filtres ###############
  ta.raw <- merge(dbio, info.transect[,c("St","Geomorpho","N_Impact")])
  # ôter les observations N=0
  ta.raw <- ta.raw[ta.raw$N > 0,]
  
  wf <- paste("T",AS,"inv",sep="_") # colonne du filtre
  ta.raw <- filtreTable(ta.raw, wf)

  # Filtre les espèces au besoin
  ta.raw <- filtreTaxo(ta.raw, action=taxoF.incl, taxtype=taxoF.utaxo, taxnom=taxoF.nom)
  
  ### 2. ######################################################################
  ### Nombre d'espèces par Campagne/GéomorphologieOuSite/Groupe_Taxonomique ###
  e1 <- new.env() # créer environnment pour assembler les tableaux au fur et à mesure
  
  rsfunk <- function(grtax) { 
    # Richesse spécifique par STATION:
    ff <- c("Campagne","Geomorpho","St")
    if(aj.impact) ff[3:4] <- c("N_Impact","St")
    if(aggr=="T") ff <- c(ff,"T") # Rajouter le transect aux facteurs si spécifié

    tb.1 <- aggregate(list("unique.tax"=ta.raw[,grtax]),
                      as.list(ta.raw[,ff]), unique)
    tb.1$St.RS <- sapply(tb.1$unique.tax,length)
    tb.1 <- tb.1[,-(grep("unique.tax",names(tb.1)))]

    # **Moyenne** de la richesse spécifique par Géomorphologie
    if(aggr == "geom") {
      ff <- ff[ff != "St"]
      tb.2 <- aggregate(list("Moy.RS"=tb.1$St.RS),
                        as.list(tb.1[,ff]), mean)
      tb.2.sd <- aggregate(list("ET.RS"=tb.1$St.RS),
                           as.list(tb.1[,ff]), sd)
      tb.all <- merge(tb.2, tb.2.sd, by=ff)
      names(tb.all) <- c(ff, paste(c("Moy.RS","ET.RS"),grtax,sep="."))
    }else{ tb.all <- tb.1; names(tb.all) <- c(ff, paste("St.RS",grtax,sep=".")) }    
    
    if(grtax != "G_Sp") e1$rs.all <- merge(e1$rs.all, tb.all, by=ff)
    
    return(tb.all)
    }
  
  e1$rs.all <- rsfunk("G_Sp")
  dmm <- lapply(c("Genre","Famille","S_Grp2","Grp2"), rsfunk)

  return(e1$rs.all)  
}

########################################
## Tableau synthèse: Nombre d'espèces ##
########################################

inv.sprich.tbl <- function(AS="A",grtax="Grp2",save=FALSE, filtre=FALSE) {

  ### 1. #############################
  ### Appliquer filtres ###############
  ta.raw <- merge(dbio, info.transect[,c("St","Geomorpho")])
  wf <- paste("T",AS,"inv",sep="_") # colonne du filtre
  ta.raw <- filtreTable(ta.raw, wf)
  
  # Filtre les espèces au besoin
  ta.raw <- filtreTaxo(ta.raw, action=taxoF.incl, taxtype=taxoF.utaxo, taxnom=taxoF.nom)
  
  ### 2. ###########################################################
  ### Nombre d'espèces par Campagne/GéomorphologieOuSite/Groupe_Taxonomique ###
  ff <- c("Campagne","Geomorpho","St")
  
  tb.1 <- aggregate(list("unique.sp"=ta.raw[,"G_Sp"]),
                  as.list(ta.raw[,c("Campagne","Geomorpho",grtax)]), unique)
  tb.1$n.sp <- sapply(tb.1$unique.sp,length)

  # species richness by Campagne/Geomorphologie
  tb.2 <- aggregate(list("unique.sp"=ta.raw$G_Sp),
                  as.list(ta.raw[,c("Campagne","Geomorpho")]), unique)
  tb.2$sp.all <- sapply(tb.2$unique.sp,length)
  tb.all <- merge(tb.1[,c("Campagne","Geomorpho",grtax,"n.sp")],
                tb.2[,c("Campagne","Geomorpho","sp.all")],by=c("Campagne","Geomorpho"))
  tb.all$ratio.sp <- tb.all$n.sp/tb.all$sp.all
  
  if(save) {
    # rajouter colonnes
    tb.all <- merge(info.transect.INV.geo, tb.all, by=c("Campagne","Geomorpho"), sort=FALSE)
    # définir info filtre pour nom de fichier
    ftag <- c("_Filtre_Absent_",paste("_Filtre_",AS,"_",sep=""))[filtre + 1]
    taxotag <- taxotagFunk()
    write.csv(tb.all,file=paste(tabl.dir,"Inv_NumEspeceParGeomorph_",grtax, ftag,taxotag, 
                            Sys.Date(),".csv",sep=""),row.names=FALSE)}
  return(tb.all)
  }

sprich.by.aggrtaxo <- function(AS="A", grtax="Grp2", filtre=TRUE, save=FALSE) {

  ### 1. #############################
  ### Appliquer filtres ###############
  wf <- paste("T",AS,"inv",sep="_") # colonne du filtre
  ta.raw <- filtreTable(dbio, wf)

  # Filtre les espèces au besoin
  ta.raw <- filtreTaxo(ta.raw, action=taxoF.incl, taxtype=taxoF.utaxo, taxnom=taxoF.nom)

  # Créer colonnes de présence
  ta.raw$presence <- ifelse(ta.raw$N == 0, 0, 1)
  ta.raw$aggrTax <- ta.raw[,grtax]
  
  ### 2. ###########################################################
  ### Nombre d'espèces par Campagne/Transect/Groupe_Taxonomique ###
  ff <- c("Campagne","St","T")

  tb.1 <- aggregate(list("unique.sp"=ta.raw$presence),
                  as.list(ta.raw[,c(ff,"aggrTax")]), sum)

  tb.all <- cast(tb.1, Campagne*St*T ~ aggrTax, value="unique.sp", fun=sum)
  
  if(save) {
    # définir info filtre pour nom de fichier
    ftag <- c("_Filtre_Absent_",paste("_Filtre_",AS,"_",sep=""))[filtre + 1]
    taxotag <- taxotagFunk()
    write.csv(tb.all,file=paste(tabl.dir,"Inv_NumEspeceParAggrTaxon",grtax, ftag,taxotag, 
                            Sys.Date(),".csv",sep=""),row.names=FALSE)}
  return(tb.all)
  }
