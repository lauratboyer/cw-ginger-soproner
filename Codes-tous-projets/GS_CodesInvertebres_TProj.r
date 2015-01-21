############################################################
############################################################
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

############################################################
############################################################

## General version of inv.dens.tbl
INV.dens.gnrl <- function(fspat=fspat.defaut, ftemp=ftempo.defaut,
                          agtaxo=agtaxo.defaut,
                          filt.camp="X", smpl.unit="St",
                          wZeroSt=FALSE, wZeroTransect=TRUE, save=FALSE) {

    departFunk() # message de depart
    on.exit(EM())

    ### 1. ########################################################
    ### Appliquer filtres campagnes et taxonomiques ###############
    # Filtre campagnes -- défini dans les arguments de la fonction
    wf <- paste("T",filt.camp,"inv",sep="_") # formatter nom du filtre
    ta.rawF <- filtreTable(dbio, wf)
    # Filtre espèces -- défini dans les options générales
    ta.rawF <- filtreTaxo(ta.rawF, action=taxoF.incl, taxtype=taxoF.utaxo, taxnom=taxoF.nom)

    if(!wZeroTransect) ta.rawF <- ta.rawF[ta.rawF$D > 0,]

    ## Densité par transect:
    ## Nombre d'individus/Aire du transect -> N/(Largeur*Longueur.Transect)
    ## Calculée lors du formattage initial de dbio

    ### 1. #################################
    ### Densité par transect ###############
    ## Définir les facteurs d'aggrégation à partir des arguments donnés
    ff <- c(ftemp, fspat, agtaxo)

    ## 3a. ########################################
    ### Densité de chaque grtax sur le transect ###
    ff.transect <- unique(c(ff, "St", "T"))
    dens.par.t <- aggregate(list("D"=ta.rawF$D), as.list(ta.rawF[,ff.transect]), sum)

    # on continue si l'aggrégation spatiale n'est pas
    # au niveau du transect
    if(!("T" %in% fspat)) {

      # on ôte "T" des facteurs d'aggrégation et on calcule
      # la densité moyenne/ET des transects par stations
      # (note: par défaut les transects avec densité de zéro sont inclus)
      ff.St <- ff.transect[-length(ff.transect)]
      tb.St <- aggr.multi(list(dens.par.t$D,
                        as.list(dens.par.t[,ff.St]), mean.sd))

      ### Si l'aggrégation spatiale est plus élevée que la station,
      ### on fait la moyenne des valeurs par station pour fspat
      ### Densité moyenne(SD) by Campagne/St/Groupe taxo ###############
      ### Faire moyenne sur stations par aggrégration spatiale fspat
      if(!("St" %in% fspat)) {

        if(wZeroSt) {
          ### Rajouter les zéros sur les stations/campagnes où l'espèce
          ### n'a pas été observée sur l'aggrégation spatiale
          fact.vals <- sapply(ff.St, function(fct) unique(tb.St[,fct]))
          all.fact.df <- do.call(expand.grid, fact.vals)

          tb.St.1 <- merge(tb.St, all.fact.df, all.y=TRUE)
          tb.St.1$Moy[is.na(tb.St.1$Moy)] <- tb.St.1$ET[is.na(tb.St.1$Moy)] <- 0
          # Oter les combinaisons de Campagne/St non-échantillonées
          # (vu que pas toutes les stations ont été échantl à chaque campagne)
          St.Camp <<- unique(dbio[,ff.St[ff.St != agtaxo]])
          tb.St <- merge(tb.St.1, St.Camp)
          t1 <<- tb.St
    }

        tb.all <- aggr.multi(list(list("dens.moy"=tb.St$Moy),
                        as.list(tb.St[,ff]), mean.sd))
        } else { tb.all <- tb.St }

      } else {
        tb.all <- dens.par.t # si smpl.unit == T, on continue avec le tableau transect seulement
        names(tb.all)[names(tb.all)=="D"] <- "dens" # renomme colonne densité, pas de moyenne
        }

    # Ordonner colonnes:
    ord.ff <- c(ftemp, fspat)
    tb.all <- tb.all[do.call(order, tb.all[,ord.ff]),]

    ### 5. ###############################
    ### Sauvegarde des tableaux ##########
    if(save) {
    ftag <- paste("_Filtre-",filt.camp,"_",sep="")
    fact.tag <- paste(fspat,ftemp,sep="-",collapse="-")
    taxotag <- taxotagFunk()

    write.csv(tb.all,file=paste(tabl.dir,"Inv_DensitePar_", smpl.unit,
                       ftag, agtaxo, taxotag,"_",
                            Sys.Date(),".csv",sep=""),row.names=FALSE)
  }

  return(tb.all)
}

#################################################
## Tableau Indices de Biodiversité par station ##
#################################################

inv.biodiv <- function(filt.camp="X", qunit="T", save=FALSE) {
    departFunk() # message de depart
    on.exit(EM())

    # définition des facteurs d'aggrégation
    mf <- c("Campagne","Geomorpho","St","T")
    if(qunit%in%fact.expl) mf <-  c("Campagne","Geomorpho",qunit,"St","T")

    # ôte les abondances = 0
    dbn <- dbio[dbio$N > 0,]

    # Filtre sur les stations
    wf <- paste("T",filt.camp,"inv",sep="_") # colonne du filtre
    dbn <- filtreTable(dbn, wf)

    # Filtre sur les espèces au besoin
    dbn <- filtreTaxo(dbn, action=taxoF.incl, taxtype=taxoF.utaxo, taxnom=taxoF.nom)

    # Indice de Shannon #
    # nombre d'individus observés sur le transect
    N.St <- aggregate(list(N.Tot=dbn$N), as.list(dbn[,mf]),sum)

    # nombre d'individus *par espèce* observés sur le transect
    N.taxon <- aggregate(list(N.Tax=dbn$N), as.list(dbn[,c(mf,"G_Sp")]),sum)

    # calcul final de l'indice de Shannon
    shannon.df <- merge(N.taxon,N.St,by=mf)
    shannon.df$Nratio <- shannon.df$N.Tax/shannon.df$N.Tot
    shannon.df$logNratio <- log(shannon.df$N.Tax/shannon.df$N.Tot)
    shannon.df$NratioMult <- shannon.df$Nratio * shannon.df$logNratio
    shannon.final <- aggregate(list("H"=shannon.df$NratioMult), as.list(shannon.df[,mf]),sum)
    shannon.final$H <- -shannon.final$H # index de Shannon Index par Campagne/Station

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
    print("Attention certaines valeurs NA ne sont pas causées par N=1 ou S=1 sur un transect")}

    if(qunit!="T") {

      mf <- mf[1:which(mf==qunit)] # on ôte le transect des facteurs d'aggrégation pour faire les stats

      # Moyenne et écart type des indices par station/Campagne:
      all.bio.mean <- aggregate(list("Moy"=all.bio[,c("H","J","d")]),
                            as.list(all.bio[,mf]), mean, na.rm=TRUE)
      all.bio.sd <- aggregate(list("ET"=all.bio[,c("H","J","d")]),
                          as.list(all.bio[,mf]), sd, na.rm=TRUE)
      stat.cols <- c("Moy.H","Moy.J","Moy.d","ET.H","ET.J","ET.d")

      # Assembler + arrondir
      all.bioFN <- merge(all.bio.mean, all.bio.sd)
      all.bioFN[,stat.cols] <- round(all.bioFN[,stat.cols],5)

      } else { all.bioFN <- all.bio }

    # Rajouter colonne géomorpho:
#    all.bioFN <- merge(unique(info.transect[,c("St","Geomorpho")]),all.bioFN)
#  all.bioFN <- merge(info.transect.INV, all.bioFN, by=c("Campagne","St","Geomorpho"))

  # Rajout de la richesse taxonomique:
  rsdf <- inv.RichSpecifique(qunit="T",AS="Absent")
  all.bioFN <- merge(all.bioFN, rsdf)

  # Ordonner le tableau par campagne
#  if(qunit == "T") {
#  all.bioFN <- with(all.bioFN, all.bioFN[order(Annee, Mois, Geomorpho, N_Impact, St, T),])
#  } else {
#    all.bioFN <- with(all.bioFN, all.bioFN[order(Annee, Mois, Geomorpho, N_Impact, St),])}

  if(save) {

    ftag <- paste("_Filtre-",filt.camp,"_",sep="")
    taxotag <- taxotagFunk()
      write.csv(all.bioFN,file=paste(tabl.dir,"Inv_IndexBiodivPar",
                          qunit,ftag,taxotag,
                        Sys.Date(),".csv",sep=""),row.names=FALSE) }

  invisible(all.bioFN)
}


INV.biodiv.gnrl <- function(AS="X", fspat="Geomorpho", save=FALSE) {
    # Campagnes "A"nnuelles ou "S"emestrielles

    departFunk() # message de depart
    on.exit(EM())
    all.bio <- inv.biodiv(qunit="St") # requiert tableau all.bio (~ 1 sec)
    all.bio <- merge(all.bio, unique(info.transect[,c("St",fspat)]))

    # Définir facteurs d'aggrégation
    ff <- c("Campagne",fspat)

    ### 1. #############################
    ### Appliquer filtre ###############
  wf <- paste("T",AS,"inv",sep="_") # colonne du filtre
  dd.filt <- filtreTable(all.bio, wf)

  ### 2. #################################
  ### Définir fonction d'aggrégation #####
    print("Pour le moment les valeurs J et d -> inf sont ignorees")
    dd.geo <- aggregate(dd.filt[,c("Moy.H","Moy.J","Moy.d")],
                      as.list(dd.filt[,ff]),mean,na.rm=TRUE)
    dd.geo.SE <- aggregate(dd.filt[,c("Moy.H","Moy.J","Moy.d")],
                           as.list(dd.filt[,ff]),sd,na.rm=TRUE)
    names(dd.geo.SE)[grep("Moy",names(dd.geo.SE))] <-
    sub("Moy","ET",names(dd.geo.SE)[grep("Moy",names(dd.geo.SE))])
    dd.geo.both <- merge(dd.geo, dd.geo.SE, by=ff)

    # Certaines geomorphologies/campagne sont echantillonees sur des
    # mois differents, donc on selecte le mois le plus recent
    IT.tmp <- aggregate(list("Mois"=info.transect.INV$Mois),
                        as.list(info.transect.INV[,c("Annee","Mission","Campagne","Geomorpho")]),
                                lastval)

    # Rajouter info complementaire Mission/Mois/etc
    dd.geo.both <- merge(IT.tmp,dd.geo.both,by=ff[1:2])
    dd.geo.both <- dd.geo.both[,c("Annee","Mission","Mois",ff,
                                  "Moy.H","ET.H","Moy.J","ET.J","Moy.d","ET.d")]

    # Rajouter info richesse spécifique
      rsdf <- inv.RichSpecifique(fspat=ff[-1], AS=AS)
      dd.geo.both <- merge(dd.geo.both, rsdf)
      dd.geo.both <- dd.geo.both[order(dd.geo.both$Annee, dd.geo.both$Mois,
                                       dd.geo.both$Geomorpho),]


    ##### Sauvegarde #####

    if(save) {
      # format nom de fichier
      ftag <- paste("_Filtre-",AS,"_",sep="")
      gtag <- gsub("N_","",paste(ff[-1],collapse="-"))
      taxotag <- taxotagFunk()
      write.csv(dd.geo.both,file=paste(tabl.dir,"Inv_IndexBiodiv_",gtag,ftag,taxotag,
                         Sys.Date(),".csv",sep=""),row.names=FALSE) }

    invisible(dd.geo.both)
}

######################################################################
######################################################################
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
    # Ordonnés de grand à petit
    ff <- c("Campagne","Geomorpho","St","T")
    # on rajoute facteur explicatif 'fspat' (ou qunit) entre Geomorpho et St au besoin
    if(!is.na(fspat[1]) | (qunit %in% fact.expl)) {
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

    # **Moyenne** de la richesse spécifique par facteurs > St
    if(!is.na(fspat[1])) {

      ff <- ff[1:which(ff==last(fspat))] # on ôte les unités spatiales plus petites que fspat
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
