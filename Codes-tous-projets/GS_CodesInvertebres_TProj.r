## GS_CodesInvertebres_TProj.r
##
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: January 22, 2015

############################################################
## General version of inv.dens.tbl
INV.dens.gnrl <- function(fspat=fspat.defaut, ftemp=ftempo.defaut,
                          agtaxo=agtaxo.defaut, par.transect=FALSE,
                          filt.camp="X",
                          wZeroSt=FALSE, wZeroT=TRUE, save=FALSE,
                          silent=FALSE) {

    if(!silent) {departFunk() # message de depart
    on.exit(EM())}

    ### 1. ########################################################
    ### Appliquer filtres campagnes et taxonomiques ###############
    # Filtre campagnes -- défini dans les arguments de la fonction
    wf <- paste("T",filt.camp,"inv",sep="_") # formatter nom du filtre
    if(!silent) start.timer()
    ta.rawF <- filtreTable(dbio, wf)
    if(!silent) stop.timer()

    # Filtre espèces -- défini dans les options générales
    ta.rawF <- filtreTaxo(ta.rawF, action=taxoF.incl,
                          taxtype=taxoF.utaxo, taxnom=taxoF.nom, silent=silent)
    if(!wZeroT) ta.rawF <- ta.rawF[ta.rawF$D > 0,]

    ## Densité par transect:
    ## Nombre d'individus/Aire du transect -> N/(Largeur*Longueur.Transect)
    ## Calculée lors du formattage initial de dbio

    ### 1. #################################
    ### Densité par transect ###############
    ## Définir les facteurs d'aggrégation à partir des arguments donnés
    ff <- c(ftemp, fspat, agtaxo)

    ## 3a. ########################################
    ### Densité de chaque grtax sur le transect ###
    ff.transect <- unique(c(ff, "Campagne","St", "T"))

    # using this solution for now
    # see also:
    # gb.args <- lapply(ff.transect, as.symbol)
    # dens.par.t <- ta.rawF %>% group_by_(gb.args) %>% summarize(D=sum(D))
    dens.par.t <- ta.rawF %>% s_group_by(ff.transect) %>% summarize(D=sum(D))

    # on continue si l'aggrégation spatiale n'est pas
    # au niveau du transect
    if(!par.transect) {

      # on ôte "T" des facteurs d'aggrégation et on calcule
      # la densité moyenne/ET des transects par stations
      # (note: par défaut les transects avec densité de zéro sont inclus)
      ff.St <- ff.transect[-length(ff.transect)]
      tb.St <- aggr.multi(list(dens.par.t$D,
                        as.list(dens.par.t[,ff.St]), mean.sd))

        if(wZeroSt) {
          ### Rajouter les zéros sur les stations/campagnes où l'espèce
          ### n'a pas été observée sur l'aggrégation spatiale
          fact.vals <- sapply(ff.St, function(fct) unique(tb.St[,fct]))
          all.fact.df <- do.call(expand.grid, fact.vals)
          tb.St.1 <- merge(tb.St, all.fact.df, all.y=TRUE)
          tb.St.1$Moy[is.na(tb.St.1$Moy)] <- 0

          # Oter les combinaisons de Campagne/St non-échantillonées
          # (vu que pas toutes les stations ont été échantl à chaque campagne)
          St.Camp <- unique(dbio[,ff.St[ff.St != agtaxo]])
          tb.St <- merge(tb.St.1, St.Camp)
          t1 <- tb.St
    }

     ### Si l'aggrégation spatiale est plus élevée que la station,
      ### on fait la moyenne des valeurs par station pour fspat
      ### Densité moyenne(SD) by Campagne/St/Groupe taxo ###############
      ### Faire moyenne sur stations par aggrégration spatiale fspat

      if(!("St" %in% fspat)) {
        tb.all <- aggr.multi(list(list("dens.moy"=tb.St$Moy),
                        as.list(tb.St[,ff]), mean.sd))
        } else { tb.all <- tb.St }

      } else {
        tb.all <- dens.par.t # on continue avec le tableau transect seulement
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
    fact.tag <- paste(fact.tag, ifelse(par.transect,"-T",""), sep="")
    taxotag <- taxotagFunk()

    write.csv(tb.all,file=paste(tabl.dir,"Inv_Densite",
                       ftag, agtaxo, "_", taxotag, fact.tag,"_",
                            Sys.Date(),".csv",sep=""),row.names=FALSE)
  }

  invisible(tb.all)
}

#################################################
## Tableau Indices de Biodiversité par station ##
#################################################
INV.biodiv.gnrl <- function(ftemp="Campagne", fspat="St",
                            par.transect=FALSE, unit.base="T",
                            filt.camp="X", save=FALSE) {
    departFunk() # message de depart
    on.exit(EM())

    # définition des facteurs d'aggrégation
    # par défaut pas d'aggrégation, calcul sur Campagne x St x T
    mf <- unique(c(ftemp, "Campagne", fspat, "St", unit.base))
    # ignore les abondances = 0 dans les calculs de biodiv
    dbn <- dbio[dbio$N > 0,]

    # Filtre sur les stations
    wf <- paste("T",filt.camp,"inv",sep="_") # colonne du filtre
    dbn <- filtreTable(dbn, wf)

    # Filtre sur les espèces au besoin
    dbn <- filtreTaxo(dbn, action=taxoF.incl, taxtype=taxoF.utaxo, taxnom=taxoF.nom)

    # Indice de Shannon #
    # nombre d'individus observés sur le transect
    N.St <- aggregate(list(N.Tot=dbn$N), as.list(dbn[,mf]), sum)

    # nombre d'individus *par espèce* observés sur le transect
    N.taxon <- aggregate(list(N.Tax=dbn$N), as.list(dbn[,c(mf,"G_Sp")]), sum)

    # calcul final de l'indice de Shannon
    shannon.df <- merge(N.taxon,N.St,by=mf)
    shannon.df$Nratio <- shannon.df$N.Tax/shannon.df$N.Tot
    shannon.df$logNratio <- log(shannon.df$N.Tax/shannon.df$N.Tot)
    shannon.df$NratioMult <- shannon.df$Nratio * shannon.df$logNratio
    shannon.final <- aggregate(list("H"=shannon.df$NratioMult), as.list(shannon.df[,mf]),sum)
    shannon.final$H <- -shannon.final$H # index de Shannon Index par Campagne/Station

    # nombre d'espèces uniques observées sur le transect
    bio.St <- aggregate(list(bd=dbn$G_Sp), as.list(dbn[,mf]), count)

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

    if(!par.transect) {

      # on calcule les stats sur les facteurs d'aggrégation en gardant
      # les transects comme unité de base
      mf <- c(ftemp, fspat)
      # Moyenne et écart type des indices par station/Campagne:
      all.bioFN <- aggr.multi(list(list("Moy"=all.bio[,c("H","J","d")]),
                            as.list(all.bio[,mf]), mean.sd))
      } else {

        all.bioFN <- all.bio }

    # Rajout de la richesse taxonomique:
    rsdf <- INV.RS.gnrl(ftemp,fspat,par.transect,unit.base,filt.camp=filt.camp)
    all.bioFN <- merge(all.bioFN, rsdf)

    if(save) {

      ftag <- paste("_Filtre-",filt.camp,"_",sep="")
      fact.tag <- paste(fspat,ftemp,sep="-",collapse="-")
      fact.tag <- paste(fact.tag, ifelse(par.transect,"-T",""), sep="")
      taxotag <- taxotagFunk()
      write.csv(all.bioFN,file=paste(tabl.dir,"Inv_IndexBiodiv",
                            ftag,taxotag,fact.tag,"_",
                            Sys.Date(),".csv",sep=""),row.names=FALSE) }

  invisible(all.bioFN)
}

######################################################################
######################################################################

INV.RS.gnrl <- function(ftemp="Campagne", fspat="St",
                            par.transect=FALSE, unit.base="T",
                            filt.camp="X", save=FALSE) {

    departFunk() # message de depart
    on.exit(EM())

    # Échelles taxonomiques sur lesquelles on calcule la RS
    niv.taxo <- c("G_Sp","Genre","Famille","S_Groupe","Groupe")

    ### 1. #############################
    ### Appliquer filtres ##############
    # on commence par ôter les observations N=0
    ta.raw <- dbio[dbio$N > 0,]

    wf <- paste("T",filt.camp,"inv",sep="_") # colonne du filtre
    ta.raw <- filtreTable(ta.raw, wf)

    # Filtre les espèces au besoin
    ta.raw <- filtreTaxo(ta.raw, action=taxoF.incl, taxtype=taxoF.utaxo, taxnom=taxoF.nom)

    # Définition des facteurs d'aggrégation temporels/spatiaux
    ff <- unique(c(ftemp, "Campagne", fspat, "St",  unit.base))

    ### 2. ######################################################################
    ### RS par niveau taxonomique par facteurs sur l'unité de base (St ou T)
    tb.all <- aggregate(list(RS=ta.raw[,niv.taxo]), as.list(ta.raw[,ff]), count)

    ### 3. Statistiques de la richesse spécifique sur aggrégation spatiale ou temporelle
    ### seulement si par.transect = FALSE
    if(!par.transect) {

      # on commence par la richesse totale sur toute l'aggrégation spécifiée,
      # tous unit.base confondus
      ff <- c(ftemp, fspat)

      # Richesse totale sur toute l'aggrégation
      # (donc on utilise ta.raw vu qu'on refait le calcul sur toute
      # l'aggrégation spécifiée)
      tb.tot <- aggregate(list(RS.Tot=ta.raw[,niv.taxo]),
                              as.list(ta.raw[,ff]), count)

      # Statistiques de la RS sur l'aggrégation calculée sur l'unité de base
      colnoms <- paste("RS",niv.taxo,sep=".")
      tb.stat <- aggr.multi(list(list(tb.all[,colnoms]),
                              as.list(tb.all[,ff]), mean.sd))

      tb.all <- merge(tb.tot, tb.stat)
    }

    invisible(tb.all)
}

################################################################
################################################################
INV.stats <- function(form="Geomorpho",
                      agtaxo="Groupe",
                      agtaxo.val="Crustaces", ftemp.val="A_2013",
                      fspat.val,
                      type.model="aov") {

  # aller chercher les données par transect
  if(!exists("dat.stat.inv")) {
  dat.stat.inv <<- INV.dens.gnrl(fspat=c(facteurs.spatio,"St"),
                      ftemp=c(facteurs.tempo, "Campagne"),
                             agtaxo=agtaxo, par.transect=TRUE)
}
print(agtaxo.val)
  dstat.now <- dat.stat.inv[dat.stat.inv[,agtaxo]%in% agtaxo.val,]
  dstat.now <- dstat.now[dstat.now$Campagne %in% ftemp.val,]
print(head(dstat.now))
#  formule <- as.formula(sprintf("log(dens+0.01) ~ %s + %s", paste(ftemp, collapse="+"), paste(fspat, collapse="+")))
 formule <- as.formula(sprintf("log(dens+0.01) ~ %s", form))

  obj <- do.call(type.model, list(formula=formule, data=dstat.now))

  #boxplot(residuals(obj) ~ Campagne, data=dat.stat)

  obj
}

