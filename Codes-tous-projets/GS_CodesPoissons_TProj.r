# Ginger/Soproner
# Code pour calcul des densités et indices de biodiversité pour les poissons
# Time-stamp: <2015-05-04 07:37:37 Laura>

# Fichiers des données de comptage = dpoissons, produit dans prep.analyse()
########################################################
########################################################
# Densité
POIS.dens.gnrl <- function(fspat=fspat.defaut, ftemp=ftempo.defaut,
                           agtaxo="G_Sp",
                           par.transect=FALSE,
                           wZeroT=TRUE,
                           filt.camp="X", silent=FALSE, save=FALSE) {

  ################################
    if(!silent) {
        departFunk() # message de depart
        on.exit(EM()) }
  ################################

  ### 1. ########################################################
  ### Appliquer filtres
### Filtre campagne:
  d1 <- filtreTable(dpoissons, filt.camp) # filtre années seulement pour l'instant

### Filtre espèces -- défini dans les options générales
  ds.calc <- filtreTaxo(d1, action=taxoF.incl,
                        taxtype=taxoF.utaxo, taxnom=taxoF.nom, silent=silent)
  if(!wZeroT) ds.calc <- ds.calc[ds.calc$N>0,] # ôter zéros

  # Calcul densité/biomasse par espèce/transect/campagne
  # unité de base pour le calcul est "St" par défaut
  ff <- unique(c(ftemp,"Campagne",fspat,"St","T", agtaxo)) # variables d'aggrégation grande -> petite

  LT <- ds.calc$Long.Transct[1]

  # calcul de la densité/biomasse sur unité spatiale aggrégée et par espèce
  # on dérive la distance moyenne pondérée
  # ... c-a-d la distance moyenne de chaque observation pondérée
                                        # par le nombre d'individus observés, i.e. sum(dm*N_i)/sum(N_i)

  ds.calc$dm.int <- ds.calc$N * (0.5*(ds.calc$D1+ds.calc$D2)+0.5)
  # Calcul de la biomasse et de la taille moyenne
  ds.calc$bio.int <- ds.calc$N * ds.calc$a * ds.calc$L^ds.calc$b
  ds.calc$TM <- ds.calc$N*ds.calc$L # numérateur taille moyenne

  # ds.df <- aggregate(ds.calc[,c("dm.int","N")],as.list(ds.calc[,ff]),sum)
  # num.bio <- aggregate(ds.calc[,c("bio.int","TM")], as.list(ds.calc[,ff]), sum)
  # replaced two lines above with single dplyr call for speed
  ds.df <- ds.calc %>% s_group_by(ff) %>% summarize(dm.int=sum(dm.int), N=sum(N), bio.int=sum(bio.int), TM=sum(TM)) %>% data.frame
  # boom!
  ds.df$dm <- ds.df$dm.int/ds.df$N # valeur finale, moyenne pondérée
  ds.df$dens <- ds.df$N/(2*ds.df$dm*LT) # Densité
  ds.df$biomasse <- ds.df$bio.int/(2*ds.df$dm*LT) # Biomasse

  # Taille moyenne pondérée sur le nombre d'individus
  ds.df$taille.moy <- ds.df$TM/ds.df$N

  # Biomasse/Dens -> zero si allN=0 mais laisser la taille.moy a NA
  # pour que les moyennes de tailles soient calculees sur les individus
  # observés seulements
  ds.df[ds.df$N==0,c("dens","biomasse")] <- 0

  if(!par.transect) {
  # Calcul des biomasses/densités moyennes par unité spatiale 'fspat'
  # (par St par défaut)
  ff <- c(ftemp, fspat, agtaxo)

  pois.metr <- c("dens","N","biomasse","taille.moy")
  t1.moy <- aggr.multi(list(ds.df[,pois.metr],
                      as.list(ds.df[,ff]), tot.mean.sd))
  t1.moy <- t1.moy[,!(names(t1.moy)=="taille.moy.Tot")] # ôte la somme de la taille moyenne
  ds.df <- t1.moy
}

  # wZeroSt de invertébrés... à adapter pour les poissons
  #      if(wZeroSt) {
          ### Rajouter les zéros sur les stations/campagnes où l'espèce
          ### n'a pas été observée sur l'aggrégation spatiale
  #        fact.vals <- sapply(ff.St, function(fct) unique(tb.St[,fct]))
  #        all.fact.df <- do.call(expand.grid, fact.vals)
  #        tb.St.1 <- merge(tb.St, all.fact.df, all.y=TRUE)
  #        tb.St.1$Moy[is.na(tb.St.1$Moy)] <- 0

         # Oter les combinaisons de Campagne/St non-échantillonées
          # (vu que pas toutes les stations ont été échantl à chaque campagne)
   #       St.Camp <- unique(dbio[,ff.St[ff.St != agtaxo]])
   #       tb.St <- merge(tb.St.1, St.Camp)
   #       t1 <- tb.St
   # }


                                        # Richesse spécifique sur l'unité spatiale par défaut (St)
  ff.RS <- ff[!(ff == agtaxo)]
  t1.RS <- aggregate(list(RS.aggr.taxo=ds.df[,agtaxo]), as.list(ds.df[,ff.RS]), count)
  BDtable <- merge(ds.df, t1.RS)

    if(save) {

      ftag <- paste("_Filtre-",filt.camp,"_",sep="")
      fact.tag <- paste(fspat,ftemp,sep="-",collapse="-")
      fact.tag <- paste(fact.tag, ifelse(par.transect,"-T",""), sep="")
      taxotag <- taxotagFunk()

    write.csv(BDtable,file=paste(tabl.dir,"Poissons_DensBio",
                       ftag, agtaxo, taxotag, "_", fact.tag, "_",
                            Sys.Date(),".csv",sep=""),row.names=FALSE)
  }

  invisible(BDtable)
}
