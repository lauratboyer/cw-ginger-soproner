# Ginger/Soproner
# Code pour calcul des densités et indices de biodiversité pour les poissons
# Time-stamp: <2015-01-23 07:25:04 Laura>

# Fichiers des données de comptage = dpoissons, produit dans prep.analyse()
########################################################
########################################################
# Densité
POIS.dens.gnrl <- function(fspat=fspat.defaut, ftemp=ftempo.defaut,
                           par.transect=FALSE,
                           agtaxo="G_Sp", wZeroT=FALSE) {
# formerly BD.by.sp()
print("changer nom argument zeros")
  ################################
  departFunk() # message de depart
  on.exit(EM())
  ################################

  # Calcul densité/biomasse par espèce/transect/campagne
  # unité de base pour le calcul est "St" par défaut
  ff <- unique(c(ftemp,"Campagne",fspat,"St","T", agtaxo)) # variables d'aggrégation grande -> petite

  ds.calc <- dpoissons # transfert à l'objet ds.calc
  if(!wZeroT) ds.calc <- ds.calc[ds.calc$N>0,]

  LT <- ds.calc$Long.Transct[1]
print(sum(is.na(ds.calc$Long.Transct)))
  if(diff(range(ds.calc$Long.Transct))!=0) print("Attention plus d'une longueur de transect par projet")

  # calcul de la densité/biomasse sur unité spatiale aggrégée et par espèce
  # on dérive la distance moyenne pondérée
  # ... c-a-d la distance moyenne de chaque observation pondérée
  # par le nombre d'individus observés, i.e. sum(dm*N_i)/sum(N_i)
  ds.calc$dm.int <- ds.calc$N * (0.5*(ds.calc$D1+ds.calc$D2)+0.5)
  ds.df <- aggregate(ds.calc[,c("dm.int","N")],as.list(ds.calc[,ff]),sum)
  ds.df$dm <- ds.df$dm.int/ds.df$N # valeur finale, moyenne pondérée
  ds.df$dens <- ds.df$N/(2*ds.df$dm*LT) # Densité

  # Calcul de la biomasse et de la taille moyenne
  ds.calc$bio.int <- ds.calc$N * ds.calc$a * ds.calc$L^ds.calc$b
  ds.calc$TM <- ds.calc$N*ds.calc$L # numérateur taille moyenne
  num.bio <- aggregate(ds.calc[,c("bio.int","TM")], as.list(ds.calc[,ff]), sum)
  ds.df$biomasse <- num.bio$bio.int/(2*ds.df$dm*LT) # Biomasse

  # Taille moyenne pondérée sur le nombre d'individus
  ds.df$taille.moy <- num.bio$TM/ds.df$N

  # Biomasse/Dens -> zero si allN=0 mais laisser la taille.moy a NA
  # pour que les moyennes de tailles soient calculees sur les individus
  # observés seulements
  ds.df[ds.df$N==0,c("dens","biomasse")] <- 0
  message("Refaire statistiques par station pour projets avec num transect > 1")

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

  # Richesse spécifique sur l'unité spatiale par défaut (St)
  ff.RS <- ff[!(ff == agtaxo)]
  t1.RS <- aggregate(list(RS.aggr.taxo=ds.df[,agtaxo]), as.list(ds.df[,ff.RS]), count)
  BDtable <- merge(ds.df, t1.RS)

  invisible(BDtable)
}
