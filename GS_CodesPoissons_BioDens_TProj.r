# Ginger/Soproner
# Code pour calcul des densités et indices de biodiversité pour les poissons
# Time-stamp: <2015-01-13 15:07:32 Laura>
setwd(dossier.R)

# Fichiers des données de comptage = dpoissons, produit dans prep.analyse()
########################################
############# Tableaux #################

#########################################################
# reformattage du tableau des données de comptage initial
poissons.tableau.brut <- function(save=FALSE, abond.nulles=FALSE) {

  ################################
  departFunk() # message de depart
  on.exit(EM())
  ################################

  # Calcul de densité et de biomasse par observation/espèce/transect/campagne
  # Formule de distance sampling ajustée pour les observations individuelles
  # D.obs = N/(2xdm.obsxL); dm.obs= 0.5*(d1+d2) + 0.5
  # Bio.obs = N*a*T^b/(2xdm.obsxL)
  fc <- c("Campagne","An","Mois","St","T","Code_SP")
  ds.calc <- dpoissons[,c(fc,"Long.Transct","N","L","D1","D2","a","b")]
  start.timer()

  # Densité/biomasse par *observation*
  # D1 et D2 sont les distances entre la ligne du transect et l'individu observé
  # Lorsqu'il y a un banc de poissons D1 est la distance de l'individu
  # ... le plus proche, D2 celle de l'indiv le plus loin
  # Sinon D1 = D2
  # Donc 0.5*(D1 + D2) est la distance moyenne, à laquelle on rajoute 0.5 (facteur de correction?)
  # L'aire du transect pour l'espèce i est la largeur x la longueur; largeur = 2 x distance moyenne
  # ... la densité est donc N / aire
  # ... et la biomasse: biomasse totale / aire
  ds.calc$dm.obs <-  0.5*(ds.calc$D1+ds.calc$D2)+0.5 # distance moyenne observation
  ds.calc$dens.obs <- ds.calc$N/(2*ds.calc$dm.obs*ds.calc$Long.Transct) # densité.observation
  ds.calc$bio.obs <- ds.calc$N*ds.calc$a*(ds.calc$L^ds.calc$b)/(2*ds.calc$dm.obs*ds.calc$Long.Transct)

  # ** Assemblage ** Tableau données brutes **
  # Densité/biomasse par observation/espèce/transect/année/campagne
  tb.1 <- ds.calc[,c("Campagne","An","Mois","St","T","Code_SP","N",
                     "dens.obs","bio.obs","D1","D2","L","a","b")]

  # Joint à:
  # Charactéristiques du transect:
  # année/campagne/transect/plongeur/zone impact/...
  # ... géomorphologie/pêche cpue/pêche effort/taux de sédimentation/...
  # Variables non-disponible: cpus, effort_ha, profondeur moyenne (Z)
  # (ce merge est le plus long)
  tb.1 <- merge(tb.1,
                info.transect[,c("An","Mois","St","T","Geomorpho","N_Impact")],
                by=c("An","Mois","St","T"))

  # + Charactéristiques de l'espèce: code_espèce/famille/genre/espèce/...
  # ... Etat commercial/groupe trophique
  tb.2 <- merge(tb.1, bioeco.all[,c("Code_SP","Famille","Genre","Espece","G_Sp",
                                    "Peche","Cible","GTlabel","moblabel")],
                by="Code_SP",sort=FALSE)


  tb.2 <- tb.2[order(tb.2$Campagne,tb.2$St, tb.2$Code_SP),] # re-ordonner

  # selectionner colonnes a conserver
  tb.2 <- tb.2[,c("An","Campagne","Geomorpho","N_Impact","St","T",
                  "Code_SP","Famille","Genre","Espece","G_Sp",
                  "N","L","D1","D2",
                  "Peche","Cible","a","b","dens.obs","bio.obs",
                  "GTlabel","moblabel")]

  # renommer colonnes
  names(tb.2) <- c("An","Campagne","Geomorpho","Zone.Impact","St","T",
                  "Code_SP","Famille","Genre","Espece","G_Sp","N","L","D1","D2",
                  "Peche","Cible","Coeff.a","Coeff.b","dens.obs","bio.obs",
                   "Groupe.Trophique","Groupe.Mobil")

  if(!abond.nulles) {
    message("On conserve seulement les abondances positives pour diminuer la taille du fichier")
    tb.2 <- tb.2[tb.2$N>0,]
  }
  # garde seulement observations N>0 pour réduire la taille du fichier .csv

  # Filtrer les espèces au besoin
  dbn <- filtreTaxo(tb.2, action=taxoF.incl, taxtype=taxoF.utaxo, taxnom=taxoF.nom)

  if(save) {
    write.csv(dbn,file=paste(tabl.dir, "GS_Poissons_TableauDonneesBrutes_",Sys.Date(),".csv",sep=""),
  row.names=FALSE) }

  invisible(dbn)
}


#########################################################
BioDens.sp.poissons <- function(fspat="St", incl.zero=FALSE,
                                par.espece=TRUE, variable.esp) { # formerly BD.by.sp()

  ################################
  departFunk() # message de depart
  on.exit(EM())
  ################################

  # Calcul densite/biomasse par espece/transect/campagne
  # unité de base pour le calcul est "St" par défaut
  ff <- unique(c("Campagne",fspat,"St","T","Code_SP")) # variables d'aggrégation grande -> petite
  if(!missing(variable.esp)) {
    ff <- c(ff[-length(ff)], variable.esp, "Code_SP")
    if(!(variable.esp %in% c("GTlabel","moblabel","Peche","Cible"))) {
      stop("La valeur de variable.esp doit être: GTlabel, moblable, Peche ou Cible") }
    }

  ds.calc <- dpoissons # transfert à l'objet ds.calc
  if(!incl.zero) ds.calc <- ds.calc[ds.calc$N>0,]
  LT <- ds.calc$Long.Transct[1]
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

  # Calcul des biomasses/densités moyennes par unité spatiale 'fspat'
  # (par St par défaut)
  tot.mean.sd <- function(x,...) round(c(Tot=sum(x,...), Moy=mean(x,...), ET=sd(x,...)),4)
  ff <- ff[1:which(ff==fspat[length(fspat)])]
  if(!missing(variable.esp)) ff <- c(ff,variable.esp)
  if(par.espece) ff <- c(ff,"Code_SP")
  pois.metr <- c("dens","N","biomasse","taille.moy")
  t1.moy <- aggregate(ds.df[,pois.metr],
                      as.list(ds.df[,ff]), tot.mean.sd, na.rm=TRUE)
  # ... reconversion à data.frame
  t2 <- do.call("data.frame", t1.moy[,pois.metr])
  t1.moy <- data.frame(t1.moy[,1:length(pois.metr)], t2)

  t1.moy <- t1.moy[,!(names(t1.moy)=="taille.moy.Tot")] # ôte la somme de la taille moyenne
  # Richesse spécifique sur l'unité spatiale par défaut (St)
  t1.moy$RS <- aggregate(ds.df$Code_SP, as.list(ds.df[,ff]), count)$x

  # ordonner par campagne/sites
#  ds.df <- ds.df[order(ds.df$St, ds.df$Campagne),]

  # rajouter colonnes infos additionelles
#  facteurs <- c("An","St","Zone.Impact","Geomorpho","Code_SP",
#                "Peche","Groupe.Trophique","Groupe.Mobil","Cible")
#  ttble <- unique(ds.calc[,c("Campagne",facteurs)])
#  BDtable.i <- merge(ds.df,ttble,c("Campagne","St","Code_SP"))
#  BDtable <- merge(BDtable.i,bioeco.all[,c("Code_SP","Famille","Genre","G_Sp","fmlabel")],
#                   by="Code_SP")

 # return(BDtable)

  return(t1.moy)
}
