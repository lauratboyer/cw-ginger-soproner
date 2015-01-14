
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

      # Rajout des colonnes additionelles
      all.bioFN <- merge(all.bio.mean, all.bio.sd)
      all.bioFN[,stat.cols] <- round(all.bioFN[,stat.cols],5)

      } else { all.bioFN <- all.bio }

    # Rajouter colonne géomorpho:
#    all.bioFN <- merge(unique(info.transect[,c("St","Geomorpho")]),all.bioFN)

#  all.bioFN <- merge(info.transect.INV, all.bioFN, by=c("Campagne","St","Geomorpho"))

  # Rajout de la richesse taxonomique:
  rsdf <- inv.RichSpecifique(qunit="T",fspat=qunit,AS="Absent")
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
