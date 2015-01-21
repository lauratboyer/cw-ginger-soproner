#info type de coraux
# Nouveau format LIT
# Corail_All -> General
# Corail_Forme -> Forme
# Corail_Acro -> Acroporidae
# Corail_Sensi -> Sensibilite
# CODE_DET -> All
# Misc: ôté 'Coraux mous' de cat2keep$Acroporidae car absent du nouvel index.LIT
    type.corail <- c("General","Acroporidae",
                     "Forme","Sensibilite","All")
    cat2keep <- list(c("Coraux","Corail mort","Coraux mous","Algues",
                       "Abiotique","Autre faune"),
                     c("Acroporidae","Non-Acroporidae"),
                     c("Corail branchu","Corail tabulaire",
                       "Corail massif","Corail encroutant",
                       "Corail foliaire","Corail sub-Massif","Corail digite"),
                     c("Corail sensible 1","Corail sensible 3"),
                     c("Macro-Algues","Assemblage d'Algues"))

names(cat2keep) <- type.corail
print("on ne tourne pas catégorie Sensibilite pour l'instant")
type.corail <- type.corail[!(type.corail == "Sensibilite")]
cat2keep <- cat2keep[type.corail]

# formerly TB.lit
LIT.tableau.brut <- function(save=FALSE,AS="X") {

    departFunk() # message de depart
    on.exit(EM())

    # Creer tableau donnees brutes pour analyses subsequentes
    # Appliquer filtres (ref: LIT.doc)
    DL <- dLIT
    LIT.transect.info <- unique(DL[,c("Campagne","St","T")])

    # Filtrer par campagne
    if(AS %in% c("A","S")) DL <- DL[grepl(AS, DL$Campagne),]

    # 2. Rajouter infos sur les coraux
    DL <- data.frame(index.LIT[DL$Code_LIT, type.corail], DL)

    # 3. Calculer couverture moyenne par Campagne/St/Transect/TypeDeCoraux
    ccfunk <- function(wc) {

      # défini les niveaux de 'wc' existants dans la base
      nivfact <- unique(DL[,wc])
      nivfact <- nivfact[nivfact %in% cat2keep[[wc]]]

        dl.i <- aggregate(list("PC"=DL$X.),
                   as.list(DL[,c("Campagne","St","T", wc)]),sum)

        dl.ii <- reshape(dl.i, timevar=wc,
                         idvar=c("Campagne","St","T"),direction="wide")
      dl.ii[is.na(dl.ii)] <- 0


        names(dl.ii) <- gsub("PC.","",names(dl.ii))

        dl.ii <- dl.ii[,c("Campagne","St", "T",
                          nivfact)]
      return(list(niv=nivfact, dat=dl.ii))
        }

    # Combiner les tableaux en 1
    dd <- lapply(type.corail, ccfunk) # boucle sur tous les types de coraux
    dd.i <- merge(dd[[1]]$dat, dd[[2]]$dat, by=c("Campagne","St", "T"))
    dd.i <- merge(dd.i, dd[[3]]$dat, by=c("Campagne","St", "T"))
    dd.i <- merge(dd.i, dd[[4]]$dat, by=c("Campagne","St", "T"))
    #dd.i <- merge(dd.i, dd[[5]], by=c("Campagne","St", "T")) # pour corail sensible

    # 4. Rajouter colonnes infos additionelles
    # Geomorphologies
    dd.i2 <- merge(unique(info.transect[,c("St",facteurs.spatio)]),dd.i,by="St")
    dd.i2 <- merge(LIT.transect.info, dd.i2, by=c("Campagne","St","T"))

    # 5. Reordonner selon instructions
    niv.all <- unlist(lapply(dd, function(x) x$niv)) # niveaux existants dans la base

    dd.i2 <- dd.i2[,c("Campagne","St","T", facteurs.spatio, niv.all)]

    # 6. Option d'appliquer un filtre sur le tableau final
    # AS a la valeur "A" ou "S"
      if(AS %in% c("A","S")) {
          wf <- paste("T",AS,"LIT",sep="_") # colonne du filtre
          dd.i2 <- filtreTable(dd.i2, wf) }

    # nametag pour filtre au besoin
    ftag <- ifelse(AS %in% c("A","S"),sprintf("%s_",AS),"")

    if(save) write.csv(dd.i2, file=paste(tabl.dir,"GS_LIT_TableauBrut_",ftag,
                             Sys.Date(),".csv",sep=""),
                             row.names=FALSE)

    attr(dd.i2, "AS") <- AS
    attr(dd.i2, "Projet") <- projet
    invisible(dd.i2) }

#################################################################


# Tableau moyenne/SE suivant categories dans "S_Corail_All",
# par geomorphologie pour une année spécifique
# formerly LIT.ts1()
LIT.couvrt.gnrl <- function(yy=2013, fspat="Geomorpho", ff="General",
                            AS="A", save=FALSE) {
  # set AS to "A" or "S" based on campagne type

    ################################
    departFunk() # message de depart
    on.exit(EM())
    ################################

    # vérifier si LIT.brut existe déjà pour la bonne campagne, sinon retourner
    if(!exists("LIT.brut")) {
      LIT.brut <<- LIT.tableau.brut(AS=AS)
      }else{
        if((attr(LIT.brut,"AS")!=AS)|(attr(LIT.brut, "Projet")!=projet))
          LIT.brut <<- LIT.tableau.brut(AS=AS) }

    # sélectionner l'année:
    dl <- LIT.brut[grepl(yy, LIT.brut$Campagne),]
    if(nrow(dl)==0) stop("\n\n---> Pas de données restantes pour l'année/filtre-campagne spécifiés")

    # (filtre par campagne appliqué dans LIT.tableau.brut())

    # Groupes de coraux
    gC <- coraux.fig[[ff]]
    gC <- gC[gC %in% names(dl)]

    # Moyenne par station
    t1 <- aggregate(list(dl[,gC]),as.list(dl[,c(fspat,"St")]),mean,na.rm=TRUE)

    # Moyenne/SE par type corail par facteur 'fspat'
    if(length(fspat)==1) { fc.list <- list(t1[,fspat])
                           names(fc.list) <- fspat
                       }else{ fc.list <- as.list(t1[,fspat]) }

    # redéfinir fonction pour arrondir
    mean.rnd <- function(x, ...) round(mean(x,na.rm=TRUE,...),4)
    stand.err.rnd <- function(x, ...) round(stand.err(x,...),4)
    dmean <- aggregate(list("Moy"=t1[,make.names(gC)]),fc.list,mean.rnd)
    dse <- aggregate(list("SE"=t1[,make.names(gC)]),fc.list,stand.err.rnd)

    # Assemblage du tableau
    d.all <- merge(dmean,dse)
    d.all <- data.frame(Annee=yy, d.all)

  if(save) write.csv(d.all,
                     file=paste(tabl.dir,"GS_LIT_Tableau1A_",
                       paste(fspat, collapse="-"), "_", ff,"_",AS,"_",
                       yy,"_",Sys.Date(),".csv",sep=""),row.names=FALSE)

  invisible(d.all)
}
