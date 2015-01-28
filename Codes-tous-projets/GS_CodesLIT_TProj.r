#info type de coraux
# Nouveau format LIT
# Corail_All -> General
# Corail_Forme -> Forme
# Corail_Acro -> Acroporidae
# Corail_Sensi -> Sensibilite
# CODE_DET -> All
# Misc: ôté 'Coraux mous' de cat2keep$Acroporidae car absent du nouvel index.LIT
    type.corail <- c("General","Acroporidae",
                     "Forme","All") #"Sensibilite",
    cat2keep <- list(c("Coraux","Coraline","Corail mort","Coraux mous","Algues",
                       "Abiotique","Autre faune"),
                     c("Acroporidae","Non-Acroporidae"),
                     c("Corail branchu","Corail tabulaire",
                       "Corail massif","Corail encroutant",
                       "Corail foliaire","Corail sub-Massif","Corail digite"),
                     c("Corail sensible 1","Corail sensible 2"),
                     c("Macro-Algues","Assemblage d'Algues"))

names(cat2keep) <- type.corail
cat2keep <- cat2keep[!(names(cat2keep)=="Sensibilite")]

# formerly TB.lit
LIT.tableau.brut <- function(save=FALSE,filt.camp="X",type.db="LIT") {


  departFunk() # message de depart
  on.exit(EM())

  # Creer tableau donnees brutes pour analyses subsequentes
  # Sélectionner tableau de données LIT ou quadrat
    if(type.db=="LIT") { DL <- dLIT
                         smpl <- "T"}
    if(type.db=="Quadrat") { DL <- dQuad
                             smpl <- "Quadrat" }

    # Filtrer par campagne
    # Option d'appliquer un filtre sur le tableau final
    # filt.camp a la valeur "A" ou "S"
      if(filt.camp %in% c("A","S")) {
          DL <- filtreTable(DL, filt.camp) }
    LIT.transect.info <- unique(DL[,c("Campagne","St",smpl)])

    # 2. Rajouter infos sur les coraux
    DL <- data.frame(index.LIT[DL$Code_LIT, type.corail], DL)

    # 3. Calculer couverture moyenne par Campagne/St/Transect/TypeDeCoraux
    ccfunk <- function(wc) {

      # défini les niveaux de 'wc' existants dans la base
      nivfact <- unique(DL[,wc])
      nivfact <- nivfact[nivfact %in% cat2keep[[wc]]]

        dl.i <- aggregate(list("PC"=DL$X.),
                   as.list(DL[,c("Campagne","St",smpl,wc)]),sum)

        dl.ii <- reshape(dl.i, timevar=wc,
                         idvar=c("Campagne","St",smpl),direction="wide")
      dl.ii[is.na(dl.ii)] <- 0


        names(dl.ii) <- gsub("PC.","",names(dl.ii))

        dl.ii <- dl.ii[,c("Campagne","St",smpl,
                          nivfact)]
      return(list(niv=nivfact, dat=dl.ii))
        }

    # Combiner les tableaux en 1
    dd <- lapply(type.corail, ccfunk) # boucle sur tous les types de coraux
    dd.i <- merge(dd[[1]]$dat, dd[[2]]$dat, by=c("Campagne","St", smpl))
    dd.i <- merge(dd.i, dd[[3]]$dat, by=c("Campagne","St", smpl))
    dd.i <- merge(dd.i, dd[[4]]$dat, by=c("Campagne","St", smpl))
#    dd.i <- merge(dd.i, dd[[5]]$dat, by=c("Campagne","St", smpl))

    # 4. Rajouter colonnes infos additionelles
    # Geomorphologies
    dd.i2 <- merge(unique(info.transect[,c("St",facteurs.spatio, facteurs.tempo)]),dd.i)
    dd.i2 <- merge(LIT.transect.info, dd.i2)

    # 5. Reordonner selon instructions
    niv.all <- unlist(lapply(dd, function(x) x$niv)) # niveaux existants dans la base

    dd.i2 <- dd.i2[,c("Campagne","St",smpl, facteurs.spatio, facteurs.tempo, niv.all)]

    # nametag pour filtre au besoin
    ftag <- ifelse(filt.camp %in% c("A","S"),sprintf("%s_",filt.camp),"")

    if(save) {
      write.csv(dd.i2, file=paste(tabl.dir,"GS_",type.db,"_TableauBrut_",ftag,
                             Sys.Date(),".csv",sep=""),
                             row.names=FALSE)
    }
    attr(dd.i2, "AS") <- filt.camp
    attr(dd.i2, "Projet") <- projet
    attr(dd.i2, "type.db") <- type.db
    invisible(dd.i2) }


#################################################################


# Tableau moyenne/SE suivant categories dans "S_Corail_All",
# par geomorphologie pour une année spécifique
# formerly LIT.ts1()
LIT.couvrt.gnrl <- function(ftemp=ftempo.defaut, fspat=fspat.defaut,
                            LIT.cat="General",type.db="LIT",
                            filt.camp="X", save=FALSE) {
  # set AS to "A" or "S" based on campagne type

    ################################
    departFunk() # message de depart
    on.exit(EM())
    ################################

    # vérifier si LIT.brut existe déjà pour la bonne campagne, sinon retourner
    # (filtre par campagne appliqué dans LIT.tableau.brut())
    if(!exists("LIT.brut")) {
      LIT.brut <<- LIT.tableau.brut(filt.camp=filt.camp, type.db=type.db)
      }else{
        if((attr(LIT.brut,"AS")!=filt.camp)|
           (attr(LIT.brut, "Projet")!=projet)|
           (attr(LIT.brut, "type.db")!=type.db))
          LIT.brut <<- LIT.tableau.brut(filt.camp=filt.camp, type.db=type.db) }

    t1 <- LIT.brut # transfert des données à objet t1

    # Creer tableau donnees brutes pour analyses subsequentes
    # Appliquer filtres (ref: LIT.doc)
    if(type.db=="LIT") smpl <- "T"
    if(type.db=="Quadrat") smpl <- "Quadrat"

    # Groupes de coraux
    gC <- coraux.fig[[LIT.cat]]
    gC <- gC[gC %in% names(t1)] # conserver ceux présent dans la DB

    # Moyenne par station
    #t1 <- aggregate(list(dl[,gC]),as.list(dl[,c(fspat,"St")]),mean,na.rm=TRUE)
    message("ôté moyenne par station")

    # Moyenne/ET par type corail par facteur 'fspat'
    fc.list <- as.list(t1[,c(ftemp, fspat)])
    d.all <- aggr.multi(list(t1[,gC],fc.list,mean.sd))

  if(save) write.csv(d.all,
                     file=paste(tabl.dir,type.db,"_MoySD-couvrt_",
                       paste(ftemp, fspat, sep="-",collapse="-"),
                       "_", LIT.cat,"_",filt.camp,"_",
                       Sys.Date(),".csv",sep=""),row.names=FALSE)

  invisible(d.all)
}

# wrap pour les deux fonctions LIT en Quad
Quad.tableau.brut <- function(...) LIT.tableau.brut(..., type.db="Quadrat")
Quad.couvrt.gnrl <- function(...) LIT.couvrt.gnrl(..., type.db="Quadrat")
