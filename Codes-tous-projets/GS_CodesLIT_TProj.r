#info type de coraux
# Nouveau format LIT
# Corail_All -> General
# Corail_Forme -> Forme
# Corail_Acro -> Acroporidae
# Corail_Sensi -> Sensibilite
# CODE_DET -> All
# Misc: ôté 'Coraux mous' de cat2keep$Acroporidae car absent du nouvel index.LIT
    type.corail <- c("General","Acroporidae",
                     "Forme","Genre","All") #"Sensibilite",
    cat2keep <- list(General=c("Coraux","Coraline","Corail mort","Coraux mous","Algues",
                       "Abiotique","Autre faune", "Coraux blanchis"),
                    Acroporidae=c("Acroporidae","Non-Acroporidae"),
                     Forme=c("Corail branchu","Corail tabulaire",
                       "Corail massif","Corail encroutant",
                       "Corail foliaire","Corail sub-Massif","Corail digite"),
                     Genre=c("Montipora","Acropora","Pavona"),
                     All=c("Macro-Algues","Assemblage d'Algues"))

# Catégories sensibilité (ôtées car non-utilisées):
#c("Corail sensible 1","Corail sensible 2"),
#names(cat2keep) <- type.corail
#cat2keep <- cat2keep[!(names(cat2keep)=="Sensibilite")]

# formerly TB.lit
LIT.tableau.brut.OLD <- function(save=FALSE,filt.camp="X",type.db="LIT",
                             wZeroT=TRUE, wZeroSt=FALSE, silent=FALSE) {


    if(!silent) {
  departFunk() # message de depart
  on.exit(EM())
}

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

  # Rajouter infos sur les coraux
    DL <- data.frame(index.LIT[DL$Code_LIT, type.corail], DL)

  # 2. Filter les zeros en fonction de wZeroT et wZeroSt
  # DL contient deja les zeros pour toutes les combinaisons de campagnes/transects/Code_LIT

  if(wZeroSt) { # si wZeroSt = TRUE, ote les zeros si l'espece est observee sur la station mais pas tous les transects
      # independemment de la valeur de wZeroT
      DL.ZeroSt <- DL %>% group_by(Campagne, St, Code_LIT) %>% summarize(status.X=all(X.==0))
      DL %<>% inner_join(DL.ZeroSt) %>% filter( ((X.>0) | status.X))
  } else {
      # si wZeroSt ET wZeroT = FALSE, on ote les zeros partout
      if(!wZeroT) DL <- filter(DL, X.>0)
  }

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
    dd.i <- merge(dd.i, dd[[5]]$dat, by=c("Campagne","St", smpl))
    # ... à rajouter au besoin si une sixième catégorie est nécéssaire (e.g. Sensibilite)
    #dd.i <- merge(dd.i, dd[[6]]$dat, by=c("Campagne","St", smpl))

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
  attr(dd.i2, "wZeroT") <- wZeroT
    attr(dd.i2, "wZeroSt") <- wZeroSt
    invisible(dd.i2) }


#################################################################

# Tableau moyenne/SE suivant categories dans "S_Corail_All",
# par geomorphologie pour une année spécifique
# formerly LIT.ts1()
LIT.couvrt.gnrl <- function(ftemp=ftempo.defaut, fspat=fspat.defaut,
                            wZeroT=TRUE, wZeroSt=FALSE, par.transect="pas implemente",
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
        LIT.brut <<- LIT.tableau.brut(filt.camp=filt.camp, type.db=type.db,
                                      wZeroT=wZeroT, wZeroSt=wZeroSt)
      }else{
        if((attr(LIT.brut,"AS")!=filt.camp)|
           (attr(LIT.brut, "Projet")!=projet)|
           (attr(LIT.brut, "type.db")!=type.db)|
           (attr(LIT.brut, "wZeroT")!=wZeroT)|
           (attr(LIT.brut, "wZeroSt")!=wZeroSt))
            LIT.brut <<- LIT.tableau.brut(filt.camp=filt.camp, type.db=type.db)
    }

    t1 <- LIT.brut # transfert des données à objet t1

    # Creer tableau donnees brutes pour analyses subsequentes
    # Appliquer filtres (ref: LIT.doc)
    if(type.db=="LIT") smpl <- "T"
    if(type.db=="Quadrat") smpl <- "Quadrat"

    # Groupes de coraux
    gC <- cat2keep[[LIT.cat]]
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



###########################################################################
###########################################################################
###########################################################################
## exploring alternative way of processing transect summaries by LIT category
## derniere version utilisee
LIT.tableau.brut <- function(save=FALSE,filt.camp="X",type.db="LIT",
                     wZeroT=TRUE, wZeroSt=FALSE,
                     frmt.ret="wide", silent=FALSE) {

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

  # Rajouter infos sur les coraux
    DL <- data.frame(index.LIT[DL$Code_LIT, type.corail], DL)

  # 2. Filter les zeros en fonction de wZeroT et wZeroSt
  # DL contient deja les zeros pour toutes les combinaisons de campagnes/transects/Code_LIT

  if(wZeroSt) { # si wZeroSt = TRUE, ote les zeros si l'espece est observee sur la station mais pas tous les transects
      # independemment de la valeur de wZeroT
      DL.ZeroSt <- DL %>% group_by(Campagne, St, Code_LIT) %>% summarize(status.X=all(X.==0))
      DL %<>% inner_join(DL.ZeroSt) %>% filter( ((X.>0) | status.X))
  } else {
      # si wZeroSt ET wZeroT = FALSE, on ote les zeros partout
      if(!wZeroT) DL <- filter(DL, X.>0)
  }

    DL$smpl.unit <- DL[,smpl]
    subs <- DL[,c("Id","Code_LIT",names(cat2keep),"Campagne","St","smpl.unit","X.")]
    indx.X <- DL[, c("Id","Code_LIT","X.")] # couverture
    d1 <- melt(subs, id.vars=c("Campagne","St","smpl.unit","Id","Code_LIT", "X."))
    cat2keep.id <- paste(rep(names(cat2keep), sapply(cat2keep, length)), unlist(cat2keep))
    d2 <- d1 %>% filter(paste(variable, value) %in% cat2keep.id) %>%
        group_by(Campagne, St, smpl.unit, Id, variable, value) %>%
            summarize(sum.val=sum(X., na.rm=TRUE)) %>% data.frame
    d2$sum.val[is.na(d2$sum.val)] <- 0
    d2.wide <- dcast(d2, Campagne + St + smpl.unit + Id ~ value, value.var="sum.val", fill=0)
    d2[,smpl] <- d2$smpl.unit
    d2.wide[,smpl] <- d2.wide$smpl.unit

    # 4. Rajouter colonnes infos additionelles
                                        # Geomorphologies
    indx.fields <- unique(c("Campagne","St","Id",facteurs.spatio, facteurs.tempo))
    dd.i2 <- merge(unique(info.transect[,indx.fields]), d2.wide)
    dd.i2 <- merge(LIT.transect.info, dd.i2) # only keep transect initially present in dataset
                                             # ... could be removed I think
    if(frmt.ret == "long") {
        dd.i2.long <- inner_join(unique(info.transect[,indx.fields]), d2)
        names(dd.i2.long)[names(dd.i2.long)=="variable"] <- "LIT.cat"
        names(dd.i2.long)[names(dd.i2.long)=="value"] <- "LIT.lev"
        names(dd.i2.long)[names(dd.i2.long)=="sum.val"] <- "pcouv"

    }

    # 5. Reordonner selon instructions
    niv.all <- unique(d2$value) # niveaux existants dans la base
    dd.i2 <- dd.i2[,c(indx.fields, smpl, niv.all)]

    if(save) {
        # nametag pour filtre au besoin
        ftag <- ifelse(filt.camp %in% c("A","S"),sprintf("%s_",filt.camp),"")
        write.csv(dd.i2, file=paste(tabl.dir,"GS_",type.db,"_TableauBrut_",ftag,
                             Sys.Date(),".csv",sep=""),
                             row.names=FALSE)
    }

    if(frmt.ret == "long") dd.i2 <- dd.i2.long
    attr(dd.i2, "AS") <- filt.camp
    attr(dd.i2, "Projet") <- projet
    attr(dd.i2, "type.db") <- type.db
    attr(dd.i2, "wZeroT") <- wZeroT
    attr(dd.i2, "wZeroSt") <- wZeroSt

    invisible(dd.i2)
}
