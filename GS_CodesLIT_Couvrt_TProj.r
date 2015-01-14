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
LIT.tableau.brut <- function(save=FALSE,AS="pas de filtre") {

    departFunk() # message de depart
    on.exit(EM())

    # Creer tableau donnees brutes pour analyses subsequentes
    # Appliquer filtres (ref: LIT.doc)
    # 1. AQCQ == NON
    DL <- dLIT[dLIT$AQCQ == "NON",]
    LIT.transect.info <- unique(DL[,c("Campagne","Date","St",
                                            "Obs","T","AQCQ")])
    LIT.transect.info$Camp.ID <- paste(LIT.transect.info$St,
                                   LIT.transect.info$Campagne,sep="_")

    # 2. Rajouter
    DL <- data.frame(index.LIT[DL$Code_LIT, type.corail], DL)
    #DL <- merge(index.LIT[,c("Code_LIT",type.corail)],DL,by="Code_LIT") #above faster

    # 3. Calculer couverture moyenne par Campagne/St/Transect/TypeDeCoraux

    ccfunk <- function(wc) {
      print(paste(wc, paste(cat2keep[[wc]],collapse=", "),sep=": "))
        dl.i <- aggregate(list("PC"=DL$X.),
                   as.list(DL[,c("Campagne","St","T", wc)]),sum)

        dl.ii <- reshape(dl.i, timevar=wc,
                         idvar=c("Campagne","St","T"),direction="wide")

        names(dl.ii) <- gsub("PC.","",names(dl.ii))
        dl.ii <- dl.ii[,c("Campagne","St", "T",
                          cat2keep[[wc]])]

        }

    # Combiner les tableaux en 1
    dd <- lapply(type.corail, ccfunk) # boucle sur tous les types de coraux
    dd.i <- merge(dd[[1]], dd[[2]], by=c("Campagne","St", "T"))
    dd.i <- merge(dd.i, dd[[3]], by=c("Campagne","St", "T"))
    dd.i <- merge(dd.i, dd[[4]], by=c("Campagne","St", "T"))
    #dd.i <- merge(dd.i, dd[[5]], by=c("Campagne","St", "T")) # pour corail sensible


    # 4. Rajouter colonnes infos additionelles
    # Geomorphologies
    dd.i2 <- merge(unique(info.transect[,c("St","Geomorpho","N_Impact")]),dd.i,by="St")
    dd.i2 <- merge(LIT.transect.info, dd.i2, by=c("Campagne","St","T"))

    # 5. Reordonner selon instructions
    dd.i2 <- dd.i2[,c("Campagne","Date","St","Camp.ID",
                      "Geomorpho","Obs","T","AQCQ","N_Impact",
                   unlist(cat2keep))]

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
    invisible(dd.i2) }

#################################################################


# Tableau moyenne/SE suivant categories dans "S_Corail_All",
# par geomorphologie, en 2011
# formerly LIT.ts1()
LIT.resume <- function(yy=2011, fspat="Geomorpho", ff="General", AS="A", save=FALSE) {
  # set AS to "A" or "S" based on campagne type

    ################################
    departFunk() # message de depart
    on.exit(EM())
    ################################

    # vérifier si LIT.brut existe déjà pour la bonne campagne, sinon refaire
    if(!exists("LIT.brut")) {
      LIT.brut <<- LIT.tableau.brut(AS=AS)
      }else{
        if(attr(LIT.brut,"AS")!=AS) LIT.brut <<- LIT.tableau.brut(AS=AS) }

    dl <- LIT.brut[LIT.brut$Campagne %in% paste(AS,yy,sep="_"),]

    # Groupes de coraux
    gC <- coraux.fig[[ff]]

  # Moyenne/SE par type corail par facteur 'fspat'
  if(length(fspat)==1) { fc.list <- list(dl[,fspat])
                       }else{ fc.list <- as.list(dl[,fspat]) }

    # redéfinir fonction pour arrondir
    mean.rnd <- function(x, ...) round(mean(x,...),4)
    stand.err.rnd <- function(x, ...) round(stand.err(x,...),4)
    dmean <- aggregate(list("Moy"=dl[,gC]),fc.list,mean.rnd)
    dse <- aggregate(list("SE"=dl[,gC]),fc.list,stand.err.rnd)

  # Assemblage du tableau
  d.all <- merge(dmean,dse)
  if(save) write.csv(d.all,
                     file=paste(tabl.dir,"GS_LIT_Tableau1A_",ff,"_",AS,"_",
                       yy,"_",Sys.Date(),".csv",sep=""),row.names=FALSE)

  invisible(d.all)
}
