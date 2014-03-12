# Ginger/Soproner
# Code pour analyses des donnees LIT
# Time-stamp: <2014-03-12 12:51:01 Laura>
try.wd <- try(setwd(dossier.R),silent=TRUE)
if(class(try.wd)=="try-error") {
    message("\nCommencez par charger le fichier GS_KNS_MotherCode.r
 pour que dossier.R soit defini \n") }
fig.dir <- paste(dossier.R,"//Graphiques//",sep='')
tabl.dir <- paste(dossier.R,"//Tableaux//",sep='')

 # Extraction tableaux des bases de donnees
if(!exists("data.read")) source("GS_ExtractionDonnees.r")

# formerly TB.lit
LIT.tableau.brut <- function(save=FALSE,AS="pas de filtre") {

    departFunk() # message de depart
    on.exit(EM())

    # Creer tableau donnees brutes pour analyses subsequentes
    # Appliquer filtres (ref: LIT.doc)
    # 1. AQCQ == NON
    DL <- data.LIT[data.LIT$AQCQ == "NON",]
    LIT.transect.info <- unique(DL[,c("Campagne","Date","St",
                                            "Obs","T","AQCQ")])
    LIT.transect.info$Camp.ID <- paste(LIT.transect.info$St,
                                   LIT.transect.info$Campagne,sep="_")

    # 2. Rajouter info type de coraux
    type.corail <- c("S_Corail_All","S_Corail_Acro",
                     "S_Corail_Forme","S_Corail_Sensi","CODE_DET")
    cat2keep <- list(c("Coraux","Coraux morts","Algues",
                       "Abiotique","Autre faune"),
                     c("Coraux mous","Acroporidae","Non-acroporidae"),
                     c("Corail branchu","Corail tabulaire",
                       "Corail massif","Corail encroutant",
                       "Corail foliaire","Corail submassif","Corail digite"),
                     c("Coraline","Corail sensible1","Corail sensible3"),
                     c("Macro-algues","Assemblage d'algues"))


    DL <- merge(index.LIT[,c("Code_LIT",type.corail)],DL,by="Code_LIT")

    # 3. Calculer couverture moyenne par Campagne/St/Transect/TypeDeCoraux
    ccfunk <- function(wc) {
        dl.i <- aggregate(list("PC"=DL$X.),
                    as.list(DL[,c("Campagne","St","T", wc)]),sum)
        dl.ii <- reshape(dl.i, timevar=wc,
                         idvar=c("Campagne","St","T"),direction="wide")
        names(dl.ii) <- gsub("PC.","",names(dl.ii))
        dl.ii <- dl.ii[,c("Campagne","St", "T",
                          cat2keep[[which(type.corail==wc)]])]
        }

    # Combiner les tableaux en 1
    dd <- lapply(type.corail, ccfunk) # boucle sur tous les types de coraux
    dd.i <- merge(dd[[1]], dd[[2]], by=c("Campagne","St", "T"))
    dd.i <- merge(dd.i, dd[[3]], by=c("Campagne","St", "T"))
    dd.i <- merge(dd.i, dd[[4]], by=c("Campagne","St", "T"))
    dd.i <- merge(dd.i, dd[[5]], by=c("Campagne","St", "T"))

    # 4. Rajouter colonnes infos additionelles
    # Geomorphologies
    dd.i2 <- merge(info.transect[,c("St","Geomorpho","N_Impact")],dd.i,by="St")
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

    invisible(dd.i2) }


# Tableau moyenne/SE suivant categories dans "S_Corail_All",
# par geomorphologie, en 2011
# formerly LIT.ts1()
LIT.resume <- function(yy=2011, ff="Coraux_Gen", AS="A", save=FALSE) {
  # set AS to "A" or "S" based on campagne type

    ################################
    departFunk() # message de depart
    on.exit(EM())
    ################################

    if(!exists("LIT.brut")) LIT.brut <<- LIT.tableau.brut(AS=AS)
    dl <- LIT.brut[LIT.brut$Campagne %in% paste(AS,yy,sep="_"),]

    # Groupes de coraux
    gC <- coraux.fig[[ff]]

  # Moyenne/SE par type corail par geomorphologie
  fc.list <- list("Geomorpho"=dl$Geomorpho)
  dmean <- aggregate(list("Moy"=dl[,gC]),fc.list,mean)
  dse <- aggregate(list("SE"=dl[,gC]),fc.list,stand.err)

  # Assemblage du tableau
  d.all <- merge(dmean,dse,by="Geomorpho")
  if(save) write.csv(d.all,
                     file=paste(tabl.dir,"GS_LIT_Tableau1A_",ff,"_",AS,"_",
                       yy,"_",Sys.Date(),".csv",sep=""),row.names=FALSE)

  invisible(d.all)
}

LIT.bp1 <- function(yy=2011, ff2="Coraux_Gen", AS="A") {

    ################################
    departFunk() # message de depart
    on.exit(EM())
    ################################
    tb1.lit <- LIT.resume(yy=yy, ff=ff2, AS=AS)
    hist.funk <- function(cat) {

    par(family="serif",omi=c(0.25,0,0,0))

    xx <- tb1.lit[,make.names(paste("Moy",cat))]
    xx.se <- tb1.lit[,make.names(paste("SE",cat))]
    se.up <- t(xx+xx.se); se.down <- t(xx-xx.se)
    yy <- tb1.lit$Geomorpho

    colv <- c("grey25","grey45","grey65")[1:length(cat)]
    bp=barplot(t(xx), ylab=expression(paste("Couverture moyenne (","%" %+-% SE,")")),
            ylim=c(0,1.25*max(se.up)),beside=T,col=colv,
           las=1, border=NA,cex.lab=1.2)

    ya <- par("usr")[4]/diff(par("plt")[4:3])
    text(apply(bp,2,mean), rep(0.075*ya,5), gsub(" "," \n ",yy),xpd=NA)
    abline(h=0)
    mtext("Geomorphologie",side=1,outer=T,cex=1.2)


    dmm <- sapply(1:length(se.up),function(i) arrows(bp[i],se.down[i],bp[i],se.up[i],
                                           code=3,angle=90,length=0.1))
    legend("topleft",cat,fill=colv,bty="n",horiz=TRUE, inset=c(0,-0.1),
           x.intersp=0.5, cex=1.2, xpd=NA)

    dev.copy2pdf(file=paste(fig.dir,"GS_LIT_Hist1_",yy,"_",make.names(cat),".pdf",sep=""))
  }

  if(ff2=="Coraux_Gen") {
      catall <- list(c("Coraux","Coraux morts","Coraux mous"),
                     c("Algues","Abiotique","Autre faune"))
      } else { catall <- coraux.fig[ff2] }

  dmm <- sapply(catall, hist.funk)
}


LIT.ts1 <- function(AS="A") { # AS = "A" pour annuelles, "S" pour semestrielles

  ################################
  departFunk() # message de depart
  on.exit(EM())
  ################################

  cnow <- coraux.fig$TS_All # categorie a utiliser

  # Stations echantillonees depuis 2006, a part pour "Recif barriere externe"
  # (distinction deja prise en compte dans le tableau filtre)
  wf <- paste("T",AS,"LIT",sep="_") # colonne du filtre
  dd.filt <- filtreTable(LIT.brut, wf)

  # Calculs moyenne/SE
  val.mean <- aggregate(list("Moy"=dd.filt[,cnow]),
                        dd.filt[,c("Campagne","Geomorpho","N_Impact")],mean)
  val.se <- aggregate(list("SE"=dd.filt[,cnow]),
                      dd.filt[,c("Campagne","Geomorpho","N_Impact")],stand.err)
  val.all <- merge(val.mean, val.se, by=c("Campagne","Geomorpho","N_Impact"))

  # Calcul de la difference de couverture par Campagne/Geomorpho
  # s'assurer qu'il y a une rangee pour toutes les combinaisons Geomorpho/Impact
  tmpl <- expand.grid("Campagne"=unique(val.all$Campagne),
                      "Geomorpho"=unique(val.all$Geomorpho),
                      "N_Impact"=unique(val.all$N_Impact))
  val.all2 <- merge(val.all, tmpl, by=c("Campagne","Geomorpho", "N_Impact"), all.y=TRUE)

  # et aussi que les rangees sont dans l'ordre impact 1, impact 0
  val.all2 <- val.all2[order(val.all2$Campagne, val.all2$Geomorpho,
                                 -val.all2$N_Impact),]

  # Calculer la difference par stratum
  val.impact <- aggregate(list("Diff"=val.all2[,make.names(paste("Moy",cnow))]),
                          by=list("Campagne"=val.all2$Campagne,
                            "Geomorpho"=val.all2$Geomorpho), diff)


  VM.split <- split(val.all, list(val.all$Geomorpho))

  fig.funk <- function(wgeo,wmorph) {
      par(omi=rep(0,4), family="serif")

      dnow <- VM.split[[wgeo]] # selectionne la geomorphologie pour le graphique courant
      # identifie les colonnes contentant la moyenne et le SE:
      wcol <- which(names(dnow) %in% make.names(paste(c("Moy","SE"),wmorph, sep=".")))
      dnow <- dnow[,c(1:3,wcol),]
      dnow$annee <- as.numeric(gsub("([AS_])","",dnow$Campagne)) # converti campagne -> annee
      dnow$annee[grepl("S",dnow$Campagne)] <- dnow$annee[grepl("S",dnow$Campagne)]-0.5

      # extraire SE
      dnow$seU <- dnow[,grepl("Moy",names(dnow))]+dnow[,grepl("SE",names(dnow))]+0.01
      dnow$seL <- dnow[,grepl("Moy",names(dnow))]-dnow[,grepl("SE",names(dnow))]-0.01
      dnow <- dnow[order(dnow$annee),]

      yl <- range(c(dnow[,c("seU","seL")])) # limites de l'axe Y pour zones impact et non-impact
      sdnow <- split(dnow, dnow$N_Impact) # diviser le data frame selon impact/pas-impact

      # lignes principales
      plot(sdnow[[1]]$annee, sdnow[[1]][,grepl("Moy",names(dnow))],
         xlab="Campagne", ylab=expression(paste("Couverture moyenne (","%" %+-% SE,")")),
         type="b",ylim=yl, cex.lab=1.2, xaxt="n")
      axis(1,at=sdnow[[1]]$annee, labels=sdnow[[1]]$Campagne) #axe des X
      mtext(paste(wgeo," - ", wmorph,":", sep=""), cex=1.3, adj=0, line=0.25)

      if(length(sdnow)>1) {
          lines(sdnow[[2]]$annee, sdnow[[2]][,grepl("Moy",names(dnow))],
                type="b", pch=2, lty=2) }

      # rajouter les lignes SE
      dmm <- sapply(1:nrow(dnow),function(i) arrows(dnow[i,"annee"],dnow[i,"seU"],
                                                    dnow[i,"annee"],dnow[i,"seL"],
                                                    code=3,angle=90,length=0.1))

      # rajouter legende
      if(length(sdnow)>1) {
          legend("topright",legend=c("Zone de reference","Zone d'impact"),
                 lty=1:2,pch=1:2,bty="n",xpd=NA, inset=c(0,-0.175), cex=1.2)
          } else {
              legend("topright",
                     legend=ifelse(sdnow[[1]]$N_Impact[1] == 0, "Zone de reference","Zone d'impact"),
                     lty=1,pch=1,bty="n",xpd=NA, inset=c(0,-0.15), cex=1.2) }

      # sortir la figure en pdf
      dev.copy2pdf(file=paste(fig.dir,"GS_LIT_SerieTempTable_",AS,"_GeoMorphImpact_",wmorph,
                   "_",wgeo,"_",Sys.Date(),".pdf",sep=""))

}

  # double boucle sur la fonction fig.funk par geomorphologie X type de corail
  geo.id <- unique(LIT.brut$Geomorpho)
  dmm <- sapply(geo.id, function(gg)
                sapply(cnow,function(cc) fig.funk(wgeo=gg, wmorph=cc)))

  # sauvegarde tableaux
  write.csv(val.all, paste(tabl.dir,"GS_LIT_SerieTempTable_",AS,"_GeoMorphImpact_",Sys.Date(),".csv",sep=""),
            row.names=FALSE)
  write.csv(val.impact,
            paste(tabl.dir,"GS_LIT_SerieTempTable_",AS,"_GeoMorphImpact_DiffCouv_",Sys.Date(),".csv",sep=""),
            row.names=FALSE)

  invisible(list(val.all, val.impact))
}

LIT.ts2 <- function(AS="A") { # AS = "A" pour annuelles, "S" pour semestrielles

    ################################
    departFunk() # message de depart
    on.exit(EM())
    ################################

    check.dev.size(8.5, 7) # ouvrir une fenetre graphique de la bonne grandeur au besoin
    cnow <- coraux.fig$TS_All

    # Stations echantillonees depuis 2006, a part pour "Recif barriere externe"
    # (distinction deja pris en compte dans le tableau filtre)
    wf <- paste("T",AS,"LIT",sep="_") # colonne du filtre
    dd.filt <- filtreTable(LIT.brut, wf)

    # calcul moyenne/SE
    val.mean <- aggregate(list("Moy"=dd.filt[,cnow]),
                        dd.filt[,c("Campagne","Geomorpho","N_Impact","St")],mean)
    val.se <- aggregate(list("SE"=dd.filt[,cnow]),
                      dd.filt[,c("Campagne","Geomorpho","N_Impact","St")],stand.err)
    val.all <- merge(val.mean, val.se, by=c("Campagne","Geomorpho","N_Impact","St"))

    # diviser le tableau par geomorpho/impact
    VM.split <- split(val.all, list(val.all$Geomorpho,val.all$N_Impact))

    fig.funk <- function(wgeo,wmorph,wimpact) {
        par(family="serif", omi=c(0,0,0,1))

        dnow <- VM.split[[paste(wgeo,wimpact,sep=".")]]

        if(nrow(dnow)>0){
            # id des colonnes avec moyenne et SE
            wcol <- which(names(dnow) %in% make.names(paste(c("Moy","SE"),wmorph, sep=".")))
            dnow <- dnow[,c(1:4,wcol),]

            dnow$annee <- as.numeric(gsub("([AS_])","",dnow$Campagne))
            dnow$annee[grepl("S",dnow$Campagne)] <- dnow$annee[grepl("S",dnow$Campagne)]-0.5
            dnow$seU <- dnow[,grepl("Moy",names(dnow))]+dnow[,grepl("SE",names(dnow))]+0.01
            dnow$seL <- dnow[,grepl("Moy",names(dnow))]-dnow[,grepl("SE",names(dnow))]-0.01
            dnow <- dnow[order(dnow$annee),]

            yl <- range(c(dnow[,c("seU","seL")])) # limites axe Y pour impact / pas impact
            if(max(yl)<10) yl <- c(0,20) # si valeur max trop petite imposer axe Y entre 0 et 20

            sdnow <- split(dnow, dnow$St) #diviser tableau entre impact / pas d'impact

            # lignes principales
            colv=rep(1:length(sdnow), nrow(sdnow[[1]]))
            plot(sdnow[[1]]$annee, sdnow[[1]][,grepl("Moy",names(dnow))],
                 xlab="Campagne", ylab=expression(paste("Couverture moyenne (","%" %+-% SE,")")),
                 type="n",ylim=yl, cex.lab=1.2, xaxt="n", las=1)
            axis(1,at=sdnow[[1]]$annee, labels=sdnow[[1]]$Campagne) # axe X

            # rajouter lignes SE
            dmm <- sapply(1:nrow(dnow),function(i) arrows(dnow[i,"annee"],dnow[i,"seU"],
                                                  dnow[i,"annee"],dnow[i,"seL"],
                                                  code=3,angle=90,length=0.1,
                                                  col=colv[i]))

            sapply(1:length(sdnow), function(i)
                   lines(sdnow[[i]]$annee, sdnow[[i]][,grepl("Moy",names(dnow))], type="b",
                         col=i,pch=14+i))

            #titre
            mtext(paste(wgeo," - ", wmorph,":", sep=""), cex=1.1, adj=0, line=0.25)

            #legende
            legend("topright",legend=names(sdnow),lty=1,col=colv[1:length(sdnow)],
                   pch=15:20,bty="n",xpd=NA, inset=c(-0.21,0), cex=1.1)

            #sauvegarde pdf
            dev.copy2pdf(file=paste(fig.dir,"GS_LIT_SerieTemp_",AS,"_GeoMorphImpact_bySt_",wmorph,
                   "_",wgeo,"_",wimpact, "_",Sys.Date(),".pdf",sep=""))
  }
  }
    # lancer la boucle sur la fonction fig.funk par impact X geomorpho X type corail
    geo.id <- unique(LIT.brut$Geomorpho)
    dmm <- sapply(c(0,1), function(wi) sapply(geo.id, function(gg)
                                              sapply(cnow,function(cc)
                                                     fig.funk(wgeo=gg, wmorph=cc,wimpact=wi))))

    #sauvegarder tableaux
    write.csv(val.all,
              paste(tabl.dir,"GS_LIT_SerieTempTable_",AS,"_GeoMorphImpact_bySt_",Sys.Date(),".csv",sep=""),
              row.names=FALSE)
    invisible(val.all)
}

