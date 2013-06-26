# Ginger/Soproner
# Time-stamp: <2013-06-26 15:41:46 Laura>
#library(XLConnect) # pas pour le moment

setwd("/Users/Laura/Documents/Projects/Soproner_Noumea")
# Extraction tableaux des bases de données
if(!exists("data.read")) source("GS_ExtractionDonnees.r")

# Définition lien fichiers pour sauvegarder graphiques/tableaux
fig.dir <- "/Users/Laura/Documents/Projects/Soproner_Noumea/Graphiques/"
tabl.dir <- "/Users/Laura/Documents/Projects/Soproner_Noumea/Tableaux/"

library(plotrix)

# Fichiers des données de comptage = data.poissons

########################################
############# TABLEAUX #################
########################################

tableau.brut <- function(save=FALSE) {

  # Calcul de densité et de biomasse par observation/espèce/transect/campagne
  # Formule de distance sampling ajustée pour les observations individuelles
  # D.obs = N/(2xdm.obsxL); dm.obs= 0.5*(d1+d2) + 0.5
  # Bio.obs = N*a*T^b/(2xdm.obsxL)

  # par défaut mais changer une fois que la colonne est ajoutée dans poissons.xls
  LT <- 50 # Longueur du transect
  fc <- c("Campagne","St","Code.SP")

  ds.calc <- data.poissons[,c(fc,"OBSERVATEUR","NOMBRE","TAILLE","D1","D2")]
  names(ds.calc)[names(ds.calc)=="NOMBRE"] <- "N" # renommer colonne NOMBRE
  print(paste(nrow(ds.calc[!(ds.calc$Code.SP %in% bioeco$Code.SP),]),
              "rangées sans info bioeco"))
  print(unique(ds.calc[!(ds.calc$Code.SP %in% bioeco$Code.SP),"Code.SP"]))
  ds.calc <- merge(ds.calc, bioeco, by="Code.SP")

  # Densité/biomasse par *observation*
  ds.calc$dm.obs <-  0.5*(ds.calc$D1+ds.calc$D2)+0.5 # dm.observation
  ds.calc$dens.obs <- ds.calc$N/(2*ds.calc$dm.obs*LT) # densité.observation
  ds.calc$bio.obs <- ds.calc$N *ds.calc$a*(ds.calc$TAILLE^ds.calc$b)/(2*ds.calc$dm.obs*LT)

  # ** Assemblage ** PRODUIT 1 ** Tableau données brutes **
  # Densité/biomasse par observation/espèce/transect/année/campagne
  tb.1 <- ds.calc[,c("Campagne","St","Code.SP","OBSERVATEUR","N",
                     "dens.obs","bio.obs","D1","D2","TAILLE","a","b")]

  # Joint à:
  # Charactéristiques du transect: année/campagne/transect/plongeur/zone impact/...
  # ... géomorphologie/pêche cpue/pêche effort/taux de sédimentation/...
  # profondeur moyenne (Z),
  print("rajouter *taux sédimentation* au tableau brut poissons")

  tb.1 <- merge(tb.1,
                info.transect[,c("St","Geomorpho","N_Impact","cpus","Effort_ha","Z")],
                by="St")
  tb.1$Tx.Sed <- NA

  # + Données lit
  print("rajouter données LIT au tableau brut poissons")
  tb.1$Donnees.LIT <- NA

  # + Charactéristiques de l'espèce: code_espèce/famille/genre/espèce/...
  # ... état commercial/groupe trophique
  tb.2 <- merge(tb.1, bioeco.all[,c("Code.SP","Famille","Genre","Espece",
                                    "Peche","Cible","GTlabel","moblabel")],
                by="Code.SP",sort=FALSE)

  tb.2$An <- apply(sapply(2006:2011, function(x) x*grepl(x, tb.2$Campagne)),1,sum)
  tb.2$Camp <- gsub("(\\w)(\\w*)","\\U\\1",tb.2$Campagne,perl=TRUE)

  tb.2$Camp[tb.2$Camp == "S"] <- "semestrielle"
  tb.2$Camp[tb.2$Camp == "A"] <- "annuelle"

  tb.2 <- tb.2[order(tb.2$Camp,tb.2$St, tb.2$Code.SP),]

  tb.2 <- tb.2[,c("An","Camp","Campagne","St","OBSERVATEUR","N_Impact","Geomorpho",
                  "cpus","Effort_ha","Tx.Sed","Z",
                  "Donnees.LIT","Code.SP","Famille","Genre","Espece","N","TAILLE","D1","D2",
                  "Peche","Cible","a","b","dens.obs","bio.obs","GTlabel","moblabel")]

  names(tb.2) <- c("An","Campagne","Campagne.ID","St","Plongeur","Zone.Impact","Geomorpho",
                  "CPUS","Effort","Tx.Sed","Profondeur",
                  "Donnees.LIT","Code.SP","Famille","Genre","Espece","N","TAILLE","D1","D2",
                  "Peche","Cible","Coeff.a","Coeff.b","dens.obs","bio.obs",
                   "Groupe.Trophique","Groupe.Mobil")

  tb.2.0 <<- tb.2[tb.2$N>0,]
  # garde seulement observations N>0 pour réduire la taille du fichier .xls
  # déclare tb.2 dans l'environnement global (pour fct write.xls)

  if(save) {
  # write.csv(t2b,file=paste(tabl.dir,
  # "GS_Poissons_TableauDonnéesBrutes_",Sys.Date(),".csv",sep=""),
  # row.names=FALSE)
  write.xls(tb.2.0,file=paste(tabl.dir,
  "GS_Poissons_TableauDonnéesBrutes_",Sys.Date(),".xls",sep=""))
  }

  return(tb.2) }

BD.by.sp <- function() {

  # Calculate density/biomass by *species/transect/campagne* (not by observation)(??)
  fc <- c("Campagne","St","Code.SP")
  LT <- 50 # Longueur du transect

  ds.calc$dm.int <- ds.calc$N * (0.5*(ds.calc$D1+ds.calc$D2)+0.5)
  ds.list <- list("Campagne"=ds.calc$Campagne.ID,"St"=ds.calc$St,"Code.SP"=ds.calc$Code.SP)
  num.dm <- aggregate(ds.calc$dm.int,ds.list,sum)
  denom.dm <- aggregate(ds.calc$N,ds.list,sum)
  ds.df <- data.frame(num.dm[,fc], "dm"=num.dm$x/denom.dm$x)
  ds.df$dens <- denom.dm$x/(2*ds.df$dm*LT) # Densité
  ds.df$allN <- denom.dm$x # nombre individus total espèce/site/campagne

  # Calcul de la biomasse:
  ds.calc$bio.int <- ds.calc$N * ds.calc$Coeff.a * ds.calc$TAILLE^ds.calc$Coeff.b
  num.bio <- aggregate(ds.calc$bio.int, ds.list, sum)
  ds.df$biomasse <- num.bio$x/(2*ds.df$dm*LT) # Biomasse

  # Calcul de la taille moyenne:
  ds.calc$TM <- ds.calc$N*ds.calc$TAILLE
  ds.df$taille.moy <- aggregate(ds.calc$TM,ds.list,sum)[,"x"]/ds.df$allN

  # Set biomasse/density to 0 if allN = 0, leave taille.moy at NA
  # (we want taille.moy of observed individuals)
  ds.df[ds.df$allN==0,c("dens","biomasse")] <- 0

  # ordonner par campagne/sites
  ds.df <- ds.df[order(ds.df$St, ds.df$Campagne),]

}



poissons.ts1 <- function(save=FALSE) {

  # PRODUIT 2:
  # Tableau synthèse 1:
  # moyenne et écart type pour: densité (d), biomasse (b),
  # richesse spécifique (rs), taille moyenne (tm)
  # calculé sur:
  # toutes espèces confondues (tot), commerciales (com), carnivore (car),
  # herbivore (her), piscivore (pis), planctonophage (pla),
  # sédentaire (sed), territoriale (ter), mobile (mo), très mobile (tmo)
  # cible pêche nouvelle-calédonie, ou non

  dfh <- BDtable;
  dfh$Seen <- 0+ifelse(BDtable$allN==0,0,1) # colonne pour présence/absence
  ds.tmpl <- unique(BDtable[,c("St","Campagne")])
  ds.tmpl <- ds.tmpl[order(ds.tmpl$St, ds.tmpl$Campagne),]

  aggr.funk <- function(wm="mean",facteur="total") {
    # wm can be richness (sum), mean, sd
    # facteur all, Type_Comptage (commercial), Groupe_Trophique, Mobilite

    # liste de facteurs
    flist <- list("Campagne"=dfh[,"Campagne"],"St"=dfh[,"St"])
    if(facteur!="total") flist[[facteur]] <- dfh[,facteur] # rajoute facteur explicatif si specifie

    # aggrégation par facteur
    if(wm=="richness") { frez <- aggregate(list("sp.rich"=dfh[,"Seen"]), flist, sum)
                         } else {
                           frez <- aggregate(dfh[,c("biomasse","dens","taille.moy")],
                                             flist, wm,na.rm=TRUE)}

    if(facteur=="total"){
      nn <- ncol(frez)
      names(frez) <- c("Campagne","St",paste(names(frez)[3:nn],".Tot",sep=""))}

    # changement d'orientation du tableau, si besoin
    if(facteur != "total") {
      frez <- reshape(frez, timevar=facteur, idvar=c("Campagne","St"),direction="wide")}

    # fondre dans tableau cadre avec tous les facteurs
    frez.lab <- names(frez)
    frez.lab <- frez.lab[!(frez.lab %in% c("Campagne","St","biomasse","dens","taille.moy"))]



    frez.b <- merge(ds.tmpl,frez,by=c("St","Campagne"),all=TRUE)
    frez.c <- as.data.frame(frez.b[,frez.lab]); names(frez.c) <- frez.lab

    # Rajouter etiquette mean/sd
    if(wm %in% c("mean","sd")){
      names(frez.c) <- paste(ifelse(wm=="sd","ET","Moy"),
                             names(frez.c),sep=".")}

    return(frez.c)
    }

  # Pour assembler le tableau:
  vtype <- c("mean","sd","richness") # différencie richness vu qu'il n'y pas de moyenne/ET à calculer
  # Définir les facteurs d'aggrégation des espèces
  vcat <- c("total","Peche","Groupe.Trophique","Groupe.Mobil","fmlabel", "Cible")
  mean.all <- do.call(data.frame,lapply(vcat, function(i) aggr.funk(facteur=i)))
  sd.all <- do.call(data.frame,lapply(vcat, function(i) aggr.funk(wm="sd",facteur=i)))
  richness.all <- do.call(data.frame,lapply(vcat, function(i) aggr.funk(wm="richness",facteur=i)))
  df.all <- data.frame(ds.tmpl[,c("Campagne","St")],mean.all,sd.all,richness.all)
  df.all <- merge(info.transect[,c("St","Geomorpho","N_Impact","cpus","Effort_ha")],df.all,by="St")
  df.all <<- df.all # declare globally for write.xls function

  if(save) write.xls(df.all,file=paste(tabl.dir,"GS_Poissons_TS1_",Sys.Date(),".xls",sep=""))
  return(df.all)
}

# PRODUIT 3
# Tableau synthèse 2:
# pour chaque espèce, moyenne de densité pour chaque année et totale,
# pour chaque station et totale

poissons.ts2 <- function(save=FALSE) {

  ds.df <- BDtable

  # Filtrer les stations pour garder les stations répétées sur toutes les campagnes annuelles
  wfiltre <- "T_A_poissons"; print(paste("stations filtrées par",wfiltre))
  fdf <- data.frame(filtre.Camp[,wfiltre]); names(fdf) <- "key"
  st.keep <- merge(fdf,index.Camp)
  ds.df <- merge(st.keep, ds.df, by=c("Campagne","St"))

  # Format
  ds.df$Year <- apply(sapply(2006:2011, function(x) x*grepl(x, ds.df$Campagne)),1,sum)
  ds.df$Code.SP <- as.numeric(ds.df$Code.SP)
  ds.df <- ds.df[order(ds.df$Code.SP),]

  # Toutes années/stations confondues
  alltog <- aggregate(list("TtStAn"=ds.df$dens),
                       list("Code.SP"=ds.df$Code.SP),mean,na.omit=TRUE)

  # Par année/toutes stations confondues
  alltog.yr <- aggregate(list("TTSt"=ds.df$dens),
                         as.list(ds.df[,c("Year","Code.SP")]),mean,na.omit=TRUE)

  alltyr.h <- reshape(alltog.yr, timevar="Year",idvar="Code.SP",direction="wide")
  alltyr.h <- alltyr.h[,sort(names(alltyr.h))]

  # Par station/toutes années
  allyears <- aggregate(list("TtAn"=ds.df$dens),
                         as.list(ds.df[,c("St","Code.SP")]),mean, na.omit=TRUE)
  ay.horiz <- reshape(allyears, timevar="St", idvar="Code.SP",direction="wide")
  ay.horiz <- ay.horiz[,sort(names(ay.horiz))]
  # Par station/année
  byyear <- aggregate(list("Dens"=ds.df$dens),
                       as.list(ds.df[,c("St","Year","Code.SP")]),mean, na.omit=TRUE)
  by.horiz <- reshape(byyear, timevar="St", idvar=c("Code.SP","Year"),direction="wide")
  by.h.yr <- reshape(by.horiz, timevar="Year", idvar="Code.SP",direction="wide")
  by.h.yr <- by.h.yr[,sort(names(by.h.yr))]

  df <- do.call(data.frame,list(alltog,alltyr.h,ay.horiz, by.h.yr))

  # ote colonnes Code.SP redondantes
  df <- df[,!(names(df) %in% paste("Code.SP.",1:3,sep=""))]
  names(df) <- gsub("Dens.","",names(df)) # Ote "Dens." des noms de colonne

  df.i <<- merge(bioeco.all[,c("Code.SP","Famille","Genre","Espece")],df,by="Code.SP")

  if(save) write.xls(df.i, file=paste(tabl.dir,"GS_Poissons_TS2_",Sys.Date(),".xls",sep=""))
  return(df)
}

########################################
############ GRAPHIQUES ################
########################################

#################################
########### PRODUIT 3 ###########
#################################
poissons.p3 <- function(quel.graph="all",save=FALSE) {

  # Info par type de graphique:
  tgraph <- c("Densite","Biomasse","Richesse specifique")
  gr.lab <- c("Moy.dens","Moy.biomasse","sp.rich")
  yx.lab <- c(expression(paste("Densite (individus/",m^2,")",sep="")),
              expression(paste("Biomasse (g/",m^2,")",sep="")),
              "Richesse specifique")

# Filtrer comme spécifié dans Tableau ...
  wfiltre <- "T_A_poissons"; print(paste("stations filtrées par",wfiltre))
  fdf <- data.frame(filtre.Camp[,wfiltre]); names(fdf) <- "key"
  st.keep <- merge(fdf,index.Camp)
  dd <- merge(st.keep, TS1, by=c("Campagne","St"))
  dd$annee <- as.numeric(gsub("A_","",dd$Campagne))

  graph.funk <- function(wt="Densité",wgroup="Tot") {

    fig.counter <<- 0 # reset fig.counter
    par(mai=c(0.2,0.15,0.2,0.1),omi=c(0.5,0.75,0.5,0.1),
        mfrow=c(3,2), tcl=0.5, family="serif")
    # Aggrégation par campagne, toutes stations confondues
    wc=paste(gr.lab[tgraph==wt],wgroup,sep=".")

    # Par Campagne -> ôtée pour le moment
    # df.dens <- aggregate(dd[,wc], list("Campagne"=dd$Campagne), mean,na.rm=TRUE)
    # Par Campagne/Géométrie:
    df.dens <- aggregate(list("mean"=dd[,wc]), list("Annee"=dd$annee,
                                       "Geom"=dd$Geom), mean,na.rm=TRUE)

    df.dens.SD <- aggregate(list("sd"=dd[,wc]), list("Annee"=dd$annee,
                                          "Geom"=dd$Geom), sd,na.rm=TRUE)
    df.split.geom <- split(df.dens,df.dens$Geom)

    # Par géométrie:
    df.dens.geo.1 <- aggregate(list("mean"=dd[,wc]), list("Annee"=dd$annee,
                                           "Impact"=dd$N_Impact,
                                           "Geom"=dd$Geom), mean,na.rm=TRUE)
    df.dens.geo.2 <- aggregate(list("sd"=dd[,wc]), list("Annee"=dd$annee,
                                           "Impact"=dd$N_Impact,
                                           "Geom"=dd$Geom), sd,na.rm=TRUE)
    df.dens.geo <- merge(df.dens.geo.1,df.dens.geo.2,by=c("Annee","Impact","Geom"))

     # when there is only 1 observateur sd is calculated as NA
    df.dens.geo$sd[is.na(df.dens.geo$sd)] <- 0 # set NA sd to zero

    # ylim for all plots below
    yl <- c(0,max(df.dens.geo$mean+df.dens.geo$sd,na.rm=TRUE))
    df.split.impact <- split(df.dens.geo, list(df.dens.geo$Impact,df.dens.geo$Geom))

    sub.funk <- function(wgeom) {

      fig.counter <<- fig.counter + 1
      onames <- paste(c(0,1),wgeom, sep=".")

      # Moyenne
      plot(2006:2011, df.split.geom[[wgeom]]$mean,type="b",lwd=3,ylim=yl,
       xlab="",ylab=yx.lab[tgraph==wt],xaxt="n",yaxt="n",las=1)
      axis(1,labels=ifelse(fig.counter%in%c(4,5),TRUE,FALSE),cex.axis=1.5)
      axis(2,labels=ifelse(fig.counter%in%c(1,3,5),TRUE,FALSE),cex.axis=1.5,las=1)

      # impact=0
      # Barres d'erreurs (écart type)
      add.impact <- function(oo) {
        dn <- df.split.impact[[paste(oo-1,wgeom,sep=".")]]

        yvect=2006:2011
        sapply(1:length(yvect),function(x)
               arrows(yvect[x], dn$mean[x]-dn$sd[x]-0.005,yvect[x],
                      dn$mean[x]+dn$sd[x]+0.005,code=3,angle=90,
                      col="grey50",length=0.1))
        lines(2006:2011,dn$mean,lwd=2,lty=oo+1)}

      sapply(1:2,add.impact)
      mtext(wgeom)
      lines(2006:2011, df.split.geom[[wgeom]]$mean,lwd=3,type="b")
    }


    dmm=sapply(names(df.split.geom), sub.funk)
    mtext(paste(wt,":",graph.key[graph.key$df.id==wgroup,"titre.fig"]),
          line=1,outer=T,cex=1.2)

    mtext(yx.lab[tgraph==wt],side=2,outer=TRUE,line=3)

    plot(1,1,type="n",ann=FALSE,axes=FALSE)
    # Légende:
    legend("left",c("Moyenne","Temoin","Impact"),xpd=NA,
           lty=1:3,lwd=c(3,rep(2,length(df.dens.geo))),bty="n",
           title="Legende:",title.adj=0.045,cex=1.5)
    if(save){
    dev.copy2pdf(file=paste(fig.dir,"GS_PoissonsST_",paste(wt,wgroup,sep="."),"_",Sys.Date(),".pdf",                   sep="")) }
}

  if(quel.graph=="all"){
  dmm=sapply(graph.key$df.id, function(cc) sapply(tgraph, function(tt)
                                              graph.funk(wt=tt,wgroup=cc)))}
}

#######################################
########### Poissons cibles ###########
#######################################
poissons.cible.graph <- function(save=FALSE) {

  if(dev.size()[1] != 9 | round(dev.size()[2],1) != 5.3) quartz(width=9, height=5.3)

  # Pour chaque geomorphologie, faire une série de temps (+ filtrées)
  # démontrant biomasse/densité/richesse moyenne avec espèces ciblées/non-ciblées
  # sur la même fenêtre, divisée par impact
  gr.lab <- c("Moy.dens","Moy.biomasse","sp.rich")
  gr.lab.2 <- make.names(paste(paste(gr.lab, "Cible"),rep(c("Oui","Non"),3)))
  yx.lab <- c(expression(paste("Densite (individus/",m^2,")",sep="")),
              expression(paste("Biomasse (g/",m^2,")",sep="")),
              "Richesse specifique")

  # Filtre sur stations
  wfiltre <- "T_A_poissons"; print(paste("stations filtrées par",wfiltre))
  fdf <- data.frame(filtre.Camp[,wfiltre]); names(fdf) <- "key"
  st.keep <- merge(fdf,index.Camp)
  dd <- merge(st.keep, TS1, by=c("Campagne","St"))
  dd$annee <- as.numeric(gsub("A_","",dd$Campagne))


  # Extraction des colonnes cible/non-cibles
  dd <- dd[,c("annee","St","Geomorpho","N_Impact",gr.lab.2)]
  nm.dd <- unique(sub(".Non","",sub(".Oui","",names(dd)[grep("Cible",names(dd))])))

  # reshape pour mettre facteur cible/non-cible en une colonne
  dd2 <- reshape(dd, varying=names(dd)[grep("Cible",names(dd))],
                 idvar=c("annee","St","Geomorpho","N_Impact"),
                 v.names=nm.dd, times=c("Oui","Non"),
                 direction="long")
  names(dd2)[grep("sp.rich",names(dd2))] <-
    sub("sp.rich","Moy.sp.rich",names(dd2)[grep("sp.rich",names(dd2))])

  # Aggrégation des valeurs par géomorphologie/impact
  dd.imp <- lapply(c("mean","std.error"),
                     function(x) aggregate(dd2[,grep("Cible",names(dd2))],                                                 list("Annee"=dd2$annee,"Geomorpho"=dd2$Geomorpho,                                     "Impact"=dd2$N_Impact, "PecheNC"=dd2$time),x))

  names(dd.imp[[2]])[grep("Cible",names(dd.imp[[2]]))] <-
    sub("Moy","SE",names(dd.imp[[2]])[grep("Cible",names(dd.imp[[2]]))])
  dd.imp.all <- merge(dd.imp[[1]],dd.imp[[2]],
                      by=c("Annee","Geomorpho","Impact","PecheNC"))

  # Note: SE est NA s'il n'y a qu'une observation
  dd.imp.all[,grep("SE",names(dd.imp.all))][is.na(dd.imp.all[,grep("SE",names(dd.imp.all))])] <- 0
  dd.split <- split(dd.imp.all,dd.imp.all$Geomorpho)

  # Définition de la fonction graphique par géomorphologie:
  graph.geo <- function(wgeo, wtype) {

    par(mfrow=c(1,2),family="serif",mai=c(1,0.35,0.5,0.1),omi=c(0,0.75,0.5,0.1))
    dnow <- dd.split[[wgeo]]

    dnow <- dnow[,c(1,3,4,grep(sub("Moy.","",wtype),names(dnow)))]

    # Calcul des barres d'erreurs
    dnow$seU <- dnow[,grep("Moy",names(dnow))]+dnow[,grep("SE",names(dnow))]
    dnow$seL <- dnow[,grep("Moy",names(dnow))]-dnow[,grep("SE",names(dnow))]

    yl <- max(dnow$seU,na.rm=TRUE)

    # Diviser par facteurs
    dnow.spl <- split(dnow, dnow$PecheNC)

    # Get SE limits out
    selist <- lapply(c(0,1), function(x) dnow[dnow$Impact == x,c("Annee","seU","seL")])

    # Fonction graphique par impact/témoin
    lab.sf <- c("Zone temoin","Zone d'impact")
    sub.fig <- function(ii){

      dnow.fig <- dnow.spl[["Oui"]]
      plot(dnow.fig[dnow.fig$Impact==ii,c(1,grep("Moy",names(dnow.fig)))],
           type="b", las=1, xlab="Annee", ylab ="", ylim=c(0,yl),
           yaxt=ifelse(ii==0,"t","n"))
      dnow.fig <- dnow.spl[["Non"]]
      lines(dnow.fig[dnow.fig$Impact==ii,c(1,grep("Moy",names(dnow.fig)))], lty=2, type="b")
      mtext(lab.sf[ii+1],adj=0)
      draw.SE(selist[[ii+1]])

    }

    dmm <- sapply(c(0,1),sub.fig) # draw two frames

    # Rajouter la déco:
    legend("topright",c("Especes cible","Autres especes"),lty=1:2,cex=1,
           inset=c(0,-0.25),bty="n",xpd=NA)
    mtext(paste(wgeo,":",sep=""),outer=TRUE,cex=1.5, adj=0)
    mtext(yx.lab[gr.lab==wtype],side=2,outer=TRUE,line=2, cex=1.1)

    if(save) {
    dev.copy2pdf(file=paste(fig.dir,"GS_PoissonsST_ParCibleNC_",paste(wgeo,wtype,sep="."),"_",Sys.Date(),".pdf",sep=""))

  }
  }


  sapply(gr.lab[1:3], function(gg) sapply(names(dd.split), function(x) graph.geo(x, gg)))

  return(dd.split)
}


###################################################################################
##################################################################################

if(!(exists("ds.calc"))) ds.calc <- tableau.brut()

if(!(exists("BDtable"))) { # fix this
  BDtable <- BD.by.sp()
  facteurs <- c("St","Zone.Impact","Geomorpho","Code.SP",
                "Peche","Groupe.Trophique","Groupe.Mobil","Cible")
  ttble <- unique(ds.calc[,c("Campagne.ID",facteurs)])
  names(ttble)[1] <- "Campagne"
  BDtable.i <- merge(BDtable,ttble,c("Campagne","St","Code.SP"))
  BDtable <- merge(BDtable.i,bioeco.all[,c("Code.SP","fmlabel")],
                   by="Code.SP")
  }

# Séries temporelles pour toutes les espèces, et par tous les types
if(!exists("TS1")) TS1 <- poissons.ts1()
if(!exists("TS2")) TS2 <- poissons.ts2()

# Tableau index pour catégories figures
graph.key <- data.frame("df.id"=c("Tot","Comm","Non.Comm","Cible.Oui","Cible.Non",
                          "Carn","Pisc","Herb","Planct",
                          "Chaet","Pom"),
                        "titre.fig"=c("Toutes especes",
                          "Especes commerciales","Especes non-commerciales",
                          "Especes ciblees", "Especes non-ciblees",
                          "Carnivores","Piscivores","Herbivores","Planctivores",
                          "Chaetodontidae","Pomacentridae"))

############################################
############################################
# Pour re-sauver tout

save.all <- FALSE
if(save.all) {
  ds.calc <- tableau.brut(save=TRUE)
  TS1 <- poissons.ts1(save=TRUE)
  TS2 <- poissons.ts2(save=TRUE)
  poissons.p3(save=TRUE)
}


