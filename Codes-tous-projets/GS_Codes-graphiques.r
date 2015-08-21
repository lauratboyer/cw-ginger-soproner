## GS_Codes-graphiques.r
## Codes pour faire les graphiques à partir de la librairie ggplot2
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: March 15, 2015
## Time-stamp: <2015-05-29 08:17:32 Laura>

## Graphiques à produire:
## one variable: barplot
## two variables: single fig -- two colours; panel fig -- one var per panel
## add filtre campagne only
## three variables: panel fig -- two var per panel

## To-do: aes to change space between bars, legend title
message("rajouter nouvel LIT de Tom Aug 21st")
## dat.stat.inv: see function INV.stats
require(ggplot2)
fig.etiq <- c(Geomorpho="Géomorphologie", N_Impact="Niveau d'impact",
              dens="Densité moyenne",
              Campagne="Campagne",
              Saison="Saison", St="Station",
              taille.moy="Taille moyenne (cm)")

geomorpho.lab <- c("Recif barriere externe"="RBE",
                   "Recif barriere interne"="RBI",
                   "Recif frangeant"="RF",
                   "Recif reticule"="RR",
                   "Passe"="PS")

## Niveau des variables catégoriques dans l'ordre désiré pour les graphiques
xfact.levels <- list(Campagne=c("A_2006", "S_2007", "A_2007", "S_2008", "A_2008",
                     "S_2009","A_2009", "S_2010", "A_2010", "S_2011", "A_2011",
                     "S_2012", "A_2012","S_2013", "A_2013", "S_2014", "A_2014"),
                     N_Impact=c("Reference","Impact"),
                     Geomorpho=c("Recif barriere externe",
                       "Recif barriere interne",
                       "Recif frangeant",
                       "Recif reticule", "Passe"),
                     Geomorpho.abbrev=c("RBE","RBI","RF","RR","PS"),
                     Saison=c("Inter","Chaude"),
                     St=c("CLC2", "FR1", "FR2", "FR3", "IBR1", "IBR2", "IBR2B", "IBR3",
                       "IRD01", "IRD05", "IRD06", "IRD12", "IRD15", "IRD23", "IRD24",
                       "IRD28", "IRD37", "IRD41", "IRD44", "IRD53", "IRD57", "KNS01",
                       "KNS02", "KNS04", "KNS05", "KNS06", "KNS31", "KO4", "LG1B", "LG2B",
                       "LG5B", "LG6B", "LG7B", "LR1", "LR2", "LR3", "OBRK3", "OBRK4",
                       "OBRK5A", "OF1C", "OF2C", "OF3C", "OF4C", "OF9", "PD4B", "PD5B",
                       "PI2", "PK4", "PK5", "PO1", "PROC4", "PV3B",
                       "SOP1", "SOP2", "ST1", "ST2", "ST3", "ST4", "ST5", "ST6"))

# Oter les campagnes au besoin si filtre sur campagne A ou S
if(!is.na(fCampagne)) xfact.levels$Campagne <- grep(fCampagne, xfact.levels$Campagne, value=TRUE)

###### ###### ###### ###### ###### ###### ###### ###### ######
if(!exists("bio.fig")) bio.fig <<- "inv"
fig.catego <- function(x) {
    if(missing(x)) {
        stop("L'une des valeurs suivantes devrait être spécifiée en argument: LIT/lit, quadrats, poissons, INV/inv")
   }else if (!(x %in% c("LIT","lit","quadrats","poissons","INV","inv"))) {
    stop("L'une des valeurs suivantes devrait être spécifiée en argument: LIT/lit, quadrats, poissons, INV/inv")
   }

    bio.fig <<- tolower(x)

    }

###### ###### ###### ###### ###### ###### ###### ###### ######
prep.var <- function(wvar, df, filt.camp="X") {

  df <- data.frame(df)
  if(wvar=="Geomorpho.abbrev") col <- geomorpho.lab[df$Geomorpho]
  else col <- df[,wvar]
  if(!(wvar %in% c("Groupe","S_Groupe","Famille","Genre"))) {
      all.levs <- xfact.levels[[wvar]]
      if(wvar == "Campagne" & filt.camp %in% c("A","S")) all.levs <- grepv(filt.camp, all.levs)
    factor(col, levels=all.levs, ordered=TRUE)
  } else { factor(col) }
}

##############################################################
# fonction pour definir les filtres selon le format requis par dplyr
# a donner en argument a filtre et filtre2
filtre.incl <- function(a, b) { # inclusion des niveaux 'b' de la variable 'a'

    b <- paste0("'",b,"'", collapse=",") # format b pour inclusion dans c('...','...') avec guillemets simples
    sprintf("%s %%in%% c(%s)", a, b)
}
filtre.excl <- function(a, b) { # inclusion des niveaux 'b' de la variable 'a'

    b <- paste0("'",b,"'", collapse=",") # format b pour inclusion dans c('...','...') avec guillemets simples
    sprintf("!(%s %%in%% c(%s))", a, b)
}


##############################################################
############ DEFINITION FONCTIONS GRAPHIQUES #################
###### ###### ###### ###### ###### ###### ###### ###### ######
fig.1var <- function(var1="Geomorpho", var.expl="dens", filtre,
                     agtaxo="Groupe", typ.taxo="Crustaces",
                     type.fig="Pas.Panneau",
                     dat=dat.stat.inv, tous.niveaux=TRUE) {

  dat <- INV.dens.gnrl(fspat=c(facteurs.spatio,"St"),
                      ftemp=c(facteurs.tempo, "Campagne"),
                             agtaxo=agtaxo, par.transect=TRUE,
                                silent=TRUE)

  dat <- data.frame(dat)
  if(!missing(typ.taxo) & !(typ.taxo=="tous")) dat <- filter_(dat, paste(agtaxo,"%in%", typ.taxo))
  if(!missing(filtre)) dat <- filter_(dat, filtre)

  dat$var.expl <- dat[,var.expl]

  dat.plot <- dat %>% s_group_by(var1) %>%
    summarize(mean = mean(var.expl, na.rm=TRUE), sd=sd(var.expl, na.rm=TRUE))

  dat.plot$vx <- prep.var(var1, dat.plot, filt.camp=filtre.camp)
  dat.plot$v.expl <- dat.plot$mean

  if(tous.niveaux & missing(filtre)) {
  all.lev <- expand.grid(var1=levels(dat.plot$vx))
  names(all.lev) <- c("vx")
  dat.plot <- left_join(all.lev, dat.plot)
}

  # defining data and main aes/mapping
  # see also aes_string maybe
  p1 <- ggplot(data=dat.plot, aes(x=vx, y=v.expl,
                 ymin=v.expl-sd, ymax=v.expl+sd))

  # defining layers
  p1 <- p1 + geom_bar(stat="identity", width=0.8,
                      position=position_dodge(0.95),
                      fill="dodgerblue2", colour="dodgerblue4")

  # changing display parameters
  p1 <- p1 +
    xlab(fig.etiq[var1]) + ylab(fig.etiq[var.expl]) + theme_bw() +
      theme(plot.margin=unit(c(0.35,0.25,0.15,0.15),"in"),
            axis.title.x=element_text(vjust=-1), axis.title.y=element_text(vjust=2)) +
              guides(fill=guide_legend(title=var1))
p1
}


###### ###### ###### ###### ###### ###### ###### ###### ######
fig.2var <- function(var1="Geomorpho", var2="Campagne",
                     var3=NULL, var.expl="dens", filtre,
                     filtre2, filtre.camp="A",
                     agtaxo="Groupe", typ.taxo="Crustaces",
                     groupe=NULL, panneau=NULL,
                     dat=dat.stat.inv, tous.niveaux=TRUE,
                     typ.fig="barre", ...) {


    var.expl.defaut <- c(inv="dens", poissons="dens", lit="Coraux")
    var.expl <- var.expl.defaut[bio.fig]

    ##################################################
                                        # extraction des données selon le type d'organisme
    start.timer()
    if(bio.fig=="inv"){
        dat <- INV.dens.gnrl(filt.camp=filtre.camp, fspat=c(facteurs.spatio,"St"),
                      ftemp=c(facteurs.tempo, "Campagne"),
                             agtaxo=agtaxo, par.transect=TRUE,
                             silent=TRUE)
    }
stop.timer()
    if(bio.fig=="lit")  dat <- LIT.tableau.brut(filt.camp=filtre.camp)
    if(bio.fig=="poissons")  dat <- POIS.dens.gnrl(fspat=c(facteurs.spatio,"St"),
           agtaxo="Famille", ftemp=c(facteurs.tempo, "Campagne"), par.transect=TRUE, filt.camp=filtre.camp)
    dat <- data.frame(dat)
    ##################################################
    if(!missing(typ.taxo) & typ.taxo!="tous" & type.fig=="inv") dat <- filter_(dat, paste(agtaxo,"%in%", typ.taxo))
    if(!missing(filtre)) dat <- filter_(dat, filtre) # voir aussi fonction filtre.incl() et filtre.excl()
    if(!missing(filtre2)) dat <- filter_(dat, filtre2)

    dat$var.expl <- dat[,var.expl] # define response variable for graph (Y)
    if(typ.fig!="boxplot") {
        dat.plot <- dat %>% s_group_by(var1, var2, var3, panneau) %>%
            summarize(mean = mean(var.expl, na.rm=TRUE), sd=sd(var.expl, na.rm=TRUE)) %>%
                mutate(be.top=mean + sd, be.bot=max(0, mean - sd))
        dat.plot$var.expl <- dat.plot$mean

    }else{
        dat.plot <- dat }

  dat.plot$vx <- prep.var(var1, dat.plot, filt.camp=filtre.camp)
  dat.plot$vy <- prep.var(var2, dat.plot, filt.camp=filtre.camp)
  if(!missing(var3)) dat.plot$vz <- prep.var(var3, dat.plot, filt.camp=filtre.camp)
  if(!is.null(panneau)) dat.plot$panneau <- prep.var(panneau, dat.plot, filt.camp=filtre.camp)

  if(tous.niveaux & missing(filtre)) {
  all.lev <- expand.grid(var1=levels(dat.plot$vx), var2=levels(dat.plot$vy))
  names(all.lev) <- c("vx","vy")
  dat.plot <- left_join(all.lev, dat.plot)
}

  # objet avec les données pour le graph
  raw.dat <<- dat
  fig.dat <<- dat.plot

  # defining data and main aes/mapping
  # see also aes_string maybe
  p0 <- ggplot(data=dat.plot, aes(x=vx, y=var.expl))

    pdod <- position_dodge(width=0.4)

  if(typ.fig=="barre") {
    p1 <- p0 + aes(fill=vy, ymin=var.expl-sd, ymax=var.expl+sd) + scale_fill_discrete("")
    p1 <- p1 + geom_bar(stat="identity", width=0.5,
                      position=position_dodge(0.95), ...)
} else if(typ.fig=="boxplot"){

    yl <- quantile(dat.plot$var.expl, 0.975)
    p1 <- p0 + ylim(0,yl) + geom_boxplot(aes(fill=vy), colour="grey50", width=1, outlier.colour="royalblue3", outlier.size=1.5, alpha=0.75) + facet_grid(~vy)
}else{
    p1 <- p0 + aes(group=vy, ymin=be.bot, ymax=be.top, colour=vy) + geom_errorbar(width=0.5, position=pdod) + geom_line(position=pdod, lwd=1, ...) + geom_point(position=pdod)}


  if(!is.null(panneau)) {
    if(missing(var3)) {
      p1 <- p1  + facet_wrap(~ panneau)#, scales="free_x")
    } else {
   #   p1 <- p1 + aes(x=vx, y=v.expl, fill=vy) #+ facet_wrap(~ vx)
      p1 <- p1 + facet_wrap(~ vz, scales="free_x")
    }
}
  # changing display parameters
  p1 <- p1 +
    xlab(fig.etiq[var1]) + ylab(fig.etiq[var.expl]) + theme_bw() +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"in"),
            axis.title.x=element_text(vjust=-1),
            axis.title.y=element_text(vjust=2)) +
              guides(fill=guide_legend(title=var2),
                     col=guide_legend(title=var2))

#  p1 + geom_point() #position=pd)
p1
}



fig.3var <- function(var1="Geomorpho", var2="Campagne", var3="N_Impact",
                     var.expl="dens", typ.taxo="Mollusques", type.fig="Panneau", dat=dat.stat.inv) {

  dat <- dat[dat$Groupe %in% typ.taxo,]
  dat.plot <- aggregate(list(dens=dat[,var.expl]), as.list(dat[,c(var1,var2,var3)]), mean)

  vx <- as.factor(dat.plot[,var1])
  vy <- factor(dat.plot[,var2],
                   levels=xfact.levels[[var2]])
  vz <- as.factor(dat.plot[,var3])
  dat.plot$vx <- vx
  dat.plot$vy <- vy
  dat.plot$vz <- vz
  dat.plot$v.expl <- dat.plot[,var.expl]

  # defining data and main aes/mapping
  p1 <- ggplot(data=dat.plot, aes(x=vx, y=v.expl,
                 fill=vy))

  # defining layers
  p1 <- p1 + geom_bar(stat="identity", width=0.4,
                      position= position_dodge(width=0.5))

  # changing display parameters
  p1 <- p1 +
    xlab(fig.etiq[var1]) + ylab(fig.etiq[var.expl]) + theme_bw() +
      theme(plot.margin=unit(c(0.25,0.25,0.35,0.35),"in"),
            axis.title.x=element_text(vjust=-1), axis.title.y=element_text(vjust=2)) +
              guides(fill=guide_legend(title=var2))

  if(type.fig=="Panneau") {
    p1 <- p1 + facet_wrap(~ vz, scales="free_x")
#    p1 <- p1 + facet_wrap(~ vz) }
  }

#  p1 + geom_point() #position=pd)

  p1
}

###################################
# Raccourcis pour fonctions ggplot:
#fig.2lign <- function(...) fig.2var(..., typ.fig="ligne")
verti.x.val <- function(angl=90) theme(axis.text.x = element_text(angle = angl))
#no.x.val <- theme(axis.text.x=element_blank())

#fig.2var() + geom_errorbar()
