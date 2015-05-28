## GS_Codes-graphiques.r
## Codes pour faire les graphiques à partir de la librairie ggplot2
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: March 15, 2015
## Time-stamp: <2015-05-28 12:10:45 Laura>

## Graphiques à produire:
## one variable: barplot
## two variables: single fig -- two colours; panel fig -- one var per panel
## three variables: panel fig -- two var per panel

## To-do: aes to change space between bars, legend title

## dat.stat.inv: see function INV.stats
require(ggplot2)
fig.etiq <- c(Geomorpho="Géomorphologie", N_Impact="Niveau d'impact",
              dens="Densité moyenne",
              Campagne="Campagne",
              Saison="Saison", St="Station")

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
                       "PI2", "PK4", "PK5", "PO1", "PROC04", "PROC4", "PV3B", "SOP01",
                       "SOP1", "SOP2", "ST1", "ST2", "ST3", "ST4", "ST5", "ST6"))

# Oter les campagnes au besoin si filtre sur campagne A ou S
if(!is.na(fCampagne)) xfact.levels$Campagne <- grep(fCampagne, xfact.levels$Campagne, value=TRUE)

###### ###### ###### ###### ###### ###### ###### ###### ######
prep.var <- function(wvar, df) {

  df <- data.frame(df)
  if(wvar=="Geomorpho.abbrev") col <- geomorpho.lab[df$Geomorpho]
  else col <- df[,wvar]
  if(!(wvar %in% c("Groupe","S_Groupe","Famille","Genre"))) {
    factor(col, levels=xfact.levels[[wvar]], ordered=TRUE)
  } else { factor(col) }
}

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

  dat.plot$vx <- prep.var(var1, dat.plot)
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
      theme(plot.margin=unit(c(0.25,0.25,0.35,0.35),"in"),
            axis.title.x=element_text(vjust=-1), axis.title.y=element_text(vjust=2)) +
              guides(fill=guide_legend(title=var1))
p1
}

###### ###### ###### ###### ###### ###### ###### ###### ######
fig.2var <- function(var1="Geomorpho", var2="Campagne",
                     var3=NULL, var.expl="dens", filtre,
                     filtre2,
                     agtaxo="Groupe", typ.taxo="Crustaces",
                     type.fig="Pas.Panneau",
                     dat=dat.stat.inv, tous.niveaux=TRUE) {

  dat <- INV.dens.gnrl(fspat=c(facteurs.spatio,"St"),
                      ftemp=c(facteurs.tempo, "Campagne"),
                             agtaxo=agtaxo, par.transect=TRUE,
                                silent=TRUE)

  dat <- data.frame(dat)
  if(!missing(typ.taxo) & typ.taxo!="tous") dat <- filter_(dat, paste(agtaxo,"%in%", typ.taxo))
  if(!missing(filtre)) dat <- filter_(dat, filtre)
  if(!missing(filtre2)) dat <- filter_(dat, filtre2)

  dat$var.expl <- dat[,var.expl]

  dat.plot <- dat %>% s_group_by(var1, var2, var3) %>%
    summarize(mean = mean(var.expl, na.rm=TRUE), sd=sd(var.expl, na.rm=TRUE))

  dat.plot$vx <- prep.var(var1, dat.plot)
  dat.plot$vy <- prep.var(var2, dat.plot)
  if(!missing(var3)) dat.plot$vz <- prep.var(var3, dat.plot)
  dat.plot$v.expl <- dat.plot$mean

  if(tous.niveaux & missing(filtre)) {
  all.lev <- expand.grid(var1=levels(dat.plot$vx), var2=levels(dat.plot$vy))
  names(all.lev) <- c("vx","vy")
  dat.plot <- left_join(all.lev, dat.plot)
}



  # defining data and main aes/mapping
  # see also aes_string maybe
  p1 <- ggplot(data=dat.plot, aes(x=vx, y=v.expl,
                 ymin=v.expl-sd, ymax=v.expl+sd,
                 fill=vy))

  #p1 <- ggplot(data=dat.plot, aes(x=vy, y=v.expl,
   #              ymin=v.expl-sd, ymax=v.expl+sd,
    #             fill=vy))

  # defining layers
  p1 <- p1 + geom_bar(stat="identity", width=0.5,
                      position=position_dodge(0.95))

  if(type.fig=="Panneau") {
    if(missing(var3)) {
      p1 <- p1 + aes(x=vy, y=v.expl, fill=vy) + facet_wrap(~ vx)#, scales="free_x")
    } else {
      p1 <- p1 + aes(x=vx, y=v.expl, fill=vy) #+ facet_wrap(~ vx)
      p1 <- p1 + facet_wrap(~ vz, scales="free_x")
    }
}
  # changing display parameters
  p1 <- p1 +
    xlab(fig.etiq[var1]) + ylab(fig.etiq[var.expl]) + theme_bw() +
      theme(plot.margin=unit(c(0.25,0.25,0.35,0.35),"in"),
            axis.title.x=element_text(vjust=-1), axis.title.y=element_text(vjust=2)) +
              guides(fill=guide_legend(title=var2))

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
verti.x.val <- theme(axis.text.x = element_text(angle = 90))
no.x.val <- theme(axis.text.x=element_blank())
