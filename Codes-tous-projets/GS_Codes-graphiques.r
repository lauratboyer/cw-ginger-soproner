## GS_Codes-graphiques.r
## Codes pour faire les graphiques à partir de la librairie ggplot2
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: March 15, 2015
## Time-stamp: <2015-05-29 08:17:32 Laura>

# missing 2014 for quadrats? fig.2var(typ.fig="boxplot", panneau="N_Impact")

## dat.stat.inv: see function INV.stats
require(ggplot2)

##################### Section THEME ########################
theme.GS <- theme_bw() + theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"in"),
                               axis.title.x=element_text(vjust=-1),
                               axis.title.y=element_text(vjust=2),
                               legend.key=element_blank(), strip.text=element_text(colour="white", size=12.5, vjust=0.35),
                               strip.background=element_rect(fill="slategray", colour=NA))

##################### Section COULEUR ######################
palette.de.couleurs <- "custom.colpal.a" # a definir entre guillemets
custom.colpal.a <- function(n=12) {
    colvect <- c("royalblue3","deepskyblue1","gold","orange1","indianred1","firebrick2","#AC2020")
    colorRampPalette(colvect)(n)
}
custom.colpal.b <- function(n=12) { # other option going from neutral to red
     colvect <- c("wheat3","wheat2","orange1","indianred1","firebrick2","#AC20202")
    colorRampPalette(colvect)(n)
}
custom.colpal.c <- function(n=12) { # other option going from neutral to red
    # <-- modifier 'colvect' pour avoir une nouvelle palette
    colvect <- c("wheat3","gold","seagreen2","turquoise3", "dodgerblue","navy")
    colorRampPalette(colvect)(n)
}
custom.colpal.d <- function(n=12, which.brewer.pal="Spectral") { # other option using RColorBrewer
    colvect <- brewer.pal(11, which.brewer.pal)
    # voir display.brewer.all() pour options de palettes toutes faites Brewer, et remplacer l'argument which.brewer.pal
    # e.g.
    # dev.new(); display.brewer.all(); brewer.pal(7,"Set1")
}
custom.colpal <- get(palette.de.couleurs) # ne pas changer: enregistre la palette a utiliser sous bas dans le code

##################### Section ETIQUETTES ######################
fig.etiq <- c(Geomorpho="Géomorphologie", Geomorpho.abbrev="Géomorphologie", N_Impact="Niveau d'impact",
              dens="Densité moyenne",
              RS="Richesse spécifique",
              taille.moy="Taille moyenne",
              biomasse="Biomasse",
              pcouv="Couverture en %",
              Campagne="Campagne",
              Saison="Saison", St="Station",
              taille.moy="Taille moyenne (cm)",
              Coraux="Coraux (couverture en %)",
              Abiotique="Abiotique (couverture en %)",
              Coraline="Coraline (couverture en %)",
              Acroporidae="Acroporidae (couverture en %)",
              Non.Acroporidae="Non-Acroporidae (couverture en %)",
              Corail.digite="Corail digité (couverture en %)",
              Corail.branchu="Corail branchu (couverture en %)",
              Corail.massif="Corail massif (couverture en %)",
              Corail.encroutant="Corail encroutant (couverture en %)",
              Corail.foliaire="Corail foliaire (couverture en %)",
              Corail.sub.Massif="Corail sub-Massif (couverture en %)",
              Corail.tabulaire="Corail tabulaire (couverture en %)")

geomorpho.lab <- c("Recif barriere externe"="RBE",
                   "Recif barriere interne"="RBI",
                   "Recif frangeant"="RF",
                   "Recif reticule"="RR",
                   "Passe"="P")

## Niveau des variables catégoriques dans l'ordre désiré pour les graphiques
xfact.levels <- list(Campagne=c("A_2006", "S_2007", "A_2007", "S_2008", "A_2008",
                     "S_2009","A_2009", "S_2010", "A_2010", "S_2011", "A_2011",
                         "S_2012", "A_2012","S_2013", "A_2013", "S_2014", "A_2014",
                                "S_2015","A_2015"),
                     N_Impact=c("Reference","Impact"),
                     Geomorpho=c("Recif barriere externe",
                       "Recif barriere interne",
                       "Recif frangeant",
                       "Recif reticule", "Passe"),
                     Geomorpho.abbrev=c("RBE","RBI","RF","RR","P"),
                     Saison=c("Froide","Inter","Chaude"),
                     St=c("CLC2", "FR1", "FR2", "FR3", "IBR1", "IBR2", "IBR2B", "IBR3",
                       "IRD01", "IRD05", "IRD06", "IRD12", "IRD15", "IRD23", "IRD24",
                       "IRD28", "IRD37", "IRD41", "IRD44", "IRD53", "IRD57", "KNS01",
                       "KNS02", "KNS04", "KNS05", "KNS06", "KNS31", "KO4", "LG1B", "LG2B",
                       "LG5B", "LG6B", "LG7B", "LR1", "LR2", "LR3", "OBRK3", "OBRK4",
                       "OBRK5A", "OF1C", "OF2C", "OF3C", "OF4C", "OF9", "PD4B", "PD5B",
                       "PI2", "PK4", "PK5", "PO1", "PROC4", "PV3B",
                         "SOP1", "SOP2", "ST1", "ST2", "ST3", "ST4", "ST5", "ST6", "CAGES", "REF"))
xfact.levels$LIT.lev <- lapply(index.LIT, function(x) na.omit(unique(x)))

# avec accents pour étiquette légende et axes seulements
xfact.levels.accents <- xfact.levels
xfact.levels.accents$N_Impact <- c(Reference="Référence", Impact="Impact")
xfact.levels.accents$Geomorpho <- factor(c("Récif barrière externe", "Récif barrière interne",
                                           "Récif frangeant", "Récif réticulé", "Passe"),
                                         levels=c("Récif barrière externe", "Récif barrière interne",
                                    "Récif frangeant", "Récif réticulé", "Passe"), ordered=TRUE)
names(xfact.levels.accents$Geomorpho) <- xfact.levels$Geomorpho

# Oter les campagnes au besoin si filtre sur campagne A ou S
if(!is.na(fCampagne)) xfact.levels$Campagne <- grep(fCampagne, xfact.levels$Campagne, value=TRUE)

########################################################################
########################## Section CODE DE SUPPORT #####################
########################################################################
if(!exists("bio.fig")) bio.fig <<- "inv"
fig.catego <- function(x) {
    if(missing(x)) {
        stop("L'une des valeurs suivantes devrait être spécifiée en argument: LIT/lit, quadrats, poissons, INV/inv")
   }else if (!(x %in% c("LIT","lit","Quad","quadrats","poissons","INV","inv"))) {
    stop("L'une des valeurs suivantes devrait être spécifiée en argument: LIT/lit, quadrats, poissons, INV/inv")
   }
    bio.fig <<- tolower(x)
}
########################################################################
## prep.var() formatte les variables et attribue les niveaux selon xfact.levels ci-haut
prep.var <- function(wvar, df, filt.camp="X", agLIT=NULL) {

    df <- data.frame(df)
    if(wvar=="Geomorpho.abbrev") col <- geomorpho.lab[df$Geomorpho]
  else col <- df[,wvar]
  if(!(wvar %in% c("Groupe","S_Groupe","Famille","Genre"))) {
      all.levs <- xfact.levels[[wvar]]
      if(wvar == "LIT.lev") all.levs <- all.levs[[agLIT]] # niveaux pour la categorie LIT specifiee seulement
      if(wvar == "Campagne" & filt.camp %in% c("A","S")) all.levs <- grepv(filt.camp, all.levs)
      if(wvar=="Geomorpho.abbrev") col <- geomorpho.lab[df$Geomorpho]
      factor(col, levels=all.levs, ordered=TRUE)
  } else { factor(col) }
}

##############################################################
# Barres d'erreur superieure: moyenne + ecart type
sd.top <- function(x) mean(x, na.rm=TRUE) + sd(x, na.rm=TRUE)
# Barres d'erreur inferieure: moyenne - ecart type, ou zero
sd.bot <- function(x) {
    rval <- mean(x, na.rm=TRUE) - sd(x, na.rm=TRUE)
    sapply(rval, max, 0) # mettre a 0 si valeur negative
}
##############################################################
# Fonctions raccourcies pour definir les filtres selon le format requis par dplyr
# a donner en argument a filtre et filtre2
# (peuvent s'utiliser directement dans la console pour voir la syntaxe requise)
filtre.incl <- function(a, b) { # inclusion des niveaux 'b' de la variable 'a'

    b <- paste0("'",b,"'", collapse=",") # format b pour inclusion dans c('...','...') avec guillemets simples
    sprintf("%s %%in%% c(%s)", a, b)
}
filtre.excl <- function(a, b) { # inclusion des niveaux 'b' de la variable 'a'

    b <- paste0("'",b,"'", collapse=",") # format b pour inclusion dans c('...','...') avec guillemets simples
    sprintf("!(%s %%in%% c(%s))", a, b)
}

##############################################################
############ Section CODES GRAPHIQUES ########################
##############################################################

# fig.2var() est la fonction principale utilisee pour faire les graphiques
fig.2var <- function(var1="Geomorpho", var2="Campagne",
                     var.expl, filtre,
                     filtre2, filtre.camp="X",
                     agLIT="General", agtaxo="Groupe", typ.taxo="Crustaces",
                     groupe=NULL, panneau=NULL,
                     dat, tous.niveaux=TRUE,
                     wZeroT=TRUE,
                     typ.fig="barre", panel.scale="free_x", ...) {


    # Définition de la variable de réponse
    var.expl.defaut <- c(inv="dens", poissons="dens", LIT="pcouv", lit="pcouv", Quad="pcouv")
    if(missing(var.expl)) var.expl <- var.expl.defaut[bio.fig] # si non spécifiée, utilise la valeur par défaut

    ##################################################
    # Extraction des données selon le type d'organisme
    ntbl <- paste0("figdat.", bio.fig, ifelse(bio.fig=="inv", paste0(".",agtaxo), ""),
                   ifelse(wZeroT, ".w0", "")) # re-creer le nom de l'objet contenant les donnes necessaires
    if(missing(dat)) {
        message(sprintf("Tableau utilisé: %s", ntbl))
        dat <- get(ntbl) %>% data.frame
    } #else { dat %<>% data.frame}
    #################################################
    # pour les LIT/Quad

    if(bio.fig %in% c("LIT","Quad","Quadrat")) dat <- filter(dat, LIT.cat == agLIT)

    ##################################################
    # filtre les donnees par campagne
    if(bio.fig=="inv") filtre.camp <- paste("T",filtre.camp,"inv",sep="_") # formatter nom du filtre pour invertebres
    dat <- filtreTable(dat, filtre.camp)

    ##################################################
    # Filtre sur les données selon les arguments
    if(!missing(typ.taxo) & typ.taxo[1]!="tous" & bio.fig == "inv") dat <- filter_(dat, paste0(agtaxo," %in% c(", paste0("'",typ.taxo,"'",collapse=", "), ")")) # invertébrés
    if(!missing(filtre)) dat <- filter_(dat, filtre) # voir aussi fonction filtre.incl() et filtre.excl()
    if(!missing(filtre2)) dat <- filter_(dat, filtre2)
    

    ##################################################
    var.expl <- gsub("[ -]", ".", var.expl)
    dat$var.expl <- dat[,var.expl] # define response variable for graph (Y)
    dat.plot <- dat

    dat.plot$vx <- prep.var(var1, dat.plot, filt.camp=filtre.camp, agLIT=agLIT) # variable en X
    dat.plot$vy <- prep.var(var2, dat.plot, filt.camp=filtre.camp, agLIT=agLIT) # variable explicative en Y
    if(!is.null(panneau)) dat.plot$panneau <- prep.var(panneau, dat.plot, filt.camp=filtre.camp, agLIT=agLIT)
    # rajouter rangees manquantes si tous.niveaux=TRUE
    if(tous.niveaux & missing(filtre)) {
        if(is.null(panneau)) all.lev <- expand.grid(vx=levels(dat.plot$vx), vy=levels(dat.plot$vy))
        else all.lev <- expand.grid(vx=levels(dat.plot$vx), vy=levels(dat.plot$vy), panneau=levels(dat.plot$panneau))
        dat.plot <- left_join(all.lev, dat.plot)
        }

    dat.plot <<- dat.plot # rajouter dat.plot dans l'environnement global pour acces au besoin

    # defining data and main aes/mapping
                                        # (could also use aes_string...?)

    if(var1 %in% c("Geomorpho","N_Impact")) {
        dat.plot$vx <- xfact.levels.accents[[var1]][dat.plot$vx]
    }
    if(var2 %in% c("Geomorpho","N_Impact")) { # ajuster la legende pour variables avec accents
        vy.labs <- xfact.levels.accents[[var2]][sort(unique(dat.plot$vy[!is.na(dat.plot$var.expl)]))]
    }else{vy.labs <- sort(unique(dat.plot$vy[!is.na(dat.plot$var.expl)])) }
    colv <- custom.colpal(length(vy.labs))
    p0 <- ggplot(data=dat.plot) +
        scale_fill_manual(labels=vy.labs, values=colv) +
            scale_colour_manual(labels=vy.labs, values=colv)
    pdod <- position_dodge(width=0.95)

    # type de graphique 1
    # barres = geom_bar
    if(typ.fig=="barre") {
      p1 <- p0 + aes(x=vx, y=var.expl, fill=vy)  +
          stat_summary(fun.ymin=sd.bot, fun.ymax=sd.top, geom="linerange", position=position_dodge(0.95)) +
              stat_summary(fun.y=mean, geom="bar", width=0.5, position=position_dodge(0.95))


    # type de graphique 2
    # fonction geom_boxplot
  } else if(typ.fig=="boxplot") {
                                        # calcul des parametres de BP via boxplot car plus rapide que ggplot
      if(is.null(panneau)) {
          bp.raw <- boxplot(as.numeric(var.expl) ~ vx*vy, dat.plot, plot=FALSE)
      }else{
          bp.raw <- boxplot(as.numeric(var.expl) ~ vx*vy*panneau, dat.plot, plot=FALSE)
      }
      bp.df <- t(bp.raw$stats) %>% data.frame
      names(bp.df) <- c("ymin","low","mid","top","ymax")
      bp.df$fact <- bp.raw$names
      lab.string <- ifelse(is.null(panneau),"(.*)\\.(.*)","(.*)\\.(.*)\\.(.*)")
      bp.df$vx <- gsub(lab.string,"\\1",bp.df$fact)
      bp.df$vy <- gsub(lab.string,"\\2",bp.df$fact)
      if(!is.null(panneau)) bp.df$panneau <- gsub(lab.string,"\\3",bp.df$fact)

      bp.df <- filter(bp.df, !is.na(mid))
      yl <- quantile(dat.plot$var.expl, 1, na.rm=TRUE)
      p1 <- ggplot(data=bp.df, aes(x=vx, ymin=ymin, lower=low, middle=mid, upper=top, ymax=ymax, fill=vy)) +
          scale_fill_manual(labels=vy.labs, values=colv) +
              ylim(0,yl) + geom_boxplot(stat="identity",
                                        color="grey", size=0.25,
                                        position=position_dodge(0.5), width=0.5,
                                        outlier.colour="royalblue3", outlier.size=1.5, alpha=0.95)

    # type de graphique 3
    #
    } else if (typ.fig=="ligne") {
      if(var1 != var2) p1 <- p0 + aes(x=vx, y=var.expl, group=vy, colour=vy) 
      if(var1 == var2) p1 <- p0 + aes(x=vx, y=var.expl, group=1, colour=vy) 
          p1 <- p1 +
          stat_summary(fun.y=mean, geom="line", position=position_dodge(0.2), size=0.6) + 
          stat_summary(fun.y=mean, fun.ymin=sd.bot, fun.ymax=sd.top, geom="pointrange", alpha=0.75, position=position_dodge(0.2))
              #stat_summary(fun.y=mean, geom="line", position=position_dodge(0.2), size=0.6)


        # type de graphique 4
        }else{

                                        #p1 <- p0 + aes(group=vy, ymin=sd.bot, ymax=sd.top, colour=vy) +
            #geom_errorbar(width=0.5, position=pdod) +
                                        #   geom_line(position=pdod, lwd=1, ...) + geom_point(position=pdod)
            message("\n--------------------
Erreur! Spécifiez typ.fig = 'barre', 'boxplot', ou 'ligne'"); stop()
        }


    if(!is.null(panneau)) p1 <- p1  + facet_wrap(~ panneau, scales=panel.scale)

    # Paramètres visuels (axes, etc.)
    ############################################
    # fonction theme(...) et, pour la légende, guides(...)
    p1 <- p1 +
        xlab(fig.etiq[var1]) + ylab(fig.etiq[var.expl]) + theme.GS  +
            guides(fill=guide_legend(title=var2), col=guide_legend(title=var2))

    # Fini!!!
    #########################################
    p1 # objet graphique ggplot
}


######################################################################
######################################################################
# Fonctions raccourcis pour fig.2var par type de données
inv.fig <- function(...) { bio.fig <<- "inv"; fig.2var(...)}
poissons.fig <- function(...) {bio.fig <<- "poissons"; fig.2var(...)}
LIT.fig <- function(...) {bio.fig <<- "LIT"; fig.2var(...)}
Quad.fig <- function(...) {bio.fig <<- "Quad"; fig.2var(...)}


######################################################################
# Raccourcis pour fonctions ggplot:
verti.x.val <- function(angl=90, size=1, ...) theme(axis.text.x = element_text(angle = angl, size=size, ...))
no.x.val <- theme(axis.text.x=element_blank())

######################################################################
# Raccourci pour sauvegarder directement dans le dossier du projet
ggsave.proj <- function(...) {
    cdir <- getwd()
    setwd(fig.dir)
    ggsave(...)
    setwd(cdir)
}
