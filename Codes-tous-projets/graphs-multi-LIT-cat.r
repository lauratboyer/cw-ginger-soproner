#####################
# VERSION 1

ab <- LIT.fig(var1="Campagne",panneau="LIT.lev", agLIT="General", filtre=filtre.incl("LIT.lev",c("Coraux", "Corail mort", "Coraux mous"))) +
    no.x.val + xlab("") + theme(legend.position="none")
ab2 <- LIT.fig(var1="Campagne",panneau="LIT.lev", agLIT="Acroporidae", filtre=filtre.incl("LIT.lev",c("Acroporidae", "Non-Acroporidae"))) +
    no.x.val + xlab("") + theme(legend.position="none")
ab3 <- LIT.fig(var1="Campagne",panneau="LIT.lev", agLIT="All", filtre=filtre.incl("LIT.lev","Macro-Algues")) + no.x.val + xlab("") + ylab("") + theme(legend.position="none")


# tu auras peut-etre besoin des librairies grid et gridExtra
vplayout <- function(x,y) viewport(layout.pos.row=x, layout.pos.col=y)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,3))) # layout de 2 rangees par 3 colonnes
print(ab, vp=vplayout(1,1:3)); print(ab2, vp=vplayout(2,1:2)); print(ab3, vp=vplayout(2,3))

############################################################################################
############################################################################################
############################################################################################
                                        # VERSION 2 (un peu plus long mais axe individuel pour chaque)
                                        # fais un graph pour chaque categorie individuelle
                                        # fonction raccourci pour pas avoir a repeter tout le temps
# verti.x.val defini l'angle des etiquettes, et la position par rapport a l'axe des x... legend.position="none" ote la legende
arg.extra <- function() verti.x.val(angl=45,size=8,vjust=0.5) + theme(legend.position="none", axis.title=element_blank())

general1 <- LIT.fig(var1="Campagne",panneau="LIT.lev", agLIT="General", filtre=filtre.incl("LIT.lev","Coraux")) + arg.extra()
general2 <- LIT.fig(var1="Campagne",panneau="LIT.lev", agLIT="General", filtre=filtre.incl("LIT.lev","Corail mort")) + arg.extra()
general3 <- LIT.fig(var1="Campagne",panneau="LIT.lev", agLIT="General", filtre=filtre.incl("LIT.lev","Coraux mous")) + arg.extra()
acro1 <- LIT.fig(var1="Campagne",panneau="LIT.lev", agLIT="Acroporidae", filtre=filtre.incl("LIT.lev","Acroporidae")) + arg.extra()
acro2 <- LIT.fig(var1="Campagne",panneau="LIT.lev", agLIT="Acroporidae", filtre=filtre.incl("LIT.lev","Non-Acroporidae")) + arg.extra()
all1 <- LIT.fig(var1="Campagne",panneau="LIT.lev", agLIT="All", filtre=filtre.incl("LIT.lev","Macro-Algues")) + arg.extra()

# ... et on place individuellement chaque objet a la position qu'on veut
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,3))) # layout de 2 rangees par 3 colonnes
print(general1, vp=vplayout(1,1)) # rangee 1, colonne 1
print(general2, vp=vplayout(1,2)) # rangee 1, colonne 2
print(general3, vp=vplayout(1,3))
print(acro1, vp=vplayout(2,1)) # rangee 2, colonne 1
print(acro2, vp=vplayout(2,2)) # rangee 2, colonne 2
print(all1, vp=vplayout(2,3)) # rangee 2, colonne 3


# pour rajouter un titre
#grid.text("Titre si voulu", y=unit(1, "npc") - unit(0.5, "lines"), gp=gpar(col="navy"))

# pour sauvegarder, utiliser dev.copy2pdf (pas ggsave)
dev.copy2pdf(file="test-LIT.pdf")

############################################################################################
############################################################################################
############################################################################################
vplayout <- function(x,y) viewport(layout.pos.row=x, layout.pos.col=y)
arg.extra <- function() verti.x.val(angl=45,size=8,vjust=0.5) + theme(legend.position="none", axis.title=element_blank())

RBE_general1 <- LIT.fig(filtre.camp="A", filtre2=filtre.incl("Geomorpho", "Recif reticule"), var1="Campagne", var2="St", panneau="LIT.lev", agLIT="General", filtre=filtre.incl("LIT.lev","Coraux"), typ.fig="ligne") + arg.extra()+ ylim(0,60)
RBE_general2 <- LIT.fig(filtre.camp="A", filtre2=filtre.incl("Geomorpho", "Recif reticule"), var1="Campagne",var2="St", panneau="LIT.lev", agLIT="General", filtre=filtre.incl("LIT.lev","Corail mort"), typ.fig="ligne") + arg.extra()+ ylim(0,60)
RBE_general3 <- LIT.fig(filtre.camp="A", filtre2=filtre.incl("Geomorpho", "Recif reticule"), var1="Campagne",var2="St", panneau="LIT.lev", agLIT="General", filtre=filtre.incl("LIT.lev","Coraux mous"), typ.fig="ligne") + arg.extra()+ ylim(0,60)
RBE_acro1 <- LIT.fig(filtre.camp="A", filtre2=filtre.incl("Geomorpho", "Recif reticule"), var1="Campagne",var2="St", panneau="LIT.lev", agLIT="Acroporidae", filtre=filtre.incl("LIT.lev","Acroporidae"), typ.fig="ligne") + arg.extra()+ ylim(0,60)
RBE_acro2 <- LIT.fig(filtre.camp="A", filtre2=filtre.incl("Geomorpho", "Recif reticule"), var1="Campagne",var2="St", panneau="LIT.lev", agLIT="Acroporidae", filtre=filtre.incl("LIT.lev","Non-Acroporidae"), typ.fig="ligne") + arg.extra()+ ylim(0,60)
RBE_all1 <- LIT.fig(filtre.camp="A", filtre2=filtre.incl("Geomorpho", "Recif reticule"), var1="Campagne",var2="St", panneau="LIT.lev", agLIT="All", filtre=filtre.incl("LIT.lev","Macro-Algues"), typ.fig="ligne") + arg.extra()+ ylim(0,60)

# ... et on place individuellement chaque objet a la position qu'on veut
grid.newpage()
pushViewport(viewport(layout=grid.layout(3,3))) # layout de 2 rangees par 3 colonnes
print(RBE_general1, vp=vplayout(1,1)) # rangee 1, colonne 1
print(RBE_general2, vp=vplayout(1,2)) # rangee 1, colonne 2
print(RBE_general3 + theme(legend.position="none"), vp=vplayout(1,3))
print(RBE_acro1, vp=vplayout(2,1)) # rangee 2, colonne 1
print(RBE_acro2, vp=vplayout(2,2)) # rangee 2, colonne 2
print(RBE_all1, vp=vplayout(2,3)) # rangee 2, colonne 3

grid_arrange_shared_legend <- function(..., ylab="Couverture (%)", legpos="right") {
  
  plots <- list(...) # graphiques ggplot faits au prealable
  g <- ggplotGrob(plots[[1]] + theme(legend.position=legpos, legend.title=element_blank()))$grobs 
  # on refait le premier graph en attrapant l'objet legende
  legend <- g[sapply(g, "[[", "name") == "guide-box"][[1]] # on trouve le grob contenant le nom 'guide-box'
 
  # (grob: 'gr'aphical 'ob'ject, nom des composantes d'un plot ggplot)
  lwidth <- sum(legend$width) # largeur de la legende
  # fait le plot avec la fonction grid.arrange
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 2, # une colonne pour les graphs principaux (assembles par arrangeGrob), une colonne pour la legende
    widths = unit.c(unit(1, "npc") - lwidth, lwidth),
    left=textGrob(ylab, rot = 90, vjust = 1))
}