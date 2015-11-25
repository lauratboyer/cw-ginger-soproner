## GS_Selection-Donnees-Par-Projet.r
## Import des tableaux d'échantillonages à partir du projet
## de base, e.g. KNS
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: November 25, 2014
## Time-stamp: <2015-05-29 14:42:44 Laura>

fCampagne <- NA # filtre appliqué sur les campagnes? modifié par selection.projet()

nom.projets <- function() {
  message("ADECAL_TOUHO\nKNS_KONIAMBO\nSLN_THIO\nSLN_KAALA\nSLN_NEPOUI\nSMSP_BORENDI\n")
}
selection.projet <- function(x="KNS_KONIAMBO", filtre.Campagne) {
  if(!(x %in% c("ADECAL_TOUHO","KNS_KONIAMBO","SLN_THIO",
              "SLN_KAALA","SLN_NEPOUI","SMSP_BORENDI"))) {
    stop("\n\nLe nom du projet devrait être l'un des suivants:
\nADECAL_TOUHO\nKNS_KONIAMBO\nSLN_THIO\nSLN_KAALA\nSLN_NEPOUI\nSMSP_BORENDI")

  }
  projet <<- x
  fcamp <- !missing(filtre.Campagne)
  message(sprintf("Projet choisi pour l'analyse: %s",x))

  # définir tableaux filtrés

  # rajouter dossier avec le nom du projet si non-existant
  # dans dossier.R
  if(!file.exists(x)) {
    dir.create(paste(dossier.R, x, sep='//'))
    dir.create(paste(dossier.R, x, "Graphiques", sep='//'))
    dir.create(paste(dossier.R, x, "Tableaux", sep='//'))
  }
  # Définition lien fichiers pour sauvegarder graphiques/tableaux
  fig.dir <<- paste(dossier.R, x, "Graphiques/", sep='/')
  tabl.dir <<- paste(dossier.R, x, "Tableaux/", sep='/')

  filt.funk <- function(dat) {

    d1 <- filter(dat, Projet == x)
    # rajoute filtre sur campagne si spécifié
    if(fcamp) {
      d1 <- filter(d1, grepl(filtre.Campagne,Campagne))
      assign("fCampagne", filtre.Campagne, .GlobalEnv)
      } else { assign("fCampagne", NA, .GlobalEnv) }
    as.data.frame(d1) # switch back to data.frame
  }

  info.transect <<- filt.funk(info.transect.TProj)
  dbio <<- filt.funk(dbio.TProj)
  dpoissons <<- filt.funk(dpoissons.TProj)
  dLIT <<- filt.funk(data.LIT.TProj)
  dQuad <<- filt.funk(dQuad.TProj)
  # et on re-source le code graphique pour changer les valeurs
  # des niveaux de campagne, au besoin
  source.with.encoding("GS_Codes-graphiques.r",encoding="UTF-8")
  do.fig.tables() # faire tableaux de bases pour graphiques
}


do.fig.tables <- function() {

    message("");
    message("Création de tableaux pré-formattés pour faire les graphiques")
    message("------------------------------------------------------------")

    # Invertebres:
    inv.taxo <- c("Groupe","S_Groupe","Famille","Genre")
    message("");
    if(nrow(dbio) > 0) {
    message("Invertébrés:")
    dmm <- lapply(inv.taxo, function(i) {
                      obj <- INV.dens.gnrl(filt.camp="X", fspat=c(facteurs.spatio,"St"),
                                       ftemp=c(facteurs.tempo, "Campagne"),
                                       agtaxo=i, par.transect=TRUE,
                                           silent=TRUE);
                      objname <- paste0("figdat.inv.",i,".w0");
                      message(objname);
                      assign(objname, obj, .GlobalEnv)
                  })
    # Tableaux sans zeros
    dmm <- lapply(inv.taxo, function(i) {
                      obj <- INV.dens.gnrl(filt.camp="X", fspat=c(facteurs.spatio,"St"),
                                       ftemp=c(facteurs.tempo, "Campagne"),
                                           agtaxo=i, par.transect=TRUE, wZeroT=FALSE, silent=TRUE);
                       objname <- paste0("figdat.inv.",i);
                      message(objname);
                      assign(objname, obj, .GlobalEnv)
                  })
} else { message("Pas de données invertébrés") }
    if(nrow(dpoissons)>0) {
                                        # Poissons
    message("");     message("Poissons:")
    figdat.poissons.w0 <<- POIS.dens.gnrl(fspat=c(facteurs.spatio,"St"),
                                          agtaxo="Famille", ftemp=c(facteurs.tempo, "Campagne"),
                                          par.transect=TRUE, filt.camp="X", silent=TRUE)
    message("figdat.poissons.w0")
    figdat.poissons <<- POIS.dens.gnrl(fspat=c(facteurs.spatio,"St"),
                                       agtaxo="Famille", ftemp=c(facteurs.tempo, "Campagne"),
                                       wZeroT=FALSE, par.transect=TRUE, filt.camp="X", silent=TRUE)
    message("figdat.poissons")
}else { message("Pas de données pour les poissons") }
    if(nrow(dLIT)>0) {
                                        # LIT
    message("")
    message("LIT:")
    figdat.LIT.w0 <<- LIT.tableau.brut(filt.camp="X", silent=TRUE, frmt.ret="long")
    message("figdat.LIT.w0")
    figdat.LIT <<- LIT.tableau.brut(filt.camp="X", wZeroT=FALSE, silent=TRUE, frmt.ret="long") # sans zeros
    message("figdat.LIT")
    message("")
} else { message("Pas de données LIT")}

    if(nrow(dQuad)>0) {
    message("Quadrats:")
    figdat.Quad.w0 <<- Quad.tableau.brut(filt.camp="X", silent=TRUE, frmt.ret="long")
    message("figdat.Quad.w0")
    figdat.Quad <<- Quad.tableau.brut(filt.camp="X", wZeroT=FALSE, silent=TRUE, frmt.ret="long") # sans zeros
    message("figdat.Quad")
} else { message("Pas de données quadrats")}
}
