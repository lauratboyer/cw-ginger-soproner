## GS_Selection-Donnees-Par-Projet.r
## Import des tableaux d'échantillonages à partir du projet
## de base, e.g. KNS
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: November 25, 2014
## Time-stamp: <2015-05-28 12:12:02 Laura>

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

  if(x=="ADECAL_TOUHO") message("Pas de données INV disponibles")
}


