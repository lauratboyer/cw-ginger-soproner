## GS_Selection-Donnees-Par-Projet.r
## Import des tableaux d'échantillonages à partir du projet
## de base, e.g. KNS
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: November 25, 2014
## Time-stamp: <2015-01-20 15:37:59 Laura>

nom.projets <- function() {
  message("KNS_KONIAMBO\nSLN_THIO\nSLN_KAALA\nSLN_NEPOUI\nSMSP_BORENDI")
}
selection.projet <- function(x="KNS_KONIAMBO") {
  if(!(x %in% c("ADECAL_TOUHO","KNS_KONIAMBO","SLN_THIO",
              "SLN_KAALA","SLN_NEPOUI","SMSP_BORENDI"))) {
    stop("\n\nLe nom du projet devrait être l'un des suivants:
\nKNS_KONIAMBO\nSLN_THIO\nSLN_KAALA\nSLN_NEPOUI\nSMSP_BORENDI")

  }
  projet <<- x
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

  info.transect <<- info.transect.TProj[info.transect.TProj$Projet == x,]
  dbio <<- dbio.TProj[dbio.TProj$Projet == x,]
  dpoissons <<- dpoissons.TProj[dpoissons.TProj$Projet == x,]
  dLIT <<- data.LIT.TProj[data.LIT.TProj$Projet == x,]


}


