## GS_Selection-Donnees-Par-Projet.r
## Import des tableaux d'échantillonages à partir du projet
## de base, e.g. KNS
## -------------------------------------------------------
## Author: Laura Tremblay-Boyer (l.boyer@fisheries.ubc.ca)
## Written on: November 25, 2014
## Time-stamp: <2014-12-11 07:59:18 Laura>

dossier.data <- "~/Projects/cw-ginger-soproner/Bases projet"
cw <- getwd()
setwd(dossier.data)
proj.fspat <- read.csv("Facteurs_spatiaux.csv", encoding="latin1")
proj.inv <- read.csv("Invertebres.csv")
proj.inv <- read.delim2("Invertebres-utf16.txt", encoding="utf-16")
proj.poissons <- read.csv("Poissons.csv", encoding="latin1")
# extraire nom du projet de la colonne ID
proj.inv$Projet <- gsub("([A-Z]*_[A-Z]*)_.*","\\1",proj.inv$ID)
proj.poissons$Projet <- gsub("([A-Z]*_[A-Z]*)_.*","\\1",proj.poissons$Id)

nom.projet <- function(x="KNS_KONIAMBO") {
  if(!(x %in% c("ADECAL_TOUHO","KNS_KONIAMBO","SLN_THIO",
              "SLN_KAALA","SLN_NEPOUI","SMSP_BORENDI"))) {
    stop("Le nom du projet devrait être: KNS_KONIAMBO,SLN_THIO,
              SLN_KAALA,SLN_NEPOUI,SMSP_BORENDI")
  }
  projet.main <<- x
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
  fig.dir <<- paste(dossier.R, x, "Graphiques", sep='')
  tabl.dir <<- paste(dossier.R, x, "Tableaux", sep='')

}

setwd(cw)
