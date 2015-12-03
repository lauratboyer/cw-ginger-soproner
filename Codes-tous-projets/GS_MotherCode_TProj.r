# Analyses des données KNS (Ginger/Soproner)
# Auteur: Laura Tremblay-Boyer, contact: l.boyer@fisheries.ubc.ca
# Time-stamp: <2015-12-01 15:03:05 lauratb>

################################################################
###### Définition des variables principales pour l'analyse #####
## Cette partie du code est à MODIFIER MANUELLEMENT au besoin ##
### *** N'oubliez pas d'enregistrer tous les changements ***####
Sys.setlocale("LC_ALL","fr_FR.UTF-8") # encodage pour les accents
message("à vérifier options pour fenêtre extérieure sur R Studio")
# Dossier contenant les codes et où les sorties vont être sauvegardées
dossier.R <- getwd() # le 'working directory',
# ... ou sinon mettre le nom du dossier désiré, e.g. C:/Documents/Codes_R
# Dossier de sauvegarde des fichiers .csv des bases de données
dossier.DB <- paste(getwd(), "/DBs tous projets", sep="")
setwd(dossier.R)
############################################
## Reformatter les tableaux ou charger ceux deja formattés? TRUE ou FALSE
## Pour importer les nouvelles données (mises sous .csv dans le dossier DBs tous projets), mettre à TRUE
refaire.tableaux <- FALSE

############################################
## Variables spatiales et temporelles dans les tableaux
## 'Facteurs_...' à mettre disponible pour les analyses
facteurs.spatio <- c("Geomorpho","N_Impact","Cote","Lieu")
facteurs.tempo <- c("Période.BACI","Saison")
facteurs.taxo <- list(INV=c("Groupe","S_Groupe","Famille","Genre","G_Sp"),
                       poissons=c("Famille","Genre","G_Sp","GTlabel","moblabel","Peche","Cible"),
                      LIT=c("General","Forme","Acroporidae","Sensibilite","Famille","Genre"))
facteurs.var.expl <- list(INV="dens",
                          poissons=c("dens","biomasse","RS","taille.moy"),
                          LIT="pcouv")
attr(facteurs.var.expl,"Note") <- "Le premier élement est utilisé par défaut dans fig.2var()"

############################################
## Valeurs par défaut pour l'aggrégation taxonomique,
## l'aggrégation spatiale, et l'aggrégation temporelle
## (i.e. les fonctions de calcul vont utiliser ces
## valeurs à moins qu'une valeur différente soit directement
## spécifiée dans l'argument de la fonction)
agtaxo.defaut <- "Groupe"
fspat.defaut <- "St"
ftempo.defaut <- "Campagne"

################################################################
### Contenu du code:
### 1. Filtre années
### 2. Filtre espèces
### 3. Définition de l'emplacement des fichiers
################################################################

############################
### 1. Filtre sur ANNEES ###
############################
filtre.annees <- 2006:2015 # indiquer quelles années à inclure dans l'analyse - tableau filtre généré automatiquement
                            # exemples de format:
                            # c(2006,2007,2009,2012), 2006:2012, seq(2006,2012,by=2)

#############################
### 2. Filtre sur ESPECES ###
#############################
## a. Filtre famille par défaut (spécifié dans Filtre-taxo_Famille.csv)
filtre.famille <- TRUE

#############################
## b. Général
filtre.sur.especes <- FALSE # pour inclure un filtre sur especes: filtre.sur.especes <- TRUE
if(filtre.sur.especes) {

    ## Modifiez les valeurs pour le filtre sur espèces ici!!! Voir ***
    ## Ces valeurs sont prises en compte seulement lorsque:
    ## filtre.sur.especes = TRUE
    ## Voir aussi la fonction def.filtre.especes() pour définir ces valeurs
    ## dans la console, import.filtre.taxo() pour importer les valeurs à
    ## filtrer directement d'un fichier .csv et la fonction voir.filtre.taxo()
    ## pour voir les valeurs présentement enregistrées

    ## *** Valeurs à modifier, les trois lignes suivantes ***
       taxoF.incl <<- "inclure" # Inclure ou exclure le niveau taxonomique donnée?
       taxoF.utaxo <<- "Groupe" # Niveau taxonomique ciblé
       taxoF.nom <<- c("Crustaces","Mollusques","Echinodermes") # Nom des membres à inclure ou exclure

       }else{
    ## Valeurs automatiques (ne pas changer) qui indique à la fonction
    ## filtreTaxo() de ne pas filtre les espèces
       taxoF.incl <<- "inclure"
       taxoF.utaxo <<- "Groupe"
       taxoF.nom <<- "Tous" }

######################################################
######################################################
####### FIN DE LA PARTIE DU CODE MODIFIABLE ##########
### (en théorie il n'y a rien à modifier plus bas) ###
######################################################
######################################################
charger.codes <- function() {

options('stringsAsFactors'=FALSE) # option générale ôte les colonnes de type "factor" lors de l'import des données

# Crée le dossier Data dans dossier.R s'ils n'existent pas
# Ce dossier contient une copie des données utilisées pour faire
# l'analyse
dir.create("Data",FALSE)
dossier.donnees <<- paste(dossier.R,"/Data/",sep="")

# Importe les fonctions nécéssaires au lancement des analyses
source.with.encoding <<- function (path, encoding, echo = getOption("verbose"), print.eval = echo,
    max.deparse.length = 150, chdir = FALSE)
{
    con = file(path, open = "r", encoding = encoding)
    on.exit(close(con))
    source(con, echo = echo, print.eval = print.eval, max.deparse.length = max.deparse.length,
        chdir = chdir)
}

source.with.encoding("GS_Codes-utils.r",encoding="UTF-8") # lance les codes LIT
source.with.encoding("GS_Selection-Donnees-Par-Projet.r",encoding="UTF-8")
source.with.encoding("GS_MiseAJourDB_TProj.r",encoding="UTF-8") # contient la fonction import.tableaux()
source.with.encoding("GS_ExtractionDonnees_TProj.r",encoding="UTF-8") # contient la fonction prep.analyse()
source.with.encoding("GS_CodesInvertebres_TProj.r",encoding="UTF-8") # contient la fonction Run.INV.biodiv()
source.with.encoding("GS_CodesPoissons_TProj.r",encoding="UTF-8") # lance les codes poissons
source.with.encoding("GS_CodesLIT_TProj.r",encoding="UTF-8") # lance les codes LIT

# si l'objet "data.read" n'existe pas, ou data.read existe
# mais a la valeur "FALSE"
# (ré)extraire et (re)formatter les données
if(!(exists("data.read"))) { import.tableaux()
                         } else { if(!data.read) import.tableaux() }
# (par défaut import.tableaux() lance la fonction prep.analyse()
# ... une fois les tableaux importés)
# ... et on charge les codes graphiques une fois les autres tableaux importés
# pour faciliter la création des légendes de graphs
source.with.encoding("GS_Codes-graphiques.r",encoding="UTF-8") # lance les codes LIT
}

library(Hmisc)
library(tools)
lp <- try(library(reshape2)) # load package reshape
lp2 <- try(library(Rcpp)) # load package Rcpp
lp3 <- try(library(dplyr)) # load package dplyr
if(class(lp)=="try-error") {
      install.packages("reshape2") # installe reshape si requis
library(reshape2)}
if(class(lp2)=="try-error") {
      install.packages("Rcpp") # installe reshape si requis
library(Rcpp)}
if(class(lp3)=="try-error") {
      install.packages("dplyr") # installe reshape si requis
library(dplyr)}
library(grid) # installé automatiquement par dplyr

charger.db <- function() {
aa <- try(load(file="Objets-DB-Codes-R.Rdata"))
if(class(aa)=="try-error") message("\nLe fichier R 'Objects-DB-Codes-R.Rdata' n'existe pas; \nlancer d'abord GS_MotherCode_TProj.r avec refaire.image = TRUE")
}

refaire.image <- TRUE
# (Pas encore implémenté)
# Lorsque refaire.image <- TRUE on sauvegarde un fichier R
# avec tous les tableaux/fonctions formattés une fois le chargement
# terminé. En faisant refaire.image <- FALSE dans une session future,
# on charge directement ce tableau dans la nouvelle session R (et donc
# on peut sauver du temps vu qu'on a pas à refaire tout l'import/
# nettoyage, etc.

if(refaire.image) {
charger.codes()
#save.image(file="Objets-DB-Codes-R.Rdata")
} else {
charger.db()
}
