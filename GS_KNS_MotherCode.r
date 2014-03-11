# Analyses des données KNS (Ginger/Soproner)
# Auteur: Laura Tremblay-Boyer, contact: l.boyer@fisheries.ubc.ca
# Time-stamp: <2014-03-12 09:11:37 Laura>

################################################################
###### Définition des variables principales pour l'analyse #####
## Cette partie du code est à MODIFIER MANUELLEMENT au besoin ##
### *** N'oubliez pas d'enregistrer tous les changements ***####

################################################################
### Contenu du code:
### 1. Filtre années
### 2. Filtre espèces
### 3. Sorties à produire (invertébrés, poissons, LIT)
### 4. Définition de l'emplacement des fichiers
################################################################

# Utilisateur pour spécifier l'emplacement des dossiers
usernow <- "Laura" # ou "Antoine", ou modifiez directement
                   # l'emplacement des dossiers plus bas (section 4)

############################
### 1. Filtre sur ANNEES ###
############################
filtre.annees <- 2006:2013 # indiquer quelles années à inclure dans l'analyse - tableau filtre généré automatiquement
                           # exemples de format:
                           # c(2006,2007,2009,2012), 2006:2012, seq(2006,2012,by=2)

#############################
### 2. Filtre sur ESPECES ###
#############################

filtre.sur.especes <- FALSE # pour inclure un filtre sur especes: filtre.sur.especes <- TRUE
if(filtre.sur.especes) {

    ## Modifiez les valeurs pour le filtre sur espèces ici!!!
    ## Ces valeurs sont prises en compte seulement lorsque:
    ## filtre.sur.especes = TRUE
    ## Voir aussi la fonction def.filtre.especes() pour définir ces valeurs
    ## dans la console, import.filtre.taxo() pour importer les valeurs à
    ## filtrer directement d'un fichier .csv et la fonction voir.filtre.taxo()
    ## pour voir les valeurs présentement enregistrées

       taxoF.incl <- "inclure" # Inclure ou exclure le niveau taxonomique donnée?
       taxoF.utaxo <- "Groupe" # Niveau taxonomique ciblé
       taxoF.nom <- c("Crustaces","Mollusques","Echinodermes") # Nom des membres à inclure ou exclure

       }else{
    ## Valeurs automatiques (ne pas changer) qui indique à la fonction
    ## filtreTaxo() de ne pas filtre les espèces
       taxoF.incl <- "inclure"
       taxoF.utaxo <- "Groupe"
       taxoF.nom <- "Tous" }

#############################
### 3. Sorties à produire ###
#############################

sorties.INV <- FALSE # pour produire les sorties des invertébrés: sorties.INV <- TRUE
sorties.LIT <- FALSE # pour produire les sorties pour le LIT: sorties.LIT <- TRUE
sorties.poissons <- FALSE # pour produire les sorties des poissons: sorties.POISSONS <- TRUE

#############################
## 4. Définition DOSSIERS ###
#############################

# Dossier où les analyses sont faites, et qui contient le code R et les graphiques/tableaux produits
# Attention de rajouter un "/" à la fin, selon le format "C:/.../dossier.R/"
dossier.R <- ifelse(usernow=="Antoine", "C:/dossier.R/","/Users/Laura/Projects/cw-ginger-soproner/")

# Dossier de sauvegarde Dropbox
dossier.DB <- ifelse(usernow=="Antoine","C:/dossier.DB/",
                     "/Users/Laura/Dropbox/KNS_GINGER-SOPRONER/DB_Dernieres_Versions/")


######################################################
######################################################
####### FIN DE LA PARTIE DU CODE MODIFIABLE ##########
### (en théorie il n'y a rien à modifier plus bas) ###
######################################################
######################################################

setwd(dossier.R) # défini le working directory de R
options('stringsAsFactors'=FALSE) # option générale ôte les colonnes de type "factor" lors de l'import des données

# Crée les dossiers Data/Tableaux/Graphiques dans dossier.R s'ils n'existent pas
ndossiers <- c("Data","Tableaux","Graphiques")
dossier.in <- file.info(ndossiers)$isdir
if(sum(is.na(dossier.in))>0) sapply(ndossiers[is.na(dossier.in)], function(x) dir.create(x))
fig.dir <- paste(dossier.R,"/Graphiques/",sep="")
tabl.dir <- paste(dossier.R,"/Tableaux/",sep="")
dossier.donnees <- paste(dossier.R,"/Data/",sep="")

# Importe les fonctions nécéssaires au lancement des analyses
source("GS_MiseAJourDB.r") # contient la fonction import.tableaux()
source("GS_ExtractionDonnees.r") # contient la fonction prep.analyse()
source("GS_CodesInvertebres_Launch.r") # contient la fonction Run.INV.biodiv()
source("GS_CodesInvertebres_Densite.r")
source("GS_CodesInvertebres_Biodiv.r")
source("GS_CodesPoissons_Launch.r") # lance les codes poissons
source("GS_CodesPoissons_BioDens.r") # contient les codes densités
                                     # et biodiversités pour les poissons
source("GS_CodesLIT_Launch.r") # lance les codes LIT
source("GS_CodesLIT_Couvrt.r") # contient les codes pour statistiques LIT


# si l'objet "data.read" n'existe pas, ou data.read existe
# mais a la valeur "FALSE"
# (ré)extraire et (re)formatter les données
if(!(exists("data.read"))) { import.tableaux()
                         } else { if(!data.read) import.tableaux() }
# (par défaut import.tableaux() lance la fonction prep.analyse()
# ... une fois les tableaux importés de la DB DropBox)

########################
# Lancement des analyses
# Pour les INV, fonctions Run.INV.biodiv() et Run.INV.densite()
# Pour les poissons, fonction Run.poissons.all()
# Pour le LIT, fonction Run.LIT.all()
if(sorties.INV) { Run.INV.biodiv(); Run.INV.densite() }
if(sorties.poissons) Run.poissons.all()
if(sorties.LIT) Run.LIT.all()



