# Analyses des données KNS (Ginger/Soproner)
# Auteur: Laura Tremblay-Boyer, contact: l.boyer@fisheries.ubc.ca
# Time-stamp: <2013-01-08 16:34:48 Laura>

# Utilisateur pour spécifier l'emplacement des dossiers
usernow <- "Laura"

# Mettre cette option à TRUE si vous voulez produire *toutes* les analyses, ou FALSE si vous
# voulez ajuster certaines fonctions manuellement, ou tourner seulement certains des codes
tourner.tout <- FALSE

# Définition des variables principales pour l'analyse
filtre.annees <- 2006:2012 # indiquer quelles années à inclure dans l'analyse - tableau filtre généré automatiquement
                           # exemples de format: c(2006,2007,2009,2012), 2006:2012, seq(2006,2012,by=2)
filtre.sur.especes <- TRUE # pour inclure un filtre sur especes: filtre.sur.especes <- TRUE
if(filtre.sur.especes) {
  
       taxoF.incl <- "inclure"
       taxoF.utaxo <- "Grp2"
       taxoF.nom <- c("Crustaces","Mollusques","Echinodermes")
       }else{
         taxoF.incl <- "inclure"; taxoF.utaxo <- "Grp2"; taxoF.nom <- "tous" }

sorties.INV <- FALSE # pour produire les sorties des invertébrés: sorties.INV <- TRUE
sorties.LIT <- FALSE # pour produire les sorties pour le LIT: sorties.LIT <- TRUE
sorties.POISSONS <- FALSE # pour produire les sorties des poissons: sorties.POISSONS <- TRUE

# Dossier où les analyses sont faites, et qui contient le code R et les graphiques/tableaux produits
dossier.R <- ifelse(usernow=="Antoine", "C:/dossier.R/","/Users/Laura/Documents/Projects/Soproner_Noumea/Code-R/")
#dossier.R <- "C:/Users/Utilisateur/Documents/Laura/KNS/dossier.R/" # pour ordi de Sara T-B

# Dossier de sauvegarde Dropbox
dossier.DB <- ifelse(usernow=="Antoine","C:/dossier.DB/",
                     "/Users/Laura/Dropbox/KNS_GINGER-SOPRONER/DB_Dernieres_Versions/")
#dossier.DB <- "C:/Users/Utilisateur/Documents/Laura/KNS/dossier.DB/"

options('stringsAsFactors'=FALSE)
######################################################
######################################################
setwd(dossier.R) # défini le working directory de R

# Crée les dossiers Data/Tableaux/Graphiques dans dossier.R s'ils n'existent pas
ndossiers <- c("Data","Tableaux","Graphiques")
dossier.in <- file.info(ndossiers)$isdir
if(sum(is.na(dossier.in))>0) sapply(ndossiers[is.na(dossier.in)], function(x) dir.create(x))
fig.dir <- paste(dossier.R,"/Graphiques/",sep="")
tabl.dir <- paste(dossier.R,"/Tableaux/",sep="")
dossier.donnees <- paste(dossier.R,"/Data/",sep="")

# Importe les fonctions nécéssaires ? l'analyse:
source("GS_MiseAJourDB.r") # contient la fonction import.tableaux()
source("GS_ExtractionDonnees.r") # contient la fonction prep.analyse()
source("GS_CodesInvertebres_Launch.r") # contient la fonction Run.INV.biodiv()
source("GS_CodesInvertebres_Densite.r")
source("GS_CodesInvertebres_Biodiv.r")

# si l'objet "data.read" n'existe pas, ou data.read existe, mais a la valeur FALSE,
# (ré)extraire et (re)formatter les données
if(!(exists("data.read"))) import.tableaux()
if(exists("data.read")) if(!data.read) import.tableaux()
# (par défaut import.tableaux() lance la fonction prep.analyse()
# ... une fois les tableaux importés de la DB DropBox)

# Définir le filtre sur especes, seulement si l'objet "filtre.sur.especes" n'est pas "FALSE"
#!!!! lancer manuellement filtre.especes("oui") ou filtre.especes() !!!!!!!
#if(!filtre.sur.especes) { filtre.especes() }else{ filtre.especes("oui") }
#filtre.especes() # temporaire en attendant de réparer filtre.especes()

# Lancer les analyses
if(sorties.INV) { Run.INV.biodiv(); Run.INV.densite() }



