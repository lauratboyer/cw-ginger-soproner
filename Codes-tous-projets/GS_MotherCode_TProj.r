# Analyses des données KNS (Ginger/Soproner)
# Auteur: Laura Tremblay-Boyer, contact: l.boyer@fisheries.ubc.ca
# Time-stamp: <2015-01-22 15:59:15 Laura>

################################################################
###### Définition des variables principales pour l'analyse #####
## Cette partie du code est à MODIFIER MANUELLEMENT au besoin ##
### *** N'oubliez pas d'enregistrer tous les changements ***####

# Dossier contenant les codes et où les sorties vont être sauvegardées
dossier.R <- getwd() # le 'working directory',
# ... ou sinon mettre le nom du dossier désiré, e.g. C:/Documents/Codes_R
# Dossier de sauvegarde des fichiers .csv des bases de données
dossier.DB <- paste(getwd(), "/DBs tous projets", sep="")
setwd(dossier.R)

############################################
## Variables spatiales et temporelles dans les tableaux
## 'Facteurs_...' à mettre disponible pour les analyses
facteurs.spatio <- c("Geomorpho","N_Impact","Cote","Lieu")
facteurs.tempo <- c("Période.BACI","Saison")

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
### 3. Sorties à produire (invertébrés, poissons, LIT)
### 4. Définition de l'emplacement des fichiers
################################################################

# Utilisateur pour spécifier l'emplacement des dossiers
usernow <- "Laura" # ou "Antoine", ou modifiez directement
                   # l'emplacement des dossiers plus bas (section 4)
if(usernow != "Laura") Sys.setlocale("LC_ALL","French") # console R pour les accents

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


######################################################
######################################################
####### FIN DE LA PARTIE DU CODE MODIFIABLE ##########
### (en théorie il n'y a rien à modifier plus bas) ###
######################################################
######################################################
options('stringsAsFactors'=FALSE) # option générale ôte les colonnes de type "factor" lors de l'import des données

# Crée les dossiers Data/Tableaux/Graphiques dans dossier.R s'ils n'existent pas
dir.create("Data",FALSE)
dossier.donnees <- paste(dossier.R,"/Data/",sep="")

# Importe les fonctions nécéssaires au lancement des analyses
source.with.encoding <- function (path, encoding, echo = getOption("verbose"), print.eval = echo,
    max.deparse.length = 150, chdir = FALSE)
{
    con = file(path, open = "r", encoding = encoding)
    on.exit(close(con))
    source(con, echo = echo, print.eval = print.eval, max.deparse.length = max.deparse.length,
        chdir = chdir)
}

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
# ... une fois les tableaux importés de la DB DropBox)

########################
# Lancement des analyses
# Pour les INV, fonctions Run.INV.biodiv() et Run.INV.densite()
# Pour les poissons, fonction Run.poissons.all()
# Pour le LIT, fonction Run.LIT.all()
if(sorties.INV) { Run.INV.biodiv(); Run.INV.densite() }
if(sorties.poissons) Run.poissons.all()
if(sorties.LIT) Run.LIT.all()



