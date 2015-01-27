# Analyses des données KNS (Ginger/Soproner)
# Auteur: Laura Tremblay-Boyer, contact: l.boyer@fisheries.ubc.ca
# Time-stamp: <2015-01-26 08:36:57 Laura>

# Sujet: Ce code importe les bases de données multi-projets
# et nettoie les champs au besoin (notamment en ôtant les accents)

# Note pour Laura: sauvegarder les csv en mode Windows sous Excel/Mac
# Locale R = FR-UTF-8; encodage pour l'import = iso-8859-1
#Sys.setlocale("LC_ALL", "fr_FR.UTF-8") # changement de locale

###############################################################
###############################################################
# Noms des fichiers
fc.inv <- "Invertebres.csv" # données brutes invertébrés
fc.LIT <- "LIT.csv" # données brutes LIT
fc.typoLIT <- "Typo_LIT.csv"  # infos typologie LIT
fc.quad <- "Quadrats.csv" # données brutes quadrats
fc.poissons <- "Poissons.csv" # données brutes poissons
fc.bioeco <- "Bioeco.csv" # infos sur les espèces de poissons
fc.transect <- "Facteurs_spatiaux.csv" # facteurs explicatifs spatiaux
fc.fct.temprl <- "Facteurs_temporels.csv" # facteurs explicatifs temporels


# types de tableaux
type.tbl <- c("inv","data.LIT","typo.LIT","data.quad","poissons","bioeco","transect",
              "fact.temprl")

###############################################################
###############################################################
# Définition de la fonction qui importe la dernière version des tableaux:

"import.tableaux" <- function(prep.tbl=TRUE) {

  # noms des objets correspondants sous R
  # attention ne pas changer ces noms car ils sont utilisés dans plusieurs fonctions
  obj.names <- c("data.inv","data.LIT",
                 "index.LIT","data.quad","data.poissons","data.bioeco",
                 "data.info.transect","data.info.temprl")

  # objet synthèse pour les noms de fichiers
  obj.files <- c(fc.inv, fc.LIT, fc.typoLIT, fc.quad, fc.poissons, fc.bioeco,
                 fc.transect, fc.fct.temprl)
  names(obj.files) <- type.tbl

  dmm <- sapply(obj.files, function(on)
                file.copy(paste(dossier.DB,on,sep="//"),
                          paste(dossier.donnees,on,sep="//")))

  DBinf.now <- lapply(type.tbl,function(x)
                    file.info(paste(dossier.DB,obj.files[x],sep="/"))$mtime)
  names(DBinf.now) <- obj.names

  # Imprimer avertissement au sujet des guillements:
  # s'il y a des guillements non-fermés dans les
  # fichiers excels l'import du tableau ne sera pas correct --
  # e.g. "pomme" est ok, mais pas pomme"
  message("Attention: s'assurer que tous les guillemets sont fermés dans les fichiers Excel")

  ############################################################
  ############################################################
  # Chargement des tableaux sous R pour débuter les analyses
  # Importer feuilles Excel et nommer les objets comme défini dans "obj.names"
  # (objets importés dans l'environnment global .GlobalEnv pour qu'ils soient accessibles partout
  check.o <- obj.names %in% ls() # verifier s'il y a des tableaux deja chargés
  if(sum(check.o)>0) rm(list=obj.names[check.o]) # oter les objets au besoin

  # fonction pour tester si les rangees sont vides:
  is.empty <- function(x) all(c(is.na(x) | is.null(x) | x == ""))

  # traduire les codes d'accent selon le type d'encodage
  # assume encodage est "latin1", spécifié en conséquence durant l'import
  # avec read.csv()
  tradfunk <<- function(x) {

      # conversion au cas où mal fait sous .csv (dû au dec = ",")
      if(class(x) == "character") x <- type.convert(x, as.is=TRUE, dec=",")

      if(class(x) == "character") { # ... et on continu si x reste un charactère
      enc <- unique(Encoding(x))
      if(length(enc) > 1 & (!("latin1" %in% enc))) {
          warning("Attention encodage non-déclaré latin1") }

      # remplacement des charactères (sur Mac)
      x <- gsub("é|è", "e", x)
      x <- gsub("ï", "i", x)

      # et pourquoi pas on nettoie les espaces vides avant et après
      x <- gsub("^\\s+","",gsub("\\s+$","",x))

      # + on ajuste les champs minuscules + majuscule au début
      capitalize <- function(x) gsub('(\\w)([\\w|\\s]*)','\\U\\1\\L\\2',x,perl=TRUE)
      x <- capitalize(tolower(x))
      if(all(grepl(".*_.*_.*_.*",x))) x <- toupper(x) # ... à part pour colonne Id
      return(x)
      } else {return(x)}}


  # extraire des fichiers csv et assigner aux noms d'objets (voir obj.names)
  dmm <- sapply(1:length(obj.names),function(x) {
      # vérifier séparateur de colonnes
      cs <- scan(paste(dossier.donnees, obj.files[x], sep=""),
                 "character",quiet=TRUE)
      seprt <- suppressWarnings(ifelse(sum(grepl(",",cs)) > sum(grepl(";",cs)), ",",";"))
      # importer le tableau csv avec le séparateur identifié ci-haut
      objnow <- read.csv(paste(dossier.donnees, obj.files[x], sep=""),
                         sep=seprt, dec=".", fileEncoding="iso-8859-1") # point décimal "." (voir aussi tradfunk())
      objnow <- objnow[!apply(objnow,1,is.empty),]; # ôter les rangées vides
      objnow <- data.frame(lapply(objnow, tradfunk))
      assign(obj.names[x], objnow, .GlobalEnv)}) # créer l'objet dans l'envir global

  print("Tableaux importés + dernière date de modification des fichiers .csv")

  print(DBinf.now[obj.names %in% ls(.GlobalEnv)])

  message("\n\nBases de données importées avec succès. Lancement du formattage des tableaux...\n\n")
  ############################################################
  ############################################################
  if(prep.tbl) prep.analyse() # lance préparation données + création tableaux annexes
}
