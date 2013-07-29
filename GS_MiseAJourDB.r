# Analyses des donnÃ©es KNS (Ginger/Soproner)
# Auteur: Laura Tremblay-Boyer, contact: l.boyer@fisheries.ubc.ca
# Time-stamp: <2013-07-29 16:55:58 Laura>

# Sujet: Ce code vérifie que les tableaux utilisÃ©s pour les analyses sont Ã  jour
# ... et dans le cas Ã©chÃ©ant modifie la base de donnÃ©es en consÃ©quence

###############################################################
###############################################################
# Noms des fichiers
fc.inv <- "Invertebres.csv" # donnÃ©es brutes invertÃ©brÃ©s
fc.bioeco <- "Bioeco.csv" # info sur les espÃ¨ces de poissons
fc.poissons <- "Poissons.csv" # donnÃ©es brutes poissons
fc.LIT <- "LIT.csv" # donnÃ©es brutes LIT
fc.typoLIT <- "Typo_LIT.csv"
fc.transect <- "Facteurs.csv" # info sur les transects
fc.Bacip <- "Periode_BACIP_Campagnes.csv" # pÃ©riodes Bacip
#fc.filtre <- "Tables des campagnes a utiliser pour analyses temporelles.csv" # tableau filtre

# types de tableaux
type.tbl <- c("inv","bioeco","poissons","data.LIT","typo.LIT","transect","Bacip")#,"filtre")

###############################################################
###############################################################
# DÃ©finition de la fonction qui importe la derniÃ¨re version des tableaux:

"import.tableaux" <- function(prep.tbl=TRUE) {

  # noms des objets correspondants sous R
  # attention ne pas changer ces noms car ils sont utilises dans plusieurs fonctions
  obj.names <- c("data.inv","data.bioeco","data.poissons","data.LIT",
                 "index.LIT","data.info.transect","pr.Bacip")

  # objet synthÃ¨se pour les noms de fichiers
  obj.files <- c(fc.inv,fc.bioeco,fc.poissons,fc.LIT,fc.typoLIT,fc.transect,fc.Bacip)#,fc.filtre)
  names(obj.files) <- type.tbl

  ############################################################
  ############################################################
  # Comparer les dates de modifications des fichiers du dossier Dropbox
  # vs le dossier pour les analyses
  DBinf.now <- lapply(type.tbl,function(x)
                    file.info(paste(dossier.DB,obj.files[x],sep=""))$mtime)
  names(DBinf.now) <- type.tbl

  # Chargement des dates de modifications prÃ©cÃ©dentes (objet DBinf)
  if(!(file.exists(paste(dossier.donnees,"FichiersDropBox_DatesModif.Rdata",sep="")))) {
    memeVersion <- rep(TRUE,length(type.tbl)) # importer tous les fichiers de DB
  } else {
    load(paste(dossier.donnees,"FichiersDropBox_DatesModif.Rdata",sep=""))
    # Comparaison aux dates prÃ©sentes
    memeVersion <- unlist(sapply(type.tbl,function(x) DBinf[[x]] != DBinf.now[[x]]))
    # importer nouvelle version du tableau si manquant dans dossier donnees
    memeVersion[is.na(memeVersion)] <- TRUE
}
  # Importer du dossier Dropbox les fichiers modifiÃ©s (au besoin)
  if(sum(memeVersion) > 0) {
    # Suppression des vieux fichiers:
    dmm <- sapply(type.tbl[memeVersion],function(x)
         suppressWarnings(file.remove(paste(dossier.donnees,obj.files[x],sep=""))))

    # Import des nouveaux fichiers
    dmm <- sapply(type.tbl[memeVersion],function(x)
                  file.copy(paste(dossier.DB,obj.files[x],sep=""),dossier.donnees))
    if(sum(dmm)!=sum(memeVersion)) print("Attention nouveaux tableaux non importÃ©s")

    # Enregistrement des nouvelles dates de modification
    DBinf <- DBinf.now
    save(DBinf,file=paste(dossier.donnees,"FichiersDropBox_DatesModif.Rdata",sep=""))
    print(paste("Nouvelle version pour les fichers", type.tbl[memeVersion])) }

  # Imprimer avertissement au sujet des guillements: s'il y a des guillements non-fermés dans les
  # fichiers excels l'import du tableau ne sera pas correct -- e.g. "pomme" est ok, mais pas pomme"
  print("Attention: s'assurer que tous les guillements sont fermés dans les fichiers Excel")

  ############################################################
  ############################################################
  # Tous les tableaux sont maintenant à jour -> chargement sous R pour débuter les analyses
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
      ecodes <- "<e9>|<e8>|(\303\251)|(\303\250)|\216|\217|\351|\350"
      x <- gsub(ecodes,"e",x) # remplace e accent aigu/grave
      x <- gsub("<ef>|(\303\257)|\357|\225","i",x) # remplace i accent trema
      x <- gsub("<a0>|\312","",x) # mystery character removed
      # conversion à iso-8859-5 pour ôter les accents (sur PC)
      x <- iconv(x, "latin1","iso-8859-5")
  } else {x}}

  # extraire des fichiers csv et assigner aux noms d'objets (voir obj.names)
  dmm <- sapply(1:length(obj.names),function(x) {
      # vérifier séparateur
      cs <- scan(paste(dossier.donnees, obj.files[x], sep=""),"character")
      seprt <- ifelse(sum(grepl(",",cs)) > sum(grepl(";",cs)), ",",";")

      objnow <- read.csv(paste(dossier.donnees, obj.files[x], sep=""),
                         sep=seprt, dec=".") # point decimal "." (voir aussi tradfunk())
      objnow <- objnow[!apply(objnow,1,is.empty),]; # oter les rangees vides

      objnow <- data.frame(lapply(objnow, tradfunk))

      assign(obj.names[x], objnow, .GlobalEnv)}) # creer l'objet dans l'envir global

  print("Tableaux importés:")
  print(DBinf.now[obj.names %in% ls(.GlobalEnv)])

  message("\n\nBases de données importées avec succès. Lancement du formattage des tableaux...\n\n")
  ############################################################
  ############################################################
  if(prep.tbl) prep.analyse() # lance préparation données + création tableaux annexes
}
