  ## Analyses des données KNS (Ginger/Soproner)
# Auteur: Laura Tremblay-Boyer, contact: l.boyer@fisheries.ubc.ca
# Time-stamp: <2013-01-21 15:04:46 Laura>

# Sujet: Formattage des tableaux de données brutes pré-analyse,
# création de tableaux annexes + fonctions de base pour l'analyse

# Put this elsewhere:
# Définition lien fichiers pour sauvegarder graphiques/tableaux
fig.dir <- paste(dossier.R,"//Graphiques//",sep='')
tabl.dir <- paste(dossier.R,"//Tableaux//",sep='')

prep.analyse <- function() {

  lp <- require(reshape) # load package reshape
  if(!lp) install.packages("reshape") # installe reshape si requis
  #######################################################################
  ###################### Formattage des tableaux ########################
  #######################################################################

  ### Fonctions utiles pour formattage
  if(!exists("capitalize")) {
  capitalize <<- function(x) {
    gsub('(\\w)(\\w*)','\\U\\1\\L\\2',tolower(x),perl=TRUE)
  }
}

  if(!exists("last")) {
  last <<- function(x) x[length(x)]
}
  
  ####################################
  ###### Information transects #######
  ####################################
  # info.transect: informations sur les transects
  info.transect$cpus <- info.transect$Prod_ha # nouvelle colonne pour le cpus
  
  # Nettoyer accents noms géométrie
  geom.key <- data.frame("Geom"=unique(info.transect$Geom),
                         "Geomorpho"=c("Recif barriere externe","Recif barriere interne",
                           "Recif reticule","Recif frangeant","Passe","Herbiers"))
  print("make geom.key robust to unique order of info.transect$Geom")
  info.transect <- merge(info.transect, geom.key, by="Geom")
  
  ########### Données LIT ############  
  ####################################
  # data.LIT: données brutes LIT
  # index.LIT: info complementaires LIT

  # ote rangées sans stations
  data.LIT <- data.LIT[data.LIT$St != "",]
  # défini noms de colonnes pour index.LIT:
  names(index.LIT) <- c("Code_LIT","CODE_DET","S_Corail_Acro","S_Corail_Forme",
                      "S_Corail_All","S_Corail_Sensi","S_Abio_Corail_All")
  # réparer accents
  index.LIT$S_Corail_Forme <- gsub("\216","e",index.LIT$S_Corail_Forme)

  ###### Information biologie/écologie poissons #######
  #####################################################
  # bioeco.all: info complémentaires sur les poissons
  
  bioeco.all <- bioeco[,c("Code_Sp","famille","genre","espece","taille_moy","commercial_2",
                        "cible_VKP","a","b","groupe_troph1","mobilite")]
  names(bioeco.all) <- c("Code_SP","Famille","Genre","Espece","Taille","Peche","Cible",
                       "Coeff_a","Coeff_b","Groupe_Trophique","Mobilite")

  # Clarifier étiquettes
  # Peche commerciale:
  bioeco.all[bioeco.all$Peche=="c","Peche"] <- "Comm"
  bioeco.all[bioeco.all$Peche=="nc","Peche"] <- "Non.Comm"
  bioeco.all[bioeco.all$Cible == 0,"Cible"] <- "Cible.Non"
  bioeco.all[bioeco.all$Cible == 1,"Cible"] <- "Cible.Oui"

  # Rajouter id pour Pomacentridae et Chaetodontidae
  bioeco.all$fmlabel <- NA
  bioeco.all$fmlabel[bioeco.all$Famille=="POMACENTRIDAE"] <- "Pom"
  bioeco.all$fmlabel[bioeco.all$Famille=="CHAETODONTIDAE"] <- "Chaet"
  
  # Rajouter info groupes trophiques:
  GT.label <- data.frame("Groupe_Trophique"=c("P","C","Z","H"),
                       "GTlabel"=c("Pisc","Carn","Planct","Herb"))
  bioeco.all <- merge(bioeco.all,GT.label,by="Groupe_Trophique")

  # Rajouter info mobilité:
  mob.label <- data.frame("Mobilite"=0:4,"moblabel"=c("Mob.0","Terr","Sed","Mob","TrMob"))
  bioeco.all <- merge(bioeco.all,mob.label,by="Mobilite",all.x=TRUE)

  # Tableau bioeco, info a + b seulement
  bioeco <- bioeco[,c("Code_Sp","a","b")]
  names(bioeco) <- c("Code_SP","a","b") # compatibilité avec data.poissons

  # Vérifier s'il y a des espèces de poissons sans valeurs a ou b:
  print(paste(length(unique(bioeco[bioeco$a=="inconnu" | bioeco$b=="inconnnu","Code_SP"])),"especes sans valeur a ou b otees de l'analyse"))
  bioeco <- bioeco[bioeco$a!="inconnu" & bioeco$b!="inconnu",]
  bioeco$a <- as.numeric(bioeco$a); bioeco$b <- as.numeric(bioeco$b)

  ###### Données de comptage Poissons ######
  ##########################################
  # data.poissons: données de comptage brutes sur les poissons
  names(data.poissons) <- c("Campagne","Date","St","Obs","Vis","Courant","Code_SP",
                          "Famille","Genre","Espece","G_Sp","N","L","D1","D2",
                          "Secteur","Categorie","Note")

  ########################
  # Nettoyage du tableau:

  # Garder seulement les valeurs de N & L déclarées numériques
  data.poissons$N <- as.numeric(data.poissons$N)
  data.poissons$L <- as.numeric(data.poissons$L)

  # Ajuster si besoin le format pour les noms d'espèces:
  data.poissons$G_Sp <- gsub('(\\w)(\\w*)','\\U\\1\\L\\2',tolower(data.poissons$G_Sp),perl=TRUE) # Lettre majuscule en premier seulement
  data.poissons$St <- toupper(data.poissons$St) # Toutes les stations en majuscules

  # Nettoyage du tableau: colonne D1/D2/N/L
  # Espèces sans Code_SP défini ôtées
  # Imprimer avertissements
  print(paste(nrow(data.poissons[data.poissons$Code_SP  %in% c("","#N/A"),]),"rangees sans Code_SP"))
  print(paste(nrow(data.poissons[is.na(data.poissons$D1),]),"valeurs D1 -> 0"))
  print(paste(nrow(data.poissons[is.na(data.poissons$D2),]),"valeurs D2 -> 0"))
  print(paste(nrow(data.poissons[is.na(data.poissons$N),]),"valeurs N -> 0"))
  print(paste(sum(is.na(data.poissons$L)),"rangees sans valeur L"))

  # Ôter rangées:
  data.poissons <- data.poissons[!(data.poissons$Code_SP %in% c("","#N/A")),]
  data.poissons$D1[is.na(data.poissons$D1)] <- 0
  data.poissons$D2[is.na(data.poissons$D2)] <- 0
  data.poissons$N[is.na(data.poissons$N)] <- 0
  data.poissons <- data.poissons[!(is.na(data.poissons$L)),] # Taille

  # Expansion du data frame pour inclure toutes les combinaisons Campagne/St/Code_SP:
  # ... donc pour les poissons les densités sont nulles sur toutes Campagnes/St où
  # ... l'espèce est non-observée
  all.comb.poissons <- expand.grid("Campagne"=unique(data.poissons$Campagne),
                                   "St"=unique(data.poissons$St),
                                   "Code_SP"=unique(data.poissons$Code_SP))
  St.by.year <- unique(data.poissons[,c("Campagne","St")]) # year/station samples
  all.comb.poissons <- merge(all.comb.poissons, St.by.year, by=c("Campagne","St"))

  # ... reunion avec données poissons pour rajouter les densité nulles
  data.poissons2 <- merge(all.comb.poissons,data.poissons,by=c("Campagne","St","Code_SP"),all.x=TRUE)
  data.poissons2$N[is.na(data.poissons2$N)] <- 0 # abondance à zero si non-observée
  data.poissons <- data.poissons2 

  
  ##########################################
  #### Données de comptage Invertébrés #####
  ##########################################
  # data.inv: données de comptage brutes pour invertébrés
  
  dbio <- data.inv[,c("Campagne","St","T","Grp2","S_Grp2","F2","G2","G_Sp","N","D","Ltrans")]
  dbio$St <- toupper(data.inv$St)

  ############################
  # nettoyage:
  # garde seulement les rangées avec transects A, B, C
  nTnd <- nrow(dbio[!(dbio$T %in% c("A","B","C")),])
  dbio <- dbio[dbio$T %in% c("A","B","C"),]; print(paste("Approx",nTnd,"rangees otees vu transect non-def"))               
  
  ############################
  # abondances nulles:
  # commencer par oter les rangées N=0 pour eliminer les stations où l'espèce n'est observée sur aucun T
  dbio <- dbio[dbio$N > 0,]
  
  # rajouter les transects N=0 *seulement* lorsque l'espèce est observée sur la ST mais pas tous les transects
  # créer identifiant unique pour chaque combinaison campagne/st/espèce observée
  dbio$uID <- paste(dbio$Campagne,dbio$St,dbio$G_Sp,sep="_")
  index.dbio <- unique(dbio[,c("uID","Campagne","St","Grp2","S_Grp2",
                               "F2","G2","G_Sp","D","Ltrans")])
  dbio.allT <- expand.grid("T"=c("A","B","C"),"uID"=unique(dbio$uID))
  dbio.tmp <- merge(dbio[,c("uID","T","N")],dbio.allT,by=c("uID","T"),all.y=TRUE)
  dbio.tmp$N[is.na(dbio.tmp$N)] <- 0 # valeurs NA remplacées par 0 pour remplir les transects manquants
  dbio <- merge(dbio.tmp, index.dbio)

  ################
  # Oter les accents des noms des invertébrés pour éviter les erreurs dûes à l'encodage
  # Si l'encodage UTF-8 est bien lu par l'ordi (les accents apparaissent correctement)    
  dbio$Grp2 <- gsub("é","e",capitalize(tolower(trim(dbio$Grp2))))
  dbio$S_Grp2 <- gsub("é","e",dbio$S_Grp2)
  dbio$S_Grp2 <- gsub("è","e",dbio$S_Grp2)
  dbio$S_Grp2 <- gsub("ï","i",dbio$S_Grp2)
  dbio$S_Grp2 <- capitalize(tolower(trim(dbio$S_Grp2)))
  #dbio$F2 <- 
  print("trouver signe pour Pterasteridae")
        
  # Si l'encodage UTF-8 n'est pas lu par l'ordi (les accents suivent un format similaire à <U+00E9>)
  # Conversion à encodage latin1 qui permet de substituer les charactères
  if("UTF-8" %in% unique(Encoding(dbio$Grp2))) {
    dbio$Grp2 <- iconv(dbio$Grp2,"UTF-8","latin1")
    dbio$S_Grp2 <- iconv(dbio$S_Grp2,"UTF-8","latin1")
    dbio$F2 <- iconv(dbio$F2,"UTF-8","latin1")
    dbio$Grp2 <- gsub("<e9>","e",dbio$Grp2) # remplace e accent aigu
    dbio$S_Grp2 <- gsub("<e9>","e",dbio$S_Grp2) # remplace e accent aigu
    dbio$S_Grp2 <- gsub("<e8>","e",dbio$S_Grp2) # remplace e accent grave
    dbio$S_Grp2 <- gsub("<ef>","i",dbio$S_Grp2) # remplace i accent trema
    dbio$F2 <- gsub("<a0>","",dbio$F2) # remplace ae
  }
  # Première lettre en majuscule, le reste en minuscules
  dbio$G_Sp <- capitalize(tolower(dbio$G_Sp))
  dbio$F2 <- capitalize(tolower(dbio$F2))
  dbio$G2 <- capitalize(tolower(dbio$G2))

  # Renommer colonnes F2 et G2      
  names(dbio)[grep("F2",names(dbio))] <- "Famille"
  names(dbio)[grep("G2",names(dbio))] <- "Genre"

  # Correction d'erreur dans la base de données (temporaire)
  # tableau index espèce/sous-groupe/groupe
  dbio[dbio$G_Sp %in% c("Paguritta sp.", "Ciliopagurus strigatus"),"S_Grp2"] <- "Decapodes"
  dbio[dbio$G_Sp == "Drupa sp.","Famille"] <- "Muricidae"
  dbio[dbio$G_Sp == "Aplysia sp.","Famille"] <- "Aplysiidae"
  print("sous groupe de ciliopagurus strigatus devrait etre decapode et non crustace?, aussi fixed pour Drupa sp. et typo dans Aplysia sp.")

  # Fix Ltrans that are NA
  ltrans.ind <- na.omit(unique(dbio[,c("Campagne","St","T","Ltrans")]))
  names(ltrans.ind)[4] <- "Ltrans.new"
  dbio <- merge(dbio, ltrans.ind, by=c("Campagne","St","T"))
  dbio <- dbio[,!(names(dbio)=="Ltrans")]
  names(dbio)[names(dbio)=="Ltrans.new"] <- "Ltrans"
   
  # rassembler en 1 rangée les observations de la même espèce sur le transect
  dbio <- aggregate(list("N"=dbio$N), as.list(dbio[,c("Campagne","St","T","Grp2","S_Grp2",
                           "Famille","Genre","G_Sp","Ltrans")]), sum)
  
  # Calculer la densité en hectares
  # longueur du transect est de 50 mètres, largeur est définie sous Ltrans
  dbio$D <- dbio$N/(50*dbio$Ltrans) * 10000
 
  ########################
  ### Tableaux annexes ###

  # clés noms des invertébrés
  InvGr.Key <- c("Algues","Ascidies","Cnidaires","Crustac?s","Echinodermes",
                 "Eponges","Mollusques","Phan?rogame","Vers")
  
  # tableau index pour la taxonomie de toutes les espèces observée:
  index.invSp <- unique(dbio[,c("G_Sp","Genre","Famille","S_Grp2","Grp2")])

  ################
  # infos sur station/mission/année
  info.transect.INV <- unique(data.inv[,c("Annee","Mois","Mission","Campagne","St")])
  info.transect.INV <- merge(info.transect.INV, info.transect,
                             by="St")[,c("Annee","Mois","Mission","Campagne","St",
                               "Geomorpho","N_Impact")]
  info.transect.INV.geo <- aggregate(info.transect.INV[,c("Mois", "Annee")],as.list(info.transect.INV[,c("Mission","Campagne","Geomorpho")]), last)
  print("si un transect/Campagne échantilloné dans 2 mois différents, garde 1 mois/année seulement dans le tableau")


  ####################################
  ###### Périodes BACIP Campagnes ####
  ####################################
  # pr.Bacip: catégorie avant/pendant/aprés pour les années
  names(pr.Bacip) <- c("Campagne","Annee","Periode_3","Periode_2")

  #######################################################
  ###### Filtres campagnes annuelles/semestrielles ######
  #######################################################
  # filtre.Camp: tableau des campagnes à utiliser pour analyses temporelles
  # modif après discussion avec Antoine 12 Octobre 2012:
  # re-créer tableau filtre à partir des années désirées pour l'analyse
  
  creerFiltre <- function(qAnnees) {

    if(length(qAnnees)==1) stop("Attention: spécifier 2 années ou plus")
    tb <- unique(dbio[,c("St","Campagne")])
    tb$echtl <- 1
    tb2 <- cast(tb, St ~ Campagne, value="echtl")

    # sélectionner les colonnes avec les années désirées pour le filtre
    # campagnes semestrielles
    keySmstrl <- c(paste("A_",qAnnees,sep=""),paste("S_",qAnnees,sep=""))
    tb3 <- tb2[,unlist(sapply(keySmstrl,function(x) grep(x,names(tb2))))]
    wStat <- tb2$St[rowSums(tb3,na.rm=TRUE)==ncol(tb3)]
    T_S_inv <- apply(expand.grid(wStat,names(tb3)),1,paste,collapse="_")

    # campagnes annuelles
    keyAnnuel <- paste("A_",qAnnees,sep="")
    tb3 <- tb2[,unlist(sapply(keyAnnuel,function(x) grep(x,names(tb2))))]
    if("data.frame" %in% class(tb3)) { wStat <- tb2$St[rowSums(tb3,na.rm=TRUE)==ncol(tb3)]
                                   }else{
                                     wStat <- tb2$St[which(na.omit(tb3==1))]}
    T_A_inv <- apply(expand.grid(wStat,names(tb3)),1,paste,collapse="_")
    
    return(list("T_S_inv"=T_S_inv, "T_A_inv"=T_A_inv))
  }

  filtre.Camp <<- creerFiltre(filtre.annees)
  
  ###################################################
  ######## Fonctions génériques #####################
  # Définition de fonctions qui seront utilisées couramment dans le code

  # 1. Applique le filtre spécifié au tableau donné en argument
  filtreTable <<- function(wtable, wfiltre) {
    if(wfiltre %in% c("T_A_inv","T_S_inv")) { # appliquer le filtre si spécifié
    print(paste("Stations filtrées par",wfiltre))
    wtable$key <- paste(wtable$St, wtable$Campagne, sep="_")
    dd.filt <- merge(data.frame("key"=filtre.Camp[[wfiltre]]),
                     wtable,by="key", drop.x="key")
  } else {
    CmpTag <- paste(filtre.annees,collapse="|")
    wCampKeep <- grep(CmpTag, unique(wtable$Campagne), value=TRUE)
    print(wCampKeep)
    wtable <- wtable[wtable$Campagne %in% wCampKeep,]
    }
  }

  # 2. Converti les noms de campagne en année (charactère -> numérique)
  as.year <<- function(x) {
    s1 <- sub("A_","",x)
    s2 <- sub("S_","",s1)
    return(as.numeric(s2)) }

  # 3. Dessine les barres d'erreurs
  draw.SE <<- function(mat, couleur=1, typel=1) {
  # Inclure matrice avec x en colonne 1 et valeur min/max en colonne 2 et 3
  mat$diff <- mat[,3]-mat[,2]
  mat <- mat[mat$diff != 0,] # garder seulement les valeurs de SE existantes
  mat[,2][mat[,2]<0] <- 0
  dmm <- sapply(1:nrow(mat),function(i) arrows(mat[i,1],mat[i,3],
                                                  mat[i,1],mat[i,2],code=3,lty=typel,
                                                  col=couleur,angle=90,length=0.1)) }

  # 4. Calcul du nombre de stations échantillonées par groupement spatial, selon le filtre
  # Groupement spatial: géomorphologie, ou géomorphologie/impact
  # défini pour les invertébrés mais pourrait être appliqué aux poisssons
  # (Retourne aussi le nom des stations)
  nst.par.gs <<- function(AS="A",impact=FALSE, allyr=FALSE) {
    # lorsque allyr=TRUE, calcule le nombre de campagnes effectuées sur chaque
    # géomorphologie +/- impact
    # sinon, calcule le nombre de stations *par campagne* sur chaque géomorpho +/- impact

    if(AS %in% c("A","S")) { dd <- filtreTable(dbio, wfiltre=paste("T",AS,"inv",sep="_"))
                             } else {
                               keySmstrl <- c(paste("A_",filtre.annees,sep=""),
                                              paste("S_",filtre.annees,sep=""))
                               dd <- dbio[dbio$Campagne %in% keySmstrl,] }
    ff <- c("Campagne","Geomorpho","N_Impact","St")
    # conserver seulement les facteurs spécifiés en arguments
    ff <- ff[c(1, 2, impact*3, (!allyr)*4)]
    dd2 <- merge(dd, info.transect, by="St")
    dd.St <- unique(dd2[,ff])

    ff <- ff[!(ff %in% c("St",ifelse(allyr, "Campagne","")))]
    dd.St2 <- aggregate(dd.St[,ifelse(allyr,"Campagne","St")],
                     lapply(ff, function(x) dd.St[,x]), length)
    names(dd.St2) <- c(ff, "N.u")
    return(list("numSt"=dd.St2, "nomsSt"=dd.St))
    }

  # Filtre les donnees analysees par groupe taxonomique lorsque specifie
  filtreTaxo <<- function(wtable,action="inclure",taxtype="Grp2",taxnom="tous") {

    if(any(taxnom %in% index.invSp[,taxtype])) {
      if(action == "inclure") {
      wtable <- wtable[wtable[,taxtype] %in% taxnom,]
    } else {
      wtable <- wtable[!(wtable[,taxtype] %in% taxnom),] }}

    if(action == "inclure" & taxtype == "Grp2" & any(taxnom %in% "tous")) {
      print("Pas de filtre sur espèces appliqué")
    } else { print(paste(c("Filtre sur espèces",capitalize(action),taxtype,
                           paste(taxnom,collapse=", ")),collapse=" :: ")) }
    
    return(wtable)
  }

  # Fonction interactive utilisée pour définir les variables du filtre sur les espèces
  # "inclure" ou "exclure" / unité taxonomique / nom
  filtre.especes <<- function(aF="tous") {

    if(aF == "tous") {
       taxoF.incl <<- "inclure"
       taxoF.utaxo <<- "Grp2"
       taxoF.nom <<- "tous"
     } else {
       print("Definition des filtres taxonomiques:")
    cat("Inclure ou exclure? ")
    taxoF.incl <<- tolower(readLines(file("stdin"),1))
    cat("Unité taxomique? (Groupe/Sous-Groupe/Famille/Genre/Espece) ")
    taxoF.utaxo <<- capitalize(tolower(readLines(file("stdin"),1)))
       if(taxoF.utaxo == "Groupe") taxoF.utaxo <<- "Grp2"
       if(taxoF.utaxo == "Sous-Groupe") taxoF.utaxo <<- "S_Grp2"
       if(taxoF.utaxo == "Espece") taxoF.utaxo <<- "G_Sp"
    cat("Nom? ")
    taxoF.nom <<- readLines(file("stdin"),1) }
    taxoF.nom <<- capitalize(tolower(trim(unlist(strsplit(taxoF.nom,",")))))
    closeAllConnections()
  }

  # Cette fonction retourne une version abbrégée des noms des groupes taxonomiques
  # inclus (ou exclus) au besoin
  taxotagFunk <<- function() {
    if(taxoF.incl=="inclure" & taxoF.utaxo == "Grp2" & any(taxoF.nom %in% "Tous")) {
      return("") # si tous les groupes sont inclus, ne pas modifier le nom de fichiers
    } else {
      taxotag <- paste(paste(abbreviate(c(capitalize(taxoF.incl), taxoF.nom)), collapse="-"),"_",sep="")
      return(taxotag) }
  }
    
  ###################################################
  ###################################################
  # Définir objets à mettre dans l'environnement global pour utilisation subséquente:
  info.transect <<- info.transect
  info.transect.INV <<- info.transect.INV
  info.transect.INV.geo <<- info.transect.INV.geo
  dbio <<- dbio
  index.invSp <<- index.invSp
  bioeco <<- bioeco
  bioeco.all <<- bioeco.all
  pr.Bacip <<- pr.Bacip
  
  # Tag TRUE when data read with no bugs
  data.read <<- TRUE
  wd.now <- dossier.R
  setwd(wd.now)
}
