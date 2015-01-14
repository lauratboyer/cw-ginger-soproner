## Analyses des données KNS (Ginger/Soproner)
# Auteur: Laura Tremblay-Boyer, contact: l.boyer@fisheries.ubc.ca
# Time-stamp: <2015-01-14 08:48:57 Laura>

# Sujet: Formattage des tableaux de données brutes pré-analyse,
# création de tableaux annexes + fonctions de base pour l'analyse

if(!exists("dossier.R")) {print("Faire \nsource(GS_KNS_MotherCode.r)\n
pour définir l'emplacement des codes et des dossiers")
                      }else{
}

prep.analyse <- function() {
  lp <- try(library(reshape)) # load package reshape
  if(class(lp)=="try-error") {
      install.packages("reshape") # installe reshape si requis
      library(reshape)}
  #######################################################################
  ###################### Formattage des tableaux ########################
  #######################################################################

  ####################################
  ###### Information transects #######
  ####################################
  # info.transect: informations sur les transects
  info.transect <- unique(data.info.transect) # extraire tableau brut
  info.transect$St <- toupper(info.transect$St)
  info.transect$cpus <- info.transect$Prod_ha # nouvelle colonne pour le cpus
  info.transect$Geomorpho <- trim(info.transect$Geom1) # oter espaces surperflus
  info.transect$Projet <- toupper(gsub("([A-Za-z]*_[A-Za-z]*)_.*","\\1",info.transect$Id))
  info.transect$An <- gsub("[A-Za-z]+_[A-Za-z]+_[0-9]+_([0-9]+)_.*","\\1",
                           info.transect$Id)
  info.transect$Mois <- gsub("[A-Za-z]*_[A-Za-z]*_([0-9]{2})_.*","\\1",
                           info.transect$Id)

  # Correction typos
  message("Info.transect ID correction typos IRB2B -> IBR2B")
  info.transect$Id <- gsub("irb2b","ibr2b",info.transect$Id)
  message("Info.transect Transect extrait du champ Id (et non de la colonne T vu quelques typos)")
  info.transect$T <- toupper(gsub(".*_([A-Za-z][0-9]{2})$","\\1",info.transect$Id))
  info.transect <- unique(info.transect)

  # Nettoyer accents noms géométrie -- à reviser vu encodage changé durant import?
  ########### Données LIT ############
  ####################################
  # data.LIT: données brutes LIT
  # index.LIT: info complementaires LIT
  data.LIT$Projet <- toupper(gsub("([A-Za-z]*_[A-Za-z]*)_.*",
                                       "\\1",data.LIT$Id))
  # ote rangées sans stations
  data.LIT <- data.LIT[data.LIT$St != "",]
  data.LIT$St <- toupper(data.LIT$St)
  data.LIT$Code_LIT<- toupper(data.LIT$Code_LIT)
  data.LIT$AQCQ <- toupper(data.LIT$AQCQ)

  # défini noms de colonnes pour index.LIT:
  # Automne 2014: Changement de format de index.LIT (?)
  # Code_DET -> All;
  names(index.LIT) <- c("Code_LIT","All","General","Forme","Acroporidae",
                        "Famille","Genre","Sensibilite")
  index.LIT$Code_LIT <- toupper(index.LIT$Code_LIT)

  message("Conversion de index.LIT/Sensibilité valeurs 1 et 2 à 'Corail sensible 1/2'")
  index.LIT$Sensibilite[index.LIT$Sensibilite == 1] <- "Corail sensible 1"
  index.LIT$Sensibilite[index.LIT$Sensibilite == 2] <- "Corail sensible 2"

 # créer catégories de substrat pour tableaux synthèses
  coraux.fig <- list("General"=c("Coraux","Corail mort","Coraux mous",
                     "Algues","Abiotique","Autre faune"),
                   "Acroporidae"=c("Acroporidae","Non-Acroporidae"),
                   "TS_All"=c("Coraux","Coraux morts","Coraux mous",
                     "Algues","Abiotique","Autre faune","Acroporidae",
                     "Non-acroporidae","Macro-Algues","Assemblage d'algues",
                     "Corail branchu","Corail tabulaire","Corail massif",
                     "Corail encroutant","Corail foliaire","Corail sub-Massif",
                     "Corail digite"))


  ###### Information biologie/écologie poissons #######
  #####################################################
  # bioeco.all: info complémentaires sur les poissons
  bioeco.all <- data.bioeco[,c("Code_Sp","famille","genre","espece","nom_total","taille_moy","commercial_2",
                        "cible_VKP","a","b","groupe_troph1","mobilite")]
  names(bioeco.all) <- c("Code_SP","Famille","Genre","Espece","G_Sp","Taille","Peche","Cible",
                       "a","b","Groupe_Trophique","Mobilite")
  message(sprintf("Doublons ôtés dans bioeco, Code_SP = %s",
                  paste(bioeco.all$Code_SP[duplicated(bioeco.all$Code_SP)],collapse=", ")))
  bioeco.all <- bioeco.all[!duplicated(bioeco.all$Code_SP),]

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

  # Ajuster si besoin le format pour les noms d'espèces:
  tl.cap <- function(x) capitalize(tolower(x))
  bioeco.all[,c("Famille","Genre","Espece","G_Sp")] <-
      sapply(bioeco.all[,c("Famille","Genre","Espece","G_Sp")], tl.cap)

  # Tableau bioeco, info a + b seulement
  bioeco <- bioeco.all[,c("Code_SP","a","b","GTlabel","moblabel","Peche","Cible")]

  # Vérifier s'il y a des espèces de poissons sans valeurs a ou b:
  abcheck <- length(unique(bioeco[bioeco$a=="Inconnu" | bioeco$b=="Inconnnu","Code_SP"]))
  if(abcheck!=0) message(sprintf("%s espèces sans valeur a ou b ôtées de l'analyse", abcheck))
  bioeco <- bioeco[bioeco$a!="Inconnu" & bioeco$b!="Inconnu",]
  bioeco$a <- as.numeric(bioeco$a); bioeco$b <- as.numeric(bioeco$b)
  rownames(bioeco) <- bioeco$Code_SP

  ####################################################################################
  #### Données de comptage Poissons ##################################################
  ####################################################################################
  # data.poissons: données de comptage brutes sur les poissons
  names(data.poissons) <- c("ID","Client","Site","Campagne","An","Mois","Date","St",
                            "Long.Transct","T","Obs","Vis","Courant","Code_SP",
                          "Famille","Genre","Espece","G_Sp","N","L","D1","D2",
                          "Secteur","Spatiaux","Tempo")
  # rajout année de la campagne
  data.poissons$Projet <- toupper(gsub("([A-Za-z]*_[A-Za-z]*)_.*",
                                       "\\1",data.poissons$ID))

  # temporary fix:
  message("Correction typo: data.poissons Kns_koniambo_11_2013_ibr3_t01, D2=Arus -> NA")
  data.poissons$D2[data.poissons$D2 == "Arus"] <- NA
  class(data.poissons$D2) <- "numeric"
  # créer tableau index espèce
  index.Poissons <<- unique(bioeco.all[,c("Code_SP","G_Sp", "Espece","Genre", "Famille")])
  index.Poissons$Groupe <<- index.Poissons$S_Groupe <<- NA

  ########################
  # Nettoyage du tableau:
  # Garder seulement les valeurs de N & L déclarées numériques
  data.poissons$N <- as.numeric(data.poissons$N)
  data.poissons$L <- as.numeric(data.poissons$L)

  # Ajuster si besoin le format pour les noms d'espèces:
  data.poissons[,c("Famille","Genre","Espece","G_Sp")] <-
      sapply(data.poissons[,c("Famille","Genre","Espece","G_Sp")], tl.cap)

  # Toutes les stations en majuscules
  data.poissons$St <- toupper(data.poissons$St)

  # Nettoyage du tableau: colonne D1/D2/N/L
  # Espèces sans Code_SP défini ôtées
  # Imprimer avertissements
  checkRangCodeSP <- nrow(data.poissons[data.poissons$Code_SP  %in% c("","#N/A"),])
  if(checkRangCodeSP != 0) message(sprintf("%s rangées sans valeur sous Code_SP", checkRangCodeSP))
  print(paste(nrow(data.poissons[is.na(data.poissons$D1),]),"valeurs D1 -> 0"))
  print(paste(nrow(data.poissons[is.na(data.poissons$D2),]),"valeurs D2 -> 0"))
  print(paste(nrow(data.poissons[is.na(data.poissons$N),]),"valeurs N -> 0"))
  print(paste(sum(is.na(data.poissons$L)),"rangees sans valeur L"))

  # Ôter rangées:
  data.poissons <- data.poissons[!(data.poissons$Code_SP %in% c("","#N/A")),]
  data.poissons$D1[is.na(data.poissons$D1)] <- 0
  data.poissons$D2[is.na(data.poissons$D2)] <- 0
  data.poissons$N[is.na(data.poissons$N)] <- 0
  #data.poissons <- data.poissons[!(is.na(data.poissons$L)),] # Taille

  # Expansion du data frame pour inclure toutes les combinaisons Campagne/St/Code_SP:
  # ... donc pour les poissons les densités sont nulles sur toutes Campagnes/St où
  # ... l'espèce est non-observée
  all.comb.poissons <- expand.grid("Campagne"=unique(data.poissons$Campagne),
                                   "St"=unique(data.poissons$St),
                                   "T"=unique(data.poissons$T),
                                   "Code_SP"=unique(data.poissons$Code_SP),
                                   stringsAsFactors=FALSE)

  # garder seulement les combinaisons distinctes de transect/station/campagne présentes dans les données:
  St.by.year <- unique(data.poissons[,c("Projet","An","Mois","Campagne","St","T","Long.Transct")]) # year/station samples
  all.comb.poissons <- merge(all.comb.poissons, St.by.year)

  # ... reunion avec données poissons pour rajouter les densité nulles
  poissons.coln <- c("Projet","An","Mois","Campagne","St","Code_SP","Long.Transct","T",
                     "N","L","D1","D2")

  data.poissons2 <- merge(all.comb.poissons,
                          data.poissons[,poissons.coln],
                          by=c("Projet","An","Mois","Campagne","St","Code_SP","Long.Transct","T"),all.x=TRUE)
  data.poissons2$N[is.na(data.poissons2$N)] <- 0 # abondance à zero si non-observée
  data.poissons3 <- merge(data.poissons2, index.Poissons[,c("Code_SP","G_Sp","Genre","Famille")], all.x=TRUE)

  # ... rajoute les coeff a et b pour calculs de biomasse
  print(sprintf("Code_SP manquants dans tableau bioeco et ôtés de l'analyse: %s",
                paste(unique(data.poissons3$Code_SP)[!(unique(data.poissons3$Code_SP) %in% bioeco$Code_SP)],
              collapse=", ")))
  data.poissons <- merge(data.poissons3, bioeco)

  # ... rajoute info.transect:
  dpois.tmp <- merge(data.poissons, unique(info.transect[,c("Projet","St","Geomorpho","N_Impact")]))
  data.poissons <- dpois.tmp

  ####################################################################################
  #### Données de comptage Invertébrés ###############################################
  ####################################################################################
  # data.inv: données de comptage brutes pour invertébrés
  data.inv$Projet <- toupper(gsub("([A-Za-z]*_[A-Za-z]*)_.*","\\1",data.inv$ID))
  dbio <- data.inv[,c("Projet","ID","Annee","Mois","Date","Campagne",
                      "St","T","Grp2","S_Grp2","F2","G2","G_Sp","N","D","Ltrans","L")]
  dbio$St <- toupper(data.inv$St)
  names(dbio) <- c("Projet","Id","Annee","Mois","Date","Campagne","St","T",
                   "Groupe","S_Groupe","Famille","Genre","G_Sp","N","D","Larg.Transct","Long.Transct")
  message(sprintf("On ôte %s rangées avec valeurs N absentes",sum(is.na(dbio$N))))
  dbio <- dbio[!is.na(dbio$N),]

  ############################
  # nettoyage:
  # garde seulement les rangées avec transects A, B, C
  nTnd <- nrow(dbio[!(dbio$T %in% c("A","B","C","T01","T02","T03","T04","T05")),])
  dbio <- dbio[dbio$T %in% c("A","B","C","T01","T02","T03","T04","T05"),];
  message(paste("Approx",nTnd,"rangées ôtées vu transect non égal à A, B, C; T01 --> T05"))

  ############################
  # abondances nulles:
  # commencer par oter les rangées N=0 pour eliminer les stations
  # où l'espèce n'est observée sur aucun T
  dbio <- dbio[dbio$N > 0,]
  # collater les rangées où la même espèce est observée plusieurs fois
  # sur le même transect
  dbio.tmp <- aggregate(list(N=dbio$N),
                        as.list(dbio[,c("Id","Projet","Campagne","St","T","Groupe","S_Groupe",
                               "Famille","Genre","G_Sp","Larg.Transct","Long.Transct")]),
                        sum)
  dbio <- dbio.tmp

  # rajouter les transects N=0 *seulement* lorsque l'espèce
  ### est observée sur la ST mais pas tous les transects
  # créer identifiant unique pour chaque combinaison campagne/st/espèce observée
  dbio$uID <- paste(dbio$Projet,dbio$Campagne,dbio$St,dbio$G_Sp,sep="_")
  # on commence par identifier les espèces identifiées au moins une fois
  ### par station/campagne:
  dbio.uID <- unique(dbio[,c("Projet","Campagne","St","G_Sp","uID")])
  # on crée un tableau avec les noms des transects échantillonés
  ### sur chaque combinaison stations x
  ### campagne (vu que certaines campagnes ont T02 et T04)
  dbio.allT <- unique(dbio[,c("Id","Projet","Campagne","St","T")])
  # on merge ces deux tableaux ensemble pour associer des densités
  ### nulles aux espèces identifiées
  # sur une station mais pas sur tous les transects
  allT.uID <- merge(dbio.allT, dbio.uID, all.x=TRUE)
  index.dbio <- unique(dbio[,c("uID","Projet","Campagne","St","Groupe","S_Groupe",
                               "Famille","Genre","G_Sp","Larg.Transct","Long.Transct")])
  dbio.tmp <- merge(dbio[,c("uID","T","N")],allT.uID,by=c("uID","T"),all.y=TRUE)
  dbio.tmp$N[is.na(dbio.tmp$N)] <- 0 # valeurs NA remplacées par 0 pour remplir les transects manquants
  dbio <- merge(dbio.tmp, index.dbio)

  ################
  # Fix Ltrans that are NA
  if((sum(is.na(dbio$Larg.Transct)) > 0 )|(sum(is.na(dbio$Long.Transct)) > 0 )) {
    message("Attention certaines largeurs ou longueurs de transect sont manquantes")
  }

  # Calculer la densité en hectares
  # Individus/hectares
  # Larg.transct et Long.transct en mètres
  dbio$D <- dbio$N/(dbio$Larg.Transct*dbio$Long.Transct) * 10000

  ########################
  ### Tableaux annexes ###

  # clés noms des invertébrés
  InvGr.Key <- c("Algues","Ascidies","Cnidaires","Crustaces","Echinodermes",
                 "Eponges","Mollusques","Phanerogame","Vers")

  # tableau index pour la taxonomie de toutes les espèces observée:
  index.invSp <- unique(dbio[,c("G_Sp","Genre","Famille","S_Groupe","Groupe")])
  # verifier s'il y a des doubles entrées dans les informations taxonomiques
  if(length(unique(index.invSp$G_Sp)<nrow(index.invSp))) {
      InvDup <- index.invSp$G_Sp[duplicated(index.invSp$G_Sp)]
      message(sprintf("!!!!! Attention erreurs d'entrées dans les champs taxonomique pour %s",
                      paste(InvDup, collapse=", ")))

      # Conserver dans le tableau la première ligne observée pour l'espèce seulement,
      # i.e. l'erreur pourrait être dans le tableau index -- la typo dans la base doit
      # être corrigée à la source pour éviter ce problème
      index.invSp <- index.invSp[!(duplicated(index.invSp$G_Sp)),] # on ôte les doublons
  }

  message("\n!! Corrections temporaires en attendant que la base INV soit nettoyée:\n Genre Plakobranchus -> Famille Plakobranchidae")
  index.invSp$Famille[index.invSp$Genre == "Plakobranchus"] <- "Plakobranchidae"
  message("!! Corrections temporaires en attendant que la base INV soit nettoyée:\n Famille Phyllidiidae -> S_Groupe Nudibranches\n")
  index.invSp$S_Groupe[index.invSp$Famille == "Phyllidiidae"] <- "Plakobranchidae"

  # maintenant que le tableau a une ligne unique par espèce on crée un index en définissant les
  # noms de rangées
  row.names(index.invSp) <- index.invSp$G_Sp

  ################
  # infos sur station/mission/année
  info.transect.INV <- unique(data.inv[,c("Annee","Mois","Mission","Campagne","St")])
  info.transect.INV$St <- toupper(info.transect.INV$St)
  info.transect.INV <- merge(info.transect.INV, info.transect,
                             by=c("St","Mois"))[,c("Annee","Mois","Mission","Campagne","St",
                               "Geomorpho","N_Impact")]
  info.transect.INV.geo <- aggregate(info.transect.INV[,c("Mois", "Annee")],
                                     as.list(info.transect.INV[,c("Mission","Campagne","Geomorpho")]), lastval)
  message("-> Lorsqu'un transect ou Campagne est échantilloné dans 2 mois différents, 1 seul mois/année est gardé dans le tableau référence")

  # rajouter les infos transects au tableau dbio principal:
  dbio.tmp <- merge(dbio, unique(info.transect[,c("Projet","St","Geomorpho","N_Impact")]))

  if(nrow(dbio.tmp)<nrow(dbio)) {
    miss.St <- unique(dbio$St)[!(unique(dbio$St) %in% info.transect$St)]
    message(sprintf("\n***********************\nStations manquantes dans tableau info.transect: %s",
                    paste(miss.St,collapse=", ")))
  }
  dbio <- dbio.tmp

  ####################################
  #### Tableau espèces universels ####
  ####################################
  index.Poissons$type <- "Poissons"
  index.invSp$type <- "Inv"
  index.tt.especes <- rbind(index.invSp, index.Poissons[,-c(1,3)])

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
    # campagnes annuelles et semestrielles
    keySmstrl <- c(paste("A_",qAnnees,sep=""),paste("S_",qAnnees,sep=""))
    tb3 <- tb2[,unlist(sapply(keySmstrl,function(x) grep(x,names(tb2))))]
    wStat <- tb2$St[rowSums(tb3,na.rm=TRUE)==ncol(tb3)]
    T_AeS_inv <- apply(expand.grid(wStat,names(tb3)),1,paste,collapse="_")

    # campagnes annuelles
    keyAnnuel <- paste("A_",qAnnees,sep="")
    tb3 <- tb2[,unlist(sapply(keyAnnuel,function(x) grep(x,names(tb2))))]
    if("data.frame" %in% class(tb3)) { wStat <- tb2$St[rowSums(tb3,na.rm=TRUE)==ncol(tb3)]
                                   }else{
                                     wStat <- tb2$St[which(na.omit(tb3==1))]}
    T_A_inv <- apply(expand.grid(wStat,names(tb3)),1,paste,collapse="_")

    # campagnes semestrielles
    keyAnnuel <- paste("S_",qAnnees,sep="")
    tb3 <- tb2[,unlist(sapply(keyAnnuel,function(x) grep(x,names(tb2))))]
    if("data.frame" %in% class(tb3)) { wStat <- tb2$St[rowSums(tb3,na.rm=TRUE)==ncol(tb3)]
                                   }else{
                                     wStat <- tb2$St[which(na.omit(tb3==1))]}
    T_S_inv <- apply(expand.grid(wStat,names(tb3)),1,paste,collapse="_")

    return(list("T_S_inv"=T_S_inv, "T_A_inv"=T_A_inv, "T_AeS_inv"=T <- T_AeS_inv))
  }

  filtre.Camp <<- creerFiltre(filtre.annees)

  ###################################################
  ######## Fonctions génériques #####################
  # Définition de fonctions qui seront utilisées couramment dans le code

  # 1. Applique le filtre spécifié au tableau donné en argument
  filtreTable <<- function(wtable, wfiltre) {
    if(wfiltre %in% c("T_A_inv","T_S_inv","T_AeS_inv")) { # appliquer le filtre si spécifié
    message(paste("Stations filtrees par",wfiltre))
    filtre.Camp <<- creerFiltre(filtre.annees) # recreer le tableau filtre

    wtable$key <- paste(wtable$St, wtable$Campagne, sep="_")
    dd.filt <- merge(data.frame("key"=filtre.Camp[[wfiltre]]),
                     wtable,by="key", drop.x="key")
    dd.filt <- dd.filt[,names(dd.filt) != "key"] #ôter colonne key
  } else {
    CmpTag <- paste(filtre.annees,collapse="|")
    wCampKeep <- grep(CmpTag, unique(wtable$Campagne), value=TRUE)
    dd.filt <- wtable[wtable$Campagne %in% wCampKeep,]
    }

    if(nrow(dd.filt)==0) warning(
           sprintf("Attention!! \n\nAucune des données ne sont sélectionnées par le filtre sur stations:
\nFiltre %s, années %s\n", wfiltre, paste(filtre.annees,collapse=", ")))
    return(dd.filt)  }

  # 2. Converti les noms de campagne en année (charactère -> numérique)
  as.year <<- function(x) as.numeric(sub("[AS]_","",x))

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
  filtreTaxo <<- function(wtable,action="inclure",
                          taxtype="Groupe",taxnom="Tous") {


      if(action == "inclure" & taxtype == "Groupe" & any(taxnom %in% "Tous")) {
      # Pas de filtre
      message("Pas de filtre sur espèces appliqué")
          } else {

      # Extraction du type d'espèces analysées maintenant
      sc <<- sys.calls() # fonctions en cours

      # Extraction de la catégorie d'espèces présente dans le filtre
      # E.g. si toutes les espèces filtrées sont des invertébrés, le filtre
      # taxonomique n'est pas appliqué aux poissons
      qtype <- unique(index.tt.especes[index.tt.especes[,taxoF.utaxo]
                                       %in% taxoF.nom,"type"])

      # Appliquer le filtre si la categorie d'espèce mentionnée dans le filtre
      # correspond aux analyses présentes (e.g. pas de filtre poissons sur les
      # invertébrés)
      if(any(grepl(tolower(qtype),tolower(sc)))) {

      # On filtre les rangees selon les especes/groupes taxo specifie
      message(paste(c("Filtre sur espèces",capitalize(action),qtype,taxtype,
                    paste(sort(taxnom),collapse=", ")),collapse=" :: "))

      if(tolower(action) == "inclure") { # inclusion des groupes X
      wtable <- wtable[wtable[,taxtype] %in% taxnom,]
      } else { # exclusion des groupes X
      wtable <- wtable[!(wtable[,taxtype] %in% taxnom),] }}}

    return(wtable)
  }

  # Fonction interactive utilisée pour définir les variables du filtre sur les espèces
  # "inclure" ou "exclure" / unité taxonomique / nom
  def.filtre.especes <<- function(aF="tous") {

    if(aF == "tous") {
       taxoF.incl <<- "inclure"
       taxoF.utaxo <<- "Groupe"
       taxoF.nom <<- "Tous"
     } else {
       message("Definition des filtres taxonomiques (pas besoin de mettre des guillemets):")
    taxoF.incl <<- tolower(readline("Inclure ou exclure? "))
    mm <- "Unité taxomique? (Groupe/Sous-Groupe/Famille/Genre/Espece) "
    taxoF.utaxo <<- capitalize(tolower(readline(mm)))
       if(taxoF.utaxo == "Sous-Groupe") taxoF.utaxo <<- "S_Groupe"
       if(taxoF.utaxo == "Espece") taxoF.utaxo <<- "G_Sp"
    taxoF.nom <<- readline("Nom? ")
    taxoF.nom <<- capitalize(tolower(trim(unlist(strsplit(taxoF.nom,",")))))
  }
}
  # Fonction qui permet de sélectionner un fichier .csv pour importer
  # sous R une liste de noms d'espèces/genre/famille
  # à utiliser dans le filtre taxonomique
  # argument "action" défini taxoF.incl
  # argument "niveau" défini le niveau taxonique, taxoF.utaxo
  # argument "titre" défini si la première rangée est le nom de la colonne dans le
  # ... fichier .csv (si non, spécifier titre = FALSE)
  import.filtre.taxo <<- function(niveau="Famille",action="inclure",titre=TRUE,sepval=";") {

      lnoms <- read.csv(file.choose(),header=titre,sep=sepval) # sélectionner le fichier dans l'ordi
      if(class(lnoms)=="data.frame") lnoms <- lnoms[,1]
      lnoms <- capitalize(tolower(trim(lnoms))) # nettoyer format des noms
      taxoF.incl <<- action
      taxoF.utaxo <<- niveau
      taxoF.nom <<- lnoms

      voir.filtre.taxo()
  }

  # Fonction qui montre la valeur présente des filtres taxonomiques
  # Tapez "voir.filtre.taxo()" dans la console
  voir.filtre.taxo <<- function() {
      fltre.now <- list(taxoF.incl, taxoF.utaxo, sort(taxoF.nom))
      names(fltre.now) <- c("taxoF.incl","taxoF.utaxo","taxoF.nom")
      print(fltre.now) }

  # Cette fonction retourne une version abbrégée des noms des groupes taxonomiques
  # inclus (ou exclus) au besoin
  taxotagFunk <<- function() {

    if(taxoF.incl=="inclure" & taxoF.utaxo == "Groupe" & any(taxoF.nom %in% "Tous")) {
      return("") # si tous les groupes sont inclus, ne pas modifier le nom de fichiers
    } else {
        if(length(taxoF.nom) <= 5) {
            taxotag <- paste("_",paste(c(capitalize(taxoF.utaxo),abbreviate(c(capitalize(taxoF.incl),
                                                                          sort(taxoF.nom)))), collapse="-"),"_",sep="")
        } else {taxotag <- paste("_",paste(c(capitalize(taxoF.utaxo),
                                         abbreviate(c(capitalize(taxoF.incl),
                                                      sort(taxoF.nom)[1:5], "etc",
                                                      paste("N=",length(taxoF.nom),sep="")))), collapse="-"),"_",sep="")}
        return(taxotag) }
  }

  # Ces fonctions sont utilisées pour indiquer le départ et la fin
  # des codes compris dans une fonction
  departFunk <<- function() {
      sc <- sys.calls()
      fname <- as.character(sc[[length(sc)-1]])[1] # extrait le nom de la fonction parent
      packageStartupMessage(sprintf("\nDepart %s()...",fname))
      try(print(unlist(mget(names(formals(fname)), envir=sys.frame(-1)))),
          silent=TRUE)
      message("#################\n")
      }
  finFunk <<- function() {
      message("\n#################")
      sc <- sys.calls()
      fname <- as.character(sc[[length(sc)-2]])[1] # extrait le nom de la fonction parent
      packageStartupMessage(sprintf("---> fin %s().",fname)) }

  ## Si erreur dans une fonction, indiquer la fonction où l'erreur se produit
  ## EM() est donnée en argument a la fonction on.exit() qui tourne automatiquement
  ## l'argument spécifié quand une fonction se termine, naturellement ou avec
  ## une erreur. Ici EM() identifie si la fonction a eu une erreur, et si c'est
  ## le cas imprime le nom de la fonction pour faciliter l'identification du bug.
  EM <<- function() {
      tb <- try(get(".Traceback",envir=baseenv()),silent=TRUE) # extraire messages d'erreur
      if(identical(paste(lastval(tb)), paste(sys.calls()[1]))) {
          message(sprintf("Erreur dans la fonction %s()",
                      paste(sys.calls()[[1]][1])))
          assign(".Traceback"[[1]],999,envir=baseenv())
          } else { finFunk()} }

  ###################################################
  ###################################################
  # Définir objets à mettre dans l'environnement global pour utilisation subséquente:
  info.transect.TProj <<- info.transect
  info.transect.INV <<- info.transect.INV
  info.transect.INV.geo <<- info.transect.INV.geo
  dpoissons.TProj <<- data.poissons
  dbio.TProj <<- dbio
  index.invSp <<- index.invSp
  bioeco <<- bioeco
  bioeco.all <<- bioeco.all
  data.LIT.TProj <<- data.LIT
  index.LIT <<- index.LIT
  coraux.fig <<- coraux.fig
  index.tt.especes <<- index.tt.especes

  # Tag TRUE when data read with no bugs
  data.read <<- TRUE
  wd.now <- dossier.R
  setwd(wd.now)

  message("\n\n***********************\nFonction prep.analyse() complétée: tableaux formattés et analyses prêtes à lancer.
Pour lancer les analyses manuellement utiliser les fonctions:\n
-- Invertébrés: Run.INV.biodiv(), Run.INV.densite() \n-- LIT: Run.LIT.all() \n-- Poissons: Run.poissons.all()\n\n***********************\n")
}

##################### Fin de la fonction prep.analyse() #####################


### Fonctions utiles pour formattage
if(!exists("capitalize")) { # rajoute une lettre majuscule au debut
  capitalize <- function(x) {
    gsub('(\\w)([\\w|\\s]*)','\\U\\1\\L\\2',x,perl=TRUE) }}

if(!exists("getmatch")) { # extrait la partie correspondante au match
      getmatch <- function(x,str2match,...) {
          unlist(regmatches(x,regexpr(str2match,x,...))) }}

# dernier element de l'objet
lastval <- function(x) x[length(x)]

# Ote les espaces au debut et apres un charactere
# e.g. "  Famille " -> "Famille"
if(!exists("trim")) { #clumsy regexp
  trim <- function(x) gsub("^\\s+","",gsub("\\s+$","",x)) }

# Vérifie que la graphic device a la bonne grandeur,
# sinon en ouvrir une nouvelle selon les dimensions ww x hh
check.dev.size <- function(ww,hh) {

    if(dev.cur()==1){ dev.new(width=ww,height=hh)
                  } else {
    ds <- dev.size()
    if(round(ds[1],2)!=round(ww,2)
       | round(ds[2],2)!=round(hh,2)) {
        dev.off(); dev.new(width=ww,height=hh)} }
}

# Standard error
stand.err <<- function(x) sd(x)/sqrt(length(x))

# Nombre de valeur unique
count <<- function(x) length(unique(x))
