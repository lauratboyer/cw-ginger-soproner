## Analyses des données KNS (Ginger/Soproner)
# Auteur: Laura Tremblay-Boyer, contact: l.boyer@fisheries.ubc.ca
# Time-stamp: <2015-05-28 08:12:33 Laura>

# Sujet: Formattage des tableaux de données brutes pré-analyse,
# création de tableaux annexes + fonctions de base pour l'analyse

if(!exists("dossier.R")) {print("Faire \nsource(GS_KNS_MotherCode.r)\n
pour définir l'emplacement des codes et des dossiers")
                      }else{
}

prep.analyse <- function(check.typo=TRUE) {

    if(refaire.tableaux) { # défini dans mother_code
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
  # Ote Id ADECAL/CAGES EXTERIEURES (cf. mail Tom H 6/5/2015)
  info.transect <- filter(info.transect, !grepl("ADECAL.*CAGES_[E|I]", Id))
  info.transect <- info.transect[!duplicated(info.transect$Id),] # ôte les doublons restants
  ID.util.Non <-   filter(info.transect, Utilise.analyse == "Non")$Id
  info.transect <- filter(info.transect, Utilise.analyse == "Oui")
  rownames(info.transect) <- info.transect$Id

  ############################################
  ########### Infos temporelles ##############
  info.transect.tempo <- unique(data.info.temprl)
  rm.all.na <- function(x) !all(is.na(x))
  info.transect.tempo <- info.transect.tempo[,sapply(info.transect.tempo, rm.all.na)]
  # pareil pour facteurs temporels
  info.transect.tempo <- filter(info.transect.tempo, !grepl("ADECAL.*CAGES_[E|I]", Id))
  info.transect.tempo <- info.transect.tempo[!duplicated(info.transect.tempo$Id),] # ôte les doublons restants
  ID.util.Non <-  unique(c(ID.util.Non, filter(info.transect.tempo, Utilise.analyse == "Non")$Id))
  info.transect.tempo <- filter(info.transect.tempo, Utilise.analyse == "Oui")
  info.spat.temp <- merge(info.transect, info.transect.tempo, by=c("Id","Client","Site"))
  rownames(info.spat.temp) <- info.spat.temp$Id

  message(sprintf("%s IDs exclus de l'analyse", length(ID.util.Non)))
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
  data.LIT <- data.LIT[toupper(data.LIT$AQCQ)=="NON",] # on garde seulement AQCQ = NON
  data.LIT <- data.LIT[data.LIT$X. > 0,] # on ôte les zéros vu qu'ils ne sont pas entrés partout
  # consolider les observations multiples du même Code_LIT sur un transect
  data.LIT <- aggregate(list(X.=data.LIT$X.),
                        data.LIT[,c("Projet","Campagne","St","Id","T","Code_LIT")], sum)

  # rajouter des abondances nulles pour les groupes observés sur certains
  # transects seulement
  # catégories identifiées au moins une fois sur la station:
  LIT.uID <- unique(data.LIT[,c("Projet","St","Code_LIT")])
  # on crée un tableau avec les noms des transects échantillonés
  ### sur chaque combinaison stations x
  ### campagne (e.g. vu que certaines campagnes ont T02 et T04, ou T09 manquant)
  LIT.allT <- unique(data.LIT[,c("Id","Projet","Campagne","St","T")])
  # on merge ces deux tableaux ensemble pour associer des densités
  ### nulles aux espèces identifiées
  # sur une station mais pas sur tous les transects
  LIT.allT.uID <- merge(LIT.allT, LIT.uID, all.x=TRUE)

  dtmp <- merge(LIT.allT.uID,data.LIT,  all.x=TRUE);
  data.LIT <- dtmp
  data.LIT[is.na(data.LIT$X.),"X."] <- 0

  # défini noms de colonnes pour index.LIT:
  # Automne 2014: Changement de format de index.LIT (?)
  # Code_DET -> All;
  names(index.LIT) <- c("Code_LIT","All","General","Forme","Acroporidae",
                        "Famille","Genre","Sensibilite")
  index.LIT$Code_LIT <- toupper(index.LIT$Code_LIT)

  message("Conversion de index.LIT/Sensibilité valeurs 1 et 2 à 'Corail sensible 1/2'")
  index.LIT$Sensibilite[index.LIT$Sensibilite == 1] <- "Corail sensible 1"
  index.LIT$Sensibilite[index.LIT$Sensibilite == 2] <- "Corail sensible 2"

  # utilise Code_LIT pour indexer le tableau
  rownames(index.LIT) <- index.LIT$Code_LIT

  ###########################################
  ## Formatter données quadrats:
  data.quad <- data.quad[,c("Id","Campagne","St","Code_LIT","Quadrat","X.")]
  data.quad$Projet <- id2proj(data.quad$Id)
  data.quad[,c("Campagne","St","Code_LIT")] <-
    sapply(data.quad[,c("Campagne","St","Code_LIT")], toupper)
  message("Quadrats: ôte rangées avec Code_LIT = Q03, pas de couverture (typo)")
  data.quad <- data.quad[!(grepl("Q0.", data.quad$Code_LIT)),]
  data.quad$X. <- as.numeric(data.quad$X)


  # on ôte les zéros vu qu'ils ne sont pas entrés partout
  data.quad <- data.quad[data.quad$X. > 0,]
  # consolider les observations multiples du même Code_LIT sur un transect
  data.quad <- aggregate(list(X.=data.quad$X.),
                        data.quad[,c("Projet","Campagne","St","Id","Quadrat","Code_LIT")], sum)

  # rajouter des abondances nulles pour les groupes observés sur certains
  # quadrats seulement
  # catégories identifiées au moins une fois par station/campagne:
  quad.uID <- unique(data.quad[,c("Projet","Campagne","St","Code_LIT")])
  # on crée un tableau avec les noms des transects échantillonés
  ### sur chaque combinaison stations x
  ### campagne (vu que certaines campagnes ont T02 et T04)
  quad.allT <- unique(data.quad[,c("Id","Projet","Campagne","St","Quadrat")])
  # on merge ces deux tableaux ensemble pour associer des densités
  ### nulles aux espèces identifiées
  # sur une station mais pas sur tous les transects
  quad.allT.uID <- merge(data.quad, quad.uID, all.x=TRUE)

  dtmp <- merge(quad.allT.uID, data.quad, all.x=TRUE)
  data.quad <- dtmp
  data.quad[is.na(data.quad$X.),"X."] <- 0

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
  names(data.poissons) <- c("Id","Client","Site","Campagne","An","Mois","Date","St",
                            "Long.Transct","T","Obs","Vis","Courant","Code_SP",
                          "Famille","Genre","Espece","G_Sp","N","L","D1","D2",
                          "Secteur","Spatiaux","Tempo")
  # rajout année de la campagne
  data.poissons$Projet <- toupper(gsub("([A-Za-z]*_[A-Za-z]*)_.*",
                                       "\\1",data.poissons$Id))

  # temporary fix:
  message("Correction typo: data.poissons Kns_koniambo_11_2013_ibr3_t01, D2=Arus -> NA")
  data.poissons$D2[data.poissons$D2 == "Arus"] <- NA
  class(data.poissons$D2) <- "numeric"
  # créer tableau index espèce
  index.Poissons <- unique(bioeco.all[,c("Code_SP","G_Sp", "Espece","Genre", "Famille")])
  index.Poissons$Groupe <- index.Poissons$S_Groupe <- NA
  rownames(index.Poissons) <- index.Poissons$Code_SP # indéxer les rangées sur le code espèce

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

                                        # remettre les valeurs de transect au cas où certaines manquaient
  unique.sans.na <- function(x) unique(na.omit(x))
  proj.ltrans <- tapply(data.poissons$Long.Transct, data.poissons$Projet, unique.sans.na)
  if(length(unlist(proj.ltrans))!=length(proj.ltrans)) stop("Longueurs de transect différentes par projet")
  if(any(is.na(data.poissons$Long.Transct))) stop("Longueurs de transects manquantes pour les poissons")

  # Expansion du data frame pour inclure toutes les combinaisons Campagne/St/Code_SP:
  # ... donc pour les poissons les densités sont nulles sur toutes Campagnes/St où
  # ... l'espèce est non-observée
  proj.poissons <- unique(data.poissons[,c("Projet","Code_SP")])
  trans.poissons <- unique(data.poissons[,c("Projet","Campagne","St","T","Id","Long.Transct")])
  pt.all <- merge(proj.poissons, trans.poissons, all=TRUE)

  # ... reunion avec données poissons pour rajouter les densité nulles
  poissons.coln <- c("Id","Projet","Campagne","St",
                     "Code_SP","Long.Transct","T","N","L","D1","D2")
  data.poissons2 <- merge(pt.all, data.poissons[,poissons.coln], all.x=TRUE)
  data.poissons2$N[is.na(data.poissons2$N)] <- 0 # abondance à zero si non-observée
  data.poissons2$D1[is.na(data.poissons2$D1)] <- 0 # ... et les autres paramètres aussi
  data.poissons2$D2[is.na(data.poissons2$D2)] <- 0 # ... pour que biomasse et densité
  data.poissons2$L[is.na(data.poissons2$L)] <- 0 # ... soit à zero

  data.poissons3 <- merge(data.poissons2, index.Poissons[,c("Code_SP","G_Sp","Genre","Famille")], all.x=TRUE)
  # ... rajoute les coeff a et b pour calculs de biomasse
  print(sprintf("Code_SP manquants dans tableau bioeco et ôtés de l'analyse: %s",
                paste(unique(data.poissons3$Code_SP)[!(unique(data.poissons3$Code_SP) %in% bioeco$Code_SP)],
              collapse=", ")))
  data.poissons <- merge(data.poissons3, bioeco)

  # rajouter la densité/biomasse par *observation*
  # D1 et D2 sont les distances entre la ligne du transect et l'individu observé
  # Lorsqu'il y a un banc de poissons D1 est la distance de l'individu
  # ... le plus proche, D2 celle de l'indiv le plus loin
  # Sinon D1 = D2
  # Donc 0.5*(D1 + D2) est la distance moyenne, à laquelle on rajoute 0.5 (facteur de correction?)
  # L'aire du transect pour l'espèce i est la largeur x la longueur; largeur = 2 x distance moyenne
  # ... la densité est donc N / aire
  # ... et la biomasse: biomasse totale / aire
  # D.obs = N/(2xdm.obsxL); dm.obs= 0.5*(d1+d2) + 0.5
  # Bio.obs = N*a*T^b/(2xdm.obsxL)
  data.poissons$dm.obs <-  0.5*(data.poissons$D1+data.poissons$D2)+0.5 # distance moyenne observation
  data.poissons$dens.obs <- data.poissons$N/(2*data.poissons$dm.obs*data.poissons$Long.Transct) # densité.observation
  data.poissons$bio.obs <- data.poissons$N*data.poissons$a*(data.poissons$L^data.poissons$b)/(2*data.poissons$dm.obs*data.poissons$Long.Transct)

  ####################################################################################
  #### Données de comptage Invertébrés ###############################################
  ####################################################################################
  # data.inv: données de comptage brutes pour invertébrés
  data.inv$Projet <- toupper(gsub("([A-Za-z]*_[A-Za-z]*)_.*","\\1",data.inv$ID))
  dbio <- data.inv[,c("Projet","ID","Date","Campagne",
                      "St","T","Grp2","S_Grp2","F2","G2","G_Sp","N","D","Ltrans","L")]
  dbio$St <- toupper(data.inv$St)
  names(dbio) <- c("Projet","Id","Date","Campagne","St","T",
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

  ##############################
  # on créé l'index des espèces invertébrés (tableau annexe)

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

  ###########################
  ## Densités nulles
  ### Rajouter les transects N=0 *seulement* lorsque l'espèce
  ### est observée sur la ST mais pas tous les transects
  # créer identifiant unique pour chaque combinaison campagne/st/espèce observée
  dbio$uID <- paste(dbio$Projet, dbio$G_Sp, sep="_")

  # on commence par identifier les espèces identifiées au moins une fois
  ### sur un projet:
  dbio.uID <- unique(dbio[,c("uID","Projet","G_Sp")])

  # longueur/largeur transect, par transect/espèces vu valeur non-uniques sur la même station
  # en S_2010
  ll.transect <- unique(dbio[,c("Projet","Campagne","St","T","G_Sp","Larg.Transct","Long.Transct")])

  # on crée un tableau avec les noms de tous les transects échantillonés
  ### sur chaque combinaison stations x
  ### campagne (vu que certaines campagnes ont T02 et T04)
  dbio.allT <- unique(dbio[,c("Id","Projet","Campagne","St","T")])

  # on merge ces deux tableaux ensemble pour associer des densités
  ### nulles aux espèces identifiées
  # sur une station mais pas sur tous les transects
  allT.uID <- merge(dbio.allT, dbio.uID, all.x=TRUE)
  dbio.tmp <- merge(dbio[,c("uID","Campagne","St","T","N")],
                    allT.uID, all.y=TRUE)
  dbio.tmp$N[is.na(dbio.tmp$N)] <- 0 # valeurs NA remplacées par 0 pour remplir les transects manquants

  ## on rajoute le reste des infos taxonomiques via le tableau index invertébrés
  dbt <- merge(dbio.tmp, index.invSp)
  dbio <- merge(dbt, ll.transect, all.x=TRUE)

  ################
  # Fix Ltrans that are NA
  if((sum(is.na(dbio$Larg.Transct)) > 0 )|(sum(is.na(dbio$Long.Transct)) > 0 )) {
    message("Attention certaines largeurs ou longueurs de transect sont manquantes")
  }

  # Calculer la densité en hectares
  # Individus/hectares
  # Larg.transct et Long.transct en mètres
  dbio$D <- dbio$N/(dbio$Larg.Transct*dbio$Long.Transct) * 10000
  dbio$D[is.na(dbio$D)] <- 0

  ################
  # Rajouter les infos transects aux tableaux de données principaux:
  add.infos <- function(x) {

      id.manq <- unique(x$Id)[!(unique(x$Id) %in% c(info.spat.temp$Id, ID.util.Non))]
    if(length(id.manq)>0) {
    message(sprintf("Id manquants dans info.transect: %s",paste(id.manq,collapse=", ")))
  }
    id2k <- x$Id[(x$Id %in% info.spat.temp$Id)]
    df <- info.spat.temp[x$Id,c("An","Mois",facteurs.spatio,facteurs.tempo)]
    df <- data.frame(x, df, row.names=NULL)
    df <- df[!is.na(df$An),] # ôte les rangées sans Id dans info.spat.temp
    # incluant Utilise.analyse == 'Non'
  }

  dbio <- add.infos(dbio) # invertébrés
  data.poissons <- add.infos(data.poissons) # poissons
  data.LIT <- add.infos(data.LIT) # LIT

  ####################################
  #### Tableau espèces universels ####
  ####################################
  index.Poissons$type <- "Poissons"
  index.invSp$type <- "Inv"
  index.tt.especes <- rbind(index.invSp, index.Poissons[,-c(1,3)])


  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###
  ### Faire rapport d'erreurs de frappe possible
  if(refaire.tableaux & check.typo) {
    typo.finder(dbio$G_Sp, tag="Invertébrés, G_Sp")
    typo.finder(dbio$Genre, tag="Invertébrés, Genre")
    typo.finder(dbio$Famille, tag="Invertébrés, Famille")
    typo.finder(dbio$S_Groupe, tag="Invertébrés, S_Groupe")
    typo.finder(dbio$Groupe, tag="Invertébrés, Groupe")
    typo.finder(data.poissons$G_Sp, tag="Poissons, G_Sp")
    typo.finder(data.poissons$Genre, tag="Poissons, Genre")
}

  ###################################################
  ###################################################
                                        # Définir objets à mettre dans l'environnement global pour utilisation subséquente:
    # import filtre famille
    if(filtre.famille) {
        message("filtre.famille=TRUE --> chargement automatique du filtre famille Filtre-taxo_Famille.csv")
        import.filtre.taxo("Filtre-taxo_Famille.csv","Famille")
    }

  if(!any(grepl("Tableaux-pour-analyses", list.dirs()))) dir.create("Tableaux-pour-analyses")
  setwd("Tableaux-pour-analyses")
  info.transect.TProj <<- info.spat.temp
  save(info.transect.TProj, file="info-transect-TProj.RData")
  dpoissons.TProj <<- data.poissons
  save(dpoissons.TProj, file="dpoissons-TProj.RData")
  dbio.TProj <<- dbio
  save(dbio.TProj, file="dbio-TProj.RData")
  index.invSp <<- index.invSp
  save(index.invSp, file="index-invSp.RData")
  index.Poissons <<- index.Poissons
  save(index.Poissons, file="index-Poissons.RData")
  bioeco <<- bioeco
  save(bioeco, file="bioeco.RData")
  bioeco.all <<- bioeco.all
  save(bioeco.all, file="bioeco-all.RData")
  data.LIT.TProj <<- data.LIT
  save(data.LIT.TProj, file="data-LIT-TProj.RData")
  index.LIT <<- index.LIT
  save(index.LIT, file="index-LIT.RData")
  dQuad.TProj <<- data.quad
  save(dQuad.TProj, file="dQuad-TProj.RData")
  coraux.fig <<- coraux.fig
  save(coraux.fig, file="coraux-fig.RData")
  index.tt.especes <<- index.tt.especes
  save(index.tt.especes, file="index-tt-especes.RData")
  setwd("../")
  message("\n\n***********************\nFonction prep.analyse() complétée: tableaux formattés et analyses prêtes à lancer.
\n***********************\nPour sélectionner un projet, utilisez selection.projet()\n
\n***********************\nPour lancer les analyses, utilisez les fonctions:\n
-- Invertébrés: INV.dens.gnrl(), INV.biodiv.gnrl() \n-- LIT: LIT.couvrt.gnrl(); Quad.couvrt.gnrl() \n-- Poissons: POIS.dens.gnrl()\n\n***********************\n")
} else { # refaire tableaux = FALSE, on charge les tableaux pre-formattes
    # import filtre famille
    if(filtre.famille) {
        message("filtre.famille=TRUE --> chargement automatique du filtre famille Filtre-taxo_Famille.csv")
        import.filtre.taxo("Filtre-taxo_Famille.csv","Famille")
    }
  setwd("Tableaux-pour-analyses")
  load(file="info-transect-TProj.RData", .GlobalEnv)
  load(file="dpoissons-TProj.RData", .GlobalEnv)
  load(file="dbio-TProj.RData", .GlobalEnv)
  load(file="index-invSp.RData", .GlobalEnv)
  load(file="index-Poissons.RData", .GlobalEnv)
  load(file="bioeco.RData", .GlobalEnv)
  load(file="bioeco-all.RData", .GlobalEnv)
  load(file="data-LIT-TProj.RData", .GlobalEnv)
  load(file="index-LIT.RData", .GlobalEnv)
  load(file="dQuad-TProj.RData", .GlobalEnv)
  load(file="coraux-fig.RData", .GlobalEnv)
  load(file="index-tt-especes.RData", .GlobalEnv)
  setwd("../")
  message("\n\n***********************\nFonction prep.analyse() complétée: tableaux pré-formattés chargés. Pour re-formatter les tableaux,
spécifiez refaire.tableaux = TRUE dans GS_mother-code_TProj.r.\n\n
\n***********************\nPour sélectionner un projet, utilisez selection.projet()\n
\n***********************\nPour lancer les analyses utilisez les fonctions:\n
-- Invertébrés: INV.dens.gnrl(), INV.biodiv.gnrl() \n-- LIT: LIT.couvrt.gnrl(); Quad.couvrt.gnrl() \n-- Poissons: POIS.dens.gnrl()\n\n***********************\n")
}

  # Tag TRUE when data read with no bugs
  data.read <<- TRUE
  wd.now <- dossier.R
  setwd(wd.now)


}

##################### Fin de la fonction prep.analyse() #####################
#############################################################################
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
stand.err <- function(x) sd(x,na.rm=TRUE)/sqrt(length(x))

# Nombre de valeurs uniques
count <- function(x) length(unique(x))

# Fonction calculant les valeurs totales de densité, la moyenne et SD par unité
val.arrondi <- 4
tot.mean.sd <- function(x,...) round(c(Tot=sum(x,...), Moy=mean(x,...), ET=sd(x,...)),val.arrondi)
mean.sd <- function(x,...) round(c(Moy=mean(x,...), ET=sd(x,...)),val.arrondi)
aggr.multi <- function(...) {

  t1 <- do.call(aggregate, args=...)
  colcl <- sapply(t1, class)
  colnom <- names(t1)[colcl=="matrix"]

  t1.sub <- lapply(colnom, function(x) {
    d1 <- as.data.frame(t1[,x]);
    if(length(colnom)>1) names(d1) <- paste(x, names(d1), sep=".");
    d1 })

  t2.sub <- do.call("data.frame", t1.sub)
  t1.fn <- data.frame(t1[,colcl!="matrix"], t1.sub)

}

# Trouve des typos possibles dans les noms d'espèces
typo.finder <- function(vect, typo.sensi=2, tag, ote.prem=TRUE) {
# vect est le vecteur de charactères où l'on veut scanner pour
 # des similarités possibles
  if(class(vect)!="character") stop("'vect' devrait être en charactères")
  if(!missing(tag)) {message("\n");message(sprintf("Rapport typo: %s", tag))}
  vals <- unique(vect)
  vals1 <- substr(vals, 1, 1) # première lettre seulement

  # comparaison des lettres en commun
  tmat <- adist(vals, vals)

  # ... on conserve les valeurs pour le triangle du bas seulement
  tmat[upper.tri(tmat, diag=TRUE)] <- 999

  if(ote.prem) {
  # si la première lettre est différente, on assume que ce n'est pas une typo
  tmat.1 <- adist(vals1, vals1)
  tmat[tmat.1 == 1] <- 999
}

  ind.typo <- which(tmat <= typo.sensi)
  if(length(ind.typo)>0) {
   mat.typo <- arrayInd(ind.typo, dim(tmat))
   d1 <- data.frame(Val=vals[mat.typo[,1]], Match=vals[mat.typo[,2]])
   print(d1)
 } else {
   message(sprintf("Pas de typos potentielles trouvées avec sensibilité de %s",
                 typo.sensi))
 }
}

# Extraire nom du projet du champ 'Id'
id2proj <- function(x) toupper(gsub("([A-Za-z]*_[A-Za-z]*)_.*", "\\1",x))

########################################################
# Fonctions dplyr:
s_group_by <- function(.data, ...) {
  eval.string.dplyr(.data, "group_by", ...) }
eval.string.dplyr <- function(.data, .fun.name, ...) {
  args <- list(...)
  args <- unlist(args)
  code <- paste0(.fun.name,"(.data,", paste0(args, collapse=","), ")")
  df <- eval(parse(text=code,srcfile=NULL))
  df
}

#######################################################
  ###### Filtres campagnes annuelles/semestrielles ######
  #######################################################
  # filtre.Camp: tableau des campagnes à utiliser pour analyses temporelles
  # modif après discussion avec Antoine 12 Octobre 2012:
  # re-créer tableau filtre à partir des années désirées pour l'analyse

  creerFiltre <- function(qAnnees) {

    if(length(qAnnees)==1) stop("Attention: spécifiez 2 années ou plus")
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

  ###################################################
  ######## Fonctions génériques #####################
  # Définition de fonctions qui seront utilisées couramment dans le code

  # 1. Applique le filtre spécifié au tableau donné en argument
  filtreTable <- function(wtable, wfiltre) {
    if((length(filtre.annees)>1) & (wfiltre %in% c("T_A_inv","T_S_inv","T_AeS_inv"))) { # appliquer le filtre si spécifié
        message(paste("Stations filtrees par",wfiltre))
        filtre.Camp <<- creerFiltre(filtre.annees) # recreer le tableau filtre
        wtable$key <- paste(wtable$St, wtable$Campagne, sep="_")
        dd.filt <- filter(wtable, key %in% filtre.Camp[[wfiltre]])
      dd.filt <- dd.filt[,names(dd.filt) != "key"] #ôter colonne key
  } else {
    if(wfiltre %in% c("A","S")) {
    CmpTag <- paste(wfiltre,filtre.annees,sep="_",collapse="|")
  } else if (wfiltre %in% c("T_A_inv","T_S_inv")) {
    wfiltre <- gsub("._(.)_.*","\\1",wfiltre)
    CmpTag <- paste(wfiltre,filtre.annees,sep="_",collapse="|")
  } else {
    CmpTag <- paste(filtre.annees,collapse="|") }
    wCampKeep <- grep(CmpTag, unique(wtable$Campagne), value=TRUE)

    dd.filt <- data.frame(filter(wtable, Campagne %in% wCampKeep))

    }

    if(nrow(dd.filt)==0) warning(
           sprintf("Attention!! \n\nAucune des données ne sont sélectionnées par le filtre sur stations:
\nFiltre %s, années %s\n", wfiltre, paste(filtre.annees,collapse=", ")))
    return(dd.filt)  }

  # 2. Converti les noms de campagne en année (charactère -> numérique)
  as.year <- function(x) as.numeric(sub("[AS]_","",x))

  # 3. Dessine les barres d'erreurs
  draw.SE <- function(mat, couleur=1, typel=1) {
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
  nst.par.gs <- function(AS="A",impact=FALSE, allyr=FALSE) {
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
  filtreTaxo <- function(wtable,action="inclure",
                          taxtype="Groupe",taxnom="Tous", silent=FALSE) {


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
          if(!silent) {
      message(paste(c("Filtre sur espèces",capitalize(action),qtype,taxtype,
                    paste(sort(taxnom),collapse=", ")),collapse=" :: "))
  }
      if(tolower(action) == "inclure") { # inclusion des groupes X
        filt.string <- sprintf("%s %%in%% c('%s')", taxtype, paste(taxnom, collapse="', '"))
        wtable <- wtable %>% filter_(filt.string)
      } else { # exclusion des groupes X
      wtable <- wtable[!(wtable[,taxtype] %in% taxnom),] }}}

    return(as.data.frame(wtable)) # switch out of dplyr df format
  }

  # Fonction interactive utilisée pour définir les variables du filtre sur les espèces
  # "inclure" ou "exclure" / unité taxonomique / nom
  def.filtre.especes <- function(taxoF.incl="tous", taxoF.utaxo, taxoF.nom) {

    if(taxoF.incl == "tous") {
       taxoF.incl <<- "inclure"
       taxoF.utaxo <<- "Groupe"
       taxoF.nom <<- "Tous"
     } else {
       if(missing(taxoF.incl)) {
       message("Definition des filtres taxonomiques (pas besoin de mettre des guillemets):")
         taxoF.incl <- tolower(readline("Inclure ou exclure? "))
       }else{ taxoF.incl <- tolower(taxoF.incl)}

       if(missing(taxoF.utaxo)) {
         mm <- "Unité taxomique? (Groupe/Sous-Groupe/Famille/Genre/Espece) "
         taxoF.utaxo <- capitalize(tolower(readline(mm)))
       }else{ taxoF.utaxo <- capitalize(tolower(taxoF.utaxo)) }

       taxoF.utaxo <- taxoF.utaxo
       if(taxoF.utaxo == "Sous-Groupe") taxoF.utaxo <- "S_Groupe"
       if(taxoF.utaxo == "Espece") taxoF.utaxo <- "G_Sp"

       if(missing(taxoF.nom)) {
         taxoF.nom <- readline("Nom? ")
         taxoF.nom <- capitalize(tolower(trim(unlist(strsplit(taxoF.nom,",")))))
       }else{ taxoF.nom <- capitalize(tolower(taxoF.nom)) }

       taxoF.incl <<- taxoF.incl
       taxoF.utaxo <<- taxoF.utaxo
       taxoF.nom <<- taxoF.nom
  }
}
  # Fonction qui permet de sélectionner un fichier .csv pour importer
  # sous R une liste de noms d'espèces/genre/famille
  # à utiliser dans le filtre taxonomique
  # argument "action" défini taxoF.incl
  # argument "niveau" défini le niveau taxonique, taxoF.utaxo
  # argument "titre" défini si la première rangée est le nom de la colonne dans le
  # ... fichier .csv (si non, spécifier titre = FALSE)
  import.filtre.taxo <- function(fichier, niveau="Famille",action="inclure",titre=FALSE,sepval=";") {

    if(missing(fichier)) fichier <- file.choose()
      lnoms <- read.csv(fichier,header=titre,sep=sepval) # sélectionner le fichier dans l'ordi
      if(class(lnoms)=="data.frame") lnoms <- lnoms[,1]
      lnoms <- capitalize(tolower(trim(lnoms))) # nettoyer format des noms
      taxoF.incl <<- action
      taxoF.utaxo <<- niveau
      taxoF.nom <<- lnoms

      voir.filtre.taxo()
  }

  # Fonction qui montre la valeur présente des filtres taxonomiques
  # Tapez "voir.filtre.taxo()" dans la console
  voir.filtre.taxo <- function() {
      fltre.now <- list(taxoF.incl, taxoF.utaxo, sort(taxoF.nom))
      names(fltre.now) <- c("taxoF.incl","taxoF.utaxo","taxoF.nom")
      print(fltre.now) }

  # Cette fonction retourne une version abbrégée des noms des groupes taxonomiques
  # inclus (ou exclus) au besoin
  taxotagFunk <- function() {

    if(taxoF.incl=="inclure" & taxoF.utaxo == "Groupe" & any(taxoF.nom %in% "Tous")) {
      return("") # si tous les groupes sont inclus, ne pas modifier le nom de fichiers
    } else {
        if(length(taxoF.nom) <= 5) {
            taxotag <- paste(paste(c(capitalize(taxoF.utaxo),abbreviate(c(capitalize(taxoF.incl),
                                                                          sort(taxoF.nom)))), collapse="-"),"_",sep="")
        } else {taxotag <- paste(paste(c(capitalize(taxoF.utaxo),
                                         abbreviate(c(capitalize(taxoF.incl),
                                                      sort(taxoF.nom)[1:5], "etc",
                                                      paste("N=",length(taxoF.nom),sep="")))), collapse="-"),"_",sep="")}
        return(taxotag) }
  }

  # Ces fonctions sont utilisées pour indiquer le départ et la fin
  # des codes compris dans une fonction
  departFunk <- function() {
      sc <- sys.calls()
      fname <- as.character(sc[[length(sc)-1]])[1] # extrait le nom de la fonction parent
      packageStartupMessage(sprintf("\nDepart %s()...",fname))
      try(print(unlist(mget(names(formals(fname)), envir=sys.frame(-1)))),
          silent=TRUE)
      message("#################\n")
      }
  finFunk <- function() {
      message("\n#################")
      sc <- sys.calls()
      fname <- as.character(sc[[length(sc)-2]])[1] # extrait le nom de la fonction parent
      packageStartupMessage(sprintf("---> fin %s().",fname)) }

  ## Si erreur dans une fonction, indiquer la fonction où l'erreur se produit
  ## EM() est donnée en argument a la fonction on.exit() qui tourne automatiquement
  ## l'argument spécifié quand une fonction se termine, naturellement ou avec
  ## une erreur. Ici EM() identifie si la fonction a eu une erreur, et si c'est
  ## le cas imprime le nom de la fonction pour faciliter l'identification du bug.
  EM <- function() {
      tb <- try(get(".Traceback",envir=baseenv()),silent=TRUE) # extraire messages d'erreur
      if(identical(paste(lastval(tb)), paste(sys.calls()[1]))) {
          message(sprintf("Erreur dans la fonction %s()",
                      paste(sys.calls()[[1]][1])))
          assign(".Traceback"[[1]],999,envir=baseenv())
          } else { finFunk()} }
