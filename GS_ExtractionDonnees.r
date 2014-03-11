## Analyses des données KNS (Ginger/Soproner)
# Auteur: Laura Tremblay-Boyer, contact: l.boyer@fisheries.ubc.ca
# Time-stamp: <2014-03-11 16:00:26 Laura>

# Sujet: Formattage des tableaux de données brutes pré-analyse,
# création de tableaux annexes + fonctions de base pour l'analyse

if(!exists("dossier.R")) {print("Faire \nsource(GS_KNS_MotherCode.r)\n
pour définir l'emplacement des codes et des dossiers")
                      }else{
# Définition lien fichiers pour sauvegarder graphiques/tableaux
fig.dir <- paste(dossier.R,"//Graphiques//",sep='')
tabl.dir <- paste(dossier.R,"//Tableaux//",sep='')}

prep.analyse <- function() {

  lp <- try(library(reshape2)) # load package reshape
  if(class(lp)=="try-error") {
      install.packages("reshape") # installe reshape si requis
      library(reshape)}
  #######################################################################
  ###################### Formattage des tableaux ########################
  #######################################################################
  ### Fonctions utiles pour formattage
  if(!exists("capitalize")) { # rajoute une lettre majuscule au debut
  capitalize <<- function(x) {
    gsub('(\\w)(\\w*)','\\U\\1\\L\\2',x,perl=TRUE) }}

  if(!exists("getmatch")) { # extrait la partie correspondante au match
      getmatch <<- function(x,str2match,...) {
          unlist(regmatches(x,regexpr(str2match,x,...))) }}

  # dernier element de l'objet
  last <<- function(x) x[length(x)]

  # ote les espaces au debut et apres un charactere
  # e.g. "  Famille " -> "Famille"
  if(!exists("trim")) { #clumsy regexp
      trim <<- function(x) gsub("^\\s+","",gsub("\\s+$","",x)) }

  # Vérifie que la graphic device a la bonne grandeur, sinon en ouvrir une nouvelle
  check.dev.size <<- function(ww,hh) {

    if(dev.cur()==1){ dev.new(width=ww,height=hh)
                  } else {
    ds <- dev.size()
    if(round(ds[1],2)!=round(ww,2)
       | round(ds[2],2)!=round(hh,2)) {
        dev.off(); dev.new(width=ww,height=hh)} }
}
  # standard error
  stand.err <<- function(x) sd(x)/sqrt(length(x))

  ####################################
  ###### Information transects #######
  ####################################
  # info.transect: informations sur les transects
  info.transect <- data.info.transect # extraire tableau brut
  info.transect$cpus <- info.transect$Prod_ha # nouvelle colonne pour le cpus
  info.transect$Geom <- trim(info.transect$Geom) # oter espaces surperflus
  info.transect$Geomorpho <- info.transect$Geom

  # Nettoyer accents noms géométrie -- à reviser vu encodage changé durant import?
  geom.key <- data.frame("Geom"=unique(info.transect$Geom))
  geom.key$Geomorpho <- geom.key$Geom # toujours nécéssaire vu accents réglés?

  ########### Données LIT ############
  ####################################
  # data.LIT: données brutes LIT
  # index.LIT: info complementaires LIT

  # ote rangées sans stations
  data.LIT <- data.LIT[data.LIT$St != "",]
  # défini noms de colonnes pour index.LIT:
  names(index.LIT) <- c("Code_LIT","CODE_DET","S_Corail_Acro","S_Corail_Forme",
                      "S_Corail_All","S_Corail_Sensi","S_Abio_Corail_All")

 # créer catégories de substrat pour tableaux synthèses
  coraux.fig <- list("Coraux_Gen"=c("Coraux","Coraux morts","Coraux mous",
                     "Algues","Abiotique","Autre faune"),
                   "Coraux_Acro"=c("Acroporidae","Non-acroporidae"),
                   "TS_All"=c("Coraux","Coraux morts","Coraux mous",
                     "Algues","Abiotique","Autre faune","Acroporidae",
                     "Non-acroporidae","Macro-algues","Assemblage d'algues",
                     "Corail branchu","Corail tabulaire","Corail massif",
                     "Corail encroutant","Corail foliaire","Corail submassif",
                     "Corail digite"))


  ###### Information biologie/écologie poissons #######
  #####################################################
  # bioeco.all: info complémentaires sur les poissons
  bioeco.all <- data.bioeco[,c("Code_Sp","famille","genre","espece","nom_total","taille_moy","commercial_2",
                        "cible_VKP","a","b","groupe_troph1","mobilite")]
  names(bioeco.all) <- c("Code_SP","Famille","Genre","Espece","G_Sp","Taille","Peche","Cible",
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

  # Ajuster si besoin le format pour les noms d'espèces:
  tl.cap <- function(x) capitalize(tolower(x))
  bioeco.all[,c("Famille","Genre","Espece","G_Sp")] <-
      sapply(bioeco.all[,c("Famille","Genre","Espece","G_Sp")], tl.cap)

  # Tableau bioeco, info a + b seulement
  bioeco <- data.bioeco[,c("Code_Sp","a","b")]
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
  # rajout année de la campagne
  data.poissons$An <- getmatch(data.poissons$Campagne,"[0-9]+")

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
                                   "Code_SP"=unique(data.poissons$Code_SP),
                                   stringsAsFactors=FALSE)
  St.by.year <- unique(data.poissons[,c("Campagne","St")]) # year/station samples
  all.comb.poissons <- merge(all.comb.poissons, St.by.year, by=c("Campagne","St"))

  # ... reunion avec données poissons pour rajouter les densité nulles
  data.poissons2 <- merge(all.comb.poissons,
                          data.poissons,by=c("Campagne","St","Code_SP"),all.x=TRUE)
  data.poissons2$N[is.na(data.poissons2$N)] <- 0 # abondance à zero si non-observée
  data.poissons <- data.poissons2

  ##########################################
  #### Données de comptage Invertébrés #####
  ##########################################
  # data.inv: données de comptage brutes pour invertébrés

  dbio <- data.inv[,c("Campagne","St","T","Grp2","S_Grp2","F2","G2","G_Sp","N","D","Ltrans")]
  dbio$St <- toupper(data.inv$St)
  names(dbio) <- c("Campagne","St","T","Groupe","S_Groupe","Famille","Genre","G_Sp","N","D","Ltrans")

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
  index.dbio <- unique(dbio[,c("uID","Campagne","St","Groupe","S_Groupe",
                               "Famille","Genre","G_Sp","D","Ltrans")])
  dbio.allT <- expand.grid("T"=c("A","B","C"),"uID"=unique(dbio$uID),
                           stringsAsFactors=FALSE)
  dbio.tmp <- merge(dbio[,c("uID","T","N")],dbio.allT,by=c("uID","T"),all.y=TRUE)
  dbio.tmp$N[is.na(dbio.tmp$N)] <- 0 # valeurs NA remplacées par 0 pour remplir les transects manquants
  dbio <- merge(dbio.tmp, index.dbio)

  ################

  # Première lettre en majuscule, le reste en minuscules
  dbio$G_Sp <- capitalize(tolower(dbio$G_Sp))
  dbio$Famille <- capitalize(tolower(dbio$Famille))
  dbio$Genre <- capitalize(tolower(dbio$Genre))

  # Correction d'erreur dans la base de données (temporaire)
  # tableau index espèce/sous-groupe/groupe
  dbio[dbio$G_Sp %in% c("Paguritta sp.", "Ciliopagurus strigatus"),"S_Groupe"] <- "Decapodes"
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
  dbio <- aggregate(list("N"=dbio$N), as.list(dbio[,c("Campagne","St","T","Groupe","S_Groupe",
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
  index.invSp <- unique(dbio[,c("G_Sp","Genre","Famille","S_Groupe","Groupe")])

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

  ####################################
  #### Tableau espèces universels ####
  ####################################
  index.Poissons <- unique(data.poissons[,c("G_Sp", "Genre", "Famille")])
  index.Poissons$Groupe <- index.Poissons$S_Groupe <- NA
  index.Poissons$type <- "Poissons"
  index.invSp$type <- "Inv"
  index.tt.especes <- rbind(index.invSp, index.Poissons)

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
    tb2 <- dcast(tb, St ~ Campagne, value="echtl")

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
      message("Pas de filtre sur especes applique")
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
      message(paste(c("Filtre sur especes",capitalize(action),qtype,taxtype,
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
  import.filtre.taxo <<- function(niveau="Famille",action="inclure",titre=TRUE) {

      lnoms <- read.csv(file.choose(),header=titre) # sélectionner le fichier dans l'ordi
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
      tb <- try(get(".Traceback",envir=baseenv())) # extraire messages d'erreur
      if(identical(paste(last(tb)), paste(sys.calls()[1]))) {
          message(sprintf("Erreur dans la fonction %s()",
                      paste(sys.calls()[[1]][1])))
          assign(".Traceback"[[1]],999,envir=baseenv())
          } else { finFunk()} }

  ###################################################
  ###################################################
  # Définir objets à mettre dans l'environnement global pour utilisation subséquente:
  info.transect <<- info.transect
  info.transect.INV <<- info.transect.INV
  info.transect.INV.geo <<- info.transect.INV.geo
  dpoissons <<- data.poissons
  dbio <<- dbio
  index.invSp <<- index.invSp
  bioeco <<- bioeco
  bioeco.all <<- bioeco.all
  data.LIT <<- data.LIT
  index.LIT <<- index.LIT
  coraux.fig <<- coraux.fig
  pr.Bacip <<- pr.Bacip
  index.tt.especes <<- index.tt.especes

  # Tag TRUE when data read with no bugs
  data.read <<- TRUE
  wd.now <- dossier.R
  setwd(wd.now)

  message("\n\nFonction prep.analyse() complétée. Tableaux formattés et analyses prêtes à lancer.
Pour lancer les analyses manuellement utiliser les fonctions:\n
Invertébrés: Run.INV.biodiv(), Run.INV.densite() \nLIT: Run.LIT.all() \nPoissons: Run.poissons.all()")
}
