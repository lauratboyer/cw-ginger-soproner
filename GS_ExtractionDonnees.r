## Analyses des donnees KNS (Ginger/Soproner)
# Auteur: Laura Tremblay-Boyer, contact: l.boyer@fisheries.ubc.ca
# Time-stamp: <2013-08-01 15:38:43 Laura>

# Sujet: Formattage des tableaux de donnees brutes pre-analyse,
# creation de tableaux annexes + fonctions de base pour l'analyse

message("Put this elsewhere:")
# Definition lien fichiers pour sauvegarder graphiques/tableaux
fig.dir <- paste(dossier.R,"//Graphiques//",sep='')
tabl.dir <- paste(dossier.R,"//Tableaux//",sep='')

prep.analyse <- function() {

  lp <- require(reshape) # load package reshape
  if(!lp) install.packages("reshape") # installe reshape si requis
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

  if(!exists("last")) { # dernier element de l'objet
  last <<- function(x) x[length(x)] }

  # ote les espaces au debut et apres un charactere
  # e.g. "  Famille " -> "Famille"
  if(!exists("trim")) { #clumsy regexp
      trim <<- function(x) gsub("^\\s+","",gsub("\\s+$","",x)) }

# V�rifie que la graphic device a la bonne grandeur, sinon en ouvrir une nouvelle
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

  # Nettoyer accents noms g�om�trie
  geom.key <- data.frame("Geom"=unique(info.transect$Geom),
                         "Geomorpho"=c("Recif barriere externe","Recif barriere interne",
                           "Recif reticule","Recif frangeant","Passe","Herbiers"))
  print("make geom.key robust to unique order of info.transect$Geom")
  info.transect <- merge(info.transect, geom.key, by="Geom")

  ########### Donn�es LIT ############
  ####################################
  # data.LIT: donn�es brutes LIT
  # index.LIT: info complementaires LIT

  # ote rang�es sans stations
  data.LIT <- data.LIT[data.LIT$St != "",]
  # d�fini noms de colonnes pour index.LIT:
  names(index.LIT) <- c("Code_LIT","CODE_DET","S_Corail_Acro","S_Corail_Forme",
                      "S_Corail_All","S_Corail_Sensi","S_Abio_Corail_All")

 # creer categories de substrat pour tableaux syntheses
  coraux.fig <- list("Coraux_Gen"=c("Coraux","Coraux morts","Coraux mous",
                     "Algues","Abiotique","Autre faune"),
                   "Coraux_Acro"=c("Acroporidae","Non-acroporidae"),
                   "TS_All"=c("Coraux","Coraux morts","Coraux mous",
                     "Algues","Abiotique","Autre faune","Acroporidae",
                     "Non-acroporidae","Macro-algues","Assemblage d'algues",
                     "Corail branchu","Corail tabulaire","Corail massif",
                     "Corail encroutant","Corail foliaire","Corail submassif",
                     "Corail digite"))


  ###### Information biologie/�cologie poissons #######
  #####################################################
  # bioeco.all: info compl�mentaires sur les poissons
  bioeco.all <- data.bioeco[,c("Code_Sp","famille","genre","espece","taille_moy","commercial_2",
                        "cible_VKP","a","b","groupe_troph1","mobilite")]
  names(bioeco.all) <- c("Code_SP","Famille","Genre","Espece","Taille","Peche","Cible",
                       "Coeff_a","Coeff_b","Groupe_Trophique","Mobilite")

  # Clarifier �tiquettes
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

  # Rajouter info mobilit�:
  mob.label <- data.frame("Mobilite"=0:4,"moblabel"=c("Mob.0","Terr","Sed","Mob","TrMob"))
  bioeco.all <- merge(bioeco.all,mob.label,by="Mobilite",all.x=TRUE)

  # Ajuster si besoin le format pour les noms d'esp�ces:
  bioeco.all[,c("Famille","Genre","Espece")] <-
      sapply(bioeco.all[,c("Famille","Genre","Espece")], capitalize)

  # Tableau bioeco, info a + b seulement
  bioeco <- data.bioeco[,c("Code_Sp","a","b")]
  names(bioeco) <- c("Code_SP","a","b") # compatibilit� avec data.poissons

  # V�rifier s'il y a des esp�ces de poissons sans valeurs a ou b:
  print(paste(length(unique(bioeco[bioeco$a=="inconnu" | bioeco$b=="inconnnu","Code_SP"])),"especes sans valeur a ou b otees de l'analyse"))
  bioeco <- bioeco[bioeco$a!="inconnu" & bioeco$b!="inconnu",]
  bioeco$a <- as.numeric(bioeco$a); bioeco$b <- as.numeric(bioeco$b)

  ###### Donn�es de comptage Poissons ######
  ##########################################
  # data.poissons: donn�es de comptage brutes sur les poissons
  names(data.poissons) <- c("Campagne","Date","St","Obs","Vis","Courant","Code_SP",
                          "Famille","Genre","Espece","G_Sp","N","L","D1","D2",
                          "Secteur","Categorie","Note")

  # rajout ann�e de la campagne
  data.poissons$An <- getmatch(data.poissons$Campagne,"[0-9]+")

  ########################
  # Nettoyage du tableau:

  # Garder seulement les valeurs de N & L d�clar�es num�riques
  data.poissons$N <- as.numeric(data.poissons$N)
  data.poissons$L <- as.numeric(data.poissons$L)

  # Ajuster si besoin le format pour les noms d'esp�ces:
  data.poissons[,c("Famille","Genre","Espece","G_Sp")] <-
      sapply(data.poissons[,c("Famille","Genre","Espece","G_Sp")], capitalize)

  # Toutes les stations en majuscules
  data.poissons$St <- toupper(data.poissons$St)

  # Nettoyage du tableau: colonne D1/D2/N/L
  # Esp�ces sans Code_SP d�fini �t�es
  # Imprimer avertissements
  print(paste(nrow(data.poissons[data.poissons$Code_SP  %in% c("","#N/A"),]),"rangees sans Code_SP"))
  print(paste(nrow(data.poissons[is.na(data.poissons$D1),]),"valeurs D1 -> 0"))
  print(paste(nrow(data.poissons[is.na(data.poissons$D2),]),"valeurs D2 -> 0"))
  print(paste(nrow(data.poissons[is.na(data.poissons$N),]),"valeurs N -> 0"))
  print(paste(sum(is.na(data.poissons$L)),"rangees sans valeur L"))

  # �ter rang�es:
  data.poissons <- data.poissons[!(data.poissons$Code_SP %in% c("","#N/A")),]
  data.poissons$D1[is.na(data.poissons$D1)] <- 0
  data.poissons$D2[is.na(data.poissons$D2)] <- 0
  data.poissons$N[is.na(data.poissons$N)] <- 0
  data.poissons <- data.poissons[!(is.na(data.poissons$L)),] # Taille

  # Expansion du data frame pour inclure toutes les combinaisons Campagne/St/Code_SP:
  # ... donc pour les poissons les densit�s sont nulles sur toutes Campagnes/St o�
  # ... l'esp�ce est non-observ�e
  all.comb.poissons <- expand.grid("Campagne"=unique(data.poissons$Campagne),
                                   "St"=unique(data.poissons$St),
                                   "Code_SP"=unique(data.poissons$Code_SP),
                                   stringsAsFactors=FALSE)
  St.by.year <- unique(data.poissons[,c("Campagne","St")]) # year/station samples
  all.comb.poissons <- merge(all.comb.poissons, St.by.year, by=c("Campagne","St"))

  # ... reunion avec donn�es poissons pour rajouter les densit� nulles
  data.poissons2 <- merge(all.comb.poissons,
                          data.poissons,by=c("Campagne","St","Code_SP"),all.x=TRUE)
  data.poissons2$N[is.na(data.poissons2$N)] <- 0 # abondance � zero si non-observ�e
  data.poissons <- data.poissons2

  ##########################################
  #### Donn�es de comptage Invert�br�s #####
  ##########################################
  # data.inv: donn�es de comptage brutes pour invert�br�s

  dbio <- data.inv[,c("Campagne","St","T","Grp2","S_Grp2","F2","G2","G_Sp","N","D","Ltrans")]
  dbio$St <- toupper(data.inv$St)
  names(dbio) <- c("Campagne","St","T","Groupe","S_Groupe","Famille","Genre","G_Sp","N","D","Ltrans")

  ############################
  # nettoyage:
  # garde seulement les rang�es avec transects A, B, C
  nTnd <- nrow(dbio[!(dbio$T %in% c("A","B","C")),])
  dbio <- dbio[dbio$T %in% c("A","B","C"),]; print(paste("Approx",nTnd,"rangees otees vu transect non-def"))

  ############################
  # abondances nulles:
  # commencer par oter les rang�es N=0 pour eliminer les stations o� l'esp�ce n'est observ�e sur aucun T
  dbio <- dbio[dbio$N > 0,]

  # rajouter les transects N=0 *seulement* lorsque l'esp�ce est observ�e sur la ST mais pas tous les transects
  # cr�er identifiant unique pour chaque combinaison campagne/st/esp�ce observ�e
  dbio$uID <- paste(dbio$Campagne,dbio$St,dbio$G_Sp,sep="_")
  index.dbio <- unique(dbio[,c("uID","Campagne","St","Groupe","S_Groupe",
                               "Famille","Genre","G_Sp","D","Ltrans")])
  dbio.allT <- expand.grid("T"=c("A","B","C"),"uID"=unique(dbio$uID),
                           stringsAsFactors=FALSE)
  dbio.tmp <- merge(dbio[,c("uID","T","N")],dbio.allT,by=c("uID","T"),all.y=TRUE)
  dbio.tmp$N[is.na(dbio.tmp$N)] <- 0 # valeurs NA remplac�es par 0 pour remplir les transects manquants
  dbio <- merge(dbio.tmp, index.dbio)

  ################
  # Oter les accents des noms des invert�br�s pour �viter les erreurs
  # d�es � l'encodage
  # Si l'encodage UTF-8 est bien lu par l'ordi
  # (les accents apparaissent correctement)
  dbio$Groupe <- gsub("�","e",capitalize(tolower(trim(dbio$Groupe))))
  dbio$S_Groupe <- gsub("�","e",dbio$S_Groupe)
  dbio$S_Groupe <- gsub("�","e",dbio$S_Groupe)
  dbio$S_Groupe <- gsub("�","i",dbio$S_Groupe)
  dbio$S_Groupe <- capitalize(tolower(trim(dbio$S_Groupe)))

  # Si l'encodage UTF-8 n'est pas lu par l'ordi (les accents suivent un format similaire � <U+00E9>)
  # Conversion � encodage latin1 qui permet de substituer les charact�res
  if("UTF-8" %in% unique(Encoding(dbio$Groupe))) {
    dbio$Groupe <- tradfunk(dbio$Groupe)
    dbio$S_Groupe <- tradfunk(dbio$S_Groupe)
    dbio$Famille <- tradfunk(dbio$Famille)
  }

  # Premi�re lettre en majuscule, le reste en minuscules
  dbio$G_Sp <- capitalize(tolower(dbio$G_Sp))
  dbio$Famille <- capitalize(tolower(dbio$Famille))
  dbio$Genre <- capitalize(tolower(dbio$Genre))

  # Correction d'erreur dans la base de donn�es (temporaire)
  # tableau index esp�ce/sous-groupe/groupe
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

  # rassembler en 1 rang�e les observations de la m�me esp�ce sur le transect
  dbio <- aggregate(list("N"=dbio$N), as.list(dbio[,c("Campagne","St","T","Groupe","S_Groupe",
                           "Famille","Genre","G_Sp","Ltrans")]), sum)

  # Calculer la densit� en hectares
  # longueur du transect est de 50 m�tres, largeur est d�finie sous Ltrans
  dbio$D <- dbio$N/(50*dbio$Ltrans) * 10000

  ########################
  ### Tableaux annexes ###

  # cl�s noms des invert�br�s
  InvGr.Key <- c("Algues","Ascidies","Cnidaires","Crustac?s","Echinodermes",
                 "Eponges","Mollusques","Phan?rogame","Vers")

  # tableau index pour la taxonomie de toutes les esp�ces observ�e:
  index.invSp <- unique(dbio[,c("G_Sp","Genre","Famille","S_Groupe","Groupe")])

  ################
  # infos sur station/mission/ann�e
  info.transect.INV <- unique(data.inv[,c("Annee","Mois","Mission","Campagne","St")])
  info.transect.INV <- merge(info.transect.INV, info.transect,
                             by="St")[,c("Annee","Mois","Mission","Campagne","St",
                               "Geomorpho","N_Impact")]
  info.transect.INV.geo <- aggregate(info.transect.INV[,c("Mois", "Annee")],as.list(info.transect.INV[,c("Mission","Campagne","Geomorpho")]), last)
  print("si un transect/Campagne �chantillon� dans 2 mois diff�rents, garde 1 mois/ann�e seulement dans le tableau")


  ####################################
  ###### P�riodes BACIP Campagnes ####
  ####################################
  # pr.Bacip: cat�gorie avant/pendant/apr�s pour les ann�es
  names(pr.Bacip) <- c("Campagne","Annee","Periode_3","Periode_2")

  #######################################################
  ###### Filtres campagnes annuelles/semestrielles ######
  #######################################################
  # filtre.Camp: tableau des campagnes � utiliser pour analyses temporelles
  # modif apr�s discussion avec Antoine 12 Octobre 2012:
  # re-cr�er tableau filtre � partir des ann�es d�sir�es pour l'analyse

  creerFiltre <- function(qAnnees) {

    if(length(qAnnees)==1) stop("Attention: sp�cifier 2 ann�es ou plus")
    tb <- unique(dbio[,c("St","Campagne")])
    tb$echtl <- 1
    tb2 <- cast(tb, St ~ Campagne, value="echtl")

    # s�lectionner les colonnes avec les ann�es d�sir�es pour le filtre
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
  ######## Fonctions g�n�riques #####################
  # D�finition de fonctions qui seront utilis�es couramment dans le code

  # 1. Applique le filtre sp�cifi� au tableau donn� en argument
  filtreTable <<- function(wtable, wfiltre) {
    if(wfiltre %in% c("T_A_inv","T_S_inv","T_AeS_inv")) { # appliquer le filtre si sp�cifi�
    message(paste("Stations filtrees par",wfiltre))
    filtre.Camp <<- creerFiltre(filtre.annees) # recreer le tableau filtre

    wtable$key <- paste(wtable$St, wtable$Campagne, sep="_")
    dd.filt <- merge(data.frame("key"=filtre.Camp[[wfiltre]]),
                     wtable,by="key", drop.x="key")
    dd.filt <- dd.filt[,names(dd.filt) != "key"] #�ter colonne key
  } else {
    CmpTag <- paste(filtre.annees,collapse="|")
    wCampKeep <- grep(CmpTag, unique(wtable$Campagne), value=TRUE)
    dd.filt <- wtable[wtable$Campagne %in% wCampKeep,]
    }

    if(nrow(dd.filt)==0) warning(
           sprintf("Attention!! \n\nAucune des donn�es ne sont s�lectionn�es par le filtre sur stations:
\nFiltre %s, ann�es %s\n", wfiltre, paste(filtre.annees,collapse=", ")))
    return(dd.filt)  }

  # 2. Converti les noms de campagne en ann�e (charact�re -> num�rique)
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


  # 4. Calcul du nombre de stations �chantillon�es par groupement spatial, selon le filtre
  # Groupement spatial: g�omorphologie, ou g�omorphologie/impact
  # d�fini pour les invert�br�s mais pourrait �tre appliqu� aux poisssons
  # (Retourne aussi le nom des stations)
  nst.par.gs <<- function(AS="A",impact=FALSE, allyr=FALSE) {
    # lorsque allyr=TRUE, calcule le nombre de campagnes effectu�es sur chaque
    # g�omorphologie +/- impact
    # sinon, calcule le nombre de stations *par campagne* sur chaque g�omorpho +/- impact

    if(AS %in% c("A","S")) { dd <- filtreTable(dbio, wfiltre=paste("T",AS,"inv",sep="_"))
                             } else {
                               keySmstrl <- c(paste("A_",filtre.annees,sep=""),
                                              paste("S_",filtre.annees,sep=""))
                               dd <- dbio[dbio$Campagne %in% keySmstrl,] }
    ff <- c("Campagne","Geomorpho","N_Impact","St")
    # conserver seulement les facteurs sp�cifi�s en arguments
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
      # On filtre les rangees selon les especes/groupes taxo specifie
      message(paste(c("Filtre sur especes",capitalize(action),taxtype,
                    paste(sort(taxnom),collapse=", ")),collapse=" :: "))

      if(action == "inclure") { # inclusion des groupes X
      wtable <- wtable[wtable[,taxtype] %in% taxnom,]
      } else { # exclusion des groupes X
      wtable <- wtable[!(wtable[,taxtype] %in% taxnom),] }}

    return(wtable)
  }

  # Fonction interactive utilis�e pour d�finir les variables du filtre sur les esp�ces
  # "inclure" ou "exclure" / unit� taxonomique / nom
  def.filtre.especes <<- function(aF="tous") {

    if(aF == "tous") {
       taxoF.incl <<- "inclure"
       taxoF.utaxo <<- "Groupe"
       taxoF.nom <<- "Tous"
     } else {
       message("Definition des filtres taxonomiques:")
    taxoF.incl <<- tolower(readline("Inclure ou exclure? "))
    mm <- "Unit� taxomique? (Groupe/Sous-Groupe/Famille/Genre/Espece) "
    taxoF.utaxo <<- capitalize(tolower(readline(mm)))
       if(taxoF.utaxo == "Sous-Groupe") taxoF.utaxo <<- "S_Groupe"
       if(taxoF.utaxo == "Espece") taxoF.utaxo <<- "G_Sp"
    taxoF.nom <<- readline("Nom? ")
    taxoF.nom <<- capitalize(tolower(trim(unlist(strsplit(taxoF.nom,",")))))
  }

  # Fonction qui permet de s�lectionner un fichier .csv pour importer
  # sous R une liste de noms d'esp�ces/genre/famille
  # � utiliser dans le filtre taxonomique
  # argument "action" d�fini taxoF.incl
  # argument "niveau" d�fini le niveau taxonique, taxoF.utaxo
  # argument "titre" d�fini si la premi�re rang�e est le nom de la colonne dans le
  # ... fichier .csv (si non, sp�cifier titre = FALSE)
  import.filtre.taxo() <<- function(action="inclure",niveau="Famille",titre=TRUE) {

      lnoms <- read.csv(file.choose(),header=titre) # s�lectionner le fichier dans l'ordi
      if(class(lnoms)=="data.frame") lnoms <- lnoms[,1]
      lnoms <- capitalize(tolower(trim(lnoms))) # nettoyer format des noms
      taxoF.incl <<- action
      taxoF.utaxo <<- niveau
      taxoF.nom <<- lnoms

      voir.filtre.taxo()
  }


  # Fonction qui montre la valeur pr�sente des filtres taxonomiques
  # Tapez "voir.filtre.taxo()" dans la console
  voir.filtre.taxo <<- function() {
      fltre.now <- list(taxoF.incl, taxoF.utaxo, sort(taxoF.nom))
      names(fltre.now) <- c("taxoF.incl","taxoF.utaxo","taxoF.nom")
      print(fltre.now) }

  # Cette fonction retourne une version abbr�g�e des noms des groupes taxonomiques
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

  # Ces fonctions sont utilis�es pour indiquer le d�part et la fin
  # des codes compris dans une fonction
  departFunk <<- function() {
      sc <- sys.calls()
      fname <- as.character(sc[[length(sc)-1]])[1] # extrait le nom de la fonction parent
      packageStartupMessage(sprintf("\nDepart %s()...",fname))
      print(unlist(mget(names(formals(fname)), envir=sys.frame(-1))))
      message("#################\n")
      }
  finFunk <<- function() {
      message("\n#################")
      sc <- sys.calls()
      fname <- as.character(sc[[length(sc)-2]])[1] # extrait le nom de la fonction parent
      packageStartupMessage(sprintf("---> fin %s().",fname)) }

  ## Si erreur dans une fonction, indiquer la fonction o� l'erreur se produit
  ## EM() est donn�e en argument a la fonction on.exit() qui tourne automatiquement
  ## l'argument sp�cifi� quand une fonction se termine, naturellement ou avec
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
  # D�finir objets � mettre dans l'environnement global pour utilisation subs�quente:
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

  # Tag TRUE when data read with no bugs
  data.read <<- TRUE
  wd.now <- dossier.R
  setwd(wd.now)

  message("\n\nFonction prep.analyse() compl�t�e. Tableaux formatt�s et analyses pr�tes � lancer.
Pour lancer les analyses manuellement utiliser les fonctions:\n
Invert�br�s: Run.INV.biodiv(), Run.INV.densite() \nLIT: Run.LIT.all() \nPoissons: Run.poissons.all()")
}
