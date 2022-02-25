#Version R 4.0.5
#copyright Urbain Aurelie, Caillaud Sylvain & INRAE
#Last release may 14 2021

#Compilation : source('Analyse_FTIR.R') or ctrl+alt+r

#Data : File.brut (choose a file in a list), 

#Export : this code will export a csv file of your data (after a normalisation and a cleaning) at the direction /Documents/FTIR/data_output

# Install all packages and the latest version of R 
library(ggplot2)
library(dplyr)
library(reshape)
library(cowplot)
library(stringr)
library(fs)
source("C:/Users/scaillaud/Documents/FTIR/Codes R/Console/Fonctions_FTIR.R")

#################
## LES DONNEES ##
#################
# Choix du fichier .brut a traiter
FILE = file.choose()
print("File")

# Verification de l'extention du fichier, n'accepte que des .brut
if ( !(as.integer( regexpr(".brut", FILE, fixed= TRUE)) !=(-1) |
  as.integer( regexpr(".csv", FILE, fixed= TRUE)) !=(-1))) {
  stop("File entered is not a File.brut or File.csv, software can't run")
}

# Premier traitement du fichier
if(as.integer( regexpr(".brut", FILE, fixed= TRUE)) !=(-1)) {
  ###############################
  ## LES PARAMETRES EN ENTREE  ##
  ###############################
  
  # Mise en place des parametres generaux (possibilité de les changer ou de les faire rentrer par l'utilisateur avec readline())
  # defautl values : 846, 674.99, 3999.81, 800, 1820.
  # Number of points in x-axis
  NBpoints =863
  # Value at the first point
  NOdeb = 674.99
  # Value at the last point
  NOfin = 3999.81
  # First point of the window 828
  debWIN = 828
  # Last point of the window 1801
  finWIN = 1801
  
  ###########################
  ## CALCULS PRELIMINAIRES ##
  ###########################
  
  #### Stockage des donnees dans un data frame
  
  # On place la date et la populationen type character et le reste en numéric
  SpectresBruts <-  read.table(FILE,colClasses="character")
  SpectresBruts[,3:ncol(SpectresBruts)] <- apply(SpectresBruts[,3:ncol(SpectresBruts)],2,as.numeric)
  
  # Calcul du pas entre chaque point du spectre
  pas = (NOfin - NOdeb)/(NBpoints -1)
  
  # Creation d'un vecteur comprenant toutes les valeurs de NO possibles
  NO = seq(NOdeb,NOfin,pas)
  
  # Valeurs de nombre d'onde les plus proches de la fenêtre
  NOdeb = NO[which.min(abs(NO-debWIN))]
  NOfin = NO[which.min(abs(NO-finWIN))]
  
  # Index dans le fichier source, prendre en compte les trois premieres colonnes (date, nom et la repetition).
  IndexDeb = which.min(abs(NO-debWIN))+3
  IndexFin = which.min(abs(NO-finWIN))+3
  
  # On recentre les NO autour des valeurs demandées
  NO = NO[which(NO==NOdeb):which(NO==NOfin)] #va de 829 à 1801, length = 253
  
  # Calcul du nombre de NO conserves dans le spectre reduit
  NbNO = length(NO)
  
  # Donnees tronquees en fonction de la fenetre d'etude demandees
  SpectresRed = SpectresBruts[,c(1:3,IndexDeb:IndexFin)]
  
  #########################
  ## PROGRAMME PRINCIPAL ##
  #########################
  
  #Calcul des dimensions de la matrice des spectres dans la fenetre d'etude
  Dim = dim(SpectresRed)
  
  # Pour chaque spectre
  for (i in c(1: Dim[1])) {
    # Calcul de la ligne de base et soustraction de cette base au spectre
    SB = LigneBase (SpectresRed[i,],NbNO,NO)
    
    # Normalisation des donnees
    Norm=Normalise(SB,pas)
    if(i==1)SpectreNorm =Norm
    if(i>1)SpectreNorm=rbind(SpectreNorm,Norm)}
  
  
  # Reorganistion des spectres normes (ajout des noms et classement de la matrice par alphabetique)
  SpectreNorm <-  cbind(SpectresRed[,1:3],SpectreNorm)
  # On passe NO en valeurs entière pour faciliter les lectures
  colnames(SpectreNorm) <- c("Date", "Population", "Echantillon",as.integer(NO))
  
  #On rajoute la colonne Sample
  SpectreNorm <- cbind(SpectreNorm, paste0(SpectreNorm$Population, "-", SpectreNorm$Echantillon))
  
  SpectreNorm <- SpectreNorm[c(1:3, ncol(SpectreNorm), 4:(ncol(SpectreNorm)-1))]
  
  colnames(SpectreNorm)[4] <- "Sample"
  
  SpectreNorm$Population <- as.factor(SpectreNorm$Population)
}

if (as.integer( regexpr(".csv", FILE, fixed= TRUE)) !=(-1)) {
  #Importation de la table
  SpectreNorm <- read.csv2(FILE, header = TRUE, check.names = FALSE)
  colnames(SpectreNorm)[5:ncol(SpectreNorm)] <- as.integer(colnames(SpectreNorm)[5:ncol(SpectreNorm)])
  SpectreNorm$Population <- as.factor(SpectreNorm$Population)
}

#################################################
## LES GRAPHES D'ECHANTILLONS AVEC LA MOYENNE  ##
##                Echantillonage               ##
#################################################

# l est le nombre de populations différentes dans les donnees
l = length(levels(SpectreNorm$Population))

#Création d'un nouveau dataframe qui va accueillir les moyennes et les variances par population
# Echantillon = 0 permet de reconnaitre les moyennes
# Echantillon = 0.5 permet de reconnaitre les variances
SN <-  data.frame(Date = rep(SpectreNorm[1,1], times = l*2), 
                  Population = factor(rep(levels(SpectreNorm$Population), each = 2)), 
                  Echantillon = rep(c(0,0.5), times = l))

SN <- cbind(SN, Sample = paste0(SN$Population, "-", SN$Echantillon)) #On crée un nouveau df

#On modifie les samples des moy et var pour les reconnaitre plus facilement ensuite
SN$Sample[SN$Echantillon == 0] <- as.character(SN$Population[SN$Echantillon == 0])
SN$Sample[SN$Echantillon == 0.5] <- paste(SN$Population[SN$Echantillon == 0.5],"var", sep = "-")

#Pour chaque nombre d'onde
for (i in 5:ncol(SpectreNorm)){
  new <- c()
  # Calcul de la moyenne et de la variance par population
  for (pop in levels(SN$Population)){
    new <- rbind(new, 
                 mean(filter(SpectreNorm, Population == pop)[,i]),
                 var(filter(SpectreNorm, Population == pop)[,i]))
  }
  # On place moyenne et variance dans SN
  SN[,ncol(SN)+1] = new
}
colnames(SN) <-  colnames(SpectreNorm)

# On ajoute les spectres des individus au moy et var dans SN
SN <- rbind(SN, SpectreNorm)

#passage en format long, plus simple pour la manipulation et les graphiques
SN <- reshape::melt(SN, id.vars = c("Date", "Population", "Echantillon", "Sample"))
colnames(SN)[5:6] <- c("NO", "Absorbance")


#Passage en facteur pour pouvoir gérer les répétitions du format long
SN$NO <- as.factor(SN$NO)
SN$Sample <- as.factor(SN$Sample)

#Graphiques qui sortent dans Plots
ncourbe = 3 # Nombre d'échantillon par graphe (3 est un bon compromis)

#On parcourt les populations
for (pop in levels(SN$Population)){
  # Nombre d'échantillon d'une population
  #nbechantillon = max(SN$Echantillon[SN$Population == pop])
  nbechantillon = 10
  # Pour ne pas avoir de pb de range (pour s'arréter)
  if ((nbechantillon + 1) %% 5 == 0) nbechantillon = nbechantillon - 1
  
  #Liste qui détermine quels échantillons sont ensembles
  decoupe = seq(1, nbechantillon, by = ncourbe) 
  
  for (i in decoupe){
    #Quand on arrive au bout des echantillons
    if (i == decoupe[length(decoupe)]){
      titre = paste("Spectre Infrarouge des echantillons", pop, i,"à", nbechantillon, sep = " ")
    }
    # Sinon
    else {
      titre = paste("Spectre Infrarouge des echantillons", pop, i,"à", i+ncourbe-1, sep = " ")
    }
    
    # Sortie graphique
    graphe = ggplot(filter(SN, Population == pop & Echantillon %in% c(0, i:(i+ncourbe-1)))) +
      aes(x = reorder(NO,desc(NO)), y = Absorbance, group = Sample, color = Sample) + 
      geom_line() +
      labs(title = titre,
           x = "Nombre d'ondes (1/cm)",
           y = "Absorbance normalisée")+
      scale_x_discrete(breaks = as.integer(NO[seq(1,length(NO),by =15)]))
    print(graphe)
    
    # Sauvergade sous Jpg
    # ggsave(paste(titre, "jpg", sep ="."))
  }
}

# Suppression des mauvaises données
SN <- Suppression(SN)

#################################################
##      ENREGISTREMENT DES DONNEES ELAGUEES    ##
#################################################

#### Creation d'un fichier csv des donnees elaguees et normalisée
SpectreOut <- reshape(SN, idvar = c("Date", "Population", "Echantillon","Sample"), timevar = "NO", direction = "wide")
colnames(SpectreOut)[5:dim(SpectreOut)[2]] <- levels(SN$NO)
SpectreOut <- filter(SpectreOut, Echantillon != 0 & Echantillon !=0.5)

# Exportation dans /DocumentFTIR/data_output du fichier
SplitFile = path_split(FILE)[[1]]
NameFile = strsplit(SplitFile[length(SplitFile)], ".brut")[[1]]
name = sprintf("%s_%i_n&c.csv",NameFile,NbNO) #n&c = normed and cleared
dest = substr(FILE,1, nchar(FILE) - nchar(NameFile) - 5)

Sauvegarde()


#################################################
##         Choix du traitement de donnees      ##
#################################################

#cf "Fonctions_FTIR.R"
Choix()


