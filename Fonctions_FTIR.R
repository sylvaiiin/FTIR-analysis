###################
## LES FONCTIONS ##
###################

################################################################################
# Fonction qui calcule la ligne de base en effectuant une regression lineaire sur 10 points
# (extremites du spectre) et soustrait cette ligne de base aux donnees

LigneBase = function (Spectre,NbNOred,NOred) {
  
  # Creation de vecteurs comportant les valeurs necessaires au calcul de la regression
  #soit les 5 premiers et les 5 derniers points
  PtsReg = 5;
  Y = cbind(Spectre[,4:(PtsReg+3)],Spectre[1,(NbNOred-PtsReg+1):NbNOred+3]);
  X = c(NOred[1:PtsReg],NOred[(NbNOred-PtsReg+1):NbNOred]);
  
  # Transformation des donnees pour obtenir les memes dimensions
  # application pour fonction lm
  Y =t(Y);
  X=as.matrix(X);
  Reg = lm(Y~X)$coefficients;
  
  # Calcul de la base a partir des resultats de la regression
  Base = Reg[1]+(Reg[2]*NOred);
  
  # Retrait de la base au spectre
  SB = Spectre[1,4:(NbNOred+3)] - Base; 
  
  return (SB)}

#################################################################################
# Permet une normalisation des donnees apres le retrait de la ligne de base

Normalise = function(SB,pas) {
  
  # Elimination des valeurs negatives
  neg = which(SB<0);
  for (i in neg) SB[i]=0;
  
  # Normalisation des donnees  
  Norm=SB/(sum(SB)*pas);
  Norm = signif(Norm,digits=6)
  
  return (Norm) }

#################################################################################
# Demande à l'utilisateur quelles donnéesil veut supprimer
Suppression <- function(SN){
  
  Supr <- readline(prompt='Quel échantillon voulez-vous éliminer ? ("q" pour quitter): ')
  Supr <- tolower(Supr)
  if (Supr == "q"){
    
    # Il faut recalculer les moyennes et les variances
    for (no in levels(SN$NO)){
      for (pop in levels(SN$Population)){
        SN$Absorbance[SN$NO == no & SN$Population == pop & SN$Echantillon == 0] <-  mean(filter(SN, NO == no & SN$Population == pop & ! (Echantillon %in% c(0, 0.5)))$Absorbance)
        SN$Absorbance[SN$NO == no & SN$Population == pop & SN$Echantillon == 0.5] <-  var(filter(SN, NO == no & SN$Population == pop & ! (Echantillon %in% c(0, 0.5)))$Absorbance)
      }
    } 
  }
  
  else if (Supr %in% tolower(levels(SN$Sample)) & 
           !(Supr %in% tolower(levels(SN$Population)))){
    SN <- filter(SN, tolower(SN$Sample) != Supr)#Enlever de SN
    SN <- Suppression(SN)
  }
  
  else{
    cat("Il faut rentrer un nom d'échantillon ou aucun !!!")
    SN <- Suppression(SN)
  }
  return(SN)
}

Sauvegarde <- function(){
  demande <- readline(paste0("Voulez-vous enregistrer les données élaguées en .csv/n",
                             "! Attention à ne pas overrun d'anciens fichiers ! (Oui ou Non) : "))
  
  if (tolower(demande) == "oui"){
    write.csv2(SpectreOut, file = paste0(dest,name),row.names = FALSE ) 
  }
  
  else if (tolower(demande) != "non"){
    Sauvergarde()
  }
}
#################################################################################
# Demande à l'utilisateur quels sorties graphiques il veut voir
Choix <- function(){
  choix <- readline(paste0("####################################\n",
                           "#     1. Graph with all samples    #\n",
                           "#       2.Means of population      #\n",
                           "#     3.Means of  one population   #\n",
                           "#         4.Mean comparison        #\n",
                           "####################################\n",
                           "Your choice (1,2,3,q to quit): "))
  
  
  if (choix == "q") {
    return(cat("Au revoir"))
    }
  
  
  else if (choix == 1){
    # GRAPHE ALL SAMPLES
    titre = paste("Spectre Infrarouge de tous les échantillons")
    
    # Absorbance en fonction du NO (dans l'ordre décroissant)
    # On groupe et colore par sample
    graph_allsamples = ggplot(SN, aes(x = reorder(NO,desc(NO)), y = Absorbance, group = factor(Sample))) +
      geom_line(aes(color=factor(Sample))) +
      labs(title = titre,
           x = "Nombre d'ondes (1/cm)",
           y = "Absorbance normalisée")+
      scale_x_discrete(breaks = as.integer(NO[seq(1,length(NO),by =15)]))
    
    print(graph_allsamples)
    #ggsave(paste(titre, "jpg", sep ="."))
  }
  
  else if (choix == 2){
    #GRAPHE MOYENNE DES POPULATIONS

    titre = paste0("Spectre Infrarouge Moyen des différentes populations")
    
    # Absorbance en fonction du NO (dans l'ordre décroissant)
    # On ne prend que les moyennes
    graphe_moy = ggplot(filter(SN, Echantillon == 0)) +
      aes(x = reorder(NO,desc(NO)), y = Absorbance, group = Population, color = Population) + 
      geom_line() +
      labs(title = titre,
           x = "Nombre d'ondes (1/cm)",
           y = "Absorbance normalisée")+
      scale_x_discrete(breaks = as.integer(NO[seq(1,length(NO),by =15)]))
    
    print(graphe_moy)
    # Sauvergade sous Jpg
    #ggsave(paste(titre, "jpg", sep ="."))
  }
  
  else if (choix == 3){
    
    pop = readline("Quelle population voulait vous voir ? : ")
    
    titre = paste0("Spectre Infrarouge Moyen de ", pop)
    
    if (tolower(pop) %in% tolower(levels(SN$Population))){
      
      graphe_mean= ggplot(filter(SN, Echantillon == 0 & Population == pop)) +
      aes(x = reorder(NO,desc(NO)), y = Absorbance, group = Population) + 
      geom_line() +
      labs(title = titre,
           x = "Nombre d'ondes (1/cm)",
           y = "Absorbance normalisée")+
      scale_x_discrete(breaks = as.integer(NO[seq(1,length(NO),by =15)]))
    
    print(graphe_mean)
    }
    
    else{
      cat("Veuillez rentrer le nom d'une population existante !!! ")
    }
    
  }
  else if (choix == 4){
    #COMPARAISON DE MOYENNE DE 2 POPULATIONS (GRAPHE + T-TEST)
    Ttest()
    
  }
  
  else cat("Veuillez rentrer 1,2,3 ou q !")
  
  Choix()
}

Ttest <- function(){
  pop1 <- readline(prompt = "Population de référence : ")
  
  if (!(tolower(pop1) %in% tolower(levels(SN$Population)))) {
    cat("Veuillez rentrer un nom de population qui existe !")
    Ttest()
  }
  
  #Calcul du t-test
  # Attention les données ne suivent pas forcément un loi normale centrée réduite (il faut valider par la biologie)
  else{
    # Degre de liberte pour le controle
    ddlCTRL <-  table(filter(SN, NO == NO[1])$Population)[pop1] - 3
    
    # Effectif du controle
    effCTRL <-  table(filter(SN, NO == NO[1])$Population)[pop1] - 2 #On prend le nb d'indiv sur la première longueur d'onde et on enlève mean et var
    
    # Moyenne de la population de référence
    MeanCTRL <-  select(filter(SN, Sample == pop1),
                        Absorbance)
    nrow(MeanCTRL)
    
    # Variance de la population de référence
    VarCTRL <-  select(filter(SN, Sample == paste(pop1,"var", sep = "-")),
                       Absorbance)
    nrow(VarCTRL)
    
    for (pop2 in (levels(SN$Population)[levels(SN$Population) != pop1])){
      #creation du data-frame pour les valeurs des t-tests
      #On le renouvelle pour chaque comparaison
      dft <- data.frame (NO = as.integer(levels(SN$NO)),
                         MeanCTRL = MeanCTRL,
                         VarCTRL = VarCTRL)
      
      
      # du degre de liberte
      ddlMUT <- table(filter(SN, NO == NO[1])$Population)[pop2] - 3
      
      # de l'effectif
      effMUT <- table(filter(SN, NO == NO[1])$Population)[pop2] - 2
      
      # de la moyenne par colonne
      MeanMUT <- select(filter(SN, Sample == pop2),
                        Absorbance)
      dft <- cbind(dft, MeanMUT)
      
      # de la variance par colonne
      VarMUT <- select(filter(SN, Sample == paste(pop2, "var", sep = "-")),
                       Absorbance)
      dft <- cbind(dft, VarMUT)
      
      # Calcul de la variance empirique S-carre
      S <- ((ddlCTRL*VarCTRL)+(ddlMUT*VarMUT))/(ddlCTRL+ddlMUT)
      dft <- cbind(dft, S)
      
      # Calcul de la statistique de test (Student)
      T <- (MeanCTRL-MeanMUT)/sqrt((S/effCTRL)+(S/effMUT))
      dft <- cbind(dft, T)
      
      colnames(dft) <- c("NO", "MeanCTRL", "VarCTRL", "MeanMUT", "VarMUT", "S","T")
      
      #Correction des valeurs de T en fonction des valeurs de S
      #si la variance empirique est nulle le test est nul
      #dft$T[dft$S == 0] <- 0
      
      # Calcul du quantile de Stutent pour un test alpha=5%
      Alpha = 0.05
      Quantil <- qt((Alpha/2),(effMUT + effCTRL)-2)
      
      # Tests multiples car Student calculee point par point 
      #et non sur l'ensemble du spectre, doit subir une correction !!
      Quantil2 <-  qt(((Alpha/2)/ Dim[2]),(effMUT + effCTRL)-2);
      
      
      titre = "Spectres Infrarouges Moyens"
      Mean = ggplot(filter(SN, Echantillon == 0 & (Population == pop1 | Population == pop2))) +
        aes(x = reorder(NO,desc(NO)), y = Absorbance, group = Population, color = Population) + 
        geom_line() +
        labs(title = titre,
             x = "Nombre d'ondes (1/cm)",
             y = "Absorbance normalisée")+
        theme(legend.position="top",
              legend.title = element_blank(),
              plot.title = element_text(hjust = 0.5))+
        scale_x_discrete(breaks = as.integer(NO[seq(1,length(NO),by =15)])) # Pour avoir un axe x lisible
      
      titre = "T-test "
      T_test = ggplot(dft) +
        aes(x = reorder(NO,desc(NO)), y = T, group = 1) + 
        geom_line() +
        labs(title = titre,
             x = "Nombre d'ondes (1/cm)",
             y = "t-value")+
        geom_hline(yintercept = 2*Quantil, colour = "red") + # On rajoute les lignes des quantils
        geom_hline(yintercept = 2*abs(Quantil), colour = "red") +
        theme(plot.title = element_text(hjust = 0.5))+
        scale_x_discrete(breaks = as.integer(NO[seq(1,length(NO),by =15)])) # Pour avoir un axe x lisible
      
      # On aligne les 2 graphes
      ploot <- plot_grid(Mean, T_test, ncol = 1, align = "v")
      print(ploot)
      save_plot(paste0(pop1, "_vs_", pop2, ".jpg"), ploot)
    }
  }
}
