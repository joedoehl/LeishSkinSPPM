#Generates all tables in which data is going to be combined
XYCoxList1 <- list();
XYCoxMod1 <- c();
XYCoxAIC1 <- c();
XYCauList1 <- list();
XYCauMod1 <- c();
XYCauAIC1 <- c();
XYMatList1 <- list();
XYMatMod1 <- c();
XYMatAIC1 <- c();
XYTomList1 <- list();
XYTomMod1 <- c();
XYTomAIC1 <- c();
XYVarList1 <- list();
XYVarMod1 <- c();
XYVarAIC1 <- c();

source("zipFastener.R");

# Applying different model approaches to data to find best fit
for (ww in 1:csvCount) {
  print(ww); #to know in which cycle the code is
  
  #Log-Gaussian Cox Process Model
  XYCox <- kppm(listXY[[paste0("RAG", ww)]] ~ polynom(x, y, 1), "LGCP", method = 'clik2'); #Model generation
  tryCatch({ #In case a model is invalid it will through an error in the "step" function
    XYCox <- step(XYCox, trace = 0); wm <<- 0},
    error = function(e) { wm <<- 1 });
  if (wm == 0) {
    XYCox <- step(XYCox, trace = 0); #Stepwise best model identification by Akaike’s Information Criterion (AIC)
    XYCoxList1 <- list.append(XYCoxList1, XYCox); #Stores models in a list
    XYCoxMod1 <- c(XYCoxMod1, as.character(XYCox$po$trend)); #Extracts model type after AIC
    XYCoxAIC1 <- c(XYCoxAIC1, tail(XYCox$anova$AIC, n=1)); #Extracts AIC value of last step for final model
  } else if (wm == 1) {
    XYCoxList1 <- list.append(XYCoxList1, "invalid"); #In case model is invalid
    XYCoxMod1 <- c(XYCoxMod1, "invalid"); #In case model is invalid
    XYCoxAIC1 <- c(XYCoxAIC1, "invalid"); #In case model is invalid
  }

  #Cluster Models
  
  #Cauchy Cluster Process Model
  XYCau <- kppm(listXY[[paste0("RAG", ww)]] ~ polynom(x, y, 1), "Cauchy", method = 'clik2'); #Model generation
  tryCatch({ #In case a model is invalid it will through an error in the "step" function
    XYCau <- step(XYCau, trace = 0); wm <<- 0},
    error = function(e) { wm <<- 1 });
  if (wm == 0) {
    XYCau <- step(XYCau, trace = 0); #Stepwise best model identification by Akaike’s Information Criterion (AIC)
    XYCauList1 <- list.append(XYCauList1, XYCau); #Stores models in a list
    XYCauMod1 <- c(XYCauMod1, as.character(XYCau$po$trend)); #Extracts model type after AIC
    XYCauAIC1 <- c(XYCauAIC1, tail(XYCau$anova$AIC, n=1)); #Extracts AIC value of last step for final model
  } else if (wm == 1) {
    XYCauList1 <- list.append(XYCauList1, "invalid"); #In case model is invalid
    XYCauMod1 <- c(XYCauMod1, "invalid"); #In case model is invalid
    XYCauAIC1 <- c(XYCauAIC1, "invalid"); #In case model is invalid
  }

  #Matern Cluster Process Model
  XYMat <- kppm(listXY[[paste0("RAG", ww)]] ~ polynom(x, y, 1), "MatClust", method = 'clik2'); #Model generation
  tryCatch({ #In case a model is invalid it will through an error in the "step" function
    XYMat <- step(XYMat, trace = 0); wm <<- 0},
    error = function(e) { wm <<- 1 });
  if (wm == 0) {
    XYMat <- step(XYMat, trace = 0); #Stepwise best model identification by Akaike’s Information Criterion (AIC)
    XYMatList1 <- list.append(XYMatList1, XYMat); #Stores models in a list
    XYMatMod1 <- c(XYMatMod1, as.character(XYMat$po$trend)); #Extracts model type after AIC
    XYMatAIC1 <- c(XYMatAIC1, tail(XYMat$anova$AIC, n=1)); #Extracts AIC value of last step for final model
  } else if (wm == 1) {
    XYMatList1 <- list.append(XYMatList1, "invalid"); #In case model is invalid
    XYMatMod1 <- c(XYMatMod1, "invalid"); #In case model is invalid
    XYMatAIC1 <- c(XYMatAIC1, "invalid"); #In case model is invalid
  }

  #Thomas Process Model
  XYTom <- kppm(listXY[[paste0("RAG", ww)]] ~ polynom(x, y, 1), "Thomas", method = 'clik2'); #Model generation
  tryCatch({ #In case a model is invalid it will through an error in the "step" function
    XYTom <- step(XYTom, trace = 0); wm <<- 0},
    error = function(e) { wm <<- 1 });
  if (wm == 0) {
    XYTom <- step(XYTom, trace = 0); #Stepwise best model identification by Akaike’s Information Criterion (AIC)
    XYTomList1 <- list.append(XYTomList1, XYTom); #Stores models in a list
    XYTomMod1 <- c(XYTomMod1, as.character(XYTom$po$trend)); #Extracts model type after AIC
    XYTomAIC1 <- c(XYTomAIC1, tail(XYTom$anova$AIC, n=1)); #Extracts AIC value of last step for final model
  } else if (wm == 1) {
    XYTomList1 <- list.append(XYTomList1, "invalid"); #In case model is invalid
    XYTomMod1 <- c(XYTomMod1, "invalid"); #In case model is invalid
    XYTomAIC1 <- c(XYTomAIC1, "invalid"); #In case model is invalid
  }

  #Variance Gamma Cluster Process Model
  XYVar <- kppm(listXY[[paste0("RAG", ww)]] ~ polynom(x, y, 1), "VarGamma", method = 'clik2'); #Model generation
  tryCatch({ #In case a model is invalid it will through an error in the "step" function
    XYVar <- step(XYVar, trace = 0); wm <<- 0},
    error = function(e) { wm <<- 1 });
  if (wm == 0) {
    XYVar <- step(XYVar, trace = 0); #Stepwise best model identification by Akaike’s Information Criterion (AIC)
    XYVarList1 <- list.append(XYVarList1, XYVar); #Stores models in a list
    XYVarMod1 <- c(XYVarMod1, as.character(XYVar$po$trend)); #Extracts model type after AIC
    XYVarAIC1 <- c(XYVarAIC1, tail(XYVar$anova$AIC, n=1)); #Extracts AIC value of last step for final model
  } else if (wm == 1) {
    XYVarList1 <- list.append(XYVarList1, "invalid"); #In case model is invalid
    XYVarMod1 <- c(XYVarMod1, "invalid"); #In case model is invalid
    XYVarAIC1 <- c(XYVarAIC1, "invalid"); #In case model is invalid
  }
}
XYModels1 <- cbind.data.frame("Log-Gaussian Cox Process" = XYCoxMod1, "Cauchy Cluster Process" = XYCauMod1, "Matern Cluster Process" = XYMatMod1,
                              "Thomas Cluster Process" = XYTomMod1, "Variance Gamma Cluster Process" = XYVarMod1);
XYAICs1 <- cbind.data.frame("Log-Gaussian Cox Process" = XYCoxAIC1, "Cauchy Cluster Process" = XYCauAIC1, "Matern Cluster Process" = XYMatAIC1,
                            "Thomas Cluster Process" = XYTomAIC1, "Variance Gamma Cluster Process" = XYVarAIC1);

XYModOutput1 <- zipFastener(XYModels1, XYAICs1,2); #Combines both data frames by adding columns in alternation; index 1 by row, index 2 by column

names(XYModOutput1)[seq(2, ncol(XYModOutput1), 2)] <- "AIC";

rownames(XYModels1) <- MouseIDs;
rownames(XYAICs1) <- MouseIDs;
rownames(XYModOutput1) <- MouseIDs;

#Saves the output data compilation file
write.csv(XYModels1, "Models 1.csv");
write.csv(XYAICs1, "Model AICs 1.csv");
write.csv(XYModOutput1, "Model Data 1.csv");

names(XYCoxList1) <- MouseIDs;
names(XYCauList1) <- MouseIDs;
names(XYMatList1) <- MouseIDs;
names(XYTomList1) <- MouseIDs;
names(XYVarList1) <- MouseIDs;

rm(XYCox, XYCau, XYMat, XYTom, XYVar, XYCoxMod1, XYCauMod1, XYMatMod1, XYTomMod1, XYVarMod1, XYCoxAIC1, XYCauAIC1, XYMatAIC1,
   XYTomAIC1, XYVarAIC1);

# Model test with graphic representation ----------------------------------------------------------------------------------------------------------------------------------

#Generates all tables in which data is going to be combined

XYCoxEnvInhList1 <- list();
XYCauEnvInhList1 <- list();
XYMatEnvInhList1 <- list();
XYTomEnvInhList1 <- list();
XYVarEnvInhList1 <- list();

# Generate graphs to the inhomogeneous models

for (ww in 1:csvCount) {
  print(ww); #to know in which cycle the code is
  
  #Log-Gaussian Cox Process Model
  wm <- XYCoxList1[[paste0("RAG", ww)]] == "invalid";
  if (wm[1] == TRUE) {
    XYCoxEnvInhList1 <- list.append(XYCoxEnvInhList1, "invalid");
  } else {
    XYCoxEnv <- envelope.kppm(XYCoxList1[[paste0("RAG", ww)]], Linhom, simulate = expression(rpoispp(DDlist[[paste0("RAG", ww)]]$estimate)), use.theory = TRUE, nsim = 39, global = TRUE);
    XYCoxEnvInhList1 <- list.append(XYCoxEnvInhList1, XYCoxEnv);
    #To export the resulting plot
    Graph.Save(XYCoxEnv, ww, 1, "XYCoxInh", "Fitted Model");
  }
  
  #Cluster Models
  
  #Cuachy Process Model
  wm <- XYCauList1[[paste0("RAG", ww)]] == "invalid";
  if (wm[1] == TRUE) {
    XYCauEnvInhList1 <- list.append(XYCauEnvInhList1, "invalid");
  } else {
    XYCauEnv <- envelope.kppm(XYCauList1[[paste0("RAG", ww)]], Linhom, simulate = expression(rpoispp(DDlist[[paste0("RAG", ww)]]$estimate)), use.theory = TRUE, nsim = 39, global = TRUE);
    XYCauEnvInhList1 <- list.append(XYCauEnvInhList1, XYCauEnv);
    #To export the resulting plot
    Graph.Save(XYCauEnv, ww, 1, "XYCauInh", "Fitted Model");
  }
  
  #Matern Cluster Process Model
  wm <- XYMatList1[[paste0("RAG", ww)]] == "invalid";
  if (wm[1] == TRUE) {
    XYMatEnvInhList1 <- list.append(XYMatEnvInhList1, "invalid");
  } else {
    XYMatEnv <- envelope.kppm(XYMatList1[[paste0("RAG", ww)]], Linhom, simulate = expression(rpoispp(DDlist[[paste0("RAG", ww)]]$estimate)), use.theory = TRUE, nsim = 39, global = TRUE);
    XYMatEnvInhList1 <- list.append(XYMatEnvInhList1, XYMatEnv);
    #To export the resulting plot
    Graph.Save(XYMatEnv, ww, 1, "XYMatInh", "Fitted Model");
  }
  
  #Thomas Process Model
  wm <- XYTomList1[[paste0("RAG", ww)]] == "invalid";
  if (wm[1] == TRUE) {
    XYTomEnvInhList1 <- list.append(XYTomEnvInhList1, "invalid");
  } else {
    XYTomEnv <- envelope.kppm(XYTomList1[[paste0("RAG", ww)]], Linhom, simulate = expression(rpoispp(DDlist[[paste0("RAG", ww)]]$estimate)), use.theory = TRUE, nsim = 39, global = TRUE);
    XYTomEnvInhList1 <- list.append(XYTomEnvInhList1, XYTomEnv);
    #To export the resulting plot
    Graph.Save(XYTomEnv, ww, 1, "XYTomInh", "Fitted Model");
  }
  
  #Variance Gamma Cluster Process Model
  wm <- XYVarList1[[paste0("RAG", ww)]] == "invalid";
  if (wm[1] == TRUE) {
    XYVarEnvInhList1 <- list.append(XYVarEnvInhList1, "invalid");
  } else {
    XYVarEnv <- envelope.kppm(XYVarList1[[paste0("RAG", ww)]], Linhom, simulate = expression(rpoispp(DDlist[[paste0("RAG", ww)]]$estimate)), use.theory = TRUE, nsim = 39, global = TRUE);
    XYVarEnvInhList1 <- list.append(XYVarEnvInhList1, XYVarEnv);
    #To export the resulting plot
    Graph.Save(XYVarEnv, ww, 1, "XYVarInh", "Fitted Model");
  }
}

names(XYCoxEnvInhList1) <- MouseIDs;
names(XYCauEnvInhList1) <- MouseIDs;
names(XYMatEnvInhList1) <- MouseIDs;
names(XYTomEnvInhList1) <- MouseIDs;
names(XYVarEnvInhList1) <- MouseIDs;

rm(XYCoxEnv, XYCauEnv, XYMatEnv, XYTomEnv, XYVarEnv);

#----------------------------------------------------------------------------------------------------------------------------------

#Generates all tables in which data is going to be combined

XYCoxEnvScaList1 <- list();
XYCauEnvScaList1 <- list();
XYMatEnvScaList1 <- list();
XYTomEnvScaList1 <- list();
XYVarEnvScaList1 <- list();

# Generate graphs to the locally scaled models

for (ww in 1:csvCount) {
  print(ww); #to know in which cycle the code is
  
  #Log-Gaussian Cox Process Model
  wm <- XYCoxList1[[paste0("RAG", ww)]] == "invalid";
  if (wm[1] == TRUE) {
    XYCoxEnvScaList1 <- list.append(XYCoxEnvScaList1, "invalid");
  } else {
    XYCoxEnv <- envelope.kppm(XYCoxList1[[paste0("RAG", ww)]], Lscaled, simulate = expression(rpoispp(DDlist[[paste0("RAG", ww)]]$estimate)), use.theory = TRUE, nsim = 39, global = TRUE);
    XYCoxEnvScaList1 <- list.append(XYCoxEnvScaList1, XYCoxEnv);
    #To export the resulting plot
    Graph.Save(XYCoxEnv, ww, 1, "XYCoxSca", "Fitted Model");
  }
  
  #Cluster Models
  
  #Cuachy Process Model
  wm <- XYCauList1[[paste0("RAG", ww)]] == "invalid";
  if (wm[1] == TRUE) {
    XYCauEnvScaList1 <- list.append(XYCauEnvScaList1, "invalid");
  } else {
    XYCauEnv <- envelope.kppm(XYCauList1[[paste0("RAG", ww)]], Lscaled, simulate = expression(rpoispp(DDlist[[paste0("RAG", ww)]]$estimate)), use.theory = TRUE, nsim = 39, global = TRUE);
    XYCauEnvScaList1 <- list.append(XYCauEnvScaList1, XYCauEnv);
    #To export the resulting plot
    Graph.Save(XYCauEnv, ww, 1, "XYCauSca", "Fitted Model");
  }
  
  #Matern Cluster Process Model
  wm <- XYMatList1[[paste0("RAG", ww)]] == "invalid";
  if (wm[1] == TRUE) {
    XYMatEnvScaList1 <- list.append(XYMatEnvScaList1, "invalid");
  } else {
    XYMatEnv <- envelope.kppm(XYMatList1[[paste0("RAG", ww)]], Lscaled, simulate = expression(rpoispp(DDlist[[paste0("RAG", ww)]]$estimate)), use.theory = TRUE, nsim = 39, global = TRUE);
    XYMatEnvScaList1 <- list.append(XYMatEnvScaList1, XYMatEnv);
    #To export the resulting plot
    Graph.Save(XYMatEnv, ww, 1, "XYMatSca", "Fitted Model");
  }
  
  #Thomas Process Model
  wm <- XYTomList1[[paste0("RAG", ww)]] == "invalid";
  if (wm[1] == TRUE) {
    XYTomEnvScaList1 <- list.append(XYTomEnvScaList1, "invalid");
  } else {
    XYTomEnv <- envelope.kppm(XYTomList1[[paste0("RAG", ww)]], Lscaled, simulate = expression(rpoispp(DDlist[[paste0("RAG", ww)]]$estimate)), use.theory = TRUE, nsim = 39, global = TRUE);
    XYTomEnvScaList1 <- list.append(XYTomEnvScaList1, XYTomEnv);
    #To export the resulting plot
    Graph.Save(XYTomEnv, ww, 1, "XYTomSca", "Fitted Model");
  }
  
  #Variance Gamma Cluster Process Model
  wm <- XYVarList1[[paste0("RAG", ww)]] == "invalid";
  if (wm[1] == TRUE) {
    XYVarEnvScaList1 <- list.append(XYVarEnvScaList1, "invalid");
  } else {
    XYVarEnv <- envelope.kppm(XYVarList1[[paste0("RAG", ww)]], Lscaled, simulate = expression(rpoispp(DDlist[[paste0("RAG", ww)]]$estimate)), use.theory = TRUE, nsim = 39, global = TRUE);
    XYVarEnvScaList1 <- list.append(XYVarEnvScaList1, XYVarEnv);
    #To export the resulting plot
    Graph.Save(XYVarEnv, ww, 1, "XYVarSca", "Fitted Model");
  }
}

names(XYCoxEnvScaList1) <- MouseIDs;
names(XYCauEnvScaList1) <- MouseIDs;
names(XYMatEnvScaList1) <- MouseIDs;
names(XYTomEnvScaList1) <- MouseIDs;
names(XYVarEnvScaList1) <- MouseIDs;

rm(XYCoxEnv, XYCauEnv, XYMatEnv, XYTomEnv, XYVarEnv);
