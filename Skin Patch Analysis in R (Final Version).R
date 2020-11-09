# TITLE ################################################################################################################################################################################################################################################################################################################

# Title:        Leishmania donovani macroscopic skin patch analysis

# Script Name:  Skin Patch Analysis in R.R
# Copyright:    Johannes S. P. Doehl (NIAID, USA) & Paul M. Kaye (University of York, UK)

# Purpose:      Exploration of the Leishmania parasite landscape in the mammalian host skin based on ImageJ data output
# Author:       Johannes S. P. Doehl
# Date:         May 2020

# Publication:

########################################################################################################################################################################################################################################################################################################################

BiocManager::install() #updates packages

update.packages(ask = FALSE, repos = 'http://cran.rstudio.org'); #updates packages; needs the project to be closed to be executed

# 1. LIBRARIES to load##################################################################################################################################################################################################################################################################################################

#This overcomes the problem of Rstudio not detecting Rtools, when R packages are installed anew or called from library
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

#Package to be opened before executing the script

#Packages improving R performance and versatility
# library(dplyr); #required for efficient merge of data frames in alternate column fashion
library(tidyverse); #Improves graphics and coding abilities; contains several base packages: ggplot2, dplyr, tidyr, readr, purrr, tibble, stingr, forcast*
library(rlist); #enhances manipulation capabilities of lists*
# library(units); #needs to be installed when R is first installed; required to assign units to vectors
library(scales); #improves graph building
# library(formula.tools);#required for efficient conversion of extracted formulas to characters
# library(gt); #improves table building capabilities
# library(gtsummary); #is an expansion of the 'gt' package
# library(webshot); #needed to save gt tables
# library(glue); #add on for table building
# library(stringr); #allow extraction of numbers from number/letter string

#Statistic packages
library(e1071); #required for skewness testing*
library(stats); #calls up all standard stats tests already supplied in R
# library(dunn.test); #required for Dunn's post hoc test*
# library(SuppDists); #adds 10 distributions to the already available once to R*
# library(car); #required for Levene's Test & Fligner's Test*
library(ggpubr); #required for ANOVA function aov()*
library(rms.gof); #required for the root-mean-square alternative to chi-square

#Packages required for spatial point pattern analysis
library(spatstat); #required for spatial point pattern analysis*

#These library are required by spatstat
library(deldir); #Calculates the Delaunay triangulation and the Dirichlet or Voronoi tessellation*
library(abind); #This is a generalization of 'cbind' and 'rbind'*
library(tensor); #*
library(polyclip); #Performs polygon clipping operations*
library(goftest); #Classical Goodness-of-Fit Tests for Univariate Distributions*

#Additional packages to support versatility in spatstat
library(sp); #Classes and methods for spatial data*
library(sparr); #Provides functions to estimate kernel-smoothed spatial and spatio-temporal densities and relative risk functions, and perform subsequent inference*
library(graphics); #improves graphics capabilites*
# library(DCluster); #A set of functions for the detection of spatial clusters of disease using count data*
# library(ecespa); #Improves the functionality of some wrappers*
# library(splines); #Required for smoothing in alternative model fitting (p 340)*
# library(ppmlasso); #Required for the lasso model (p 339)*
# library(lgcp); #Required for Log-Gaussian Cox Process (LGCP) (p 471)*
library(RandomFields); #Required for wide range of Gaussian Random Field models (p 455)*
library(spatstat.utils); #Utilities for spatstat to enhance capabilities

# *these packages require prior installation

# 2. VARIABLES##########################################################################################################################################################################################################################################################################################################

DenCut <- 70; #This is the detection threshold in greyscale for the patch detection in ImageJ, which is subtracted from the Mean Density of each patch
ImWid <- 4; #Indicates the column in Table 2, where the original Image Width in µm is stored
ImHei <- 5; #Indicates the column in Table 2, where the original Image Height in µm is stored
ConvSma <- 21.47444; #ImageJ saved coordinates of Skin Outline in pixel counts; this is the converter to µm for a half sized image
ConvLar <- 10.78725; #ImageJ saved coordinates of Skin Outline in pixel counts; this is the converter to µm for a full sized image
ConvSca1 <- 1000; #Converts the image dimensions from µm to mm
ConvSca2 <- 1000000; #Covnerts dimensinos from µm^2 to mm^2
ConvSca3 <- 100000000; #Covnerts dimensinos from µm^2 to cm^2
ImageRes <- 1200; #Resultion for skin images
GraphRes <- 300; #Resolution for Graphs
GridX <- 5; #determines the horizontal grid size for the qadrat test and quadrat counting
GridY <- 6; #determines the vertical grid size for the qadrat test and quadrat counting
maxDis <- 3.9; #Marks the maximum distance allow between dots in nearest-neighbour calculations
StatX2 <- 3; #Sets grid for correlation-stationary test
StatY2 <- 4; #Sets grid for correlation-stationary test
StatTest <- 5; #Number of repeats of correlation-stationary test


# 3. SCRIPT#############################################################################################################################################################################################################################################################################################################

# 3.1. Data extraction, manipulation and preparation for statistical analysis===========================================================================================================================================================================================================================================

#3.1.1. Data Location
BiblioPath1 <- "E:/Work/Ongoing Work/Joe (Work for Paul Kaye)/Written Article/4 - Paper with Paul on skin parasite distribution/Data/Stereo/Data/OutputWhole(Results)"; #!!!!!! PATH MUST BE CHANGED TO YOUR DATA LOCATION !!!!!!
# BiblioPath1 <- "D:/Lab Work/Ongoing Work/Joe (Work for Paul Kaye)/Written Article/4 - Paper with Paul on skin parasite distribution/Data/Stereo/Data/OutputWhole(Results)"; #Alternative data location

csvCount <- length(list.files(path = BiblioPath1, pattern = "csv", recursive = FALSE, full.names = TRUE)); #determines number of sample data by counting data files from imageJ

options(scipen=999); #prevents scientific notation

# 3.1.2. Data extraction, resorting and preparing for use---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Generates all tables in which data is going to be combined
DataList <- list(); #Creates a list in which the raw data is stored of all samples
MouseIDs <- c(); #Creates a vector, in which the samples names are going to be loaded
AllPatchCounts <- c(); #Creates a vector, in which the Patch Counts will be loaded
TotalSkinArea <- c(); #Creates a vector, in which all Total Skin Areas (in µm^2) are loaded
TotalPatchArea <- c(); #Creates a vector, in which all Total Patch Areas (in µm^2) are loaded
AllAreas <- c(); #Creates a vector, in which all individual Patch Areas (in µm^2) are loaded from all samples in a single column
AreaSummary <- data.frame(); #Creates a data frame, in which all the Area Summary data will be loaded
AreasList <- list(); #Creates a list, in which all individual Patch Areas are loaded by sample per column
AllDensities <- c(); #Creates a vector, in which all individual Mean Patch Densities are loaded from all samples in a single column
DensityList <- list(); #Creates a list, in which all individual Mean Patch Densities are loaded by sample per column
AllNormDens <- c(); #Creates a vector, in which all individual Mean Patch normalized Densities are loaded from all samples in a single column
NormDensList <- list(); #Creates a list, in which all individual Mean normalized Patch Densities are loaded by sample per column
AllMinDia <- c(); #Creates a vector, for shortest diameter of patch based Feret's diameter calculation
MinDiaList <- list(); #Creates a list, for shortest diameter of patch based Feret's diameter calculation
MinDiaSummary <- data.frame(); #Summary statistics for shortest Feret's diameter
AllMaxDia <- c(); #Creates a vector, for longest diameter of patch based Feret's diameter calculation
MaxDiaList <- list(); #Creates a list, for longest diameter of patch based Feret's diameter calculation
MaxDiaSummary <- data.frame(); #Summary statistics for longest Feret's diameter
ArDenList <- list(); #Creates a list of all Area & Density Marks
XMYMList <- list(); #Creates a list of all X & Y coordinates
TableXY <- list(); #Table for saving X & Y coordinates
TableXYIDs <- c(); #IDs for TableXY
DupCheck1 <- c(); #Creates a vector, in which all results of the duplicate coordinates check is stored
DupCheck2 <- c(); #Creates a vector, in which all results of the duplicate coordinates check is stored

#Opens spread sheets and extract relevant data columns
for (aa in 1:csvCount) {
  
  #Load data file (.csv)
  RAG <- read.csv(paste0(BiblioPath1,"/Rag",aa,"_Results.csv")); #Loads the .csv file into R and saves it into a variable
  DataList <- list.append(DataList, RAG); #List of raw data
  MouseIDs <- c(MouseIDs, paste0("RAG", aa)); #Adds sample name to a character vector
  
  #Data extraction by column
  AREA <- RAG$Area; #extracts Area column as a numeric vector that contains the patch area data (in µm^2)
  DENSITY <- RAG$Mean; #extracts Area column as a numeric vector that contains the mean patch density data (in grayscale [0-255])
  XM <- RAG$XM; #extracts XM column as a numeric vectorand places it into a variable; data is in µm
  YM <- RAG$YM; #extracts YM column as a numeric vector and places it into a variable; data is in µm
  MINDIA <- RAG$Minor; #extracts the column of the shortest Feret's diameter
  MAXDIA <- RAG$Major; #extracts the column of the longest Feret's diameter
  
  #Data manipulation
  #Patch Counts
  PatchCounts <- length(AREA)-2; #Give me the number of patches detected in the skin by removing measurements for total patch are and total skin area
  AllPatchCounts <- c(AllPatchCounts, PatchCounts); #Adds patch count to a vector
  
  #Patch Area
  SkinArea <- AREA[AllPatchCounts[aa]+2]; #Gets total Skin Area from Area vector (in µm^2)
  TotalSkinArea <- c(TotalSkinArea, SkinArea); #Adds total skin are to a vector
  
  PatchArea <- AREA[AllPatchCounts[aa]+1]; #Gets total Patch Area from Area vector (in µm^2)
  TotalPatchArea <- c(TotalPatchArea, PatchArea); #Adds total skin are to a vector
  
  AREA <- AREA[1:AllPatchCounts[aa]]; #Removes the last two rows, which are total patch and skin area
  AllAreas <- c(AllAreas, AREA);
  
  ArSum <- summary(AREA); #Calculates the median of the Area vector
  AreaSummary <- rbind(AreaSummary, ArSum); #Collects all summary data into a data frame
  rm(ArSum);
  
  AreasList <- list.append(AreasList, AREA); #Creates a list of all area variables, which will be converted into a data frame to build a table of different variable lengths
  
  #Patch Density
  DENSITY <- DENSITY[1:AllPatchCounts[aa]]; #Removes the last two rows, which are total patch and skin area
  AllDensities <- c(AllDensities, DENSITY);
  DensityList <- list.append(DensityList, DENSITY); #Creates a list of all density variables, which will be converted into a data frame to build a table of different variable lengths
  
  NormDens <- DENSITY-DenCut; #Converts all individual mean patch densities, by subtracting the threshold, rendering only the value over the threshold
  AllNormDens <- c(AllNormDens, NormDens);
  NormDensList <- list.append(NormDensList, NormDens); #Creates a list of all normalized density variables, which will be converted into a data frame to build a table of different variable lengths
  
  #Patch min./max. Feret's diameter
  MINDIA <- MINDIA[1:AllPatchCounts[aa]]; #Removes the last two rows, which are total patch and skin area
  AllMinDia <- c(AllMinDia, MINDIA);
  MinDiaList <- list.append(MinDiaList, MINDIA);
  
  MinDiaSum <- summary(MINDIA); #Calculates the median of the Area vector
  MinDiaSummary <- rbind(MinDiaSummary, MinDiaSum); #Collects all summary data into a data frame
  rm(MinDiaSum);
  
  MAXDIA <- MAXDIA[1:AllPatchCounts[aa]]; #Removes the last two rows, which are total patch and skin area
  AllMaxDia <- c(AllMaxDia, MAXDIA);
  MaxDiaList <- list.append(MaxDiaList, MAXDIA);
  
  MaxDiaSum <- summary(MAXDIA); #Calculates the median of the Area vector
  MaxDiaSummary <- rbind(MaxDiaSummary, MaxDiaSum); #Collects all summary data into a data frame
  rm(MaxDiaSum);
  
  #Create Marks
  ArDen <- data.frame("Area" = AREA, "Density" = NormDens); #Creates a double mark for the spatial point pattern analysis
  ArDenList <- list.append(ArDenList, ArDen);
  
  #X and Y coordinates (in µm) of center of mass of individual patches
  XM <- XM[1:AllPatchCounts[aa]]; #Removes the last two rows, which are total patch and skin area
  YM <- YM[1:AllPatchCounts[aa]]; #Removes the last two rows, which are total patch and skin area
  XMYMList <- list.append(XMYMList, list("X" = XM, "Y" = YM));
  TableXY <- list.append(TableXY, XM, YM);
  TableXYIDs <- c(TableXYIDs,paste0("RAG", aa, "-X"), paste0("RAG", aa, "-Y"));
  
  #Checks patch coordinates for potential duplicates, which is required for spacial point pattern analysis
  CoorXY <- cbind(XM, YM); #creates a martix with XM and YM coordinates
  colnames(CoorXY) <- c("X", "Y"); #Names the columns in the data frame
  CoorXY1 <- all(duplicated(CoorXY) == FALSE); #Tests whether coordinates occur in duplicates, which is not permitted for many spatstat functions
  if (CoorXY1 == TRUE) {
    CoorXY2 <- "No Duplicates"
  } else if (CoorXY1 == FALSE) {
    CoorXY2 <- "Duplicates"
  }
  DupCheck1 <- c(DupCheck1, CoorXY1); #lists all test results in a single vector: TRUE = no duplicate coordinates found, FALSE = duplicates coordinates found
  DupCheck2 <- c(DupCheck2, CoorXY2); #lists all test results in a single vector: TRUE = no duplicate coordinates found, FALSE = duplicates coordinates found
  
  rm(RAG, AREA, DENSITY, XM, YM, NormDens, PatchCounts, SkinArea, PatchArea, ArDen, CoorXY, CoorXY1, CoorXY2);
}

MouseIDs <- factor(MouseIDs);

names(XMYMList) <- MouseIDs;

names(DataList) <- MouseIDs; #Adds sample names to the raw data tables
colnames(AreaSummary) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."); #Names the columns in the data frame
colnames(MinDiaSummary) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."); #Names the columns in the data frame
colnames(MaxDiaSummary) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."); #Names the columns in the data frame
names(ArDenList) <- MouseIDs; #Creates a list of all Area & Density Marks

#this builds a data table of the patch area data extracted from individual files of each mouse into a single table with different length columns
names(AreasList) <- MouseIDs; #assigns variable names to all list entries
AreasTable <- sapply(AreasList, '[', seq(max(sapply(AreasList, length)))); #converts List into a matrix table
AreasTable[is.na(AreasTable)] <- " "; #delets "NA"
write.csv(AreasTable, "Areas Table.csv"); #saves the collated table as csv file that can be read in excel; it will be safed in the same directory as the data

#this builds a data table of the mean patch density data extracted from individual files of each mouse into a single table with different length columns
names(DensityList) <- MouseIDs; #assigns variable names to all list entries
DensityTable <- sapply(DensityList, '[', seq(max(sapply(DensityList, length)))); #converts List into a matrix table
DensityTable[is.na(DensityTable)] <- " "; #deletes "NA"
write.csv(DensityTable, "Density Table.csv"); #saves the collated table as csv file that can be read in excel; it will be safed in the same directory as the data

#this builds a data table of the mean normalized patch density data extracted from individual files of each mouse into a single table with different length columns
names(NormDensList) <- MouseIDs; #assigns variable names to all list entries
NormDensTable <- sapply(NormDensList, '[', seq(max(sapply(NormDensList, length)))); #converts List into a matrix table
NormDensTable[is.na(NormDensTable)] <- " "; #deletes "NA"
write.csv(NormDensTable, "Norm. Density Table.csv"); #saves the collated table as csv file that can be read in excel; it will be safed in the same directory as the data

#this builds a data table of the mean normalized patch density data extracted from individual files of each mouse into a single table with different length columns
names(MinDiaList) <- MouseIDs; #assigns variable names to all list entries
MinDiaTable <- sapply(MinDiaList, '[', seq(max(sapply(MinDiaList, length)))); #converts List into a matrix table
MinDiaTable[is.na(MinDiaTable)] <- " "; #deletes "NA"
write.csv(MinDiaTable, "Minor Feret's Diameter Table.csv"); #saves the collated table as csv file that can be read in excel; it will be safed in the same directory as the data

#this builds a data table of the mean normalized patch density data extracted from individual files of each mouse into a single table with different length columns
names(MaxDiaList) <- MouseIDs; #assigns variable names to all list entries
MaxDiaTable <- sapply(MaxDiaList, '[', seq(max(sapply(MaxDiaList, length)))); #converts List into a matrix table
MaxDiaTable[is.na(MaxDiaTable)] <- " "; #deletes "NA"
write.csv(MaxDiaTable, "Major Feret's Diameter Table.csv"); #saves the collated table as csv file that can be read in excel; it will be safed in the same directory as the data

#this builds a data table of the mean normalized patch density data extracted from individual files of each mouse into a single table with different length columns
names(TableXY) <- TableXYIDs; #assigns variable names to all list entries
XMYMTable <- sapply(TableXY, '[', seq(max(sapply(TableXY, length)))); #converts List into a matrix table
XMYMTable[is.na(XMYMTable)] <- " "; #deletes "NA"
write.csv(XMYMTable, "XY coordinates Table.csv"); #saves the collated table as csv file that can be read in excel; it will be safed in the same directory as the data

rm(AreasTable, DensityTable, NormDensTable, MinDiaTable, MaxDiaTable, TableXY, TableXYIDs, XMYMTable, aa);

# 3.1.3. Peaks per patch data-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#This extracts the number of high density peaks in a parasite patch from a different set of .csv data files

#Data location
BiblioPath2 <- "E:/Work/Ongoing Work/Joe (Work for Paul Kaye)/Written Article/4 - Paper with Paul on skin parasite distribution/Data/Stereo/Data/OutputPeaks"; #!!!!!! PATH MUST BE CHANGED TO YOUR DATA LOCATION !!!!!!
# BiblioPath2 <- "D:/Lab Work/Ongoing Work/Joe (Work for Paul Kaye)/Written Article/4 - Paper with Paul on skin parasite distribution/Data/Stereo/Data/OutputPeaks"; #Alternative data location

#Generates all tables in which data is going to be combined
AllPeaks <- c(); #Creates a vector, in which all individual Density Peaks per Patch are loaded from all samples in a single column
PeaksList <- list(); #Creates a list, in which all individual Density Peaks per Patch are loaded by sample per column
ArPeaksList <- list(); #Creates a marks list (Area & Peaks)
MarksList <- list(); #Creates a marks list (Area & Norm. Density & Peaks)

for (bb in 1:csvCount) {
  
  #Load data file (.csv)
  Peaks <- read.csv(paste0(BiblioPath2,"/Rag",bb,"_Peak_Counts.csv")); #opens data sheets, NOTE: paste0("Name", j) is the same as paste("Name", j, sep="") as both omit the space between "Name" and j
  
  #Data Manipulation
  #Density Peaks per Patch Area
  PeakCounts <- Peaks$Count; #Extracts the column where peaks per patch are listed
  AllPeaks <- c(AllPeaks, PeakCounts);
  PeaksList <- list.append(PeaksList, PeakCounts);
  
  #Create Karks
  ArPeaks <- data.frame(AreasList[[paste0("RAG", bb)]], "Peaks" = PeakCounts); #Creates a double mark for the spatial point pattern analysis
  ArPeaksList <- list.append(ArPeaksList, ArPeaks);
  
  Marks <- data.frame("Areas" = AreasList[[paste0("RAG", bb)]] / ConvSca2, "Density" = NormDensList[[paste0("RAG", bb)]], "Peaks" = PeakCounts); #Creates a triple mark for the spatial point pattern analysis
  MarksList <- list.append(MarksList, Marks);
  
  rm(Peaks, PeakCounts, ArPeaks, Marks);
}

names(ArPeaksList) <- MouseIDs;
names(MarksList) <- MouseIDs;
  
#this builds a data table of the patch area data extracted from individual files of each mouse into a single table with different length columns
names(PeaksList) <- MouseIDs; #assigns variable names to all list entries
PeaksTable <- sapply(PeaksList, '[', seq(max(sapply(PeaksList, length)))); #converts List into a matrix table
PeaksTable[is.na(PeaksTable)] <- " "; #delets "NA"
write.csv(PeaksTable, "Peaks Table.csv"); #saves the collated table as csv file that can be read in excel; it will be safed in the same directory as the data

rm(PeaksTable, bb);

# 3.1.4. Additional statistics-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Statistical test of data distribution

#Checks each Area data set for normality; as soon as one data set does not complete normality, datasets are Log10 transformed and checked again for normality

#Generates all tables in which data is going to be combined
PV1Results <- c();
PV1Call <- c()
PV2Results <- c();
PV2Call <- c()

for (yy in 1:csvCount) {
  SWT <- shapiro.test(AreasList[[paste0("RAG", yy)]]); #Shapiro-Wilks test of normal data distribution
  
  pVal <- as.numeric(format(round(SWT[["p.value"]], 4), nsmall = 4));
  
  if (pVal < 0.05) {
    PV1Call <- c(PV1Call, "Non-Normal");
    if (pVal < 0.0001) {
      PV1Results <- c(PV1Results, "<0.0001");
    } else {
      PV1Results <- c(PV1Results, pVal);
    }
  } else {
    PV1Call <- c(PV1Call, "Normal");
    PV1Results <- c(PV1Results, pVal);
  }
  SWT <- shapiro.test(log10(AreasList[[paste0("RAG", yy)]])); #Shapiro-Wilks test of log-transformed distribution
  
  pVal <- as.numeric(format(round(SWT[["p.value"]], 4), nsmall = 4));
  
  if (pVal < 0.05) {
    PV2Call <- c(PV2Call, "Non-Normal");
    if (pVal < 0.0001) {
      PV2Results <- c(PV2Results, "<0.0001");
    } else {
      PV2Results <- c(PV2Results, pVal);
    }
  } else {
    PV2Call <- c(PV2Call, "Normal");
    PV2Biop <- c(PV2Biop, pVal);
  }
}
PV1Call <- factor(PV1Call);
PV2Call <- factor(PV2Call);

rm(yy, SWT, pVal);

#Assess skin parasite patch area data skewness

#Generates all tables in which data is going to be combined
AllSkewArea <- c();
AllSkewness <- c();
AllSES <- c();
AllTS <- c();
AllIoP <-  c();

for (xx in 1:csvCount) {
  SkewArea <- as.numeric(format(round(skewness(AreasList[[paste0("RAG", xx)]], na.rm = FALSE, type = 3), 2), nsmall = 2)); #excess skewness of Area data
  AllSkewArea <- c(AllSkewArea, SkewArea);
  SkewLevel <- (-1 > SkewArea | SkewArea > 1); #determines if results means data is skewed or not
  SES<- sqrt((6*AllPatchCounts[xx]*(AllPatchCounts[xx]-1)) / ((AllPatchCounts[xx]-2)*(AllPatchCounts[xx]+1)*(AllPatchCounts[xx]+3))); #calculates the standard error of the skewness (SES)
  AllSES <- c(AllSES, SES);
  TestStat <- SkewArea / SES; #generates a test statistic, which allows to tell if sample skewness is transferable to population skewness
  AllTS <- c(AllTS, TestStat);
  SkewIoP<- (-2 > TestStat | TestStat > 2); #determines if results means data is skewed or not
  
  if (SkewLevel == TRUE) {
    AllSkewness <- append(AllSkewness,"Skewed");
  } else {
    AllSkewness <- append(AllSkewness,"Not Skewed");
  }
  if (SkewIoP == TRUE) {
    AllIoP <- append(AllIoP,"Skewed");
  } else {
    AllIoP <- append(AllIoP,"Not Skewed");
  }
}
AreaSkewAna <- data.frame(AllSkewArea, AllSkewness, AllSES, AllTS, AllIoP);

rm(xx, SkewArea, SkewLevel, SES, TestStat, SkewIoP);

#Assess skin parasite patch area data skewness

#Generates all tables in which data is going to be combined
AllSkewDens <- c();
AllDensSkew <- c();
AllDensSES <- c();
AllDensTS <- c();
AllDensIoP <-  c();

for (xx in 1:csvCount) {
  SkewArea <- as.numeric(format(round(skewness(NormDensList[[paste0("RAG", xx)]], na.rm = FALSE, type = 3), 2), nsmall = 2)); #excess skewness of Area data
  AllSkewDens <- c(AllSkewDens, SkewArea);
  SkewLevel <- (-1 > SkewArea | SkewArea > 1); #determines if results means data is skewed or not
  SES<- sqrt((6*AllPatchCounts[xx]*(AllPatchCounts[xx]-1)) / ((AllPatchCounts[xx]-2)*(AllPatchCounts[xx]+1)*(AllPatchCounts[xx]+3))); #calculates the standard error of the skewness (SES)
  AllDensSES <- c(AllDensSES, SES);
  TestStat <- SkewArea / SES; #generates a test statistic, which allows to tell if sample skewness is transferable to population skewness
  AllDensTS <- c(AllDensTS, TestStat);
  SkewIoP<- (-2 > TestStat | TestStat > 2); #determines if results means data is skewed or not
  
  if (SkewLevel == TRUE) {
    AllDensSkew <- append(AllDensSkew,"Skewed");
  } else {
    AllDensSkew <- append(AllDensSkew,"Not Skewed");
  }
  if (SkewIoP == TRUE) {
    AllDensIoP <- append(AllDensIoP,"Skewed");
  } else {
    AllDensIoP <- append(AllDensIoP,"Not Skewed");
  }
}
AreaSkewAna <- data.frame(AllSkewDens, AllDensSkew, AllDensSES, AllDensTS, AllDensIoP);

rm(xx, SkewArea, SkewLevel, SES, TestStat, SkewIoP);

#Analyzes parasites / biopsy data distribution

BiblioPath4 <- "E:/Work/Ongoing Work/Joe (Work for Paul Kaye)/Written Article/4 - Paper with Paul on skin parasite distribution/Data/Stereo/Data/"; #!!!!!! PATH MUST BE CHANGED TO YOUR DATA LOCATION !!!!!!
# BiblioPath4 <- "D:/Lab Work/Ongoing Work/Joe (Work for Paul Kaye)/Written Article/4 - Paper with Paul on skin parasite distribution/Data/Stereo/Data/"; #Alternative data location

ParaBiop <- read.csv(paste0(BiblioPath4,"Parasites per Biopsy.csv")); #Loads the .csv file into R and saves it into a variable
BiopCount <- ParaBiop[,1];
ParaBiop <- ParaBiop[,c(2:length(ParaBiop))];

#Generates all tables in which data is going to be combined
ParaBiopSum <- data.frame();
PV1BiopRes <- c();
PV1Biop <- c();
PV2BiopRes <- c();
PV2Biop <- c();

for (tt in 1:csvCount) {
  ParaBiopSum <- rbind(ParaBiopSum, summary(ParaBiop[,tt])); #Generates summary statistics for parasite / biopsy
  
  SWT <- shapiro.test(ParaBiop[,tt]+1); #Shapiro-Wilks test of normal data distribution
  
  pVal <- as.numeric(format(round(SWT[["p.value"]], 4), nsmall = 4));
  
  if (pVal < 0.05) {
    PV1BiopRes <- c(PV1BiopRes, "Non-Normal");
    if (pVal < 0.0001) {
      PV1Biop <- c(PV1Biop, "<0.0001");
    } else {
      PV1Biop <- c(PV1Biop, pVal);
    }
  } else {
    PV1BiopRes <- c(PV1BiopRes, "Normal");
    PV1Biop <- c(PV1Biop, pVal);
  }
  SWT <- shapiro.test(log10(ParaBiop[,tt]+1)); #Shapiro-Wilks test of log-transformed distribution
  
  pVal <- as.numeric(format(round(SWT[["p.value"]], 4), nsmall = 4));
  
  if (pVal < 0.05) {
    PV2BiopRes <- c(PV2BiopRes, "Non-Normal");
    if (pVal < 0.0001) {
      PV2Biop <- c(PV2Biop, "<0.0001");
    } else {
      PV2Biop <- c(PV2Biop, pVal);
    }
  } else {
    PV2BiopRes <- c(PV2BiopRes, "Normal");
    PV2Biop <- c(PV2Biop, pVal);
  }
}
#Summary statistics for the parasites / biopsy data
ColSumTit <- head(names(summary(ParaBiop[,tt])), -1);
colnames(ParaBiopSum) <- ColSumTit;

#Combines Shapiro-Wilks data for Linear and Log-transformed data
SWTBiop <- data.frame(PV1Biop, PV1BiopRes, PV2Biop, PV2BiopRes)
colnames(SWTBiop) <- c("P-Value (Linear)", "Data Distribution (Linear)", "P-Value (Log)", "Data Distribution (Log)")

#Parasites / complete skin based on Median and 0.6 mm biopsy diameter
ParaPerSkin <- ParaBiopSum$Median * (TotalSkinArea/ConvSca3) / ((0.6/2)^2 * pi);
#Parasites / cm^2 based on Median and 0.6 mm biopsy diameter
ParaPercm2 <- ParaBiopSum$Median / ((0.6/2)^2 * pi);

rm(ParaBiop, BiopCount, tt, SWT, pVal, ColSumTit);

# All data save---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#This creates a data frame summarizing all the relevant skin patch data for all samples in a table
AllData <- data.frame(as.numeric(format(round(AllPatchCounts / (TotalSkinArea / ConvSca3), 1), nsmall = 1)), AllPatchCounts, 
                      AreaSummary %>% mutate(across(is.numeric, ~ round(., 1))), PV1Results, PV1Call, PV2Results, PV2Call,
                      as.numeric(format(round(AllSkewArea, 3), nsmall = 3)), AllSkewness,  as.numeric(format(round(AllSES, 3), nsmall = 3)),
                      as.numeric(format(round(AllTS, 3), nsmall = 3)), AllIoP, MinDiaSummary %>% mutate(across(is.numeric, ~ round(., 1))),
                      MaxDiaSummary %>% mutate(across(is.numeric, ~ round(., 1))), as.numeric(format(round(TotalPatchArea / ConvSca3, 3), nsmall = 3)), 
                      as.numeric(format(round(TotalSkinArea / ConvSca3, 3), nsmall = 3)), as.numeric(format(round(TotalPatchArea * 100 / TotalSkinArea, 3), nsmall = 3)),
                      as.numeric(format(round(AllSkewDens, 3), nsmall = 3)), AllDensSkew, as.numeric(format(round(AllDensSES, 3), nsmall = 3)),
                      as.numeric(format(round(AllDensTS, 3), nsmall = 3)), AllDensIoP, ParaBiopSum %>% mutate(across(is.numeric, ~ round(., 1))),
                      SWTBiop, as.numeric(format(round(ParaPercm2, 1), nsmall = 1)), ParaPerSkin);



AllDataCol <- c("Parasite Patches / cm^2", "Total Parasite Patch Counts", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.",
                "Area P-Value (Linear)", "Area Data Distribution (Linear)", "Area P-Value (Log)", "Area Data Distribution (Log)",
                "Area Data Skewness", "Area Skewness (Sample)", "S.E. of Area Skewness", "Area Test Statistic", "Area Skewness (Population)", 
                "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.",
                "Total Patch Area (cm^2)", "Total Skin Area (cm^2)", "Patch Ratio", "Density Data Skewness", "Density Skewness (Sample)",
                "S.E. of Density Skewness", "Density Test Statistic", "Density Skewness (Population)", "Min.", "1st Qu.", "Median", 
                "Mean", "3rd Qu.", "Max.", "Parasites P-Value (Linear)", "Parasites Data Distribution (Linear)", "Parasites P-Value (Log)",
                "Data Distribution (Log)",  "Parasites / cm^2", "Parasites / Skin");

colnames(AllData) <- AllDataCol;

rownames(AllData) <- MouseIDs;

#Saves the output data compilation file
write.csv(AllData, "All Data.csv");

# 3.2.	Spatial point pattern analysis==================================================================================================================================================================================================================================================================================

#Calling functions stored in other scripts required for in this part of the script
BiblioPath5 <- "E:/Work/Ongoing Work/Joe (Work for Paul Kaye)/Written Article/4 - Paper with Paul on skin parasite distribution/Data/Stereo/Skin Patch Analysis in R/Spatial Point Pattern";
# BiblioPath5 <- "D:/Lab Work/Ongoing Work/Joe (Work for Paul Kaye)/Written Article/4 - Paper with Paul on skin parasite distribution/Data/Stereo/Skin Patch Analysis in R/Spatial Point Pattern";

source(paste0(BiblioPath5, "/", "ImageSave.R"));
source(paste0(BiblioPath5, "/", "GraphSave.R"));

#Spatial Point Pattern Analysis with spatstat (based on Spatial Point Patterns Methodology and Applications in R, A. Baddeley, E. Rubak, R. Turner, 2016, 1st edition)

if (all(DupCheck1 == TRUE) == TRUE) { #This only allows the SPP to start if all patch coordinates are duplicate free

# 3.2.1.	ppp format generation-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  #Required to get the width and height of the original image from Table2
  BiblioPath3 <- "E:/Work/Ongoing Work/Joe (Work for Paul Kaye)/Written Article/4 - Paper with Paul on skin parasite distribution/Data/Stereo/Data/OutputPolygon";
  # BiblioPath3 <- "D:/Lab Work/Ongoing Work/Joe (Work for Paul Kaye)/Written Article/4 - Paper with Paul on skin parasite distribution/Data/Stereo/Data/OutputPolygon"; #Alternative data location
  
  #This will build the ppp variable that can be analysed in spatstat
  DimTable <- read.csv(paste0(BiblioPath3,"/Table2.csv")); #Opens data sheets, NOTE: paste0("Name", j) is the same as paste("Name", j, sep="") as both omit the space between "Name" and j
  ImageWidth <- DimTable[, ImWid];
  ImageHeight <- DimTable[, ImHei];
  ImageSize <- cbind("Width (µm)" = ImageWidth, "Height (µm)" = ImageHeight);
  rownames(ImageSize) <- MouseIDs;
  
  #Generates all tables in which data is going to be combined
  SPlist <- list();
  pppDup <- c();
  listXY <- list();
  listXYarea <- list();
  listXYdensity <- list();
  listXYpeak <- list();
  listXYmarks <- list();
  
  #3.2.1.1.	Skin polygon window generation
  for (cc in 1:csvCount) {
    #Defines a polygon that describes a rough outline of the mouse skin; all X and Y coordinates are comprised in a list
    DimTable <- read.csv(paste0(BiblioPath3,"/Rag",cc,"_SkinOutline.csv")); #Opens data sheets, NOTE: paste0("Name", j) is the same as paste("Name", j, sep="") as both omit the space between "Name" and j
    
    Xpoly <- rev(DimTable$X);
    Ypoly <- rev(DimTable$Y);
    
    if (cc ==1 | cc ==2 | cc ==3 | cc == 22 | cc ==23) { #These 5 images are half the number of pixels in X and Y dimension than the rest and therefore, require a different converter
      Xpoly <- Xpoly * ConvSma;
      Ypoly <- Ypoly * ConvSma;
    } else {
      Xpoly <- Xpoly * ConvLar;
      Ypoly <- Ypoly * ConvLar;
    }
    SkinPoly <- data.frame(x = Xpoly, y = Ypoly);
    
    #This generate a skin outline, which can be superimposed on skin image graphs later
    SP <- ppp(Xpoly, Ypoly, xrange = c(0, ImageWidth[cc]), yrange = c(0, ImageHeight[cc])); #includes only the patch area mark
    SP <- rescale.ppp(SP, ConvSca1, "mm"); #this adjusted the unit size of the coordinates only from µm to mm
    SPlist <- list.append(SPlist, SP);
    
    #Defines a polygon that describes a rough outline of the mouse skin; all X and Y coordinates are comprised in a list
    DimTable <- read.csv(paste0(BiblioPath3,"/Rag",cc,"_REarOutline.csv")); #Opens data sheets, NOTE: paste0("Name", j) is the same as paste("Name", j, sep="") as both omit the space between "Name" and j
    
    REarlength <- length(DimTable[["X"]]);
    if (REarlength != 0) {
      Xpoly <- DimTable$X;
      Ypoly <- DimTable$Y;
      
      if (cc ==1 | cc ==2 | cc ==3 | cc == 22 | cc ==23) { #These 5 images are half the number of pixels in X and Y dimension than the rest and therefore, require a different converter
        Xpoly <- Xpoly * ConvSma;
        Ypoly <- Ypoly * ConvSma;
      } else {
        Xpoly <- Xpoly * ConvLar;
        Ypoly <- Ypoly * ConvLar;
      }
      REarPoly <- data.frame(x = Xpoly, y = Ypoly);
    }
    
    #Defines a polygon that describes a rough outline of the mouse skin; all X and Y coordinates are comprised in a list
    DimTable <- read.csv(paste0(BiblioPath3,"/Rag",cc,"_LEarOutline.csv")); #Opens data sheets, NOTE: paste0("Name", j) is the same as paste("Name", j, sep="") as both omit the space between "Name" and j
    
    LEarlength <- length(DimTable[["X"]]);
    if (LEarlength != 0) {
      Xpoly <- DimTable$X;
      Ypoly <- DimTable$Y;
      
      if (cc ==1 | cc ==2 | cc ==3 | cc == 22 | cc ==23) { #These 5 images are half the number of pixels in X and Y dimension than the rest and therefore, require a different converter
        Xpoly <- Xpoly * ConvSma;
        Ypoly <- Ypoly * ConvSma;
      } else {
        Xpoly <- Xpoly * ConvLar;
        Ypoly <- Ypoly * ConvLar;
      }
      LEarPoly <- data.frame(x = Xpoly, y = Ypoly);
    }
    
    
    if (REarlength != 0 & LEarlength != 0) { #In case ear holes are missing
      SkinAllBoarders <- list(SkinPoly, REarPoly, LEarPoly);
    } else if (REarlength != 0 & LEarlength == 0) {
      SkinAllBoarders <- list(SkinPoly, REarPoly);
    } else if (REarlength == 0 & LEarlength != 0) {
      SkinAllBoarders <- list(SkinPoly, LEarPoly);
    } else {
      SkinAllBoarders <- list(SkinPoly);
    }
    
    #3.2.1.2.	Create ppp
    XY <- ppp(XMYMList[[paste0("RAG", cc)]]$X, XMYMList[[paste0("RAG", cc)]]$Y, poly = SkinAllBoarders, xrange = c(0, ImageSize[paste0("RAG", cc), 1]), yrange = c(0, ImageSize[paste0("RAG", cc), 2])); #includes only the patch area mark
    XY <- rescale(XY, ConvSca1, "mm"); #This adjusted the unit size of the coordinates only from µm to mm
    listXY <- list.append(listXY, XY);
    
    XYarea <- XY;
    marks(XYarea) <- (AreasList[[paste0("RAG", cc)]] / ConvSca2);
    listXYarea <- list.append(listXYarea, XYarea);
    
    XYdensity <- XY;
    marks(XYdensity) <- NormDensList[[paste0("RAG", cc)]];
    listXYdensity <- list.append(listXYdensity, XYdensity);
    
    XYpeaks <- XY;
    marks(XYpeaks) <- PeaksList[[paste0("RAG", cc)]];
    listXYpeak <- list.append(listXYpeak, XYpeaks);
    
    XYmarks <- XY;
    marks(XYmarks) <- MarksList[[paste0("RAG", cc)]];
    listXYmarks <- list.append(listXYmarks, XYmarks);
    
    #3.2.1.3.	Check for duplicate points in ppp
    DUP <- all(duplicated.ppp(XY) == FALSE); #Checks for duplicate point in point pattern
    if (DUP == TRUE) {
      DUP <- "No Doubles";
    } else if (DUP == FALSE) {
      DUP <- "Doubles";
    }
    pppDup <- c(pppDup, DUP); #Collects Output together
    
    #3.2.1.4.	Plot & image generation
    
    #To export the resulting plot
    Image.Save(XY, cc, 1, 1, "XY", "Patches by Center of Mass", cols = "red");
    # Image.Save(XYarea, cc, 1, 1, "XYarea", "Patches by Area", cols = "red", maxsize = 10);
    # Image.Save(XYdensity, cc, 1, 1, "XYdensity", "Patches by Density", cols = "red", maxsize = 10);
    # Image.Save(XYpeaks, cc, 1, 1, "XYpeak", "Patches by Peaks", cols = "red", maxsize = 10);
  }
  
  #Adds names to list and converts them to solists for grouped plotting
  names(SPlist) <- MouseIDs;
  names(listXY) <- MouseIDs;
  listXY <- as.solist(listXY);
  names(listXYarea) <- MouseIDs;
  listXYarea <- as.solist(listXYarea);
  names(listXYdensity) <- MouseIDs;
  listXYdensity <- as.solist(listXYdensity);
  names(listXYpeak) <- MouseIDs;
  listXYpeak <- as.solist(listXYpeak);
  
  rm(DimTable, SP, Xpoly, Ypoly, SkinPoly, REarPoly, LEarPoly, REarlength, LEarlength, SkinAllBoarders, XY, XYarea, XYdensity, XYpeaks, XYmarks, DUP, cc);

# PPP Data--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  #This creates a data frame summarizing all the relevant skin patch data for all samples in a table
  PPPData <- data.frame("Double Coordinates Check" = pppDup);
  
  rm(pppDup);
  
# 3.2.2. Quadrat test---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  #3.2.2. Quadrat test
  
  #Check whether the point process has a homogeneous (null hypothesis) or heterogeneous intensity
  
  #Generates all tables in which data is going to be combined
  TestRun <- c();
  QTlist <- list();
  QTp.values <- c();
  QTregcluList <- list();
  QTregP <- c();
  QTcluP <- c();
  FTlist <- list();
  FTp.values <- c();
  QTid <- c();
  HomIntRes <- c();
  HIRes <- c();
  RMSp.values <- c();
  QCp.values <- c();
  PoissonTest <- c();
  
  for (dd in 1:csvCount) {
    #Quadrat Test as described for spatstat
    # suppressWarnings(QT <- quadrat.test(listXY[[paste0("RAG", dd)]], nx = GridX, ny = GridY, alternative = "two.sided", CR = 1)); #The test will throw a warning, if to many quadrants with <5 patches or 0 patches are found
    QT <- quadrat.test(listXY[[paste0("RAG", dd)]], nx = GridX, ny = GridY, alternative = "two.sided", method = "MonteCarlo", nsim = 999); #If quadrat test has a p-value of <0.05, it diverges from a homogeneous distribution. Here we check if the alternative is "regular".
    # QTlist <- list.append(QTlist, QT); #Collects all the quadrat.test output into a single list
    QTp.values <- c(QTp.values, round(QT$p.value, 5)); #Collects all the p-values into a vector
    
    if (QTp.values[dd] < 0.05) {
      QTreg <- quadrat.test(listXY[[paste0("RAG", dd)]], nx = GridX, ny = GridY, alternative = "regular", method = "MonteCarlo", nsim = 999); #If quadrat test has a p-value of <0.05, it diverges from a homogeneous distribution. Here we check if the alternative is "regular".
      QTclu <- quadrat.test(listXY[[paste0("RAG", dd)]], nx = GridX, ny = GridY, alternative = "clustered", method = "MonteCarlo", nsim = 999); #If quadrat test has a p-value of <0.05, it diverges from a homogeneous distribution. Here we check if the alternative is "clusstered".
      QTregcluList <- list.append(QTregcluList, list("Alt. Regular" = QTreg, "Alt. Clustered" = QTclu));
      QTregP <- c(QTregP, round(QTreg$p.value, 5));
      QTcluP <- c(QTcluP, round(QTclu$p.value, 5));
      
    } else {
      QTreg = " ";
      QTclu = " ";
      QTregcluList <- list.append(QTregcluList, list("Alt. Regular" = QTreg, "Alt. Clustered" = QTclu));
      QTregP <- c(QTregP, QTreg);
      QTcluP <- c(QTcluP, QTclu);
    }
  
    #Check is counts per square are >=5 to ensure Chi-Square is valid
    QTobs <- QT$observed; #This extracts the observed values from the quadrat.test above
    TestVal <- all(QTobs >= 5);
    
    if (TestVal == FALSE) {
      #Fisher's Exact Test - alternative to the chi-square approach above, if there are many zeros in the observed data
      TestRun <- c(TestRun, "Fisher's Exact");
    
      #If squares were not counted because they were outside the window, they will be added back here with "0" counts so a GridY x GridX sized matrixed can be foremed from the vector
      if (length(QT$observed) < GridX * GridY) {
        #Quadrat Counting for Intensity Estimation
        QC <- quadratcount(listXY[[paste0("RAG", dd)]], nx = GridX, ny= GridY);
        
        #This will add the missing grid quadrants to the table so it can be converted to a matrix
        rowCou <- 0; #This is the start for the row counter
        QCdf <- as.data.frame(QC, stringsAsFactors = FALSE); #Table needs to be converted into a data frame without factorization of the string data
        
        if (nrow(QCdf) != GridX * GridY) {
          
          for (yy in 1:GridY) { #Indicates the numbe rof rows in the grid from the quadrat count
            for (xx in 1:GridX){ #Indicates the numbe rof Columns in the grid from the quadrat count
              
              rowInd <- paste0("Tile row ", yy,", col ", xx); #The tile column contains the grid quadrant locations
              
              rowNum <- which(QCdf$tile == rowInd); #This identifies the row by number
              
              vecLeng <- length(rowNum); #This measures the length of the vector
              
              rowCou <- rowCou + 1; #Increases the row count from the previous line by 1 to indectate where the missing row needs to be inserted
              
              if (vecLeng == 0) {
                
                InsertedRow <- data.frame("tile" = paste0("Tile row ", yy, ", col ", xx), "Freq" = "0", stringsAsFactors = FALSE); #Creates the row to be inserted
                
                insertRow <- function(QC, InsertedRow, r) { #This builds the function for inserting a row in data frame
                  QCdf[seq(r+1,nrow(QC)+1),] <- QC[seq(r,nrow(QC)),]
                  QCdf[r,] <- InsertedRow
                  QCdf
                }
                QCdf <- insertRow(QCdf, InsertedRow, rowCou); #Inserts the row
              }
              rowCou <- which(QCdf$tile == rowInd); #determine at which row we are after a row may have been inserted
            }
          }
          QCdf <- head(QCdf, rowCou); #removes an artifactual row from the end of the data frame
          QTobs <- matrix(as.numeric(QCdf$Freq), GridY, GridX, byrow = TRUE);
        }
        rm(rowNum, vecLeng, rowInd, InsertedRow);
      } else if (length(QT$observed) == GridX * GridY) { #In case the vector length of QT$observed corresponds to GridY * GridX
        QTobs <- matrix(QT$observed, GridY, GridX, byrow = TRUE);
      }
      
      
      #Fisher's exact test
      FT <- fisher.test(QTobs, hybrid = TRUE, simulate.p.value = TRUE, B = 1000); #Fisher's exact test; simulated P-values by Monte Carlo
      FTlist <- list.append(FTlist, FT); #Collects all the quadrat.test output into a single list
      FTp.values <- c(FTp.values, round(FT$p.value, 5)); #Collects all p-values into one vector for storage
      QTid <- cbind(QTid, paste0("RAG", dd));
      
      #Check P-values
      HomInt <- FT$p.value < 0.05; #Checks if p-values are significant or not
      
      if (HomInt == TRUE) {
        HI <- "Inhomogeneous";
      } else if(HomInt == FALSE) {
        HI <- "Homogeneous";
      }
      
    } else if (TestVal == TRUE) {
      TestRun <- c(TestRun, "Chi-Square");
      FTp.values <- c(FTp.values, " "); #in case Fisher's Exact Test was not run, the vector gets filled with nothing
      
      HomInt <- QT$p.value < 0.05; #Checks if p-values are significant or not
      if (HomInt == TRUE) {
        HI <- "Inhomogeneous";
      } else if (HomInt == FALSE) {
        HI <- "Homogeneous";
      }
    }
    
    HomIntRes <- c(HomIntRes, HomInt);
    HIRes <- c(HIRes, HI); #Creates vector with results
    
    #To export the resulting plot
    tiff(paste0("QT",dd,".tif"), width = ImageWidth[dd]/ConvSca1, height = ImageHeight[dd]/ConvSca1, units = "mm", res = ImageRes, pointsize = 6, compression = 'lzw');
    plot(listXY[[paste0("RAG", dd)]], ylim = c(ImageHeight[dd]/ConvSca1, 0), maxsize = 25, cols = "red", main = paste0("RAG", dd, " - Homogeneity Test"));
    plot(QT, ylim = c(ImageHeight[dd]/ConvSca1, 0), add = TRUE);
    dev.off();
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #Root-Mean-Square Goodness-Of-Fit Test as alternative to Chi-Square and Fisher's exact test
   
    #this is the root-mean-square approach with Monte Carlo for p-value simulation; this test suggests all skin a homogeneously distributed with patches
    QTobs <- QT$observed;
    QTexp <- round(QT$expected, 0);
    QTout <- rms.pval(QTobs, QTexp, num_sim = 1000);
    RMSp.values <- round(c(RMSp.values, QTout), 5);
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #Quadrat Counting test of Poisson distribution
    
    #This check for Poisson distribution assumes a homogeneous intensity
    qX <- quadratcount(listXY[[paste0("RAG", dd)]], nx = 16, ny = 20); #Quadrat counting with tiny grid
    tX <- factor(as.numeric(qX), levels = 0:max(qX));
    tX <- mergeLevels(tX, ">=5" = 6:length(levels(tX))-1); #Collapse rare categories into one
    tX <- table(tX); #build table
    dp <- dpois(0:4, mean(qX));
    eX <- AllData$AllPatchCounts[dd] * c(dp, 1 - sum(dp));
    X2 <- sum((tX - eX)^2/eX);
    QCp <- round(pchisq(X2, df = length(tX) - 1, lower.tail = FALSE), 5); #P-value generation
    if (QCp < 0.05) {
      TestPoi <- "Not Poisson";
    } else if (QCp > 0.05) {
      TestPoi <- "Poisson";
    }
    
    QCp.values <- c(QCp.values, QCp);
    PoissonTest <- c(PoissonTest, TestPoi);
  }
  
  homTestRes <- data.frame("Quadrat Test" = QTp.values, "Alt. Regular" = QTregP, "Alt. Clustered" = QTcluP, "Fisher's Exact Test" = FTp.values, "Accepted Test" = TestRun, "Intensity Type" = HIRes, "Root-Mean-Square GOF Test" = RMSp.values, "Quadrat Count (P-values)" = QCp.values, "Quadrat Counts for Poisson distribution" = PoissonTest);
  names(FTlist) <- QTid;
  names(QTlist) <- MouseIDs;
  QTlist <- as.anylist(QTlist);
  suppressWarnings(QTpool <- pool.quadrattest(QTlist)); #Give the p-value for the combined data, but unreliable if inhomogeneous intensities are among the samples
  
  rm(QT, QTreg, QTclu, QTobs, QC, QCdf, rowCou, TestVal, FT, HomInt, HI, QTid, QTp.values, FTp.values, TestRun, HIRes, qX, tX, dp, eX, X2, QCp, dd);
  
# PPP Data--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  #This creates a data frame summarizing all the relevant skin patch data for all samples in a table
  PPPData <- data.frame(PPPData, homTestRes);
  
  rm(homTestRes);
  
# 3.2.3. Point pattern intensity estimation-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  #Generates all tables in which data is going to be combined
  exeIntTest <- c();
  XYintCount <- c();
  XYintCouSD <- c();
  XYintArea <- c();
  XYintArSD <- c();
  XYintDens <- c();
  XYintDenSD <- c();
  XYintPeak <- c();
  XYintPeakSD <- c();
  IQCdfs <- data.frame(matrix( , nrow = GridX * GridY, ncol = 1));
  QCid <- c();
  IQCdfList <- list();
  
  for (ee in 1:csvCount) {
    if (HomIntRes[ee] == FALSE) { #If homogeneous (p>0.05)
      
      #Point Pattern Intensity Estimation
      exeIntTest <- c(exeIntTest, "Intensity");
      
      #This one gets an intensity estimate based on the point counts
      XYiCo<- round(intensity(listXY[[paste0("RAG", ee)]]), 4); #This calculation of intensity assumed homogeneity of points
      XYintCount <- c(XYintCount, XYiCo);
      XYiCoSD <- round(sqrt(XYiCo/(TotalSkinArea[ee]/ConvSca2)), 4); #The calculation of a SD of the intensity also assumes a Poisson point process
      XYintCouSD <- c(XYintCouSD, XYiCoSD);
      rm(XYiCo, XYiCoSD);
      
      #This one gets and intensity estimate based on the Patch Area mark
      XYiAr <- round(intensity(listXY[[paste0("RAG", ee)]], weights = AreasList[[paste0("RAG", ee)]] / ConvSca2), 4); #This calculation of intensity assumed homogeniety of points
      XYintArea <- c(XYintArea, XYiAr);
      XYiArSD <- round(sqrt(XYiAr/(TotalSkinArea[ee]/ConvSca2)), 4); #The calculation of a SD of the intensity also assumes a Poisson point process
      XYintArSD <- c(XYintArSD, XYiArSD);
      rm(XYiAr, XYiArSD);
      
      #This one gets and intensity estimate based on the Patch adjusted Density mark
      XYiDe <- round(intensity(listXY[[paste0("RAG", ee)]], weights = NormDensList[[paste0("RAG", ee)]]), 4); #This calculation of intensity assumed homogeniety of points
      XYintDens <- c(XYintDens, XYiDe);
      XYiDenSD <- round(sqrt(XYiDe/(TotalSkinArea[ee]/ConvSca2)), 4); #The calculation of a SD of the intensity also assumes a Poisson point process
      XYintDenSD <- c(XYintDenSD, XYiDenSD);
      rm(XYiDe, XYiDenSD);
      
      #This one gets and intensity estimate based on the Patch Peak mark
      XYiPe <- round(intensity(listXY[[paste0("RAG", ee)]], weights = PeaksList[[paste0("RAG", ee)]]), 4); #This calculation of intensity assumed homogeniety of points
      XYintPeak <- c(XYintPeak, XYiPe);
      XYiPeSD <- round(sqrt(XYiPe/(TotalSkinArea[ee]/ConvSca2)), 4); #The calculation of a SD of the intensity also assumes a Poisson point process
      XYintPeakSD <- c(XYintPeakSD, XYiPeSD);
      rm(XYiPe, XYiPeSD);
      
    } else if (HomIntRes[ee] == TRUE) { #If inhomogeneous (P<0.05)
      
      exeIntTest <- c(exeIntTest, "Quadrat Count"); #Creates the label that informs which test was performed on what sample
      
      #Quadrat Counting for Intensity Estimation
      QC <- quadratcount(listXY[[paste0("RAG", ee)]], nx = GridX, ny= GridY);
      
      #To export the resulting plot
      IQC <- intensity(QC, image = TRUE); #Calculates the intensity for each quadrant of the grid and wraps it into a plotable format
      tiff(paste0("IQC",ee,".tif"), width = ImageWidth[ee]/ConvSca1, height = ImageHeight[ee]/ConvSca1, units = "mm", res = ImageRes, pointsize = 6, compression = 'lzw');
      plot(IQC, ylim = c(ImageHeight[ee]/ConvSca1, 0), main = paste0("Estimated Intensity - RAG", ee));
      dev.off();
      
      IQC <- intensity(QC); #Calculates the intensity for each quadrant of the grid
      
      #This will add the missing grid quadrants to the table so it can be converted to a matrix
      rowCou <- 0; #This is the start for the row counter
      IQCdf <- as.data.frame(IQC, stringsAsFactors = FALSE); #Table needs to be converted into a data frame without factorization of the string data
      
      if (nrow(IQCdf) != GridX * GridY) {
        for (yy in 1:GridY) { #Indicates the numbe rof rows in the grid from the quadrat count
          for (xx in 1:GridX){ #Indicates the numbe rof Columns in the grid from the quadrat count
            
            rowInd <- paste0("Tile row ", yy,", col ", xx); #The tile column contains the grid quadrant locations
            
            zzz <- which(IQCdf$tile == rowInd); #This identifies the row by number
            
            vecLeng <- length(zzz); #This measures the length of the vector
            
            rowCou <- rowCou + 1; #Increases the row count from the previous line by 1 to indectate where the missing row needs to be inserted
            
            if (vecLeng == 0) {
              
              InsertedRow <- data.frame("tile" = paste0("Tile row ", yy, ", col ", xx), "Freq" = "0", stringsAsFactors = FALSE); #Creates the row to be inserted
              
              insertRow <- function(IQC, InsertedRow, r) { #This builds the function for inserting a row in data frame
                IQCdf[seq(r+1,nrow(IQC)+1),] <- IQC[seq(r,nrow(IQC)),]
                IQCdf[r,] <- InsertedRow
                IQCdf
              }
              IQCdf <- insertRow(IQCdf, InsertedRow, rowCou); #Inserts the row
            }
            rowCou <- which(IQCdf$tile == rowInd); #determine at which row we are after a row may have been inserted
          }
        }
        IQCdf <- head(IQCdf, rowCou); #removes an artifactual row from the end of the data frame
        rm(zzz, vecLeng, rowInd, InsertedRow);
      }
      
      IQCdfs <- data.frame(IQCdfs, round(as.numeric(IQCdf$Freq), 5));
      IQCdfR <- IQCdf$tile;
      QCid <- cbind(QCid, paste0("RAG", ee));
      
      IQCdf <- matrix(round(as.numeric(IQCdf$Freq), 5), GridY, GridX); #Converts the data frame into a GridX by GridY matrix
      IQCdfList <- list.append(IQCdfList, IQCdf);
      
      MeanI <- mean(IQC); #Average intensity across all quadrats
      XYintCount <- c(XYintCount, MeanI);
      IQCse <- sqrt(var(as.numeric(IQC)) / (length(IQC)-1)); #Standard error estimation
      XYintCouSD <- c(XYintCouSD, IQCse);
      
      XYintArea <- c(XYintArea, "");
      XYintArSD <- c(XYintArSD, "");
      XYintDens <- c(XYintDens, "");
      XYintDenSD <- c(XYintDenSD, "");
      XYintPeak <- c(XYintPeak, "");
      XYintPeakSD <- c(XYintPeakSD, "");
    }
  }
  
  QCid <- factor(QCid);
  IQCdfs <- IQCdfs[,-1];
  colnames(IQCdfs) <- QCid;
  rownames(IQCdfs) <- IQCdfR;
  names(IQCdfList) <- QCid;

  rm(QC, IQC, MeanI, IQCse, rowCou, IQCdf, IQCdfR, QCid, ee, xx, yy);
  
# PPP Data--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  #This creates a data frame summarizing all the relevant skin patch data for all samples in a table
  PPPData <- data.frame(PPPData, "Test Applied" = exeIntTest, "Est. # of Patch / mm^2 (Intensity)" = XYintCount, "S.E. of Intensity" = XYintCouSD, "Est. Patch Area / mm^2" = XYintArea, 
             "S.E. of Intensity" = XYintArSD, "Est. Density / mm^2" = XYintDens, "S.E. of Intensity" = XYintDenSD, "Est. Patch Peaks / mm^2" = XYintPeak, "S.E. of Intensity" = XYintPeakSD);
  
  rm(XYintCount, XYintCouSD, XYintArea, XYintArSD, XYintDens, XYintDenSD, XYintPeak, XYintPeakSD);

# 3.2.4. Intensity estimation by Kernal smoothing-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  #More accurate when heterogeneity is suspected
  
  #Generates all tables in which data is going to be combined
  DDlist <- list();
  DDareaList <- list();
  DDdensList <- list();
  DDpeakList <- list();
  DDmeanInt <- c();
  DDarMeanInt <- c();
  DDdenMeanInt <- c();
  DDpeakMeanInt <- c();
  DDmeanIntSE <- c();
  DDarMeanIntSE <- c();
  DDdenMeanIntSE <- c();
  DDpeakMeanIntSE <- c();
  
  for (ff in 1:csvCount) { 
    
    #Creates density image
    suppressWarnings(DD <- density(listXY[[paste0("RAG", ff)]], sigma = bw.ppl, positive=TRUE, se=TRUE)); #bw.ppl assume Poison process, while bw.diggle assume Cox process that is more clustered than Poison
    DDlist <- list.append(DDlist, DD);
    DDmean <- mean(DD$estimate);
    DDmeanInt <- c(DDmeanInt, DDmean);
    DDmeanSE <- mean(DD$SE);
    DDmeanIntSE <- c(DDmeanIntSE, DDmeanSE);
  
    suppressWarnings(DDarea <- density(listXY[[paste0("RAG", ff)]], weights = AreasList[[paste0("RAG", ff)]] / ConvSca2, sigma = bw.ppl, positive=TRUE, se=TRUE)); #bw.ppl assume Poison process, while bw.diggle assume Cox process that is more clustered than Poison
    DDareaList <- list.append(DDareaList, DDarea);
    DDmean <- mean(DDarea$estimate);
    DDarMeanInt <- c(DDarMeanInt, DDmean);
    DDmeanSE <- mean(DDarea$SE);
    DDarMeanIntSE <- c(DDarMeanIntSE, DDmeanSE);
    
    suppressWarnings(DDdens <- density(listXY[[paste0("RAG", ff)]], weights = DensityList[[paste0("RAG", ff)]], sigma = bw.ppl, positive=TRUE, se=TRUE)); #bw.ppl assume Poison process, while bw.diggle assume Cox process that is more clustered than Poison
    DDdensList <- list.append(DDdensList, DDdens);
    DDmean <- mean(DDdens$estimate);
    DDdenMeanInt <- c(DDdenMeanInt, DDmean);
    DDmeanSE <- mean(DDdens$SE);
    DDdenMeanIntSE <- c(DDdenMeanIntSE, DDmeanSE);

    suppressWarnings(DDpeak <- density(listXY[[paste0("RAG", ff)]], weights = PeaksList[[paste0("RAG", ff)]], sigma = bw.ppl, positive=TRUE, se=TRUE)); #bw.ppl assume Poison process, while bw.diggle assume Cox process that is more clustered than Poison
    DDpeakList <- list.append(DDpeakList, DDpeak);
    DDmean <- mean(DDpeak$estimate);
    DDpeakMeanInt <- c(DDpeakMeanInt, DDmean);
    DDmeanSE <- mean(DDpeak$SE);
    DDpeakMeanIntSE <- c(DDpeakMeanIntSE, DDmeanSE);
    
    #To export the resulting plot
    Image.Save(DD$estimate, ff, 1, 1, "DD", "Estimated Intensity");
    Image.Save(DD$SE, ff, 1, 1, "DDse", "S.E. of Est. Intensity");
    Image.Save(DDarea$estimate, ff, 1, 1, "DDarea", "Estimated Intensity");
    Image.Save(DDarea$SE, ff, 1, 1, "DDareaSE", "S.E. of Est. Intensity");
    Image.Save(DDdens$estimate, ff, 1, 1, "DDdensity", "Estimated Intensity");
    Image.Save(DDdens$SE, ff, 1, 1, "DDdensSE", "S.E. of Est. Intensity");
    Image.Save(DDpeak$estimate, ff, 1, 1, "DDpeaks", "Estimated Intensity");
    Image.Save(DDpeak$SE, ff, 1, 1, "DDpeakSE", "S.E. of Est. Intensity");
  }
  
  names(DDlist) <- MouseIDs;
  names(DDareaList) <- MouseIDs;
  names(DDdensList) <- MouseIDs;
  names(DDpeakList) <- MouseIDs;
  
  rm(DD, DDarea, DDdens, DDpeak, DDmean, DDmeanSE);

# PPP Data----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  #This creates a data frame summarizing all the relevant skin patch data for all samples in a table
  PPPData <- data.frame(PPPData, "Est. # of Patch / mm^2 (Intensity)" = DDmeanInt, "S.E. of Intensity" = DDmeanIntSE, "Est. Patch Area / mm^2" = DDarMeanInt, "S.E. of Intensity" = DDarMeanIntSE, 
                        "Est. Density / mm^2" = DDdenMeanInt, "S.E. of Intensity" = DDdenMeanIntSE, "Est. Patch Peaks / mm^2" = DDpeakMeanInt, "S.E. of Intensity" = DDpeakMeanIntSE);
  
  rm(DDmeanInt, DDmeanIntSE, DDarMeanInt, DDarMeanIntSE, DDdenMeanInt, DDdenMeanIntSE, DDpeakMeanInt, DDpeakMeanIntSE);
  
# 3.2.5.	Hot spot analysis------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  #this will check for hot spots and test whether they are significant via Monte Carlo
  
  #Generates all tables in which data is going to be combined
  HSlist <- list();
  HSpvalList <- list();
  
  for (gg in 1:csvCount) {
    
    #looks for hot spots by looking for elevated intensity and calculates p-values for the area; presented in plot form
    HS <- scanLRTS(listXY[[paste0("RAG", gg)]], r = 2 * bw.ppl(listXY[[paste0("RAG", gg)]])); #scan test
    HSlist <- list.append(HSlist, HS);
    HSpval <- eval.im(pchisq(HS, df = 1, lower.tail = FALSE)); #calculates p-values
    HSpvalList <- list.append(HSpvalList, HSpval);
    
    #To export the resulting plot
    tiff(paste0("Hot Spots ", gg,".tif"), width = ImageWidth[gg]/ConvSca1, height = ImageHeight[gg]/ConvSca1, units = "mm", res = ImageRes, pointsize = 6, compression = 'lzw');
    plot(HS, ylim = c(ImageHeight[gg]/ConvSca1, 0), main = paste0("RAG", gg, " - Hot Spots"));
    lines(SPlist[[gg]], ylim = c(ImageHeight[gg]/ConvSca1, 0), col = "white", lwd = 3);
    dev.off();
    
    #To export the resulting plot
    tiff(paste0("HS p-values", gg,".tif"), width = ImageWidth[gg]/ConvSca1, height = ImageHeight[gg]/ConvSca1, units = "mm", res = ImageRes, pointsize = 6, compression = 'lzw');
    plot(HSpval < 0.01, ylim = c(ImageHeight[gg]/ConvSca1, 0), col = c("white", "red"), main = paste0("RAG", gg, " - Hot Spots (Significance)"));
    lines(SPlist[[gg]], ylim = c(ImageHeight[gg]/ConvSca1, 0), lwd = 3);
    dev.off();
    
    rm(HS, HSpval);
  }
  
  names(HSlist) <- MouseIDs;
  names(HSpvalList) <- MouseIDs;
  
# 3.2.6.	High/low intensity estimation------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  #this is to check for high and low point intensity
  
  #Generates all tables in which data is going to be combined
  NNlist <- list();
  
  #Looped high/low test
  for (hh in 1:csvCount) {
    NN <- nnclean(listXY[[paste0("RAG", hh)]], k = 1, plothist = TRUE, main = paste0("RAG", hh, " - Neasrest Neighbor Distance"));
    dev.copy(tiff, paste0("NN", hh, " - Cleaning Hist.tif"), res = GraphRes, width = 100, height = 100, units = "mm", pointsize = 12, compression = 'lzw');
    dev.off();
    
    names(NN[["marks"]]) <- c("Classification", "Fitted Probability");
    NNlist <- list.append(NNlist, NN);
    
    #To export the resulting plot
    tiff(paste0("NN", hh, " High-Low (Classificaiton).tif"), width = ImageWidth[hh]*1.25/ConvSca1, height = ImageHeight[hh]/ConvSca1, units = "mm", res = ImageRes, pointsize = 6, compression = 'lzw');
    plot(subset.ppp(NN, select = "Classification"), ylim = c(ImageHeight[hh]/ConvSca1, 0), cex = 3, cols = "red", main = paste0("RAG", hh, " - Classification"));
    dev.off();
    
    #To export the resulting plot
    tiff(paste0("NN", hh, " High-Low (Fitted Probability).tif"), width = ImageWidth[hh]*1.25/ConvSca1, height = ImageHeight[hh]/ConvSca1, units = "mm", res = ImageRes, pointsize = 6, compression = 'lzw');
    plot(subset.ppp(NN, select = "Fitted Probability"), ylim = c(ImageHeight[hh]/ConvSca1, 0), cols = "red", main = paste0("RAG", hh, " - Fitted Probability"));
    dev.off();
  }
  
  names(NNlist) <- MouseIDs;  

  rm(NN);
  
# 3.2.7.	Ripley's K, Besag's L, Empty-Space, Nearest-Neighbor, J, Pair Correlation functions---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  #This assess independence of data points (correlation) and determine if data points are regularly spaced, independent or clustered
  #This Ktest assumes the point process to be stationary
  
  #Calling functions stored in other scripts required for in this part of the script
  source("NormAxis.R");
  
  #Generates all tables in which data is going to be combined
  XYscaDF <- data.frame();
  XYinhDF <- data.frame();
  KtestList <- list();
  LtestList <- list();
  FtestList <- list();
  GtestList <- list();
  JtestList <- list();
  KstestList <- list();
  LstestList <- list();
  KpcfList <- list();
  KTSElist <- list();
  LTSElist <- list();
  
  for (ii in 1:csvCount) { #The Ripley's K function takes quite a while to execute on all 23 skins
    if (HomIntRes[ii] == TRUE) {  
      #Check if correlation-stationary is true; allows distinction between inhomogeneous and scaled distribution
      #Approach by Hahn and Jensen (Chapter 16 Section 8.5 page 691)
      XYscaled <- c();
      XYinhom <- c();
      for (zzz in 1:StatTest) {
        XYsel <- sort(sample(1:(StatX2*StatY2), (StatX2*StatY2/2)));
        XYy <- quantess(listXY[[paste0("RAG", ii)]], "y", StatY2);
        XYxy <- nestsplit(listXY[[paste0("RAG", ii)]], XYy, nx = StatX2);
        XYxy$count <- c(1:(StatX2*StatY2));
        XYxy$inten <- factor(XYxy$count %in% c(XYsel), labels = c("R", "L"));
        
        tryCatch(
          {studpermu.test(XYxy, pts ~ inten, summaryfunction = Kscaled, rinterval = c(0, 2)); EM <<- 0},
          error = function(e) {EM <<- 1});
        
        if (EM == 0) {
          XYsca <- studpermu.test(XYxy, pts ~ inten, summaryfunction = Kscaled, rinterval = c(0, 2));
          XYinh <- studpermu.test(XYxy, pts ~ inten, summaryfunction = Kinhom, lambda = DDlist[[paste0("RAG", ii)]]$estimate, rinterval = c(0, 2));
          
          XYscaled <- c(XYscaled, XYsca$p.value);
          XYinhom <- c(XYinhom, XYinh$p.value);
        } else if (EM == 1) {
          XYscaled <- c(XYscaled, NA);
          XYinhom <- c(XYinhom, NA);
        }
      }
      XYscaDF <- rbind(XYscaDF, XYscaled);
      XYinhDF <- rbind(XYscaDF, XYinhom);
  
      #Graphical approach (Chapter 16 Section 8.5 page 689)
      mx <- median(coords(listXY[[paste0("RAG", ii)]])$x);
      halves <- chop.tess(Window(listXY[[paste0("RAG", ii)]]), infline(v = mx));
      
      Ks <- anylapply(split(listXY[[paste0("RAG", ii)]], halves), Kscaled);
      names(Ks) <- c("left", "right");
      Ki <- anylapply(split(listXY[[paste0("RAG", ii)]], halves), Kinhom, lambda = DDlist[[paste0("RAG", ii)]]$estimate);
      names(Ki) <- c("left", "right");
      
      #To export the resulting graphs
      Graph.Save(Ks, ii, 2, "Ks", "Correlation-Stationary test (Kscaled)");
      Graph.Save(Ki, ii, 2, "Ki", "Correlation-Stationary test (Kinhom)");
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Inhomogeneous
      #~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      #Uses the inhomogeneous distribution as a null hypothesis
      #Ripley's K function - inhomogeneous (Chapter 7 Section 10 page 242)
      #Require the point process to be correlation-stationary
      #The intensity estimation given by the DD (class" im) gives a warning about possibly fixed ranges; alternatively, the function can estimate intensity by leave-one-out kernel smoother by not giving lambda
      Kitest <- Kinhom(listXY[[paste0("RAG", ii)]], lambda = DDlist[[paste0("RAG", ii)]]$estimate); #Ripley's K function for inhomogeneous intensity
      KranInh <- envelope.ppp(listXY[[paste0("RAG", ii)]], Kinhom, simulate = expression(rpoispp(DDlist[[paste0("RAG", ii)]]$estimate)), use.theory = TRUE, nsim = 39, global = TRUE); #This will build the envelops of the theoretical Poisson
      
      #To adjust X and Y scales to same size for both plots individually
      XYlims <- Norm.Axis(Kitest, KranInh);
      
      #To export the resulting plot
      Graph.Save(Kitest, ii, 1, "Ktest (Kinhom) ", "Ripley's K Function (Kinhom)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      Graph.Save(KranInh, ii, 1, "DistTest (Kinhom) ", "Ripley's K Function (Kinhom)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      #Besag's L function - inhomogeneous (Chapter 7 Section 10 page 242)
      #Require the point process to be correlation-staionary
      Litest <- Linhom(listXY[[paste0("RAG", ii)]], lambda = DDlist[[paste0("RAG", ii)]]$estimate); #Besag's L function for inhomogeneous intensity
      LranInh <- envelope.ppp(listXY[[paste0("RAG", ii)]], Linhom, simulate = expression(rpoispp(DDlist[[paste0("RAG", ii)]]$estimate)), use.theory = TRUE, nsim = 39, global = TRUE);

      #To adjust X and Y scales to same size for both plots individually
      XYlims <- Norm.Axis(Litest, LranInh);
      
      #To export the resulting plot
      Graph.Save(Litest, ii, 1, "Ltest (Linhom) ", "Besag's L Function (Linhom)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      Graph.Save(LranInh, ii, 1, "DistTest (Linhom) ", "Besag's L Function (Linhom)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      #Empty-space F function (Chapter 8 Section 3 page 203)
      
      Fitest <- Finhom(listXY[[paste0("RAG", ii)]], lambda = DDlist[[paste0("RAG", ii)]]$estimate); 
      # FranInh <- envelope.ppp(listXY[[paste0("RAG", ii)]], Finhom, simulate = expression(rpoispp(DDlist[[paste0("RAG", ii)]]$estimate)), use.theory = TRUE, nsim = 39, global = TRUE);
      FnGtra <- function(x) {asin(sqrt(x))};
      FranInh <- envelope.ppp(listXY[[paste0("RAG", ii)]], Finhom, simulate = expression(rpoispp(DDlist[[paste0("RAG", ii)]]$estimate)), nsim = 39, nrank = 1, global = TRUE, transform = expression(FnGtra(.)));
      
      #To adjust X and Y scales to same size for both plots individually
      XYlims <- Norm.Axis(Fitest, FranInh);
      
      #To export the resulting plot
      Graph.Save(Fitest, ii, 1, "Ftest (Finhom) ", "Empty-Space Function (Finhom)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      Graph.Save(FranInh, ii, 1, "DistTest (Finhom) ", "Empty-Space Function (Finhom)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      
      #To export the resulting plot
      tiff(paste0("FTPP", ii,".tif"), width = 150, height = 150, units = "mm", res = GraphRes, pointsize = 12, compression = 'lzw');
      plot(Fitest, . ~ theo, main = paste0("RAG", ii, " - Empty-Space function (P-P plot)"));
      dev.off();
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      #Nearest-neighbor function
      
      Gitest <- Ginhom(listXY[[paste0("RAG", ii)]], lambda = DDlist[[paste0("RAG", ii)]]$estimate); 
      # GranInh <- envelope.ppp(listXY[[paste0("RAG", ii)]], Ginhom, simulate = expression(rpoispp(DDlist[[paste0("RAG", ii)]]$estimate)), use.theory = TRUE, nsim = 39, global = TRUE);
      GranInh <- envelope.ppp(listXY[[paste0("RAG", ii)]], Ginhom, simulate = expression(rpoispp(DDlist[[paste0("RAG", ii)]]$estimate)), nsim = 39, nrank = 1, global = TRUE, transform = expression(FnGtra(.)));
      
      #To adjust X and Y scales to same size for both plots individually
      XYlims <- Norm.Axis(Gitest, GranInh);
      
      #To export the resulting plot
      Graph.Save(Gitest, ii, 1, "Gtest (Ginhom) ", "Nearest-Neighbour Function (Ginhom)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      Graph.Save(GranInh, ii, 1, "DistTest (Ginhom-test1) ", "Nearest-Neighbour Function (Ginhom)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      
      tiff(paste0("GTQQ", ii,".tif"), width = 150, height = 150, units = "mm", res = GraphRes, pointsize = 12, compression = 'lzw');
      plot(QQversion(Gitest), main = paste0("RAG", ii, " - Nearest-Neighbour function (Q-Q plot)"));
      dev.off();
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      Jtest <- Jinhom(listXY[[paste0("RAG", ii)]], lambda = DDlist[[paste0("RAG", ii)]]$estimate); #J function
      
      #To export the resulting plot
      Graph.Save(Jtest, ii, 1, "Jtest (Jinhom) ", "J Function (Jinhom)");
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      #Pair correlation function
      Kipcf <- pcfinhom(listXY[[paste0("RAG", ii)]], spar = 0.5, divisor = "d", method = "c"); #Pair correlation function; helps to interpret the output of Ripley's K test
      PCFenv <- envelope.ppp(listXY[[paste0("RAG", ii)]], pcfinhom, simulate = expression(rpoispp(DDlist[[paste0("RAG", ii)]]$estimate)), divisor = "d", use.theory = TRUE, nsim = 39, global = TRUE);
      
      #To adjust X and Y scales to same size for both plots individually
      XYlims <- Norm.Axis(Kipcf, PCFenv);
      
      #To export the resulting plot
      Graph.Save(Kipcf, ii, 1, "PCF (PCFinhom) ", "Pair Correlation Function (PCFinhom)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      Graph.Save(PCFenv, ii, 1, "DistTest (PCFinhom) ", "Pair Correlation Function (PCFinhom)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      #Determines S.E. for both K and L function
      LTSE <- lohboot(listXY[[paste0("RAG", ii)]], Linhom, simulate = expression(rpoispp(DDlist[[paste0("RAG", ii)]]$estimate)), global = TRUE);
      KTSE <- eval.fv(pi * LTSE^2);
      
      #To export the resulting plot
      Graph.Save(LTSE, ii, 1, "Ltest S.E. (Linhom) 1-", "Besag's L function S.E. (Linhom)");
      Graph.Save(KTSE, ii, 1, "Ktest S.E. (Kinhom) 1-", "Ripley's K function S.E. (Kinhom)");
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Scaled
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      #Uses the scaled distribution as a null hypothesis
      #Ripley's K function
      #If correlation-staionary is not true
      Kstest <- Kscaled(listXY[[paste0("RAG", ii)]], lambda = DDlist[[paste0("RAG", ii)]]$estimate); #Ripley's K function for inhomogeneous intensity
      # KranSc <- envelope.ppp(listXY[[paste0("RAG", ii)]], Kscaled, lambda = DDlist[[paste0("RAG", ii)]]$estimate, simulate = expression(rpoispp(DDlist[[paste0("RAG", ii)]]$estimate)), use.theory = TRUE, nsim = 39, global = TRUE);
      KranSc <- envelope.ppp(listXY[[paste0("RAG", ii)]], Kscaled, simulate = expression(rpoispp(DDlist[[paste0("RAG", ii)]]$estimate)), use.theory = TRUE, nsim = 39, global = TRUE);
      
      #To adjust X and Y scales to same size for both plots individually
      XYlims <- Norm.Axis(Kstest, KranSc);
      
      #To export the resulting plot
      Graph.Save(Kstest, ii, 1, "Ktest (Kscaled) ", "Ripley's K Function (Kscaled)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      Graph.Save(KranSc, ii, 1, "DistTest (Kscaled) ", "Ripley's K Function (Kscaled)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      #Besag's L functions
      #If correlation-staionary is not true
      Lstest <- Lscaled(listXY[[paste0("RAG", ii)]], lambda = DDlist[[paste0("RAG", ii)]]$estimate); #Besag's L function for inhomogeneous intensity
      LranSc <- envelope.ppp(listXY[[paste0("RAG", ii)]], Lscaled, simulate = expression(rpoispp(DDlist[[paste0("RAG", ii)]]$estimate)), use.theory = TRUE, nsim = 39, global = TRUE);
      
      #To adjust X and Y scales to same size for both plots individually
      XYlims <- Norm.Axis(Lstest, LranSc);
      
      #To export the resulting plot
      Graph.Save(Lstest, ii, 1, "Ltest (Lscaled) ", "Besag's L Function (Lscaled)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      Graph.Save(LranSc, ii, 1, "DistTest (Lscaled-test1) ", "Besag's L Function (Lscaled)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      Kspcf <- pcf(Kstest, spar = 0.5, method = "c", savefun = TRUE); #Pair correlation function; helps to interpret the output of Ripley's K test
      
      #To export the resulting plot
      Graph.Save(Kspcf, ii, 1, "PCF (PCFscaled) ", "Pair Correlation Function (PCFscaled)");
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      Ktest <- list(Kitest, KranInh, Kstest, KranSc);
      Ltest <- list(Litest, LranInh, Lstest, KranSc);
      Ftest <- list(Fitest, FranInh);
      Gtest <- list(Gitest, GranInh);
      Kpcf <- list(Kipcf, PCFenv, Kspcf);
      
    } else if (HomIntRes[ii] == FALSE) { #Uses homogeneous distribution as null hypothesis
      #Ripley's K function & Besag's L functions require the point process to be staionary, which assumse a homogeneous intensity
      Ktest <- Kest(listXY[[paste0("RAG", ii)]]); #Ripley's K function (Chapter 7 Section 3 page 203)
      Kran <- envelope.ppp(listXY[[paste0("RAG", ii)]], Kest, nsim = 39, fix.n = TRUE, nrank = 1, global = TRUE);
      
      XYlims <- Norm.Axis(Ktest, Kran);
      
      #To export the resulting plot
      Graph.Save(Ktest, ii, 1, "Ktest (K) ", "Ripley's K Function (K)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      Graph.Save(Kran, ii, 1, "DistTest (K) ", "Ripley's K Function (K)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      Ltest <- Lest(listXY[[paste0("RAG", ii)]]); #Besag's L function (Chapter 7 Section 3 page 207)
      Lran <- envelope.ppp(listXY[[paste0("RAG", ii)]], Lest, nsim = 39, nrank = 1, global = TRUE);
      
      XYlims <- Norm.Axis(Ltest, Lran);
      
      #To export the resulting plot
      Graph.Save(Ltest, ii, 1, "Ltest (L) ", "Besag's L Function (L)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      Graph.Save(Lran, ii, 1, "DistTest (L) ", "Besag's L Function (L)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      Ftest <- Fest(listXY[[paste0("RAG", ii)]]); #Empty-space function (Chapter 8 Section 3 page 261)
      Fran <- envelope.ppp(listXY[[paste0("RAG", ii)]], Fest, nsim = 39, fix.n = TRUE, nrank = 1, global = TRUE, transform = expression(FnGtra(.)));
      
      XYlims <- Norm.Axis(Ftest, Fran);
      
      #To export the resulting plot
      Graph.Save(Ftest, ii, 1, "Ftest (F) ", "Empty-Space Function (F)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      Graph.Save(Fran, ii, 1, "DistTest (F) ", "Empty-Space Function (F)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      
      tiff(paste0("FTPP", ii,".tif"), width = 150, height = 150, units = "mm", res = GraphRes, pointsize = 12, compression = 'lzw');
      plot(Ftest, . ~ theo, main = paste0("RAG", ii, " - Empty-Space function (P-P plot)"));
      dev.off();
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      Gtest <- Gest(listXY[[paste0("RAG", ii)]]); #Nearest-neighbour function (Chapter 8 Section 3 page 261)
      Gran <- envelope.ppp(listXY[[paste0("RAG", ii)]], Gest, nsim = 39, fix.n = TRUE, nrank = 1, global = TRUE, transform = expression(FnGtra(.)));
      
      XYlims <- Norm.Axis(Gtest, Gran);
      
      #To export the resulting plot
      Graph.Save(Gtest, ii, 1, "Gtest (G) ", "Nearest-Neighbour Function (G)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      Graph.Save(Gran, ii, 1, "DistTest (G) ", "Nearest-Neighbour Function (G)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      
      tiff(paste0("GTQQ", ii,".tif"), width = 150, height = 150, units = "mm", res = GraphRes, pointsize = 12, compression = 'lzw');
      plot(QQversion(Gtest), main = paste0("RAG", ii, " - Nearest-Neighbour function (Q-Q plot)"));
      dev.off();
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      Jtest <- Jest(listXY[[paste0("RAG", ii)]]); #J finction (Chapter 8 Section 3 page 261)
      
      Graph.Save(Jtest, ii, 1, "Jtest (J) ", "J Function (J)");
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      #Pair correlation function
      Kpcf <- pcf(Khtest, divisor = "d", spar = 0.5, method = "c"); #Pair correlation function; helps to interpret the output of Ripley's K test
      PCFenv <- envelope.ppp(listXY[[paste0("RAG", ii)]], pcf, nsim = 39, fix.n = TRUE, nrank = 1, global = TRUE)
      
      #To adjust X and Y scales to same size for both plots individually
      XYlims <- Norm.Axis(Kpcf, PCFenv);
      
      #To export the resulting plot
      Graph.Save(Kpcf, ii, 1, "PCF (PCF) ", "Pair Correlation Function (PCF)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      Graph.Save(PCFenv, ii, 1, "DistTest (PCF) ", "Pair Correlation Function (PCF)", xlim = c(XYlims[1,]), ylim = c(XYlims[2,]));
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      #Determines S.E. for both K and L function
      LTSE <- lohboot(listXY[[paste0("RAG", ii)]], Lest, global = TRUE);
      KTSE <- eval.fv(pi * LTSE^2);
      
      #To export the resulting plot
      Graph.Save(LTSE, ii, 1, "Ltest S.E. (L) ", "Besag's L function S.E. (L)");
      Graph.Save(KTSE, ii, 1, "Ktest S.E. (K) ", "Ripley's K function S.E. (K)");
    }
    KtestList <- list.append(KtestList, Ktest);
    LtestList <- list.append(LtestList, Ltest);
    KpcfList <- list.append(KpcfList, Kpcf);
    LTSElist <- list.append(LTSElist, LTSE);
    KTSElist <- list.append(KTSElist, KTSE);
    FtestList <- list.append(FtestList, Ftest);
    GtestList <- list.append(GtestList, Gtest);
    JtestList <- list.append(JtestList, Jtest);
  }
  RepCount <- c();
  for (zzz in 1:StatTest) {
    RepCount <- c(RepCount, paste("Rep.", zzz));
  }
  colnames(XYscaDF) <- c(RepCount);
  colnames(XYinhDF) <- c(RepCount);
  
  rm(Ki, Ks, Kitest, Litest, Kipcf, Fitest, Gitest, Jitest, KTSE, LTSE, Kstest, Lstest, Ktest, Ltest, Kpcf, Ftest, Gtest, Jtest);
  rm(XYlim, X1Lo, X1Hi, X2Lo, X2Hi, Y1Lo, Y1Hi, Y2Lo, Y2Hi, MinX, MaxX, MinY, MaxY, LimX, LimY, zzz);

# PPP Data--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
  #This creates a data frame summarizing all the relevant skin patch data for all samples in a table
  PPPData <- data.frame(PPPData, XYscaDF, XYinhDF);
  
  rm(XYscaDF, XYinhDF);
  
# 3.2.8. Validation of Independence Assumption--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
for (jj in 1:csvCount) {
  Kres(listXY[[paste0("RAG", jj)]], correction = "best");
  Gres(listXY[[paste0("RAG", jj)]], correction = "best");
  
  #To export the resulting plot
  Graph.Save(Kres, jj, 1, "Kres", "Independence Assumtion Test (Kres)");
  Graph.Save(Gres, jj, 1, "Gres", "Independence Assumtion Test (Gres)");
}
  
# 3.2.9. Simulation of true randomness----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  #This simulation tests, whether data is truly random or not
  
  #Generates all tables in which data is going to be combined
  MADtest <- list();
  DCLFtest <- list();
  LsigList <- list();
  LprogList <- list();
    
  for(kk in 1:csvCount) {
    #This is the non-graphical equivalent
    MADhom <- mad.test(listXY[[paste0("RAG", kk)]], Lest, nsim = 99, rmax = 4, fix.n = TRUE, use.theo = TRUE, savepatterns = TRUE); #a non-graphical approach to test for true data randomness
    DCLFhom <- dclf.test(MADhom);
    MADinh <- mad.test(listXY[[paste0("RAG", kk)]], Linhom, lambda = DDlist[[paste0("RAG", kk)]]$estimate, nsim = 99, rmax = 4, fix.n = TRUE, use.theo = TRUE, savepatterns = TRUE); #a non-grafical approache to test for true data randomness
    DCLFinh <- dclf.test(MADinh);
    
    MADtest <- list.append(MADtest, list("MAD (L)" = MADhom, "MAD (Linhom)" = MADinh));
    DCLFtest <- list.append(DCLFtest, list("DCLF (L)" = DCLFhom, "DCLF (Linhom)" = DCLFinh));
    
    LsigTrahom <- dclf.sigtrace(listXY[[paste0("RAG", kk)]], Lest, nsim=19);
    Graph.Save(LsigTrahom, kk, 1, "Sig-Traces (L) ", "Significance Traces (L)");
    LsigTrainh <- dclf.sigtrace(listXY[[paste0("RAG", kk)]], Linhom, lambda = DDlist[[paste0("RAG", kk)]]$estimate, nsim=19);
    Graph.Save(LsigTrainh, kk, 1, "Sig-Traces (Linhom) ", "Significance Traces (Linhom)");
    
    Lproghom <- dclf.progress(listXY[[paste0("RAG", kk)]], Lest, nsim=19);
    Graph.Save(Lproghom, kk, 1, "Prog-Plot (L) ", "Progress Plot (L)");
    Lproginh <- dclf.progress(listXY[[paste0("RAG", kk)]], Linhom, lambda = DDlist[[paste0("RAG", kk)]]$estimate, nsim=19);
    Graph.Save(Lproginh, kk, 1, "Prog-Plot (Linhom) ", "Progress Plot (Linhom)");
    
    LsigList <- list.append(LsigList, list("Sig. Traces (L)" = LsigTrahom, "Sig. Traces (Linhom)" = LsigTrainh));
    LprogList <- list.append(LprogList, list("Progress Plot (L)" = Lproghom, "Progress Plot (Linhom)" = Lproginh));
  }
  names(MADtest) <- MouseIDs;
  names(DCLFtest) <- MouseIDs;
  names(LsigList) <- MouseIDs;
  names(LprogList) <- MouseIDs;
  
  rm(MADhom, DCLFhom, MADinh, DCLFinh ,LsigTrahom, LsigTrainh, Lproghom, Lproginh);
  
# 3.2.10.	Test of anisotropy--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  for (ll in 1:csvCount) {
    
    Khoriz <- Ksector(listXY[[paste0("RAG", ll)]], begin = -15, end = 15, units = "degrees");
    Kvert <- Ksector(listXY[[paste0("RAG", ll)]], begin = 90-15, end = 90+15, units = "degrees");
    
    #To export the resulting plot
    tiff(paste0("Anisotropy", ll, " (1).tif"), res = GraphRes, width = 150, height = 150, units = "mm", pointsize = 12, compression = 'lzw');
    plot(Khoriz, trans/theo ~ r, lty = 2, main = paste0("RAG", ll, " - Anisotropy Analysis"));
    plot(Kvert, trans/theo ~ r, add  =TRUE);
    dev.off();
    
    XY <- listXY[[paste0("RAG", ll)]]
    dK <-  function(XY, ...) {
      K1 <- Ksector(XY, ..., begin = -15, end = 15, units = "degrees");
      K2 <- Ksector(XY, ..., begin = 90-15, end = 90+15, units = "degrees");
      eval.fv(K1 - K2);
    }
    CIdK <- varblock(listXY[[paste0("RAG", ll)]], dK, nx = 3);
    
    #To export the resulting plot
    tiff(paste0("Anisotropy", ll, " (2).tif"), res = GraphRes, width = 200, height = 200, units = "mm", pointsize = 12, compression = 'lzw');
    plot(CIdK, ylab = "Horizontal - Vertical", main = paste0("RAG", ll, " - Difference"));
    dev.off();
  }
  
  rm(Khoriz, Kvert, XY, CIdK);
  
# 3.2.11. Further tests of true randomness (Hopkins-Skellon, Stienien and Dirichlet tests)------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  #This is a formal test of complete spatial randomness (CSR) based short distances
  
  #Generates all tables in which data is going to be combined
  HStestList <- list();
  HStestA <- c();
  HStestP <- c();
  
  for (nn in 1:csvCount) {
    
    HStest <- hopskel.test(listXY[[paste0("RAG", nn)]], alternative = "clustered"); #Hopkins-Skellan test to test for CSR; alternative means alternative to null hypothesis (P < 0.05 = clustered)
    HStestList <- list.append(HStestList, HStest);
    HStestA <- c(HStestA, round(HStest$statistic, 5)); #A > 1 = regular, A == 1 = random, A < 1 = clustered
    HStestP <- c(HStestP, round(HStest$p.value, 5)); #P-value for tested hypothesis
    
    stienen(listXY[[paste0("RAG", nn)]], cols = "red", main = paste0("RAG", nn, " - Stienen Diagram")); #Explorarory graph for cluster observation
    dev.copy(tiff, paste0("Stienen", nn, ".tif"), res = GraphRes, width = 175, height = 150, units = "mm", pointsize = 12, compression = 'lzw');
    dev.off();
    
    plot(flipxy.tess(rotate(dirichlet(listXY[[paste0("RAG", nn)]]))), main = paste0("RAG", nn, " - Dirichlet Tessellation")); #Explorarory graph for cluster observation
    dev.copy(tiff, paste0("Dirichlet", nn, ".tif"), res = GraphRes, width = 150, height = 150, units = "mm", pointsize = 12, compression = 'lzw');
    dev.off();
  }
  
  names(HStestList) <- MouseIDs;
  
  rm(HStest);
  
# PPP Data--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  #This creates a data frame summarizing all the relevant skin patch data for all samples in a table
  PPPData <- data.frame(PPPData, "A value for Hopkins-Skellam test" = HStestA, "P-value for alt. Hypothesis (Clustered)" = HStestP);
  
  rm(HStestA, HStestP);
  
# 3.2.11. Dao-Genton Test-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  #Dao-Genton Test to improve on Monte Carlo to eliminate conservatism
  for (oo in 1:csvCount) {
    dg.test(listXY[[paste0("RAG", oo)]], Lest, nsim = 19);
    dg.envelope(listXY[[paste0("RAG", oo)]], Lest);
    
    dg.test(listXY[[paste0("RAG", oo)]], Linhom, nsim = 19);
    dg.envelope(listXY[[paste0("RAG", oo)]], Linhom);
  }
} #This on closes the loop for duplicate measures

# PPP data save---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Saves the output data compilation file
write.csv(PPPData, "PPP Data.csv");
