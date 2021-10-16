####################### VALE INSTITUTE OF TECHNOLOGY ##########################

######################## LANDSCAPE GENOMICS TUTORIAL ##########################
############### STEP 2.5: ESTIMATED EFFECTIVE MIGRATION SURFACES ###############




### Script prepared by Jeronymo Dalapicolla ###
### Based on Petkova, D., Novembre, J. & Stephens, M. Visualizing spatial population structure with estimated effective migration surfaces. Nat Genet 48, 94–100 (2016). https://doi.org/10.1038/ng.3464
### and on the GITHUB page https://github.com/dipetkov/eems and the Supplementary Material of this paper



#### PRE-ANALYSIS #### 

##1. INPUTS FOR THIS TUTORIAL
#A. THE FILE ".VCF" CLEANED AFTER FILTERING AND WITH DELIMITED GENETIC CLUSTERS, STEP 2.
#B. THE ".TXT" FILE WITH GEOGRAPHICAL COORDINATES FOR EACH INDIVIDUAL
#C. A SHAPEFILE (.SHP) WITH THE STUDY AREA TO ESTIMATE THE MIGRATION 



##2. GOALS FOR THIS STEP
#A. VISUALIAZING SPATIAL POPULATION STRUCTURE WITH ESTIMATED EFFECTIVE MIGRATION SURFACES.




##3. CHOOSE A FOLDER FOR RUNNING THE ANALYSES. THE FILES MUST BE THERE! 
#A. IN RStudio GO TO  SESSION >> SET WORKING DIRECTORY >> CHOOSE DIRECTORY.. IN RStudio TOOL BAR OR USE THE SHORCUT CTRL+SHIFT+H



##4. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT
rm(list=ls())





##5. LOAD THE AUXILIARY FUNCTIONS  

VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}


#' Compute the dissimilarity matrix, using the "mean allele frequency"
#' imputation method, which corresponds to `bed2diffs-v2`. See the `bed2diffs`
#' README for more information about the two methods of imputing missing
#' that are implemented in `bed2diffs-v1` and `bed2diffs-v2`.
#+
bed2diffs_v2 <- function(Geno) {
  nIndiv <- nrow(Geno)
  nSites <- ncol(Geno)
  Miss <- is.na(Geno)
  ## Impute NAs with the column means (= twice the allele frequencies)
  Mean <- matrix(colMeans(Geno, na.rm = TRUE), ## a row of means
                 nrow = nIndiv, ncol = nSites, byrow = TRUE) ## a matrix with nIndiv identical rows of means
  Mean[Miss == 0] <- 0 ## Set the means that correspond to observed genotypes to 0
  Geno[Miss == 1] <- 0 ## Set the missing genotypes to 0 (used to be NA) 
  Geno <- Geno + Mean
  ## Compute similarities
  Sim <- Geno %*% t(Geno) / nSites
  SelfSim <- diag(Sim) ## self-similarities
  vector1s <- rep(1, nIndiv) ## vector of 1s
  ## This chunk generates a `diffs` matrix
  Diffs <- SelfSim %*% t(vector1s) + vector1s %*% t(SelfSim) - 2 * Sim
  Diffs
}

#' Here is a function that implements the "pairwise.complete.obs" method, 
#' which corresponds to `bed2diffs-v1`. The straightforward implementation
#' uses a double loop, so would be slow if the sample size is large.
#+
bed2diffs_v1 <- function(Geno) {
  nIndiv <- nrow(Geno)
  nSites <- ncol(Geno)
  Diffs <- matrix(0, nIndiv, nIndiv)
  
  for (i in seq(nIndiv - 1)) {
    for (j in seq(i + 1, nIndiv)) {
      x <- Geno[i, ]
      y <- Geno[j, ]
      Diffs[i, j] <- mean((x - y)^2, na.rm = TRUE)
      Diffs[j, i] <- Diffs[i, j]
    }
  }
  Diffs
}



##6. INSTALL AND LOAD THE PACKAGES
#For r2vcftools do you need install VCFTools in you computer:https://vcftools.github.io/index.html
#Basic Packages for installation:
if (!require('remotes'))      install.packages('remotes');           library('remotes')

#From Github or BiocManager:
if (!require('r2vcftools'))   remotes::install_github("nspope/r2vcftools");          library('r2vcftools')

#From CRAN R:
if (!require("adegenet"))     install.packages("adegenet");          library("adegenet")
if (!require('vcfR'))         install.packages("vcfR");              library('vcfR')
if (!require('tidyverse'))    install.packages("tidyverse");         library('tidyverse')
if (!require('raster'))       install.packages("raster");            library('raster')
if (!require('ggplot2'))      install.packages("ggplot2");           library('ggplot2')





#### ANALYSIS ---- 

#### 1. INSTALL SOFTWARES IN YOUR COMPUTER. I USED UBUNTU 20.04 ----
###A. BOOST
# Open your terminal and install the software usinf brew or pip ou snap.
# I used Homebrew. To see how to install Homebrew:
# https://www.edivaldobrito.com.br/como-instalar-o-homebrew-no-ubuntu-20-04-debian-10-e-derivados/
# In terminal:
# $ sudo brew install boost

# You need to find this software folder in your computer. In my this one is a hidden folder.
# Ask to show the hidden folders and use CRTL+F to find it.
 #/home/linuxbrew/.linuxbrew/Cellar/boost/1.76.0/


###B. EIGEN
# Open your terminal and install the software usinf brew or pip ou snap.
# $ sudo brew install eigen

# You need to find this software folder in your computer. In my this one is a hidden folder.
# Ask to show the hidden folders and use CRTL+F to find it.
#/home/linuxbrew/.linuxbrew/Cellar/eigen/3.4.0_1/


###C. EEMS
# Go to the GITHUB page https://github.com/dipetkov/eems >>
# Download all files clicking in "CODE" >> "DOWNLOAD ZIP" close to description of the github "ABOUT"; 
# Unzip the folder "eems-master" in your computer and move to any place in your computer, preferably the root (:/C or /home) 

# Open text editor (WordPad, Text, Notepad) the file "Makefile" in "/eems-master/runeems_snps/src"  
# Edit the path for Eigen and Boost folders software according to your computer:
# EIGEN_INC = /home/linuxbrew/.linuxbrew/Cellar/eigen/3.4.0_1/include/eigen3
# BOOST_LIB = /home/linuxbrew/.linuxbrew/Cellar/boost/1.76.0/lib
# BOOST_INC = /home/linuxbrew/.linuxbrew/Cellar/boost/1.76.0/include

#If you'll use the runeems_sats you need to do the same step for Makefile in "/eems-master/runeems_sats/src"  

#After you open the terminal, go the folder where Makefile is, and compile the file:
# $ cd ./eems-master/runeems_snps/src
# $ sudo make linux
#for Mac and Windows I don't know the commands

#New files will be create such as "draw.o", "eems.o", "graph.o" and others.
#Now you have the software installed!



#### 2. PREPARING THE INPUT FILES ----
###A. LONG/LAT = COORD FILE:
# two cols in a txt file saved as ".coord" NOT ".txt". First col is LONGITUDE, second LATITUDE!
# each row represent an individual in the same order than genetic input!
# the two col is separated by space not tab! Without rownames and colnames!

# Load neutral .vcf file with geographical information in the metafile
snps_neutral = vcfLink("vcf/pilocarpus_filtered_neutral.vcf", overwriteID=T)
VCFsummary(snps_neutral) #277 individuals and 5277 SNPs.
names(snps_neutral@meta) #verify col names in metafile
head(snps_neutral@meta)

#select long/lat cols
locations = snps_neutral@meta[,c(5:6)] 
head(locations)

#save the coord file
write.table(locations, "pilocarpus_filtered_neutral.coord", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ")




###B. STUDY AREA = OUTER FILE:
# two cols in a txt file saved as ".outer" NOT ".txt". First col is LONGITUDE, second LATITUDE!
# each row represent a point in the contour of the study area!
# the two col is separated by space not tab! Without rownames and colnames!

# define the contour using a shapefile:
flona = raster::shapefile("maps/shapefiles/Flona_Carajas.shp")
plot(flona)
class(flona)

#select points in the contour: 
outer_points = flona@polygons[[1]]@Polygons[[1]]@coords

#save as the same name as coord file!
write.table(outer_points, "pilocarpus_filtered_neutral.outer", quote = F, row.names = F, col.names = F, sep=" ")



###C. GENETIC DIFFERENCES IN AFS USING VCF/GENIND = DIFFS FILE:
#Load vcf:
vcf = read.vcfR("pilocarpus_filtered_neutral.vcf", verbose = TRUE)

#Convert "VCF" to "GENIND"
data = vcfR2genind(vcf)
data@tab[1:5,1:5]

# Select the genotype matrix.
Geno = data@tab

# We are using a filtered dataset but If you didn't filter you need to keep only bi-allelic alleles.
multi.loci = names(which(data@loc.n.all != 2))
multi.cols = which(grepl(paste0("^", multi.loci, "\\.\\d+$", collapse = "|"), colnames(Geno)))
if (length(multi.cols)) Geno = Geno[, - multi.cols]

dim(Geno)

#' Let's convert the matrix to 0-1 labeling. I arbitrarily choose one allele to be 
#' the "derived" allele and, for each individual, count how many copies of the derived
#'  allele it carries. This is very easy if the tab matrix is of type "codom".
stopifnot(identical(data@type, 'codom'))

#' Since the labeling does not matter for computing differences,
#' I pick the second allele to label as the "derived" allele. 
#' That is, I pick all loci whose name ends with `.1`. Some files could be ".01"
Geno = Geno[, str_detect(colnames(Geno), "\\.1$")]


# Next compute the dissimilarity matrix, using the "pairwise.complete.obs" method (v1) and the "mean allele frequency" imputation method (v2)
Diffs_v1 = bed2diffs_v1(Geno)
Diffs_v2 = bed2diffs_v2(Geno)
Diffs_v1 = round(Diffs_v1, digits = 6)
Diffs_v2 = round(Diffs_v2, digits = 6)

#' Check that the dissimilarity matrix has one positive eigenvalue and `nIndiv-1`
#' negative eigenvalues, as required by a full-rank Euclidean distance matrix.
sort(round(eigen(Diffs_v1)$values, digits = 2))
sort(round(eigen(Diffs_v2)$values, digits = 2))

#' Save the file where the condition is hold. In my case both files:
write.table(Diffs_v2, "pilocarpus_filtered_neutral.diffs", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(Diffs_v1, "pilocarpus_filtered_neutral_v1.diffs", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)



#### 3. RUNNING EEMS ----
#A. Edit .ini files. These file are in "/eems-master/runeems_snps/src". Open it in a text editor.
# See the manual to set the parameters.
# All inputs must be the same name with different types: .coord, .outer, and .diffs
# For each chain the seed number must be different
# datapath = /home/b0219.itv/Documentos/PilocarpusProject/EEMS_test/inputs/pilocarpus_filtered_neutral
# mcmcpath = /home/b0219.itv/Documentos/PilocarpusProject/EEMS_test/outputs/chain2
# nIndiv = 277
# nSites = 5277
# nDemes = 300
# diploid = true
# numMCMCIter = 10000000
# numBurnIter = 5000000
# numThinIter = 10000

#B. Add these command in the file:
# seed = 456
# mSeedsProposalS2 = 0.0005
# qSeedsProposalS2 = 0.04
# mEffctProposalS2 = 0.7
# qEffctProposalS2 = 0.02
# mrateMuProposalS2 = 0.0005


#C. Run EEMS:
# Open the terminal, using "cd" go to "/eems-master/runeems_snps/src"
# $ ./runeems_snps --params params-chain1.ini
# Depending how much data and iteration it can take hours!

#C. Tuning the models:
#After run the all the iterations, in the output open the file "eemsrun.txt" with the final results, and analyses the final acceptance proportions:

#Acceptance proportions:
# (93556/623769) = 15% for proposal type "qTileRate",		 with proposal variance "qEffctProposalS2"
# (49000/625485) = 7.8% for proposal type "qTileMove",		 with proposal variance "qSeedsProposalS2"
# (206282/625083) = 33% for proposal type "qBirthDeath"
# (610795/1876348) = 33% for proposal type "mTileRate",		 with proposal variance "mEffctProposalS2"
# (322860/1248302) = 26% for proposal type "mMeanRate",		 with proposal variance "mrateMuProposalS2"
# (629952/1875093) = 34% for proposal type "mTileMove",		 with proposal variance "mSeedsProposalS2"
# (409379/1875154) = 22% for proposal type "mBirthDeath"
# (1003379/1250766) = 80% for proposal type "degrees of freedom"

#Good values are between 20-30%, regular values between 10-20% and 30-40%.
#You can change the values in the .ini files get good values in the porportions.
#If you Increase values in .ini file, it will decrease the values of acceptance proportions
# High values in parameters >> low values in proportions.

#D. Number of Chains:
#After tuned the model test, you need run multiple chains and use the average model to build the graphics
#Authors recommend at least 8 chains. I run 10 chains.




#### 4. CHECKING CONVERGENCY AND PLOTTING RESULTS ----
###A. Install rEEMSplots
## Move the working directory in RStudio to "/eems-master/plotting/" 
## Check that the current directory contains the rEEMSplots source directory
if (file.exists("./rEEMSplots")) {
  install.packages("rEEMSplots", repos = NULL, type = "source")
} else {
  stop("Move to the directory that contains the rEEMSplots source to install the package.")
}

## Possibly change the working directory with setwd()



###B. Generate graphics
#load the new library
library (rEEMSplots)

#define where the outputs are (mcmcpath) and where to save the figures (plotpath)
mcmcpath = c("/home/b0219.itv/Documentos/PilocarpusProject/EEMS_test/outputs/chain1",
             "/home/b0219.itv/Documentos/PilocarpusProject/EEMS_test/outputs/chain2",
             "/home/b0219.itv/Documentos/PilocarpusProject/EEMS_test/outputs/chain3",
             "/home/b0219.itv/Documentos/PilocarpusProject/EEMS_test/outputs/chain4",
             "/home/b0219.itv/Documentos/PilocarpusProject/EEMS_test/outputs/chain5",
             "/home/b0219.itv/Documentos/PilocarpusProject/EEMS_test/outputs/chain6",
             "/home/b0219.itv/Documentos/PilocarpusProject/EEMS_test/outputs/chain7",
             "/home/b0219.itv/Documentos/PilocarpusProject/EEMS_test/outputs/chain8",
             "/home/b0219.itv/Documentos/PilocarpusProject/EEMS_test/outputs/chain9",
             "/home/b0219.itv/Documentos/PilocarpusProject/EEMS_test/outputs/chain10")
             
plotpath = "/home/b0219.itv/Documentos/PilocarpusProject/EEMS_test/outputs/10M_GRAPHICS"


#If you need to add other shapefile. I didn't test with raster
#set the projection:
projection_none = "+proj=longlat +datum=WGS84"
projection_mercator = "+proj=merc +datum=WGS84"

#load the shapefile and change the projection
cangas = raster::shapefile("Cangas.shp")
cangas = spTransform(cangas, CRSobj = CRS(projection_none))
plot(cangas)


#Run this script it will save the figures automatically. See documentation for details 
eems.plots(mcmcpath,
           plotpath,
           longlat = TRUE,
           plot.height = 8,
           plot.width = 7,
           res = 600, #only for png
           out.png = TRUE, # FALSE will save it in pdf
           add.grid = TRUE,
           col.grid = "gray80",
           lwd.grid = 0.25,
           add.outline = TRUE,
           col.outline = "black",
           lwd.outline = 1.5,
           add.demes = TRUE,
           col.demes = "black",
           pch.demes = 1,
           min.cex.demes = 0.5,
           max.cex.demes = 2.5,
           m.plot.xy = { plot(cangas, col = NA, add = TRUE) },
           q.plot.xy = { plot(cangas, col = NA, add = TRUE) })




###C. Reproduce manually some graphics
#After a RData is created. You can load it!
load(paste0(plotpath, "-rdist.RData"))
ls()
#> [1] "B.component" "G.component" "W.component" "xym.values"  "xyq.values"

# Redo the plotpath-rdist01.png,
## which plots observed vs fitted dissimilarities between demes
ggplot(B.component %>% filter(size > 1),
       aes(fitted, obsrvd)) +
  geom_point() +
  theme_bw()


## Reproduce plotpath-rdist02.png,
## which plots observed vs fitted dissimilarities within demes
ggplot(W.component %>% filter(size > 1),
       aes(fitted, obsrvd)) +
  geom_point()

## Reproduce plotpath-rdist03.png,
## which plots observed dissimilarities against great circle distances between demes
ggplot(W.component %>% filter(size > 1),
       aes(fitted, obsrvd)) +
  geom_point()

## Empty matrices unless you have specified additional coordinates
## at which to estimate the migration and diversity rates
## with the `xy.coords` argument to `eems.plots`
xym.values
xyq.values

##END;
