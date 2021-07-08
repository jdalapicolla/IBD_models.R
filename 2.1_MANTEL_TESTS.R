####################### VALE INSTITUTE OF TECHNOLOGY ##########################

######################## LANDSCAPE GENOMICS TUTORIAL ##########################
#####################   STEP 04: IBD BY MANTEL TEST   #########################




### Script prepared by Jeronymo Dalapicolla, Carolina S. Carvalho, Luciana C. Resende-Moreira, Jamille C. Veiga, and Rodolfo Jaffé ###



#### PRE-ANALYSIS #### 

##1. INPUTS FOR THIS TUTORIAL ----
#A. THE FILE ".VCF" CLEANED AFTER FILTERING AND WITH DELIMITED GENETIC CLUSTERS, STEP 2.

#B. THE FILE "functions_LanGen.R" WITH FUNCTIONS DESIGNED FOR THIS PIPELINE IN THE WORKING DIRECTORY. YOU CAN DOWNLOAD IN  https://github.com/jdalapicolla/ OR https://github.com/rojaff/LanGen_pipeline



##2. GOALS FOR THIS STEP:
#A. ISOLATION-BY-DISTANCE ANALYSIS USING MANTEL TEST



##3. CHOOSE A FOLDER FOR RUNNING THE ANALYSES. THE FILES MUST BE THERE! 
#A. IN RStudio GO TO  SESSION >> SET WORKING DIRECTORY >> CHOOSE DIRECTORY.. IN RStudio TOOL BAR OR USE THE SHORCUT CTRL+SHIFT+H


##4. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT:
rm(list=ls())




##5. LOAD THE FILE "functions_LanGen.R" WITH FUNCTIONS TO BE USED ON THIS STEP. MORE INFORMATION ON FUNCTIONS IN NUMBER 2.
source("functions_LanGen.R")


##6. INSTALL AND LOAD THE PACKAGES ----
#For r2vcftools do you need install VCFTools in you computer:https://vcftools.github.io/index.html
#Basic Packages for installation:
if (!require('remotes'))      install.packages('remotes');           library('remotes')
if (!require('BiocManager'))  install.packages('BiocManager');       library('BiocManager')
if (!require('pacman'))       install.packages('pacman');            library('pacman')
if (!require('devtools'))     install.packages('devtools');          library('devtools')

#From Github or BiocManager:
if (!require('r2vcftools'))   remotes::install_github("nspope/r2vcftools");          library('r2vcftools')
if (!require('LEA'))          BiocManager::install("LEA");                           library('LEA')

#From CRAN R:
if (!require('tidyverse'))    install.packages("tidyverse");         library('tidyverse') 
if (!require('raster'))       install.packages("raster");            library('raster')
if (!require('rgdal'))        install.packages("rgdal");             library('rgdal')
if (!require('geosphere'))    install.packages("geosphere");         library('geosphere')
if (!require('ggplot2'))      install.packages("ggplot2");           library('ggplot2') 
if (!require('ade4'))         install.packages("ade4");              library('ade4')
if (!require('reshape2'))     install.packages("reshape2");          library('reshape2') 
if (!require('usedist'))      install.packages("usedist");           library('usedist')
if (!require('car'))          install.packages("car");               library('car')

#Load multiple packages using the package 'pacman'. If the package is missing "p_load" will download it from CRAN. "" in packages names is not mandatory.
pacman::p_load(r2vcftools, usedist, raster, rgdal, geosphere, ade4, tidyverse, reshape2, ggplot2, car)



##7. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS ----
create_dir(c("./Results/Step04/Mantel_tests"))



##8. CREATE A PATTERN FOR GRAPHIC FOLDERS S TO SAVE THE RESULTS ----
theme_genetics = theme(axis.text=element_text(size=10, color="black"), #text in ticks axes
                       axis.title=element_text(size=12, face="bold"), #label axes
                       axis.line = element_line(colour = "black", size = 1, linetype = "solid"), #line on axes
                       axis.ticks = element_line(colour = "black", size = 1, linetype = "solid"), #line on ticks
                       axis.ticks.length = unit(.25, "cm"), #ticks length
                       axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), #space between axis and label
                       axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)), #space between axis and label
                       strip.text.x = element_text(size = 12, face="bold"), #facets label 
                       panel.grid.major = element_blank(), # remove grids
                       panel.grid.minor = element_blank(), # remove grids
                       panel.background = element_blank(), # remove background
                       panel.border = element_blank()) # remove borders)  



#### ANALYSIS ---- 


#### 1. LOAD FILES -----
#A. Project name:
project_name = "pilocarpus"

#B. Load neutral .vcf file with geographical information and genetic clusters ID, choosen in step 2. 
snps_neutral = vcfLink(paste0("vcf/", project_name, "_filtered_neutral_clusters.vcf"), overwriteID=T)
VCFsummary(snps_neutral) #277 individuals and 5268 SNPs.
names(snps_neutral@meta) #verify col names in metafile

#C. Number and name of cluster and method:
optimal_K = 4
name_clusters = c("A", "B", "C", "D")

#D.Position of samples by population/genetic clusters:
for (i in name_clusters){
  pop = which(snps_neutral@meta$POP_ID == i)
  assign(paste0("pop_POSI_", i), pop)
}




###2. GENETIC DISTANCE ----
#We calculated genetic distance by individuals or clusters in step 03. Chose one genetic distance metric for Mantel.

#A. Load distance files by individuals. We just have 4 populations, minimal 5 points for a regression.
#PCA with 95% of variance explained
PCA_95 = read.csv(paste0("Results/Step03/Distance/PCA_Distance_95Var_IND_mahalanobis_", project_name, ".csv"), row.names = 1)
PCA_95[1:5, 1:5]


#B. CREATE GENETIC 'DIST' OBJECTS FOR ALL SAMPLES:
genDIST_all = as.dist(PCA_95)


#C. CREATE GENETIC 'DIST' OBJECTS BY CLUSTERS:
#Subset By Pop A:
genDIST_popA= dist_subset(genDIST_all, snps_neutral@sample_id[pop_POSI_A])
#verify dimensions
length(genDIST_popA)

#Subset By Pop B:
genDIST_popB= dist_subset(genDIST_all, snps_neutral@sample_id[pop_POSI_B])
#verify dimensions
length(genDIST_popB)

#Subset By Pop C:
genDIST_popC= dist_subset(genDIST_all, snps_neutral@sample_id[pop_POSI_C])
#verify dimensions
length(genDIST_popC)

#Subset By Pop D:
genDIST_popD= dist_subset(genDIST_all, snps_neutral@sample_id[pop_POSI_D])
#verify dimensions
length(genDIST_popD)



###3. GEOGRAPHICAL DISTANCE ----
#You can use any geographical distance matrix in here, topographical distance, distance along rivers, etc. Here we will calculate geographic distance, based on Euclidean distance. Package geoshere works with Lat/Long coordinates. Others package may need UTM coordinates. Be careful! 

#A. Create a data frame with the geographical coordinates:
names(snps_neutral@meta)
coord_SN = snps_neutral@meta[,c(2,6:7)] # cols "sample_name", "Longitude", "Latitude" 
head(coord_SN)

#B. Set long and lat colunms
coordinates(coord_SN) = coord_SN[,c(2,3)]
projection(coord_SN) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84")

#C. Subset coordinates by populations:
coord_POPA = coord_SN[pop_POSI_A,]
coord_POPA
coord_POPB = coord_SN[pop_POSI_B,]
coord_POPB
coord_POPC = coord_SN[pop_POSI_C,]
coord_POPC
coord_POPD = coord_SN[pop_POSI_D,]
coord_POPD

#D. Create matrices of geographical distance by pop:
coord_files = list(coord_POPA, coord_POPB, coord_POPC, coord_POPD)

for (l in 1:length(coord_files)){
  coord = coord_files[[l]]
  sdist_SN  = distm(coord, fun=distGeo)
  rownames(sdist_SN) = coord$sample_name
  colnames(sdist_SN) = coord$sample_name
  sdist_SN  = sdist_SN /1000
  sdist_SN = as.dist(sdist_SN)
  assign(paste0("geoDIST_pop", name_clusters[l]), sdist_SN)
}

#E. Create matrices of geographical distance for all individuals:
sdist_SN  = distm(coord_SN, fun=distGeo)
rownames(sdist_SN) = snps_neutral@sample_id
colnames(sdist_SN) = snps_neutral@sample_id
sdist_SN  = sdist_SN /1000
geoDIST_all = as.dist(sdist_SN)


#F. Verify dimension for 'dist' objects:
length(geoDIST_all)
length(genDIST_all)

length(geoDIST_popA)
length(genDIST_popA)

length(geoDIST_popB)
length(genDIST_popB)

length(geoDIST_popC)
length(genDIST_popC)

length(geoDIST_popD)
length(genDIST_popD)



###4. MANTEL TESTS ----
#A. By Species
mantel_all = mantel.rtest(genDIST_all, geoDIST_all, 10000)

#B. By populations
mantel_popA = mantel.rtest(genDIST_popA, geoDIST_popA, 10000)
mantel_popB = mantel.rtest(genDIST_popB, geoDIST_popB, 10000)
mantel_popC = mantel.rtest(genDIST_popC, geoDIST_popC, 10000)
mantel_popD = mantel.rtest(genDIST_popD, geoDIST_popD, 10000)

#C. Compile results and save as data.frame
results_mantel = matrix(NA, 5, 2)
colnames(results_mantel) = c("r", "p-value")
rownames(results_mantel) = c("All Samples", "POPA", "POPB", "POPC", "POPD")
results_mantel

results_mantel[1,] = cbind (mantel_all$obs, mantel_all$pvalue)
results_mantel[2,] = cbind (mantel_popA$obs, mantel_popA$pvalue)
results_mantel[3,] = cbind (mantel_popB$obs, mantel_popB$pvalue)
results_mantel[4,] = cbind (mantel_popC$obs, mantel_popC$pvalue)
results_mantel[5,] = cbind (mantel_popD$obs, mantel_popD$pvalue)
results_mantel

#D. Save results
write.csv(results_mantel, paste0("./Results/Step04/Mantel_tests/Mantel_Results_", project_name, ".csv"))








#### 5. PERMUTATION TESTS BY POPULATION -----
#A. Define populations positions to subset in a list
populations = list(snps_neutral@sample_id[pop_POSI_A], snps_neutral@sample_id[pop_POSI_B], snps_neutral@sample_id[pop_POSI_C], snps_neutral@sample_id[pop_POSI_D])

#B. Define number of populations
pop_counter = length(populations)

#C. Define a object to save the results
mantel_perm = matrix(NA,pop_counter,2)

#D. Perform the permutation
for(i in 1:pop_counter) {
  samples = !(snps_neutral@sample_id %in% populations[[i]])
  genDIST_temp = dist_subset(genDIST_all, snps_neutral@sample_id[samples])
  geoDIST_temp = dist_subset(geoDIST_all, snps_neutral@sample_id[samples])
  mt = mantel.rtest(genDIST_temp, geoDIST_temp, 10000)
  mantel_perm[i, ] = cbind (mt$obs, mt$pvalue)
}

mantel_perm
colnames(mantel_perm) = c("r", "p-value")
rownames(mantel_perm) = c("POP1", "POP2", "POP3", "POP4")
mantel_perm

#E. Save the results:
write.csv(mantel_perm, paste0("./Results/Step04/Mantel_tests/Mantel_Results_Permutation_", project_name, ".csv"))







#### 6. PERFORMING ANCOVA AND TUKEY'S TESTS TO COMPARE SLOPES IN MANTEL TESTS -----
#A. Convert distance matrices to data frames
#set files in same order:
ancova_datasets = c("A", "B", "C", "D", "SPE")
gen_dist = list(genDIST_popA, genDIST_popB, genDIST_popC, genDIST_popD, genDIST_all) 
geo_dist = list(geoDIST_popA, geoDIST_popB, geoDIST_popC, geoDIST_popD, geoDIST_all)

for (k in 1:length(ancova_datasets)) {
  mta = as.matrix(gen_dist[[k]])
  mta[upper.tri(mta, diag = T)] = NA
  mta = mta %>%
    melt %>%
    na.omit %>%
    arrange(., Var1) %>%
    setNames(c("X1", "X2","GenDist"))
  
  mtb = as.matrix(geo_dist[[k]])
  mtb[upper.tri(mtb, diag = T)] = NA
  mtb = mtb %>%
    melt %>%
    na.omit %>%
    arrange(., Var1) %>%
    setNames(c("X1", "X2","GeoDist"))
  
  df = as.data.frame(cbind(mta[,3], mtb[,3])) %>%
    mutate(CLUSTER = ancova_datasets[k]) %>%
    setNames(c("GEN", "GEO", "CLUSTER")) %>%
    mutate(GEN = scale(GEN), GEO = scale(GEO))
  
  assign(paste0("ancova_pop_", ancova_datasets[k]), df)
}

df_ancova_clusters = rbind(ancova_pop_A, ancova_pop_B, ancova_pop_C, ancova_pop_D)
head(df_ancova_clusters)
tail(df_ancova_clusters)
df_ancova_species = ancova_pop_SPE #only one species in our example
head(df_ancova_species)


#B. Performing the ANCOVA. You can replace clusters by species if you have two or more species.
model1 = lm (GEN ~ GEO + CLUSTER + GEO:CLUSTER, data = df_ancova_clusters)
Anova(model1, type="II")
#Anova Table (Type II tests)

#Response: GEN
#             Sum Sq    Df F value   Pr(>F)    
#GEO             0.1     1  0.1391   0.7091    
#CLUSTER         0.0     3  0.0000   1.0000    
#GEO:CLUSTER    57.6     3 19.2942 1.79e-12 ***
#Residuals   11083.3 11137           

### Interaction is significant, so the slope among groups
### is different. 

model2 = lm (GEN ~ GEO + CLUSTER, data = df_ancova_clusters)
Anova(model2, type="II")
#Anova Table (Type II tests)

#Response: GEN
#           Sum Sq    Df F value Pr(>F)
#GEO           0.1     1  0.1385 0.7098
#CLUSTER       0.0     3  0.0000 1.0000
#Residuals 11140.9 11140   

### The category variable (Species) is not significant,
### so the intercepts among groups are not different


#### 7. PLOTING MANTEL TESTS IN GRAPHS -----
#By species
plot_mantel_species = 
ggplot(data=df_ancova_species) +
  geom_point(mapping = aes(x=GEO, y = GEN, shape=factor(CLUSTER)), size=3, color ="black", alpha=0.5) +
  scale_shape_manual(values=c(21))+
  geom_smooth(mapping = aes(x=GEO, y = GEN, linetype = factor(CLUSTER)), method='lm', color = "red", size=1) +
  theme_bw() +
  theme_genetics +
  theme(legend.text = element_text(size=10, face="italic"),
        legend.title = element_text(size=12, face="bold")) +
  xlab("Scaled Geographic Distances (Km)") + ylab("PC Distance (Broken Stick Rule)") +
  labs(color  = "Species", linetype = "Species", shape = "Species")
  

plot_mantel_species

pdf("./Results/Step04/Mantel_tests/Mantel_Plot_SPECIES.pdf", onefile=FALSE)
plot_mantel_species
dev.off()



#By clusters:
colors_pop = c('#ffff00','#ffc0cb', "#ff0000", '#0000ff') #color for the genetic cluster

plot_mantel_clusters = 
ggplot(data=df_ancova_clusters) +
  aes(x=GEO, y = GEN, shape=factor(CLUSTER), linetype = factor(CLUSTER), color = factor(CLUSTER)) +
  geom_point(size=3, alpha=0.5) +
  geom_smooth(method='lm', size=1) +
  scale_shape_manual(values=c(21,22,23,24))+
  scale_color_manual(values=colors_pop) +
  theme_bw() +
  theme_genetics +
  theme(legend.text = element_text(size=10, face="italic"),
        legend.title = element_text(size=12, face="bold")) +
  xlab("Scaled Geographic Distances (Km)") + ylab("PC Distance (Broken Stick Rule)") +
  labs(color  = "Clusters", linetype = "Clusters", shape = "Clusters")

plot_mantel_clusters

pdf("./Results/Step04/Mantel_tests/Mantel_Plot_CLUSTERS.pdf", onefile=FALSE)
plot_mantel_clusters
dev.off()


##END
