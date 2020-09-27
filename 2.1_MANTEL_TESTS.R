################ VALE INSTITUTE OF TECHNOLOGY ################
################  LANDSCAPE GENOMICS TUTORIAL ################
################ ISOLATION BY DISTANCE MODELS ################
################    STEP 01: MANTEL TESTS     ################

### Script prepared by Jeronymo Dalapicolla & Joyce Rodrigues do Prado.


###### PRE-ANALYSIS -----
##1. INFORMATION
#A. FOR ALL STEPS, PLEASE GO TO: https://github.com/jdalapicolla/
#B. BASED ON PIPELINE FROM LANGEN: https://github.com/rojaff/LanGen_pipeline
#C. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3


##2. INPUTS FOR THIS STEP:
#A. THE FILE ".VCF" CLEANED AFTER FILTERING, STEP 1.
#B. THE FILE "functions_LanGen.R" WITH FUNCTIONS DESIGNED FOR THIS PIPELINE IN THE WORKING DIRECTORY. YOU CAN DOWNLOAD IN  https://github.com/jdalapicolla/ OR https://github.com/rojaff/LanGen_pipeline


##3. GOALS FOR THIS STEP:
#A. ISOLATION-BY-DISTANCE ANALYSIS USING MANTEL TEST 

##4. CHOOSE A FOLDER FOR RUNNING THE ANALYSES. THE FILES MUST BE THERE! 
#A. IN RStudio GO TO  SESSION >> SET WORKING DIRECTORY >> CHOOSE DIRECTORY.. IN RStudio TOOL BAR OR USE THE SHORCUT CTRL+SHIFT+H

##5. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT:
rm(list=ls())


##6. LOAD THE FILE "functions_LanGen.R" WITH FUNCTIONS TO BE USED ON THIS STEP. MORE INFORMATION ON FUNCTIONS IN NUMBER 2.
source("functions_LanGen.R")


##7. INSTALL AND LOAD THE PACKAGES
#A. install the packages automatically
if("remotes" %in% rownames(installed.packages()) == FALSE){install.packages("remotes")
} else {print (paste0("'remotes' has already been installed in library"))}
if("BiocManager" %in% rownames(installed.packages()) == FALSE){install.packages("BiocManager")
} else {print (paste0("'BiocManager' has already been installed in library"))}
if("pacman" %in% rownames(installed.packages()) == FALSE){install.packages("pacman")
} else {print (paste0("'pacman' has already been installed in library"))}
if("devtools" %in% rownames(installed.packages()) == FALSE){install.packages("devtools")
} else {print (paste0("'devtools' has already been installed in library"))}

if("r2vcftools" %in% rownames(installed.packages()) == FALSE){remotes::install_github("nspope/r2vcftools")
} else {print (paste0("'r2vcftools' has already been installed in library"))}
if("geosphere" %in% rownames(installed.packages()) == FALSE){install.packages("geosphere")
} else {print (paste0("'geosphere' has already been installed in library"))}
if("raster" %in% rownames(installed.packages()) == FALSE){install.packages("raster")
} else {print (paste0("'raster' has already been installed in library"))}
if("rgdal" %in% rownames(installed.packages()) == FALSE){install.packages("rgdal")
} else {print (paste0("'rgdal' has already been installed in library"))}
if("ade4" %in% rownames(installed.packages()) == FALSE){install.packages("ade4")
} else {print (paste0("'ade4' has already been installed in library"))}
if("usedist" %in% rownames(installed.packages()) == FALSE){install.packages("usedist")
} else {print (paste0("'usedist' has already been installed in library"))}
if("reshape2" %in% rownames(installed.packages()) == FALSE){install.packages("reshape2")
} else {print (paste0("'reshape2' has already been installed in library"))}

#B. load packages multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load(r2vcftools, usedist, raster, rgdal, geosphere, ade4, reshape2)

##8. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS:
create_dir(c("./Results_IBD/Mantel_Tests"))



##### 1. Loading Files ----- 
###1.1. CHOOSE A NAME FOR THE PROJECT. MUST BE THE SAME ONE THAN OTHER STEPS:
#A. Project name:
project_name = "steerei"


###1.2. LOAD VCF FILE WITH GEOGRAPHICAL INFORMATION: 
#A. Load neutral .vcf file with geographical information:
snps_neutral = vcfLink(paste0("vcf/", project_name,"_filtered_neutral_LEA_DAPC_TESS.vcf"), overwriteID=T)
VCFsummary(snps_neutral) #19 individuals and 13971 SNPs.



##### 2. Preparing Genetic Data -----
###2.1. LOAD GENETIC DISTANCE:
#A. I used genetic distance based on PCA distance with broken stick rule

#B. Load distance file among all individuals calculate in STEP 4:
#PCA with BS rule of variance explained
PCA = read.csv(paste0("Results_Distance/PCA_Distance_BSR_IND_neutral_", project_name, ".csv"), row.names = 1)
PCA[1:5, 1:5]


###2.2. CREATE GENETIC 'DIST' OBJECTS FOR ALL SAMPLES: 
#A. All samples:
genDIST_all = as.dist(PCA)


###2.3. CREATE GENETIC 'DIST' OBJECTS BY CLUSTERS:
#A. Position of samples by populion by sNMF approach. Choose one method and change it on script:
for (i in 1:length(unique(snps_neutral@meta$PopID_snmf))){
  pop = which(snps_neutral@meta$PopID_snmf == i)
  assign(paste0("pop_SNMF_", i), pop)
}

#B. Subset By Pop 1:
genDIST_pop1= dist_subset(genDIST_all, snps_neutral@sample_id[pop_SNMF_1])
#verify dimensions
length(genDIST_pop1)

#C. Subset By Pop 2:
genDIST_pop2= dist_subset(genDIST_all, snps_neutral@sample_id[pop_SNMF_2])
#verify dimensions
length(genDIST_pop2)

#D. Subset By Pop 3:
genDIST_pop3= dist_subset(genDIST_all, snps_neutral@sample_id[pop_SNMF_3])
#verify dimensions
length(genDIST_pop3)



##### 3. Preparing Geographical Data -----
###3.1. LOAD GEOGRAPHIC INFORMATION:
#A. Create a data frame with the geographical coordenates:
coord_SN = snps_neutral@meta[,c(7,4:5)]
head(coord_SN)

#B. Set long and lat colunms
coordinates(coord_SN) = coord_SN[,c(2,3)]
projection(coord_SN) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84")

#C. Subset coord by populations:
coord_POP1 = coord_SN[pop_SNMF_1,]
coord_POP1
coord_POP2 = coord_SN[pop_SNMF_2,]
coord_POP2
coord_POP3 = coord_SN[pop_SNMF_3,]
coord_POP3


###3.2. CREATE MATRICES FOR GEOGRAPHICAL DISTANCE:
#by populations
coord_files = list(coord_POP1, coord_POP2, coord_POP3)

for (l in 1:length(coord_files)){
  coord = coord_files[[l]]
  sdist_SN  = distm(coord, fun=distGeo)
  sdist_SN  = sdist_SN /1000
  assign(paste0("sdist_POP", l), sdist_SN)
}

# all samples
sdist_SN  = distm(coord_SN, fun=distGeo)
rownames(sdist_SN) = snps_neutral@sample_id
colnames(sdist_SN) = snps_neutral@sample_id
sdist_SN
sdist_SN  = sdist_SN /1000
#verify the matrix
nrow(sdist_SN)
ncol(sdist_SN)
max(sdist_SN)
dim(sdist_SN)
class(sdist_SN)
sdist_SN


###3.3. CREATE 'DIST' OBJECTS FOR ALL SAMPLES AND BY CLUSTERS FROM MATRIXES: 
#A. All samples:
geoDIST_all = as.dist(sdist_SN)

#B. By Pop:
geoDIST_pop1 = as.dist(sdist_POP1)
geoDIST_pop2 = as.dist(sdist_POP2)
geoDIST_pop3 = as.dist(sdist_POP3)


###3.4. VERIFY DIMENSION FOR THE 'DIST' OBJECTS:
length(geoDIST_all)
length(genDIST_all)

length(geoDIST_pop1)
length(genDIST_pop1)

length(geoDIST_pop2)
length(genDIST_pop2)

length(geoDIST_pop3)
length(genDIST_pop3)



##### 4. Performing Mantel Tests -----
###4.1. PERFORM MANTEL TESTS
#A. By Species
mantel_all = mantel.rtest(genDIST_all, log(geoDIST_all), 10000)

#B. By populations
mantel_pop1 = mantel.rtest(genDIST_pop1, log(geoDIST_pop1), 10000)
mantel_pop2 = mantel.rtest(genDIST_pop2, log(geoDIST_pop2), 10000)
mantel_pop3 = mantel.rtest(genDIST_pop3, log(geoDIST_pop3), 10000)

#C. Compile results and save as data.frame
results_mantel = matrix(NA, 4, 2)
colnames(results_mantel) = c("r", "p-value")
rownames(results_mantel) = c("All Samples", "Only POP1", "Only POP2", "Only POP3")
results_mantel

results_mantel[1,] = cbind (mantel_all$obs, mantel_all$pvalue)
results_mantel[2,] = cbind (mantel_pop1$obs, mantel_pop1$pvalue)
results_mantel[3,] = cbind (mantel_pop2$obs, mantel_pop2$pvalue)
results_mantel[4,] = cbind (mantel_pop3$obs, mantel_pop3$pvalue)
results_mantel

write.csv(results_mantel, paste0("./Results_IBD/Mantel_Tests/novo_Mantel_Results_", project_name, ".csv"))


#D. Save Plot as PDF by population
#Prepare a single table with results
POP1_10 = as.data.frame(cbind(as.vector(geoDIST_pop1), as.vector(genDIST_pop1)))
POP1_10$POP = c(1)
POP2_15 = as.data.frame(cbind(as.vector(geoDIST_pop2), as.vector(genDIST_pop2)))
POP2_15$POP = c(2)
POP3_28 = as.data.frame(cbind(as.vector(geoDIST_pop3), as.vector(genDIST_pop3)))
POP3_28$POP = c(3)
#verify
df = rbind(POP1_10, POP2_15, POP3_28)
head(df)

#create the graph
graph_mantel = ggplot(data=df) +
  geom_point(mapping = aes(x=V1, y = V2, fill=factor(POP)), size=4, color = "black", shape=21) +
  scale_fill_manual(values = c("black", "grey", "white")) +
  geom_smooth(mapping = aes(x=V1, y = V2, group=POP, linetype = factor(POP)), method='lm', color = "black", size=1) +
  scale_linetype_manual(values = c("1" = 1, "2" = 2, "3" = 4))+
  theme_bw() +
  labs(x = "Geographic Distances (Km)", y = "PC Distance (Broken Stick Rule)") +
  guides(fill=guide_legend(title="Clusters"))

#verify
graph_mantel

#save graph
pdf(paste0("./Results_IBD/Mantel_Tests/Mantel_pcaBSR_", project_name,".pdf"))
graph_mantel
dev.off()


#D. Save Plot as PDF for all individuals together
#Preparing a single table
df_ALL = as.data.frame(cbind(as.vector(geoDIST_all), as.vector(genDIST_all)))
head(df_ALL)

mantel_graph_all = ggplot(data=df_ALL) +
  geom_point(mapping = aes(x=V1, y = V2), size=4, colour="black", pch=21) +
  geom_smooth(mapping = aes(x=V1, y = V2), method='lm', color = "black", size=1) +
  theme_bw() +
  labs(x = "Geographic Distances (Km)", y = "PC Distance (Broken Stick Rule)") +
  guides(fill=guide_legend(title=""))

#verify
mantel_graph_all

pdf("./Results_IBD/Mantel_Tests/Mantel_Plot_All.pdf", onefile=FALSE)
mantel_graph_all
dev.off()


###### 5. Permutation Tests by Population -----
###5.1. MANTEL TEST WITH GEOGRAPHIC DISTANCE MATRIX - PERMUTATION
#A. Define populations positions to subset in a list
populations = list(snps_neutral@sample_id[pop_TESS_1], snps_neutral@sample_id[pop_TESS_2], snps_neutral@sample_id[pop_TESS_3])

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
rownames(mantel_perm) = c("POP1", "POP2", "POP3")
mantel_perm

#E. Save the results:
write.csv(mantel_perm, paste0("./Results_IBD/Mantel_Tests/Mantel_Results_Permutation_", project_name, ".csv"))




###### 6. ANCOVA for testing IBD slopes -----
##6.1 PREPARING DATA FOR DIFFERENT POPULATIONS OR SPECIES:
#A. Population 1
sp1_gen = subset(melt(as.matrix(genDIST_pop1)), value!=0)
sp1_geo = subset(melt(as.matrix(geoDIST_pop1)), value!=0)
sp1_gen$geo = sp1_geo$value
sp1_ancova = sp1_gen[duplicated(sp1_gen$value),]
sp1_ancova = sp1_ancova[,-1]
sp1_ancova$Var2 = c("STE_1")
colnames(sp1_ancova) = c("SPE", "FST", "GEO")
sp1_ancova

#B. Population 2
sp2_gen = subset(melt(as.matrix(genDIST_pop2)), value!=0)
sp2_geo = subset(melt(as.matrix(geoDIST_pop2)), value!=0)
sp2_gen$geo = sp2_geo$value
sp2_ancova = sp2_gen[duplicated(sp2_gen$value),]
sp2_ancova = sp2_ancova[,-1]
sp2_ancova$Var2 = c("STE_2")
colnames(sp2_ancova) = c("SPE", "FST", "GEO")
sp2_ancova

#C. Population 3
sp3_gen = subset(melt(as.matrix(genDIST_pop3)), value!=0)
sp3_geo = subset(melt(as.matrix(geoDIST_pop3)), value!=0)
sp3_gen$geo = sp3_geo$value
sp3_ancova = sp3_gen[duplicated(sp3_gen$value),]
sp3_ancova = sp3_ancova[,-1]
sp3_ancova$Var2 = c("STE_3")
colnames(sp3_ancova) = c("SPE", "FST", "GEO")
sp3_ancova

ancova_pro = rbind(sp1_ancova, sp2_ancova, sp3_ancova)
rownames(ancova_pro)= NULL

#Genetic data is modeled as the dependent variable with SPE as the factor and GEO as the covariate
mod1 = aov(FST~GEO*SPE, data=ancova_pro)
summary(mod1)
capture.output(summary(mod1),file="Model1_interaction.doc")
#The summary of the results show a significant effect of GEO, SPE, and significant interaction. These results suggest that the slope of the regression between FST and GEO is different for all species.

mod2 <- aov(FST~GEO+SPE, data=ancova_pro)
summary(mod2)
capture.output(summary(mod2),file="Model2_Nointeraction.doc")
#The second model shows that SPE has a significant effect on the dependent variable which in this case can be interpreted as a significant difference in ‘intercepts’ between the regression lines of species. 

#assess if removing the interaction significantly affects the fit of the model:
anova(mod2,mod1)
capture.output(anova(mod1,mod2),file="Model1_2_Comparison.doc")
#The anova() command clearly shows that removing the interaction significantly affect the fit of the model (F=4.75, p=0.02). Therefore, we may conclude that the most parsimonious model is mod1 with more parameters.

##ANCOVA
#factor = species
#independent variable = geography
#dependent variable = fst


#Differences in intercepts are interpreted as differences in magnitude
#Differences in slopes are interpreted as differences in the rate of change. Lines should not be parallel. This means that growth is similar for both lines but one group is simply larger than the other. 
#A difference in slopes is interpreted as differences in the rate of change. In allometric studies, this means that there is a significant change in growth rates among groups.
#Slopes should be tested first, by testing for the interaction between the covariate and the factor. If slopes are significantly different between groups, then testing for different intercepts is somewhat inconsequential since it is very likely that the intercepts differ too 

##END
