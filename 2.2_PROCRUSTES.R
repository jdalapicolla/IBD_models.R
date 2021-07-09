####################### VALE INSTITUTE OF TECHNOLOGY ##########################

######################## LANDSCAPE GENOMICS TUTORIAL ##########################
#####################   STEP 05: IBD BY PROCRUSTES   #########################




### Script prepared by Jeronymo Dalapicolla, Luciana C. Resende-Moreira, Joyce R. do Prado ###



#### PRE-ANALYSIS #### 


##1. INPUTS FOR THIS TUTORIAL ----
#A. THE FILE ".VCF" CLEANED AFTER FILTERING AND WITH DELIMITED GENETIC CLUSTERS, STEP 2.

#B. THE FILE "functions_LanGen.R" WITH FUNCTIONS DESIGNED FOR THIS PIPELINE IN THE WORKING DIRECTORY. YOU CAN DOWNLOAD IN  https://github.com/jdalapicolla/ OR https://github.com/rojaff/LanGen_pipeline



##2. GOALS FOR THIS STEP:
#A. #A. ISOLATION-BY-DISTANCE ANALYSIS USING PROCRUSTES ANALYSES 



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

#From CRAN R:
if (!require('tidyverse'))    install.packages("tidyverse");         library('tidyverse') 
if (!require('raster'))       install.packages("raster");            library('raster')
if (!require('rgdal'))        install.packages("rgdal");             library('rgdal')
if (!require('vcfR'))         install.packages("vcfR");              library('vcfR')
if (!require('vegan'))        install.packages("vegan");             library('vegan')
if (!require('adegenet'))     install.packages("adegenet");          library('adegenet')
if (!require('ggplot2'))      install.packages("ggplot2");           library('ggplot2') 
if (!require('aspace'))       install.packages("aspace");            library('aspace') #R 3.5
if (!require('reshape2'))     install.packages("reshape2");          library('reshape2') 
if (!require('GISTools'))     install.packages("GISTools");          library('GISTools')
if (!require('sp'))           install.packages("sp");                library('sp')



#Load multiple packages using the package 'pacman'. If the package is missing "p_load" will download it from CRAN. "" in packages names is not mandatory.
pacman::p_load(r2vcftools, vcfR, adegenet, raster, rgdal, vegan, GISTools, tidyverse, reshape2, ggplot2, sp)



##7. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS ----
create_dir(c("./Results/Step05/Procrustes"))



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

#B. Load neutral .vcf file with geographical information:
snps_neutral = vcfLink(paste0("vcf/", project_name, "_filtered_neutral_clusters.vcf"), overwriteID=T)
VCFsummary(snps_neutral) #277 individuals and 5268 SNPs.

#C. Load the vcfR file for neutral SNPs
vcf_neutral = read.vcfR( paste0("vcf/", project_name, "_filtered_neutral_clusters.vcf"), verbose = FALSE)

#D. Convert to genind
genind_neutral = vcfR2genind(vcf_neutral)

#E. Number and name of cluster and method:
optimal_K = 4
name_clusters = c("A", "B", "C", "D")

#F.Position of samples by population/genetic clusters:
for (i in name_clusters){
  pop = which(snps_neutral@meta$POP_ID == i)
  assign(paste0("pop_POSI_", i), pop)
}

#G. Choose a extension for maps that cover all points. In this case:
long_extension = c(min(snps_neutral@meta$Longitude), max(snps_neutral@meta$Longitude))
lat_extension = c(min(snps_neutral@meta$Latitude), max(snps_neutral@meta$Latitude))
long_extension
lat_extension
ext = raster::extent(-50.7, -49.7, -6.5, -5.8)

#H. Load, crop and verify the CRS for all rasters and shapefiles:
#Landcover
landcover = raster::stack("./Inputs/rasters/NoSea_CL.tif") #multi band raster as satelite images
landcover = raster::crop(landcover, ext)
proj4string(landcover)

#landcover3
landcover2 = raster::stack("./Inputs/rasters/LC_Carajas_2018.tif")
landcover2 = raster::crop(landcover2, ext)
proj4string(landcover2)

#hillshade
hillshade = raster::stack("./Inputs/rasters/hillshade_carajas.tif")
hillshade = raster::crop(hillshade, ext)
proj4string(hillshade)

#if you need to fix the CRS base on a model, also a raster:
#raster_to_fix = spTransform(raster_to_fix, proj4string(raster_model))

#cangas
cangas = shapefile("./Inputs/shapefiles/Cangas.shp")
#if your shapefile is smaller than study area you do not need to crop. If you do it, the object will be NULL
#cangas = crop(cangas, ext) 
crs(cangas)
#if you need to fix the CRS base on a model, also a raster:
cangas = spTransform(cangas, proj4string(landcover2))
crs(cangas)








#### 2. PERFORM PROCRUSTES ANALYSES -----
#A. Perform a PCA
#choose the two first PC for view the results in a two dimension graph. You can choose different numbers of PCs to represent the variance! See STEP 2, Action #5, Operation #D for options. With more PCs less correlation among genetic and geography and the PCA pattern of two-dimension figure is not recovered. Be careful!
n_pcs = 2
input_scaled = scaleGen (genind_neutral, center = TRUE, scale = TRUE, NA.method = "mean")
pca_input = dudi.pca(input_scaled, center = TRUE, scannf = FALSE, nf = n_pcs)

#B. Lat and Long coordinates should be as matrix without row and cols names in decimal degrees 
coord_mat = as.matrix(snps_neutral@meta[, 6:7])
row.names(coord_mat) = NULL
colnames(coord_mat) = NULL

#C. Calculate the correlation between the matrices! Scale is for scaling linearly configuration Y for maximum similarity!
test_pro_1 = protest (X=coord_mat, Y=pca_input, scale=T, symmetric=F, choices = c(1:n_pcs), permutations=10000)
#view the results
test_pro_1
#p-value is "Significance"
#correlation after permutation is "Correlation in a symmetric Procrustes rotation"

#F. Save results in a csv file:
proc = cbind.data.frame (test_pro_1$ss, test_pro_1$scale, test_pro_1$signif, test_pro_1$permutations)
colnames(proc) = c("Procrustes Sum of Squares (m12 squared)", "Correlation in a symmetric Procrustes rotation", "p-value", "Number of Permutation")
rownames(proc) = project_name
proc
write.csv(proc,paste0("./Results/Step05/Procrustes/Procustes_statistics_", project_name, ".csv"))





install.packages("aspace")








#### 3. PERMUTATION FOR PROCRUSTES ANALYSES -----
#A. Define some variables: 
#list for samples position by pops
pops = list(pop_POSI_A, pop_POSI_B, pop_POSI_C, pop_POSI_D)
#vector with pop identifiers.
pop_vec = name_clusters

##B. Perform a Procrustes analysis on the total dataset first.
pro = protest(X=coord_mat,Y=pca_input,scale=T,symmetric=F,permutations=10000,choices = c(1:n_pcs))

#C. Define other variables
#asin_d(pro$rotation[2,1])
#Since the R default is to compute trigonometric functions on angular 
#measurements stored in radians,function asin_d performs the conversion 
#from degrees, reducing the need to do so a priori, outside the function.
total = pro$scale
rotationallpops = asin_d(pro$rotation[2,1])
t_prime = NULL
t_prime2 = NULL
rotation = NULL
rotation2 = NULL
signif = NULL
difference = NULL
max_diff = 0

#D. Removing one pop by loop
for(m in 1:length(pops)) {
  ## delete from coordinate matrix those related to each of the eliminated pop
  temp_mat = coord_mat[-(pops[[m]]),]
  ## delete from the matrix of mean scaled genotypes those related to each of the eliminated pop
  red = input_scaled[-(pops[[m]]),]
  ## does pca with reduced matrix and protest with the new dataset (without each eliminated pop)
  pca1 = dudi.pca(red, center = T, scale = T, scannf=F, nf = n_pcs)
  pro = protest(X=temp_mat,Y=pca1$li,scale=T,symmetric=F,permutations=10000,choices = c(1:n_pcs))
  t_prime2 = c(t_prime2,pro$scale)
  rotation2 = c(rotation2, asin_d(pro$rotation[2,1]))
  signif = c(signif, pro$signif)
  difference = c(difference,(pro$scale-total) )
  
  ## delete from pca$li matrix the results related to each of the eliminated pop
  pca_original = pca_input$li[-(pops[[m]]),]
  pro = protest(X=pca_original,Y=pca1$li,scale=T,symmetric=F,permutations=10000, choices = c(1:n_pcs))
  t_prime = c(t_prime,pro$scale)
  rotation = c(rotation, asin_d(pro$rotation[2,1]))
}

#E. Save results as figure
#xlim is number of pops
xlim=c(0, length(pops))
#order difference values to use as ylim
diff_order = sort(difference) #-0.066 to 0.11
diff_order
ylim = c((diff_order[1]-0.05),(diff_order[length(pops)]+0.05))

#set populations names and colors:
col_tess
legend_tess

#plot the results
pdf(paste0("./Results/Step05/Procrustes/Procrustes_Permutation_", project_name, ".pdf"), onefile = F)
#set margins and axes
plot.new()
par(oma=c(1,1,1,1))
par(mar=c(2,2,1,0))
plot.window(xlim=xlim, ylim=ylim)
#plot the graph
barplot (difference, ylim = ylim  , main = "Relative difference in t0", col=col)
legend (x="topright", ncol=2, pch=21, bty='n', cex= 0.8, pt.cex=2,legend = legend_tess, col = "black", pt.bg = col_tess)
abline(h=0)
dev.off()


#F. Save results as table:
#create a data frame with all relevant information and saves it in a .csv file
test_data = data.frame(t0 = total, Angle = rotationallpops,POP = pop_vec, GEO = t_prime2, DIFF = difference, ROT_GEO = rotation2, P = signif, PCA = t_prime, ROT_PCA = rotation)

write.csv(test_data, file=paste0("./Results/Step05/Procrustes/Procrustes_Permutation_table_", project_name, ".csv"))







#### 4. SAVE A MAP WITH PROCRUSTES -----
#A. Use this to get coordinates for the PCA value:
test_pro_2 = procrustes (X=coord_mat, Y=pca_input, scale=T, symmetric=F, choices = c(1:n_pcs))

#B. View the arrows
plot (test_pro_2, main=NULL, ar.col='black', xlab='', ylab='')

#C. Open pdf, choose the font letters, margins and axes.
pdf(paste0("./Results_IBD/Procrustes_Analysis/Map_TESS_Procrustes_Map_", project_name, ".pdf"), family="Helvetica", onefile = F)
#set parameters
plot.new()
# Add margins
par(oma=c(2,2,2,2))
# Add axes
plot.window(xlim= axis.x, ylim= axis.y)

#D. Plot maps for base
plotRGB(hillshade, r = 1, g = 2, b = 3, interpolate = T , ext= ext, alpha = 255, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r') #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plotRGB(landcover2,r=1, g=2, b=3, interpolate=T, ext= ext, alpha =120, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r', add=T) #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plotRGB(landcover3,r=1, g=2, b=3, interpolate=T, ext= ext, alpha =100, bgalpha =0, stretch="lin", bty="7", axes = F, yaxs='r', xaxs='r', add=T) #stretch "lin" or "hist" are contrast and alpha (0 to 255) is brightness
plot(cangas, add=TRUE, axes = F, bg = FALSE, bty ="n", xaxt='n', ann=FALSE, yaxt='n', lwd=1, col = "beige")

#E. Plot long and lat axes:
axis(1, at=seq(-50.7, -49.7, by=0.2),  cex.axis=1, pos=-6.5)
axis(2, at=seq(-6.5, -5.8, by=0.2),  cex.axis=1, las=1, pos=-50.7)

#F. Set population names and colors. SAME THAN STEP 2!
col_tess= c("yellow", 'darkslategrey', 'red', 'blue')
legend_tess = c("POP1", "POP2", "POP3", "POP4")

#G. Plot arrows
arrows(test_pro_2$Y[,1]+test_pro_2$translation[1,1], test_pro_2$Y[,2]+test_pro_2$translation[1,2], test_pro_2$X[,1]+test_pro_2$translation[1,1], test_pro_2$X[,2]+test_pro_2$translation[1,2], length=0.05, code = 0, cex=0.1, lwd=2) #code is the style of arrow, 0 = without head

#H. Plot genetic position of population using DAPC groups.Add or remove lines according to your necessity.
points(test_pro_2$Y[pop_TESS_1,1]+test_pro_2$translation[1,1],test_pro_2$Y[pop_TESS_1,2]+test_pro_2$translation[1,2], col = 'black', bg = col_tess[1], cex = 1.5, pch=21) #POP1
points(test_pro_2$Y[pop_TESS_2,1]+test_pro_2$translation[1,1],test_pro_2$Y[pop_TESS_2,2]+test_pro_2$translation[1,2], col = 'black', bg = col_tess[2], cex = 1.5, pch=21) #POP2
points(test_pro_2$Y[pop_TESS_3,1]+test_pro_2$translation[1,1],test_pro_2$Y[pop_TESS_3,2]+test_pro_2$translation[1,2], col = 'black', bg = col_tess[3], cex = 1.5, pch=21) #POP3
points(test_pro_2$Y[pop_TESS_4,1]+test_pro_2$translation[1,1],test_pro_2$Y[pop_TESS_4,2]+test_pro_2$translation[1,2], col = 'black', bg = col_tess[4], cex = 1.5, pch=21) #POP4

#I. Plot geographical position for the same population on previous operation #H
points(test_pro_2$X[pop_TESS_1,1]+test_pro_2$translation[1,1],test_pro_2$X[pop_TESS_1,2]+test_pro_2$translation[1,2], col = 'black', bg = col_tess[1], cex = 1.5, pch=24)  #POP1
points(test_pro_2$X[pop_TESS_2,1]+test_pro_2$translation[1,1],test_pro_2$X[pop_TESS_2,2]+test_pro_2$translation[1,2], col = 'black', bg = col_tess[2], cex = 1.5, pch=24)  #POP2
points(test_pro_2$X[pop_TESS_3,1]+test_pro_2$translation[1,1],test_pro_2$X[pop_TESS_3,2]+test_pro_2$translation[1,2], col = 'black', bg = col_tess[3], cex = 1.5, pch=24)  #POP3
points(test_pro_2$X[pop_TESS_4,1]+test_pro_2$translation[1,1],test_pro_2$X[pop_TESS_4,2]+test_pro_2$translation[1,2], col = 'black', bg = col_tess[4], cex = 1.5, pch=24)  #POP4

#J. Plot the legend:
legend(x=-49.85, y=-5.85, legend = legend_tess, cex = 0.8, bty="o", bg = "white", box.col="black",  col = "black", pt.bg = col_tess, pch= 22, pt.cex=1.5)

#K. Plot the scale bar:
scalebar(10, xy=c(-50.65,-6.45), type="bar", lonlat = T, below = "Km", label = c(0,5,10), adj=c(0,-1.5), cex = 0.7)

#L. Plot a north arrow:
north.arrow(xb=-49.8, yb=-6.45, len=0.015, cex.lab=0.8, lab="N", col='black', font.sub=1)

#M. Turn off pdf:
dev.off()



##END
