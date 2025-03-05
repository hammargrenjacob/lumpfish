rm(list = ls()) #clear environment
#-----PACKAGES-----####

## Packages (I just included all that we used for the r practicals)
# install packages
install.packages("radiator")
install.packages("hierfstat")
install.packages("pegas")
install.packages("adegenet")
install.packages("tidyverse")
install.packages("assignPOP")
install.packages("ggplot2")
install.packages("rlang")
install.packages("usethis")
install.packages("ade4")
# install devtools
install.packages("devtools")
# install diveRsity from github
devtools::install_github("kkeenan02/diveRsity")
devtools::install_github("thierrygosselin/radiator")
devtools::install_github('wrengels/HWxtest', 
                         subdir='pkg')
install.packages("vcfR")
install.packages("shapefiles")
install.packages("mapplots")
install.packages("tibble")
install.packages("dplyr")
install.packages("devtools")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("LEA")
install.packages("quantreg")
install.packages("forcats")
install.packages("ggthemes")
install.packages("gridGraphics")
install.packages("gridExtra")
install.packages("msa")
install.packages("phangorn")
install.packages("seqinr")
install.packages("BiocManager")
BiocManager::install("msa")
install.packages("grDevices")
install.packages("maps")
install.packages("genesis")
install.packages("readxl")

# load libraries
library("hierfstat")
library("pegas")
library("pegas")
library("adegenet")
library("diveRsity")
library("tidyverse")
library("assignPOP")
library("ggplot2")
library("rlang")
library("usethis")
library("ade4")
library("radiator")
library("vcfR")
library("shapefiles")
library("mapplots")
library("tibble")
library("dplyr")
library("devtools")
library("LEA")
library("quantreg")
library("forcats")
library("ggthemes")
library("gridGraphics")
library("gridExtra")
library("cowplot")
library("msa")
library("ape")
library("phangorn")
library("seqinr")
library("tidyverse")
library("grDevices")
library("maps")
library("readxl")









#-----READ DATA------
vcf_lump <- read.vcfR("Lumpfish_2bRAD.vcf")
metadata_lump <- read_xlsx("Lumpfish_metadata.xlsx") %>% 
  filter(used_in_2bRAD == "yes")

genind_lump_broad <- read.genepop("Lumpfish_SNP.gen", ncode = 3)
metadata_lump_broad <- read_xlsx("Lumpfish_metadata.xlsx") 


#1. 2bRAD dataset-----
##----Merging metadata to genind object----
str(genind_lump)
rownames(genind_lump@tab)  
metadata_lump <- metadata_lump[match(rownames(genind_lump@tab), 
                                     metadata_lump$sample_ID_2bRAD), ]
genind_lump@other$metadata <- metadata_lump

head(genind_lump@other$metadata$Location) #sampling site full name
head(genind_lump@other$metadata$sample_ID_SNPdata) #site and id abbreviated
head(genind_lump@other$metadata$Pop_ID_article) #Population and country

##----PCA----
?dudi.pca
x_lump <- tab(genind_lump, freq=TRUE, NA.method="mean")
pca_lump <- dudi.pca(x_lump, center=TRUE, scale=FALSE) # 3 axes
pca_lump

s.class(pca_lump$li, 
        fac = as.factor(genind_lump@other$metadata$Pop_ID_article),  
        clab = 1, 
        col = transp(funky(length(unique(genind_lump@other$metadata$Pop_ID_article)))), 
        csta = 0, 
        cpoint = 4, 
        cellipse = 1, 
        xax = 1, 
        yax = 2)

eig.perc <- 100*pca_lump$eig/sum(pca_lump$eig)
head(eig.perc) #11.052981  3.445938  2.747533 % for the first 3

##-----DAPC------

length(unique(genind_lump@other$metadata$Pop_ID_article)) #How many sampling sites were there?
kmeans_lump <- find.clusters(genind_lump, 
                             n.pca = 50, 
                             max.n = 10,
                             scale = FALSE) #4 clusters
kmeans_lump$grp %>% as.vector


dapc_lump <- dapc(genind_lump, pop = kmeans_lump$grp, n.pca = 50)

# Plotting the first 2 linear discrimination functions
scatter(dapc_lump, xax=1, yax=2, grp=dapc_lump$grp, 
        col=transp(c("deepskyblue3","goldenrod2","mediumorchid3","darkolivegreen3")),
        pch=19, bg="white",cstar = 1, cellipse = 1,clabel = 1,scree.da=FALSE,
        scree.pca=FALSE)
dapc_lump_plot <- recordPlot()

#and PCA again with the same color assignment
s.class(pca_lump$li, fac=dapc_lump$grp, clab=1,
        col=transp(c("deepskyblue3","goldenrod2","mediumorchid3","darkolivegreen3")), 
        csta=1, cpoint=2, cellipse =1, xax=1, yax=2)
pca_lump_plot <- recordPlot()

##----Map----
#Data frame with pop and cluster assignment for each individual
pop_cluster_lump <- data.frame(Population=genind_lump@other$metadata$Pop_ID_article, 
                               Group=kmeans_lump$grp)
count_cluster_per_pop <- pop_cluster_lump %>% 
                         dplyr::count(Population,Group) %>% 
                         data.frame()
head(count_cluster_per_pop)

#Location coordinates
#Ensuring that the coordinate columns contain numeric values
metadata_lump$Lat <- as.numeric(metadata_lump$Lat)
metadata_lump$Lon <- as.numeric(metadata_lump$Lon)

location_coords <- metadata_lump %>%
  dplyr::group_by(Pop_ID_article) %>%
  dplyr::summarise(
    Lat = mean(Lat, na.rm = TRUE),
    Lon = mean(Lon, na.rm = TRUE)
  ) %>%
  dplyr::rename(Population = Pop_ID_article) %>% 
  as.data.frame() 
head(location_coords)


#Merging objects based on the population column
cluster_pop_coords <- inner_join(location_coords, count_cluster_per_pop, by="Population")
head(cluster_pop_coords)

xyz <- make.xyz(cluster_pop_coords$Lon,
                cluster_pop_coords$Lat,
                cluster_pop_coords$n,
                cluster_pop_coords$Group)

#Map
shape <- read.shapefile("world_map")
str(shape)

?basemap
#1. Map size
basemap(xlim=c(-70,20), #min/max longitude
        ylim=c(40,75), #min/max latitude
        # yaxp = c(40,75,5),
        bg="white",
        xlab = "",
        ylab = "",
        frame.plot = FALSE,
        axes = FALSE) 
axis(2, at = seq(35, 85, by = 5), cex.axis = 0.6)  
axis(1, at = seq(-80, 40, by = 10), cex.axis = 0.6)



#2. Load shapefile
map <- draw.shape(shape, col="grey85")

#3. Pie charts and labels
draw.pie(xyz$x, xyz$y, xyz$z, radius=1.5,
         col=c("deepskyblue3","goldenrod2","mediumorchid3","darkolivegreen3"),
        )

text(cluster_pop_coords$Lon, cluster_pop_coords$Lat - 2,  
     labels = cluster_pop_coords$Population,
     col = "black", cex = 0.8, font = 2.5)

map_plot_lump <- recordPlot()




#2. SNP panel dataset-----
##----Merging metadata to genind object----
str(genind_lump_broad)
rownames(genind_lump_broad@tab)  
metadata_lump_broad <- 
  metadata_lump_broad[match(rownames(genind_lump_broad@tab), 
                                                 metadata_lump_broad$sample_ID_SNPdata), ]
genind_lump_broad@other$metadata <- metadata_lump_broad

head(genind_lump_broad@other$metadata$Location) #sampling site full name
head(genind_lump_broad@other$metadata$sample_ID_SNPdata) #site and id abbreviated
head(genind_lump_broad@other$metadata$Pop_ID_article) #Population and country
head(genind_lump_broad@pop)


##----PCA----
?dudi.pca
x_lump_broad <- tab(genind_lump_broad, freq=TRUE, NA.method="mean")
pca_lump_broad <- dudi.pca(x_lump_broad, center=TRUE, scale=FALSE) # 5 axes
pca_lump_broad

s.class(pca_lump_broad$li, 
        fac = as.factor(genind_lump_broad@other$metadata$Pop_ID_article),  
        clab = 1, 
        col = transp(funky(length(unique(
          genind_lump_broad@other$metadata$Pop_ID_article)))), 
        csta = 0, 
        cpoint = 4, 
        cellipse = 1, 
        xax = 1, 
        yax = 2)

s.class(pca_lump_broad$li, 
        fac = as.factor(genind_lump_broad@pop),  
        clab = 1,
        col = transp(funky(length(unique(
          genind_lump_broad@other$metadata$Pop_ID_article)))), 
        csta = 0, 
        cpoint = 4, 
        cellipse = 1, 
        xax = 1, 
        yax = 2)


unique(genind_lump_broad@pop)

eig.perc.broad <- 100*pca_lump_broad$eig/sum(pca_lump_broad$eig)
head(eig.perc.broad) #32.473924 12.555573  5.834492  2.581348 %

##-----DAPC------

length(unique(genind_lump_broad@other$metadata$Pop_ID_article)) #How many sampling sites were there?
kmeans_lump_broad <- find.clusters(genind_lump_broad, 
                                   n.pca = 50, 
                                   max.n = 44,
                                   scale = FALSE) #10 clusters
kmeans_lump_broad$grp %>% as.vector

dapc_lump_broad <- dapc(genind_lump_broad, 
                        pop = kmeans_lump_broad$grp, 
                        n.pca = 50) #4 LDFs

# Plotting the first 2 linear discrimination functions
scatter(dapc_lump_broad, xax=1, yax=2, grp=dapc_lump_broad$grp, 
        col=transp(c("lightblue", "goldenrod2", "mediumorchid3", 
                     "darkolivegreen3", "firebrick2",
                     "deepskyblue4", "chocolate2", "chartreuse4", 
                     "tomato3", "purple3")),
        pch=19, bg="white",cstar = 1, cellipse = 1,clabel = 1,scree.da=FALSE,
        scree.pca=FALSE)
dapc_lump_broad_plot <- recordPlot()

#and PCA again with the same color assignment
s.class(pca_lump_broad$li, fac=dapc_lump_broad$grp, clab=1,
        col=transp(c("lightblue", "goldenrod2", "mediumorchid3", 
                     "darkolivegreen3", "firebrick2",
                     "deepskyblue4", "chocolate2", "chartreuse4", 
                     "tomato3", "purple3")), 
        csta=1, cpoint=2, cellipse =1, xax=1, yax=2)
pca_lump_broad_plot <- recordPlot()

##----Map----
#Data frame with pop and cluster assignment for each individual
pop_cluster_lump_broad <- data.frame(Population=genind_lump_broad@other$metadata$Pop_ID_article, 
                                     Group=kmeans_lump_broad$grp)
count_cluster_per_pop_broad <- pop_cluster_lump_broad %>% 
  dplyr::count(Population,Group) %>% 
  data.frame()
head(count_cluster_per_pop_broad)

#Location coordinates
##ensuring that the coordinate columns contain numeric values
metadata_lump_broad$Lat <- as.numeric(metadata_lump_broad$Lat)
metadata_lump_broad$Lon <- as.numeric(metadata_lump_broad$Lon)

location_coords_broad <- metadata_lump_broad %>%
  dplyr::group_by(Pop_ID_article) %>%
  dplyr::summarise(
    Lat = mean(Lat, na.rm = TRUE),
    Lon = mean(Lon, na.rm = TRUE)
  ) %>%
  dplyr::rename(Population = Pop_ID_article) %>% 
  as.data.frame() 
head(location_coords_broad)


#Merging objects based on the population column
cluster_pop_coords_broad <- inner_join(location_coords_broad, count_cluster_per_pop_broad, by="Population")
head(cluster_pop_coords_broad)

xyz_broad <- make.xyz(cluster_pop_coords_broad$Lon,
                      cluster_pop_coords_broad$Lat,
                      cluster_pop_coords_broad$n,
                      cluster_pop_coords_broad$Group)

#Map
shape <- read.shapefile("world_map")

?basemap

#1. Map size
basemap(xlim=c(-45,5), #min/max longitude
        ylim=c(45,80), #min/max latitude
        # yaxp = c(40,75,5),
        bg="white",
        xlab = "",
        ylab = "",
        frame.plot = FALSE,
        axes = FALSE) 
axis(2, at = seq(35, 85, by = 5), cex.axis = 0.6)  
axis(1, at = seq(-80, 40, by = 20), cex.axis = 0.6)


#2. Load shapefile
map <- draw.shape(shape, col="grey85")

#3. Pie charts and labels
draw.pie(xyz_broad$x, xyz_broad$y, xyz_broad$z, radius=1.5,
         col=c("lightblue", "goldenrod2", "mediumorchid3", 
               "darkolivegreen3", "firebrick2",
               "deepskyblue4", "chocolate2", "chartreuse4", 
               "tomato3", "purple3"),
)

#text(cluster_pop_coords_broad$Lon, cluster_pop_coords_broad$Lat - 2,  
 #    labels = cluster_pop_coords_broad$Population,
  #   col = "black", cex = 0.8, font = 2.5)

map_plot_lump_broad <- recordPlot()


## sNMF----

##Converting to .geno format
hierfstat_lump_broad <- genind2hierfstat(genind_lump_broad)

write.struct(hierfstat_lump_broad, 
             ilab=indNames(genind_lump_broad), 
             pop = genind_lump_broad@other$metadata$Pop_ID_article, 
             fname = "lump_broad.str")

struct2geno(input.file = "lump_broad.str", 
            ploidy = 2, 
            FORMAT = 2, 
            extra.column = 2, 
            extra.row = 0) 

genofile_lump_broad <- "lump_broad.str.geno"
cic.geno_lump_broad <- read.geno(genofile_lump_broad)


#How many ancestral pops are there? (TAKES A BIT OF TIME)  
length(unique(genind_lump_broad@other$metadata$Pop_ID_article)) # Number of pops
snmf_search_lump_broad <- snmf(genofile_lump_broad, 
                               K = 2:44, #Using max number = number of pops
                               entropy = TRUE, 
                               repetitions = 2, 
                               project = "new")

# Check the plot, should converge with DAPC, in this case doesnt though
plot(snmf_search_lump_broad, col = transp("steelblue4"), pch = 19) # 10 ancestral pops

#calculating ancestry of each individual for K=10
snmf_lump_broad <- snmf(genofile_lump_broad, K=10, project="new")

#probability matrix for K=10
qmatrix_lump_broad <- Q(snmf_lump_broad, K=10)

#create a new qmatrix and add some information (ID, population, and K-means cluster assignment)
qmatrix2_lump_broad <- data.frame(ID=rownames(pop_cluster_lump_broad),
                                  pop_cluster_lump_broad,
                                  qmatrix_lump_broad)

qmatrix_new_list_lump_broad <- list() 

#for-loop to create the new columns:
for(i in 4:ncol(qmatrix2_lump_broad)){ 
  qmatrix_new_list_lump_broad[[i-3]] <- qmatrix2_lump_broad[,c(1:3,i)] %>% 
    mutate("Var"=rep(colnames(qmatrix2_lump_broad)[i],nrow(qmatrix2_lump_broad)))
  colnames(qmatrix_new_list_lump_broad[[i-3]])  <- c("ID","Population","Kmeans_cluster","Prob","Variable")
}

#transform the list into the new qmatrix
qmatrix3_lump_broad <- bind_rows(qmatrix_new_list_lump_broad)



#use a few reference individuals with certain ancestry to change the cluster names so they are consistent between the DAPC, map and barplot
for(i in c("NOR17", "ISF21", "BRE15", "BOK5", "NEW12", "BAL2_14", "BOL9", "NUU2_16", "SVAN76", "SKA50" 
)){
  max_anc_lump_broad <- filter(qmatrix3_lump_broad,ID==i) %>% filter(Prob==max(Prob))
  qmatrix3_lump_broad$Variable <- 
    replace(qmatrix3_lump_broad$Variable,
            qmatrix3_lump_broad$Variable==max_anc_lump_broad$Variable,
            max_anc_lump_broad$Kmeans_cluster)}

?replace 
qmatrix3_lump_broad%>% 
  filter(Kmeans_cluster==1)

unique(qmatrix3_lump_broad$Variable)

#plot the ancestry coefficients in a nicer-looking barplot
anc_plot_lump_broad <- ggplot(qmatrix3_lump_broad, 
                              aes(factor(ID), Prob, fill = factor(Variable))) +
  geom_col(width=1) +
  facet_grid(~fct_inorder(as.factor(Population)), switch = "x", 
             scales = "free", 
             space = "free") +
  theme_minimal() + 
  labs(x = "", title = , y = "Probability") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.05, "lines"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size=5, angle = 90)) +
  scale_fill_manual(name="Cluster",
                    values=c("lightblue", "goldenrod2", "mediumorchid3", 
                             "darkolivegreen3", "firebrick2", "deepskyblue4", 
                             "chocolate2", "chartreuse4", "tomato3", "purple3"))
anc_plot_lump_broad

grid.arrange(arrangeGrob(as_grob(dapc_lump_broad_plot),
                         as_grob(map_plot_lump_broad),
                         ncol=2,clip="on"),
             anc_plot_lump_broad,nrow = 2,heights=c(4,2))
combined_map_lump_broad <- recordPlot()


#----3. Subset: Scandinavia-----
##----Subset----
#Read and filter
metadata_lump_broad_ss1 <- read_xlsx("Lumpfish_metadata.xlsx") %>%
  filter(grepl("_(SE|NO|IS|DK|BS)$", Pop_ID_article))

# extracting metadata info
pop_ids_ss1 <- genind_lump_broad@other$metadata$Pop_ID_article

# identifying which individuals to keep
keep_pops_ss1 <- grepl("_SE|_NO|_DK|_IS|_BS", pop_ids_ss1)

# creating genind subset
genind_lump_broad_ss1 <- genind_lump_broad[keep_pops_ss1, ]
unique(genind_lump_broad_ss1@other$metadata$Pop_ID_article)

##----PCA----
?dudi.pca
x_lump_broad_ss1 <- tab(genind_lump_broad_ss1, freq=TRUE, NA.method="mean")
pca_lump_broad_ss1 <- dudi.pca(x_lump_broad_ss1, center=TRUE, scale=FALSE) # 5 axes
pca_lump_broad_ss1

s.class(pca_lump_broad_ss1$li, 
        fac = as.factor(genind_lump_broad_ss1@other$metadata$Pop_ID_article),  
        clab = 1, 
        col = transp(funky(length(unique(
          genind_lump_broad_ss1@other$metadata$Pop_ID_article)))), 
        csta = 0, 
        cpoint = 4, 
        cellipse = 1, 
        xax = 1, 
        yax = 2)

eig.perc.broad_ss1 <- 100*pca_lump_broad_ss1$eig/sum(pca_lump_broad_ss1$eig)
head(eig.perc.broad_ss1) #31.1318980 11.9283696  2.6121774  1.6230219 %

##-----DAPC------
#How many sampling sites were there?
length(unique(genind_lump_broad_ss1@other$metadata$Pop_ID_article)) 

kmeans_lump_broad_ss1 <- find.clusters(genind_lump_broad_ss1, 
                                       n.pca = 50, 
                                       max.n = 28,
                                       scale = FALSE) #7 clusters
kmeans_lump_broad_ss1$grp %>% as.vector

dapc_lump_broad_ss1 <- dapc(genind_lump_broad_ss1, 
                            pop = kmeans_lump_broad_ss1$grp, 
                            n.pca = 50) # 2 DFs retained

# Plotting the first 2 linear discrimination functions
scatter(dapc_lump_broad_ss1, xax=1, yax=2, grp=dapc_lump_broad_ss1$grp, 
        col=transp(c("lightblue", "goldenrod2", "mediumorchid3", 
                     "darkolivegreen3", "firebrick2", "deepskyblue4", 
                     "chocolate2", "chartreuse4", "tomato3", "purple3")),
        pch=19, bg="white",cstar = 1, cellipse = 1,clabel = 1,scree.da=FALSE,
        scree.pca=FALSE)
dapc_lump_broad_ss1_plot <- recordPlot()

#and PCA again with the same color assignment
s.class(pca_lump_broad_ss1$li, fac=dapc_lump_broad_ss1$grp, clab=1,
        col=transp(c("lightblue", "goldenrod2", "mediumorchid3", 
                     "darkolivegreen3", "firebrick2", "deepskyblue4", 
                     "chocolate2", "chartreuse4", "tomato3", "purple3")), 
        csta=1, cpoint=2, cellipse =1, xax=1, yax=2)
pca_lump_broad_ss1_plot <- recordPlot()

##----Map----
#Data frame with pop and cluster assignment for each individual
pop_cluster_lump_broad_ss1 <- 
  data.frame(Population=genind_lump_broad_ss1@other$metadata$Pop_ID_article, 
             Group=kmeans_lump_broad_ss1$grp)

count_cluster_per_pop_broad_ss1 <- pop_cluster_lump_broad_ss1 %>% 
  dplyr::count(Population,Group) %>% 
  data.frame()
head(count_cluster_per_pop_broad_ss1)

#Location coordinates
##ensuring that the coordinate columns contain numeric values
metadata_lump_broad_ss1$Lat <- as.numeric(metadata_lump_broad_ss1$Lat)
metadata_lump_broad_ss1$Lon <- as.numeric(metadata_lump_broad_ss1$Lon)

#calculating the average coordinates of each population
location_coords_broad_ss1 <- metadata_lump_broad_ss1 %>%
  dplyr::group_by(Pop_ID_article) %>%
  dplyr::summarise(
    Lat = mean(Lat, na.rm = TRUE),
    Lon = mean(Lon, na.rm = TRUE)
  ) %>%
  dplyr::rename(Population = Pop_ID_article) %>% 
  as.data.frame() 
head(location_coords_broad_ss1)


#Merging objects based on the population column
cluster_pop_coords_broad_ss1 <- inner_join(location_coords_broad_ss1, count_cluster_per_pop_broad_ss1, by="Population")
head(cluster_pop_coords_broad_ss1)

xyz_broad_ss1 <- make.xyz(cluster_pop_coords_broad_ss1$Lon,
                          cluster_pop_coords_broad_ss1$Lat,
                          cluster_pop_coords_broad_ss1$n,
                          cluster_pop_coords_broad_ss1$Group)

#Map
shape <- read.shapefile("world_map")

?basemap
#1. Map size
basemap(xlim=c(-5,0), #min/max longitude
        ylim=c(54,71), #min/max latitude
        # yaxp = c(40,75,5),
        bg="white",
        xlab = "",
        ylab = "",
        frame.plot = FALSE,
        axes = FALSE) 
axis(2, at = seq(35, 85, by = 5), cex.axis = 0.6)  
axis(1, at = seq(-80, 40, by = 10), cex.axis = 0.6)


#2. Load shapefile
map <- draw.shape(shape, col="grey85")

#3. Pie charts and labels
draw.pie(xyz_broad_ss1$x, xyz_broad_ss1$y, xyz_broad_ss1$z, radius=1,
         col=c("lightblue", "goldenrod2", "mediumorchid3", 
               "darkolivegreen3", "firebrick2", "deepskyblue4", 
               "chocolate2", "chartreuse4", "tomato3", "purple3"),
)

#text(cluster_pop_coords_broad_ss1$Lon, cluster_pop_coords_broad_ss1$Lat - 2,  
 #    labels = cluster_pop_coords_broad_ss1$Population,
  #   col = "black", cex = 0.8, font = 0.5)

map_plot_lump_broad_ss1 <- recordPlot()

## sNMF----

##Converting to .geno format
hierfstat_lump_broad_ss1 <- genind2hierfstat(genind_lump_broad_ss1)

write.struct(hierfstat_lump_broad_ss1, 
             ilab=indNames(genind_lump_broad_ss1), 
             pop = genind_lump_broad_ss1@other$metadata$Pop_ID_article, 
             fname = "lump_broad_ss1.str")

struct2geno(input.file = "lump_broad_ss1.str", 
            ploidy = 2, 
            FORMAT = 2, 
            extra.column = 2, 
            extra.row = 0) 

genofile_lump_broad_ss1 <- "lump_broad_ss1.str.geno"
cic.geno_lump_broad_ss1 <- read.geno(genofile_lump_broad_ss1)


#How many ancestral pops are there? (TAKES A BIT OF TIME)  
length(unique(genind_lump_broad_ss1@other$metadata$Pop_ID_article)) # Number of pops
snmf_search_lump_broad_ss1 <- snmf(genofile_lump_broad_ss1, 
                                   K = 2:29, #Using max number = number of pops
                                   entropy = TRUE, 
                                   repetitions = 2, 
                                   project = "new")


# Check the plot, should converge with DAPC, in this case doesnt though
plot(snmf_search_lump_broad_ss1, col = transp("steelblue4"), pch = 19) # 7 ancestral pops

#calculating ancestry of each individual for K=6
snmf_lump_broad_ss1 <- snmf(genofile_lump_broad_ss1, K=7, project="new")

#probability matrix for K=10
qmatrix_lump_broad_ss1 <- Q(snmf_lump_broad_ss1, K=7)

#create a new qmatrix and add some information (ID, population, and K-means cluster assignment)
qmatrix2_lump_broad_ss1 <- data.frame(ID=rownames(pop_cluster_lump_broad_ss1),
                                      pop_cluster_lump_broad_ss1,
                                      qmatrix_lump_broad_ss1)

qmatrix_new_list_lump_broad_ss1 <- list() 

#for-loop to create the new columns:
for(i in 4:ncol(qmatrix2_lump_broad_ss1)){ 
  qmatrix_new_list_lump_broad_ss1[[i-3]] <- qmatrix2_lump_broad_ss1[,c(1:3,i)] %>% 
    mutate("Var"=rep(colnames(qmatrix2_lump_broad_ss1)[i],nrow(qmatrix2_lump_broad_ss1)))
  colnames(qmatrix_new_list_lump_broad_ss1[[i-3]])  <- c("ID","Population","Kmeans_cluster","Prob","Variable")
}

#transform the list into the new qmatrix
qmatrix3_lump_broad_ss1 <- bind_rows(qmatrix_new_list_lump_broad_ss1)


#data.frame(select(metadata_lump_broad_ss1, 
#                 c(Pop_ID_article, sample_ID_SNPdata)), 
#         kmeans_lump_broad_ss1$grp)

#Which individuals best represent each cluster:
cluster_represent_lump_broad_ss1 <- apply(qmatrix2_lump_broad_ss1, 2, 
                                          function(col) 
                                            rownames(qmatrix2_lump_broad_ss1)
                                          [which.max(col)])


#use a few reference individuals with certain ancestry to change the cluster names so they are consistent between the DAPC, map and barplot
for(i in c(
  cluster_represent_lump_broad_ss1[["V1"]], 
  cluster_represent_lump_broad_ss1[["V2"]], 
  cluster_represent_lump_broad_ss1[["V3"]], 
  cluster_represent_lump_broad_ss1[["V4"]], 
  cluster_represent_lump_broad_ss1[["V5"]], 
  cluster_represent_lump_broad_ss1[["V6"]],
  cluster_represent_lump_broad_ss1[["V7"]]
)){
  max_anc_lump_broad_ss1 <- filter(qmatrix3_lump_broad_ss1,ID==i) %>% filter(Prob==max(Prob))
  qmatrix3_lump_broad_ss1$Variable <- 
    replace(qmatrix3_lump_broad_ss1$Variable,
            qmatrix3_lump_broad_ss1$Variable==max_anc_lump_broad_ss1$Variable,
            max_anc_lump_broad_ss1$Kmeans_cluster)}

#plot the ancestry coefficients in a nicer-looking barplot
anc_plot_lump_broad_ss1 <- ggplot(qmatrix3_lump_broad_ss1, 
                                  aes(factor(ID), Prob, 
                                      fill = factor(Variable))) +
  geom_col(width=1) +
  facet_grid(~fct_inorder(as.factor(Population)), switch = "x", 
             scales = "free", 
             space = "free") +
  theme_minimal() + 
  labs(x = "", title = , y = "Probability") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.05, "lines"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size=5, angle = 90)) +
  scale_fill_manual(name="Cluster",
                    values=c("lightblue", "goldenrod2", "mediumorchid3", 
                             "darkolivegreen3", "firebrick2", "deepskyblue4", 
                             "chocolate2", "chartreuse4", "tomato3", "purple3"))
anc_plot_lump_broad_ss1

#Take the 3 plots and plot them together
grid.arrange(arrangeGrob(as_grob(dapc_lump_broad_ss1_plot),
                         as_grob(map_plot_lump_broad_ss1),
                         ncol=2,clip="on"),
             anc_plot_lump_broad_ss1,nrow = 2,heights=c(4,2))
combined_map_lump_broad_ss1 <- recordPlot()

#----4. Subset: Iceland-----
##----Subset----
#Read and filter
metadata_lump_broad_ss2 <- read_xlsx("Lumpfish_metadata.xlsx") %>%
  filter(grepl("_(IS)$", Pop_ID_article))

# extracting metadata info
pop_ids_ss2 <- genind_lump_broad@other$metadata$Pop_ID_article

# identifying which individuals to keep
keep_pops_ss2 <- grepl("_IS", pop_ids_ss2)

# creating genind subset
genind_lump_broad_ss2 <- genind_lump_broad[keep_pops_ss2, ]
unique(genind_lump_broad_ss2@other$metadata$Pop_ID_article)

##----PCA----
?dudi.pca
x_lump_broad_ss2 <- tab(genind_lump_broad_ss2, freq=TRUE, NA.method="mean")
pca_lump_broad_ss2 <- dudi.pca(x_lump_broad_ss2, center=TRUE, scale=FALSE) # 5 axes
pca_lump_broad_ss2

s.class(pca_lump_broad_ss2$li, 
        fac = as.factor(genind_lump_broad_ss2@other$metadata$Pop_ID_article),  
        clab = 1, 
        col = transp(funky(length(unique(
          genind_lump_broad_ss2@other$metadata$Pop_ID_article)))), 
        csta = 0, 
        cpoint = 4, 
        cellipse = 1, 
        xax = 1, 
        yax = 2)

eig.perc.broad_ss2 <- 100*pca_lump_broad_ss2$eig/sum(pca_lump_broad_ss2$eig)
head(eig.perc.broad_ss2) #31.1318980 11.9283696  2.6121774  1.6230219 %

##-----DAPC------
#How many sampling sites were there?
length(unique(genind_lump_broad_ss2@other$metadata$Pop_ID_article)) 

kmeans_lump_broad_ss2 <- find.clusters(genind_lump_broad_ss2, 
                                       n.pca = 50, 
                                       max.n = 7,
                                       scale = FALSE) #3 clusters
kmeans_lump_broad_ss2$grp %>% as.vector

dapc_lump_broad_ss2 <- dapc(genind_lump_broad_ss2, 
                            pop = kmeans_lump_broad_ss2$grp, 
                            n.pca = 50) # 2 DFs retained

# Plotting the first 2 linear discrimination functions
scatter(dapc_lump_broad_ss2, xax=1, yax=2, grp=dapc_lump_broad_ss2$grp, 
        col=transp(c("lightblue", "goldenrod2", "mediumorchid3", 
                     "darkolivegreen3", "firebrick2", "deepskyblue4", 
                     "chocolate2", "chartreuse4", "tomato3", "purple3")),
        pch=19, bg="white",cstar = 1, cellipse = 1,clabel = 1,scree.da=FALSE,
        scree.pca=FALSE)
dapc_lump_broad_ss2_plot <- recordPlot()

#and PCA again with the same color assignment
s.class(pca_lump_broad_ss2$li, fac=dapc_lump_broad_ss2$grp, clab=1,
        col=transp(c("lightblue", "goldenrod2", "mediumorchid3", 
                     "darkolivegreen3", "firebrick2", "deepskyblue4", 
                     "chocolate2", "chartreuse4", "tomato3", "purple3")), 
        csta=1, cpoint=2, cellipse =1, xax=1, yax=2)
pca_lump_broad_ss2_plot <- recordPlot()

##----Map----
#Data frame with pop and cluster assignment for each individual
pop_cluster_lump_broad_ss2 <- 
  data.frame(Population=genind_lump_broad_ss2@other$metadata$Pop_ID_article, 
             Group=kmeans_lump_broad_ss2$grp)

count_cluster_per_pop_broad_ss2 <- pop_cluster_lump_broad_ss2 %>% 
  dplyr::count(Population,Group) %>% 
  data.frame()
head(count_cluster_per_pop_broad_ss2)

#Location coordinates
##ensuring that the coordinate columns contain numeric values
metadata_lump_broad_ss2$Lat <- as.numeric(metadata_lump_broad_ss2$Lat)
metadata_lump_broad_ss2$Lon <- as.numeric(metadata_lump_broad_ss2$Lon)

#calculating the average coordinates of each population
location_coords_broad_ss2 <- metadata_lump_broad_ss2 %>%
  dplyr::group_by(Pop_ID_article) %>%
  dplyr::summarise(
    Lat = mean(Lat, na.rm = TRUE),
    Lon = mean(Lon, na.rm = TRUE)
  ) %>%
  dplyr::rename(Population = Pop_ID_article) %>% 
  as.data.frame() 
head(location_coords_broad_ss2)


#Merging objects based on the population column
cluster_pop_coords_broad_ss2 <- inner_join(location_coords_broad_ss2, count_cluster_per_pop_broad_ss2, by="Population")
head(cluster_pop_coords_broad_ss2)

xyz_broad_ss2 <- make.xyz(cluster_pop_coords_broad_ss2$Lon,
                          cluster_pop_coords_broad_ss2$Lat,
                          cluster_pop_coords_broad_ss2$n,
                          cluster_pop_coords_broad_ss2$Group)

#Map
shape <- read.shapefile("world_map")

?basemap
#1. Map size
basemap(xlim=c(-25,-12), #min/max longitude
        ylim=c(65,65), #min/max latitude
        # yaxp = c(40,75,5),
        bg="white",
        xlab = "",
        ylab = "",
        frame.plot = FALSE,
        axes = FALSE) 
axis(2, at = seq(30, 80, by = 2), cex.axis = 0.6)  
axis(1, at = seq(-30, 10, by = 5), cex.axis = 0.6)


#2. Load shapefile
map <- draw.shape(shape, col="grey85")

#3. Pie charts and labels
draw.pie(xyz_broad_ss2$x, xyz_broad_ss2$y, xyz_broad_ss2$z, radius=0.3,
         col=c("lightblue", "goldenrod2", "mediumorchid3", 
               "darkolivegreen3", "firebrick2", "deepskyblue4", 
               "chocolate2", "chartreuse4", "tomato3", "purple3"),
)

text(cluster_pop_coords_broad_ss2$Lon, cluster_pop_coords_broad_ss2$Lat - 0.3,  
     labels = cluster_pop_coords_broad_ss2$Population,
     col = "black", cex = 0.8, font = 0.5)

map_plot_lump_broad_ss2 <- recordPlot()

## sNMF----

##Converting to .geno format
hierfstat_lump_broad_ss2 <- genind2hierfstat(genind_lump_broad_ss2)

write.struct(hierfstat_lump_broad_ss2, 
             ilab=indNames(genind_lump_broad_ss2), 
             pop = genind_lump_broad_ss2@other$metadata$Pop_ID_article, 
             fname = "lump_broad_ss2.str")

struct2geno(input.file = "lump_broad_ss2.str", 
            ploidy = 2, 
            FORMAT = 2, 
            extra.column = 2, 
            extra.row = 0) 

genofile_lump_broad_ss2 <- "lump_broad_ss2.str.geno"
cic.geno_lump_broad_ss2 <- read.geno(genofile_lump_broad_ss2)


#How many ancestral pops are there? (TAKES A BIT OF TIME)  
length(unique(genind_lump_broad_ss2@other$metadata$Pop_ID_article)) # Number of pops
snmf_search_lump_broad_ss2 <- snmf(genofile_lump_broad_ss2, 
                                   K = 2:8, #Using max number = number of pops
                                   entropy = TRUE, 
                                   repetitions = 2, 
                                   project = "new")

# Check the plot, should converge with DAPC, in this case doesnt though
plot(snmf_search_lump_broad_ss2, col = transp("steelblue4"), pch = 19) # 3 ancestral pops

#calculating ancestry of each individual for K=3
snmf_lump_broad_ss2 <- snmf(genofile_lump_broad_ss2, K=3, project="new")

#probability matrix for K=10
qmatrix_lump_broad_ss2 <- Q(snmf_lump_broad_ss2, K=3)

#create a new qmatrix and add some information (ID, population, and K-means cluster assignment)
qmatrix2_lump_broad_ss2 <- data.frame(ID=rownames(pop_cluster_lump_broad_ss2),
                                      pop_cluster_lump_broad_ss2,
                                      qmatrix_lump_broad_ss2)

qmatrix_new_list_lump_broad_ss2 <- list() 

#for-loop to create the new columns:
for(i in 4:ncol(qmatrix2_lump_broad_ss2)){ 
  qmatrix_new_list_lump_broad_ss2[[i-3]] <- qmatrix2_lump_broad_ss2[,c(1:3,i)] %>% 
    mutate("Var"=rep(colnames(qmatrix2_lump_broad_ss2)[i],nrow(qmatrix2_lump_broad_ss2)))
  colnames(qmatrix_new_list_lump_broad_ss2[[i-3]])  <- c("ID","Population","Kmeans_cluster","Prob","Variable")
}

#transform the list into the new qmatrix
qmatrix3_lump_broad_ss2 <- bind_rows(qmatrix_new_list_lump_broad_ss2)


#data.frame(select(metadata_lump_broad_ss2, 
#                 c(Pop_ID_article, sample_ID_SNPdata)), 
#         kmeans_lump_broad_ss2$grp)

#Which individuals best represent each cluster:
cluster_represent_lump_broad_ss2 <- apply(qmatrix2_lump_broad_ss2, 2, 
                                          function(col) 
                                            rownames(qmatrix2_lump_broad_ss2)
                                          [which.max(col)])


#use a few reference individuals with certain ancestry to change the cluster names so they are consistent between the DAPC, map and barplot
for(i in c(
  cluster_represent_lump_broad_ss2[["V1"]], 
  cluster_represent_lump_broad_ss2[["V2"]], 
  cluster_represent_lump_broad_ss2[["V3"]], 
  cluster_represent_lump_broad_ss2[["V4"]], 
  cluster_represent_lump_broad_ss2[["V5"]], 
  cluster_represent_lump_broad_ss2[["V6"]]
)){
  max_anc_lump_broad_ss2 <- filter(qmatrix3_lump_broad_ss2,ID==i) %>% filter(Prob==max(Prob))
  qmatrix3_lump_broad_ss2$Variable <- 
    replace(qmatrix3_lump_broad_ss2$Variable,
            qmatrix3_lump_broad_ss2$Variable==max_anc_lump_broad_ss2$Variable,
            max_anc_lump_broad_ss2$Kmeans_cluster)}

#plot the ancestry coefficients in a nicer-looking barplot
anc_plot_lump_broad_ss2 <- ggplot(qmatrix3_lump_broad_ss2, 
                                  aes(factor(ID), Prob, 
                                      fill = factor(Variable))) +
  geom_col(width=1) +
  facet_grid(~fct_inorder(as.factor(Population)), switch = "x", 
             scales = "free", 
             space = "free") +
  theme_minimal() + 
  labs(x = "", title = , y = "Probability") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.05, "lines"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size=5, angle = 90)) +
  scale_fill_manual(name="Cluster",
                    values=c("lightblue", "goldenrod2", "mediumorchid3", 
                             "darkolivegreen3", "firebrick2", "deepskyblue4", 
                             "chocolate2", "chartreuse4", "tomato3", "purple3"))
anc_plot_lump_broad_ss2

#Take the 3 plots and plot them together
grid.arrange(arrangeGrob(as_grob(dapc_lump_broad_ss2_plot),
                         as_grob(map_plot_lump_broad_ss2),
                         ncol=2,clip="on"),
             anc_plot_lump_broad_ss2,nrow = 2,heights=c(4,2))
combined_map_lump_broad_ss2 <- recordPlot()

#----5. Subset: Norway-----
##----Subset----
#Read and filter
metadata_lump_broad_ss3 <- read_xlsx("Lumpfish_metadata.xlsx") %>%
  filter(grepl("_(NO)$", Pop_ID_article))

# extracting metadata info
pop_ids_ss3 <- genind_lump_broad@other$metadata$Pop_ID_article

# identifying which individuals to keep
keep_pops_ss3 <- grepl("_NO", pop_ids_ss3)

# creating genind subset
genind_lump_broad_ss3 <- genind_lump_broad[keep_pops_ss3, ]
unique(genind_lump_broad_ss3@other$metadata$Pop_ID_article)

##----PCA----
?dudi.pca
x_lump_broad_ss3 <- tab(genind_lump_broad_ss3, freq=TRUE, NA.method="mean")
pca_lump_broad_ss3 <- dudi.pca(x_lump_broad_ss3, center=TRUE, scale=FALSE) # 5 axes
pca_lump_broad_ss3

s.class(pca_lump_broad_ss3$li, 
        fac = as.factor(genind_lump_broad_ss3@other$metadata$Pop_ID_article),  
        clab = 1, 
        col = transp(funky(length(unique(
          genind_lump_broad_ss3@other$metadata$Pop_ID_article)))), 
        csta = 0, 
        cpoint = 4, 
        cellipse = 1, 
        xax = 1, 
        yax = 2)

eig.perc.broad_ss3 <- 100*pca_lump_broad_ss3$eig/sum(pca_lump_broad_ss3$eig)
head(eig.perc.broad_ss3) #31.1318980 11.9283696  2.6121774  1.6230219 %

##-----DAPC------
#How many sampling sites were there?
length(unique(genind_lump_broad_ss3@other$metadata$Pop_ID_article)) 

kmeans_lump_broad_ss3 <- find.clusters(genind_lump_broad_ss3, 
                                       n.pca = 50, 
                                       max.n = 7,
                                       scale = FALSE) #3 clusters
kmeans_lump_broad_ss3$grp %>% as.vector

dapc_lump_broad_ss3 <- dapc(genind_lump_broad_ss3, 
                            pop = kmeans_lump_broad_ss3$grp, 
                            n.pca = 50) # 2 DFs retained

# Plotting the first 2 linear discrimination functions
scatter(dapc_lump_broad_ss3, xax=1, yax=2, grp=dapc_lump_broad_ss3$grp, 
        col=transp(c("darkolivegreen3", "firebrick2", "deepskyblue4", 
                     "chocolate2", "chartreuse4", "tomato3", "purple3")),
        pch=19, bg="white",cstar = 1, cellipse = 1,clabel = 1,scree.da=FALSE,
        scree.pca=FALSE)
dapc_lump_broad_ss3_plot <- recordPlot()

#and PCA again with the same color assignment
s.class(pca_lump_broad_ss3$li, fac=dapc_lump_broad_ss3$grp, clab=1,
        col=transp(c("lightblue", "goldenrod2", "mediumorchid3", 
                     "darkolivegreen3", "firebrick2", "deepskyblue4", 
                     "chocolate2", "chartreuse4", "tomato3", "purple3")), 
        csta=1, cpoint=2, cellipse =1, xax=1, yax=2)
pca_lump_broad_ss3_plot <- recordPlot()

##----Map----
#Data frame with pop and cluster assignment for each individual
pop_cluster_lump_broad_ss3 <- 
  data.frame(Population=genind_lump_broad_ss3@other$metadata$Pop_ID_article, 
             Group=kmeans_lump_broad_ss3$grp)

count_cluster_per_pop_broad_ss3 <- pop_cluster_lump_broad_ss3 %>% 
  dplyr::count(Population,Group) %>% 
  data.frame()
head(count_cluster_per_pop_broad_ss3)

#Location coordinates
##ensuring that the coordinate columns contain numeric values
metadata_lump_broad_ss3$Lat <- as.numeric(metadata_lump_broad_ss3$Lat)
metadata_lump_broad_ss3$Lon <- as.numeric(metadata_lump_broad_ss3$Lon)

#calculating the average coordinates of each population
location_coords_broad_ss3 <- metadata_lump_broad_ss3 %>%
  dplyr::group_by(Pop_ID_article) %>%
  dplyr::summarise(
    Lat = mean(Lat, na.rm = TRUE),
    Lon = mean(Lon, na.rm = TRUE)
  ) %>%
  dplyr::rename(Population = Pop_ID_article) %>% 
  as.data.frame() 
head(location_coords_broad_ss3)


#Merging objects based on the population column
cluster_pop_coords_broad_ss3 <- inner_join(location_coords_broad_ss3, count_cluster_per_pop_broad_ss3, by="Population")
head(cluster_pop_coords_broad_ss3)

xyz_broad_ss3 <- make.xyz(cluster_pop_coords_broad_ss3$Lon,
                          cluster_pop_coords_broad_ss3$Lat,
                          cluster_pop_coords_broad_ss3$n,
                          cluster_pop_coords_broad_ss3$Group)

#Map
shape <- read.shapefile("world_map")

?basemap
#1. Map size
basemap(xlim=c(8,8), #min/max longitude
        ylim=c(57,72), #min/max latitude
        # yaxp = c(40,75,5),
        bg="white",
        xlab = "",
        ylab = "",
        frame.plot = FALSE,
        axes = FALSE) 
axis(2, at = seq(35, 85, by = 5), cex.axis = 0.6)  
axis(1, at = seq(-80, 40, by = 10), cex.axis = 0.6)


#2. Load shapefile
map <- draw.shape(shape, col="grey85")

#3. Pie charts and labels
draw.pie(xyz_broad_ss3$x, xyz_broad_ss3$y, xyz_broad_ss3$z, radius=0.7,
         col=c("darkolivegreen3", "firebrick2", "deepskyblue4", 
               "chocolate2", "chartreuse4", "tomato3", "purple3"),
)

text(cluster_pop_coords_broad_ss3$Lon, cluster_pop_coords_broad_ss3$Lat - 0.3,  
     labels = cluster_pop_coords_broad_ss3$Population,
     col = "black", cex = 0.8, font = 0.5)

map_plot_lump_broad_ss3 <- recordPlot()

## sNMF----

##Converting to .geno format
hierfstat_lump_broad_ss3 <- genind2hierfstat(genind_lump_broad_ss3)

write.struct(hierfstat_lump_broad_ss3, 
             ilab=indNames(genind_lump_broad_ss3), 
             pop = genind_lump_broad_ss3@other$metadata$Pop_ID_article, 
             fname = "lump_broad_ss3.str")

struct2geno(input.file = "lump_broad_ss3.str", 
            ploidy = 2, 
            FORMAT = 2, 
            extra.column = 2, 
            extra.row = 0) 

genofile_lump_broad_ss3 <- "lump_broad_ss3.str.geno"
cic.geno_lump_broad_ss3 <- read.geno(genofile_lump_broad_ss3)


#How many ancestral pops are there? (TAKES A BIT OF TIME)  
length(unique(genind_lump_broad_ss3@other$metadata$Pop_ID_article)) # Number of pops
snmf_search_lump_broad_ss3 <- snmf(genofile_lump_broad_ss3, 
                                   K = 2:8, #Using max number = number of pops
                                   entropy = TRUE, 
                                   repetitions = 2, 
                                   project = "new")

# Check the plot, should converge with DAPC, in this case doesnt though
plot(snmf_search_lump_broad_ss3, col = transp("steelblue4"), pch = 19) # 3 ancestral pops!!

#calculating ancestry of each individual for K=2
snmf_lump_broad_ss3 <- snmf(genofile_lump_broad_ss3, K=3, project="new")

#probability matrix for K=10
qmatrix_lump_broad_ss3 <- Q(snmf_lump_broad_ss3, K=3)

#create a new qmatrix and add some information (ID, population, and K-means cluster assignment)
qmatrix2_lump_broad_ss3 <- data.frame(ID=rownames(pop_cluster_lump_broad_ss3),
                                      pop_cluster_lump_broad_ss3,
                                      qmatrix_lump_broad_ss3)

qmatrix_new_list_lump_broad_ss3 <- list() 

#for-loop to create the new columns:
for(i in 4:ncol(qmatrix2_lump_broad_ss3)){ 
  qmatrix_new_list_lump_broad_ss3[[i-3]] <- qmatrix2_lump_broad_ss3[,c(1:3,i)] %>% 
    mutate("Var"=rep(colnames(qmatrix2_lump_broad_ss3)[i],nrow(qmatrix2_lump_broad_ss3)))
  colnames(qmatrix_new_list_lump_broad_ss3[[i-3]])  <- c("ID","Population","Kmeans_cluster","Prob","Variable")
}

#transform the list into the new qmatrix
qmatrix3_lump_broad_ss3 <- bind_rows(qmatrix_new_list_lump_broad_ss3)


#data.frame(select(metadata_lump_broad_ss3, 
#                 c(Pop_ID_article, sample_ID_SNPdata)), 
#         kmeans_lump_broad_ss3$grp)

#Which individuals best represent each cluster:
cluster_represent_lump_broad_ss3 <- apply(qmatrix2_lump_broad_ss3, 2, 
                                          function(col) 
                                            rownames(qmatrix2_lump_broad_ss3)
                                          [which.max(col)])


#use a few reference individuals with certain ancestry to change the cluster names so they are consistent between the DAPC, map and barplot
for(i in c(
  cluster_represent_lump_broad_ss3[["V1"]], 
  cluster_represent_lump_broad_ss3[["V2"]], 
  cluster_represent_lump_broad_ss3[["V3"]], 
  cluster_represent_lump_broad_ss3[["V4"]], 
  cluster_represent_lump_broad_ss3[["V5"]], 
  cluster_represent_lump_broad_ss3[["V6"]]
)){
  max_anc_lump_broad_ss3 <- filter(qmatrix3_lump_broad_ss3,ID==i) %>% filter(Prob==max(Prob))
  qmatrix3_lump_broad_ss3$Variable <- 
    replace(qmatrix3_lump_broad_ss3$Variable,
            qmatrix3_lump_broad_ss3$Variable==max_anc_lump_broad_ss3$Variable,
            max_anc_lump_broad_ss3$Kmeans_cluster)}

#plot the ancestry coefficients in a nicer-looking barplot
anc_plot_lump_broad_ss3 <- ggplot(qmatrix3_lump_broad_ss3, 
                                  aes(factor(ID), Prob, 
                                      fill = factor(Variable))) +
  geom_col(width=1) +
  facet_grid(~fct_inorder(as.factor(Population)), switch = "x", 
             scales = "free", 
             space = "free") +
  theme_minimal() + 
  labs(x = "", title = , y = "Probability") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.05, "lines"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size=5, angle = 90)) +
  scale_fill_manual(name="Cluster",
                    values=c("darkolivegreen3", "firebrick2", "deepskyblue4", 
                             "chocolate2", "chartreuse4", "tomato3", "purple3"))
anc_plot_lump_broad_ss3

#Take the 3 plots and plot them together
grid.arrange(arrangeGrob(as_grob(dapc_lump_broad_ss3_plot),
                         as_grob(map_plot_lump_broad_ss3),
                         ncol=2,clip="on"),
             anc_plot_lump_broad_ss3,nrow = 2,heights=c(4,2))
combined_map_lump_broad_ss3 <- recordPlot()


#----Citations----

?citation
citation(package = "adegenet")
citation(package = "LEA")
citation(package = "hierfstat")
citation(package = "vcfR")
citation(package = "mapplots")
citation(package = "base")



Frichot E, Francois O (2015). “LEA: an R package for Landscape and Ecological Association
studies.” _Methods in Ecology and Evolution_.
<http://membres-timc.imag.fr/Olivier.Francois/lea.html>.

Gerritsen H (2023). _mapplots: Data Visualisation on Maps_. R package version 1.5.2,
<https://CRAN.R-project.org/package=mapplots>.

Goudet J, Jombart T (2022). _hierfstat: Estimation and Tests of Hierarchical F-Statistics_.
R package version 0.5-11, <https://CRAN.R-project.org/package=hierfstat>.

Jombart, T. (2008) adegenet: a R package for the multivariate analysis of genetic markers.
Bioinformatics 24: 1403-1405. doi: 10.1093/bioinformatics/btn129

Jombart T. and Ahmed I. (2011) adegenet 1.3-1: new tools for the analysis of genome-wide SNP
data. Bioinformatics. doi: 10.1093/bioinformatics/btr521

Knaus BJ, Grünwald NJ (2017). “VCFR: a package to manipulate and visualize variant call
format data in R.” _Molecular Ecology Resources_, *17*(1), 44-53. ISSN 757,
<https://dx.doi.org/10.1111/1755-0998.12549>.


