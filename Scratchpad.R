#2. Subset 1: -----

##----Read and Subset Data Failed attempt----
genind_lump_broad_ss1 <- read.genepop("Lumpfish_SNP.gen", ncode = 3)

#filtering metadata to retain only individuals from SE, NO, IS, DK
metadata_lump_broad_ss1 <- read_xlsx("Lumpfish_metadata.xlsx") %>%
  filter(grepl("_(SE|NO|IS|DK)$", Pop_ID_article))

# findng individuals in genind that also exist in the filtered metadata
matching_inds_ss1 <- rownames(genind_lump_broad_ss1@tab) %in% 
                     metadata_lump_broad_ss1$sample_ID_SNPdata

# creating a filtered genind object with only those individuals
genind_lump_broad_ss1@tab <- genind_lump_broad_ss1@tab[matching_inds_ss1, ]

# Align and combine
metadata_lump_broad_ss1 <- 
  metadata_lump_broad_ss1[match(rownames(genind_lump_broad_ss1@tab), 
                          metadata_lump_broad_ss1$sample_ID_SNPdata), ]
genind_lump_broad_ss1@other$metadata <- metadata_lump_broad_ss1



#----Subset: Draft-----

#Read and filter
metadata_lump_broad_ss1 <- read_xlsx("Lumpfish_metadata.xlsx") %>%
  filter(grepl("_(SE|NO|IS|DK)$", Pop_ID_article))

# extracting metadata info
pop_ids_ss1 <- genind_lump_broad@other$metadata$Pop_ID_article

# identifying which individuals to keep
keep_pops_ss1 <- grepl("_SE|_NO|_DK|_IS", pop_ids)

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
                                   max.n = 24,
                                   scale = FALSE) #12 clusters
kmeans_lump_broad_ss1$grp %>% as.vector

dapc_lump_broad_ss1 <- dapc(genind_lump_broad_ss1, 
                            pop = kmeans_lump_broad_ss1$grp, 
                            n.pca = 50) # 2 DFs retained

# Plotting the first 2 linear discrimination functions
scatter(dapc_lump_broad_ss1, xax=1, yax=2, grp=dapc_lump_broad_ss1$grp, 
        col=transp(c("dodgerblue3", "purple2", "darkolivegreen3", 
                     "firebrick2", "deepskyblue4", "chocolate2", 
                     "yellow2", "chartreuse4", "tomato3",
                     "steelblue3", "tan3", "olivedrab3", "magenta", 
                     "goldenrod2", "orangered3", "royalblue3", "peru", 
                     "mediumvioletred",  "seagreen3", "indianred3", 
                     "darkslategray3")),
        pch=19, bg="white",cstar = 1, cellipse = 1,clabel = 1,scree.da=FALSE,
        scree.pca=FALSE)
dapc_lump_broad_ss1_plot <- recordPlot()

#and PCA again with the same color assignment
s.class(pca_lump_broad_ss1$li, fac=dapc_lump_broad_ss1$grp, clab=1,
        col=transp(c("dodgerblue3", "purple2", "darkolivegreen3", 
                     "firebrick2", "deepskyblue4", "chocolate2", 
                     "yellow2", "chartreuse4", "tomato3",
                     "steelblue3", "tan3", "olivedrab3", "magenta", 
                     "goldenrod2", "orangered3", "royalblue3", "peru", 
                     "mediumvioletred",  "seagreen3", "indianred3", 
                     "darkslategray3")), 
        csta=1, cpoint=2, cellipse =1, xax=1, yax=2)
pca_lump_broad_ss1_plot <- recordPlot()

##----Map----
#Data frame with pop and cluster assignment for each individual
pop_cluster_lump_broad_ss1 <- data.frame(Population=genind_lump_broad_ss1@other$metadata$Pop_ID_article, 
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
basemap(xlim=c(0,10), #min/max longitude
        ylim=c(52,73), #min/max latitude
        # yaxp = c(40,75,5),
        bg="white",frame.plot=F,
        main = "1669 Individuals and 139 Loci") 
#axis(2, at = seq(50, 75, by = 5))  # Adjust this sequence for your preferred interval
#axis(1, at = seq(-40, 40, by = 20),)


#2. Load shapefile
map <- draw.shape(shape, col="grey85")

#3. Pie charts and labels
draw.pie(xyz_broad_ss1$x, xyz_broad_ss1$y, xyz_broad_ss1$z, radius=1,
         col=c("dodgerblue3", "purple2", "darkolivegreen3", 
               "firebrick2", "deepskyblue4", "chocolate2", 
               "yellow2", "chartreuse4", "tomato3",
               "steelblue3", "tan3", "olivedrab3", "magenta", 
               "goldenrod2", "orangered3", "royalblue3", "peru", 
               "mediumvioletred",  "seagreen3", "indianred3", 
               "darkslategray3"),
)

text(cluster_pop_coords_broad_ss1$Lon, cluster_pop_coords_broad_ss1$Lat - 2,  
     labels = cluster_pop_coords_broad_ss1$Population,
     col = "black", cex = 0.8, font = 0.5)

map_plot_lump_broad_ss1 <- recordPlot()

#2.4 sNMF broad----

#To consider genetic admixture, we can calculate the ancestry coefficients (Q) 
# of individuals. We will use a method called sparse Non-negative Matrix Factorisation (sNMF)

##--Converting to .geno format----
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

barplot(t(qmatrix_lump_broad), 
        col=c("dodgerblue3", "purple2", "darkolivegreen3", 
              "firebrick2", "deepskyblue4", "chocolate2", 
              "yellow2", "chartreuse4", "tomato3",
              "steelblue3", "tan3", "olivedrab3", "magenta", 
              "goldenrod2", "orangered3", "royalblue3", "peru", 
              "mediumvioletred",  "seagreen3", "indianred3", 
              "darkslategray3"), 
        border=NA, space=0, 
        xlab="Individuals", 
        ylab = "Ancestry")

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


#data.frame(select(metadata_lump_broad, 
 #                 c(Pop_ID_article, sample_ID_SNPdata)), 
  #         kmeans_lump_broad$grp)

#Which individuals best represent each cluster:
cluster_represent_lump_broad <- apply(qmatrix2_lump_broad, 2, 
                                      function(col) 
                                        rownames(qmatrix2_lump_broad)
                                        [which.max(col)])


#use a few reference individuals with certain ancestry to change the cluster names so they are consistent between the DAPC, map and barplot
for(i in c(
  cluster_represent_lump_broad[["V1"]], 
  cluster_represent_lump_broad[["V2"]], 
  cluster_represent_lump_broad[["V3"]], 
  cluster_represent_lump_broad[["V4"]], 
  cluster_represent_lump_broad[["V5"]], 
  cluster_represent_lump_broad[["V6"]], 
  cluster_represent_lump_broad[["V7"]], 
  cluster_represent_lump_broad[["V8"]], 
  cluster_represent_lump_broad[["V9"]], 
  cluster_represent_lump_broad[["V10"]]
)){
  max_anc_lump_broad <- filter(qmatrix3_lump_broad,ID==i) %>% filter(Prob==max(Prob))
  qmatrix3_lump_broad$Variable <- 
    replace(qmatrix3_lump_broad$Variable,
            qmatrix3_lump_broad$Variable==max_anc_lump_broad$Variable,
            max_anc_lump_broad$Kmeans_cluster)}

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
                    values=c("dodgerblue3", "goldenrod2", "mediumorchid3", 
                             "darkolivegreen3", "firebrick2", "deepskyblue4", 
                             "chocolate2", "chartreuse4", "tomato3", "purple3"))
anc_plot_lump_broad

#HMMMMM its missing cluster 1... 

#Take the 3 plots we saved previously, and plot them together
grid.arrange(arrangeGrob(as_grob(dapc_plot_),
                         as_grob(map_plot),
                         ncol=2,clip="on"),
             anc_plot,nrow = 2,heights=c(4,2))
combined_map <- recordPlot()
## sNMF ss1----

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
                                   K = 2:24, #Using max number = number of pops
                                   entropy = TRUE, 
                                   repetitions = 2, 
                                   project = "new")

# Check the plot, should converge with DAPC, in this case doesnt though
plot(snmf_search_lump_broad_ss1, col = transp("steelblue4"), pch = 19) # 6 ancestral pops

#calculating ancestry of each individual for K=6
snmf_lump_broad_ss1 <- snmf(genofile_lump_broad_ss1, K=6, project="new")

#probability matrix for K=10
qmatrix_lump_broad_ss1 <- Q(snmf_lump_broad_ss1, K=6)

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
  cluster_represent_lump_broad_ss1[["V7"]], 
  cluster_represent_lump_broad_ss1[["V8"]], 
  cluster_represent_lump_broad_ss1[["V9"]], 
  cluster_represent_lump_broad_ss1[["V10"]]
)){
  max_anc_lump_broad_ss1 <- filter(qmatrix3_lump_broad_ss1,ID==i) %>% filter(Prob==max(Prob))
  qmatrix3_lump_broad_ss1$Variable <- 
    replace(qmatrix3_lump_broad_ss1$Variable,
            qmatrix3_lump_broad_ss1$Variable==max_anc_lump_broad_ss1$Variable,
            max_anc_lump_broad_ss1$Kmeans_cluster)}

#plot the ancestry coefficients in a nicer-looking barplot
anc_plot_lump_broad_ss1 <- ggplot(qmatrix3_lump_broad_ss1, 
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
anc_plot_lump_broad_ss1

#Take the 3 plots and plot them together
grid.arrange(arrangeGrob(as_grob(dapc_lump_broad_ss1_plot),
                         as_grob(map_plot_lump_broad_ss1),
                         ncol=2,clip="on"),
             anc_plot_lump_broad_ss1,nrow = 2,heights=c(4,2))
combined_map_lump_broad_ss1 <- recordPlot()
#----3. Subset: Iceland-----
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
axis(2, at = seq(35, 85, by = 5), cex.axis = 0.6)  # Adjust this sequence for your preferred interval
axis(1, at = seq(-80, 40, by = 10), cex.axis = 0.6)


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
