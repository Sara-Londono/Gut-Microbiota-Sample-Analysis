#_______________________________________________________________________________
#     Diversity analysis 23 samples 030822      
#         Author: Sara Londoño Osorio            
#_______________________________________________________________________________

# call libraries
install.packages("rlang")  # run when there is a problem running ggplot2
library(rlang)             # run when there is a problem running ggplot2
library(ggplot2)
library(phyloseq)
library(phyloseqGraphTest)
library(BiocManager)
library(gridExtra)
library(ggsignif)  
library(RColorBrewer)
library(permute)
library(lattice)
library(vegan)
library(dplyr)
library(stringi)
library(stringr)
library(reshape2)
library(carData)
library(car) 
library(Matrix)
library(lme4)
library(mvtnorm)
library(survival)
library(MASS)
library(TH.data)
library(zoo)
library(multcomp) 
library(UpSetR)
library(circlize)
library(limma)
library(RCurl)
library(bitops)
library(ALDEx2)
library(dendextend)
library(readr)
library(tibble)
library(ComplexHeatmap) 

# nedeed
library(edgeR)  # error
BiocManager::install('edgeR')

install.packages("ggpubr")
library(ggpubr)   # to join figures 
install.packages("corrplot")
library("corrplot") # Correlation matrix

#               Import data 
# data_path:C:\Users\EQ01\Documents\Microbiomas\metagenomica\AnalisisComplete23samples\
#_______________________________________________________________________________

# 1.Read the sample metadata files
metadata <- read.table(file = "metadata23samplesQ2.tsv", header = TRUE, sep = "\t",
                       encoding = "UTF-8",na.strings = "")  

# 2.Read feature table (number of reads of each feature per sample)
feature_table <- read.table(file= "feature-table_r.txt", header = TRUE, sep= "\t")

# 3.Read taxonomy table (Features ID)
taxonomy <- read.table(file= "taxonomy.txt", header = TRUE, sep= "\t")


#                  Clean data
#_______________________________________________________________________________

# 1. 
metadata$sample.id -> rownames(metadata)  # To named the rows with the ID name
# 2.
feature_table <- feature_table %>% 
  tibble::column_to_rownames("OTU.ID") # define the row names from the sample.id column

colnames(feature_table) <- c("10isl_S9","11isl_S10","122isl_S1","123isl_S2",
                             "124isl_S3","125isl_S4","126isl_S5","127isl_S6","128isl_S7",
                             "129isl_S8","12isl_S11","130isl_S9","131isl_S10","132isl_S11",
                             "133isl_S12","134isl_S13","135isl_S14","136isl_S15",
                             "137isl_S16","138isl_S17","139isl_S18","13isl_S12",
                             "140isl_S19","141isl_S20","142isl_S21","143isl_S22",
                             "144isl_S23","145isl_S24","146isl_S25","147isl_S26",
                             "148isl_S27","149isl_S28","14isl_S13","150isl_S29",
                             "151isl_S30","15isl_S14","16isl_S15","17isl_S16",
                             "18isl_S17","19isl_S18","20isl_S19","21isl_S20",
                             "22isl_S21","23isl_S22","24isl_S23","9isl_S8")
# 3.
taxonomy <- taxonomy %>% 
  tibble::column_to_rownames("FeatureID") # define the row names from the sample.id column


#                  Create Phyloseq Object
#_______________________________________________________________________________

# Transform into matrixes feature and tax tables (sample table can be left as data frame)
feature_table <- as.matrix(feature_table)
taxonomy <- as.matrix(taxonomy)

# Use phyloseq package to link ASV count-ASV taxa-Metadata
ps_23 <- phyloseq(otu_table(feature_table, taxa_are_rows=TRUE),sample_data(metadata),
                  tax_table(taxonomy))    
ps_23
sample_names(ps_23)     # phyloseq object sample names
rank_names(ps_23)       # phyloseq object taxonomy categories
sample_variables(ps_23)  # phyloseq object metadata evaluated

# Filter: For the rest of the analyses we discard samples with low reads
ps.filsam <- prune_samples(sample_sums(ps_23)>=2500, ps_23) # To delete samples with less than 2500 individuals


#                 Data Rarefaction
#_______________________________________________________________________________

# This process does not afect a lot the alpha and beta diversity calculus. However calculus with diferential abundance are afected. 
set.seed(1) 
ps_23.rarefied <- rarefy_even_depth(ps.filsam, rngseed=1, sample.size=min(sample_sums(ps.filsam)), replace=F)
# 311 ASVs were removed because they are no longer present in any sample after random subsampling


#                 Calculate Alpha diversity index
#_______________________________________________________________________________

# Alpha index for rarefied phyloseq object
ps.alpha_div <- estimate_richness(ps_23.rarefied, split = TRUE, measure = c("Chao1","Fisher", "Shannon","Simpson"))

# Add the rest of variables to the rarefied diversity table
ps.alpha_div$sample.id <- rownames(ps.alpha_div) %>% as.factor()
ps.alpha_div$sample.id <- substring(ps.alpha_div$sample.id, 2)   # remove X from the sample.id
# Data frame with metadata and alpha diversity data 
ps.sampDiv <- sample_data(ps_23.rarefied) %>% unclass() %>% data.frame() %>% left_join(ps.alpha_div, by = "sample.id")

#Include number of ASVs in each sample as a column
ps.sampDiv$nasvs<-apply(otu_table(ps_23.rarefied),2, function(x) sum(x!=0)) #2 indica columas, y significa que mientras el conteo de otus sea mayor a cero lo pone en la nueva columna


##### import faith_pd index 
faith_alphaDiv <- read.table(file = "alpha-diversity_faithpd.tsv", header = TRUE, sep = "\t",
                             encoding = "UTF-8",na.strings = "") 
# call rows as sample ID 
faith_alphaDiv$X -> rownames(faith_alphaDiv)  
# Add new column with volunteer number
faith_alphaDiv$Volunteer <- c("11","12","1","1","8","8","9","9","10","10","12",
                              "13","13","14","14","15","16","16","17","17","2",
                              "18","18","19","19","20","20","21","21","22","22",
                              "2","23","23","3","3","4","4","5","5","6","6","7","7","11")
#Add new column with volunteer state
faith_alphaDiv$State <- c("Breastfeeding","Breastfeeding","Breastfeeding","Breastfeeding",
                          "Pregnant","Pregnant","Pregnant","Pregnant","Pregnant",
                          "Pregnant","Breastfeeding","Pregnant","Pregnant","Pregnant",
                          "Pregnant","Control","Control","Control","Control","Control",
                          "Breastfeeding","Control","Control","Control","Control",
                          "Control","Control","Control","Control","Breastfeeding",
                          "Breastfeeding","Breastfeeding","Breastfeeding","Breastfeeding",
                          "Breastfeeding","Breastfeeding","Breastfeeding","Breastfeeding",
                          "Pregnant","Pregnant","Breastfeeding","Breastfeeding","Pregnant",
                          "Pregnant","Breastfeeding")

faith_alphaDiv_less2B = faith_alphaDiv[-21,]
faith_alphaDiv_less2B = faith_alphaDiv_less2B[-31,]
faith_alphaDiv_less2B = faith_alphaDiv_less2B[-35,]
faith_alphaDiv_less2B = faith_alphaDiv_less2B[-35,]


#                  Statistical Analysis
#_______________________________________________________________________________

# Kluskal-wallis non parametric test

nonparametric_test <- kruskal.test(Shannon ~ State, data = ps.sampDiv_less2B)
nonparametric_test   # 0.8266

nonparametric_test <- kruskal.test(nasvs ~ State, data = ps.sampDiv_less2B)
nonparametric_test   # 0.8523

nonparametric_test <- kruskal.test(faith_pd ~ State, data = faith_alphaDiv)
nonparametric_test   # 0.52


nonparametric_test <- kruskal.test(Shannon ~ Volunteer, data = ps.sampDiv)
nonparametric_test   #0.006286

nonparametric_test <- kruskal.test(nasvs ~ Volunteer, data = ps.sampDiv)
nonparametric_test   #0.009078

nonparametric_test <- kruskal.test(faith_pd ~ Volunteer, data = faith_alphaDiv)
nonparametric_test   #0.007371


#                  Alpha and Beta diversity Graphics
#_______________________________________________________________________________
# To export data frame
df <- psmelt(ps_23.rarefied)
write.csv(df, "df_tesis.csv", row.names = TRUE)

ps.sampDiv_less2B = ps.sampDiv[-21,]
ps.sampDiv_less2B = ps.sampDiv_less2B[-31,]
ps.sampDiv_less2B = ps.sampDiv_less2B[-35,]
ps.sampDiv_less2B = ps.sampDiv_less2B[-35,]

# set the graphics background and axis
#ps.sampDiv_less2B have 2 observation of lactating woman missing to make the boxplots comparable
alpha.shannon <- ggplot(ps.sampDiv_less2B , aes(x=factor(State),y=Shannon, fill=factor(State))) #data, axis, fill the boxes 
alpha.notus <- ggplot(ps.sampDiv_less2B, aes(x=factor(State),y=nasvs, fill=factor(State)))
alpha.faithPD <- ggplot(faith_alphaDiv_less2B, aes(x=factor(State),y=faith_pd, fill=factor(State)))
alpha.shannon.volunteer <- ggplot(ps.sampDiv,aes(x=factor(Volunteer),y=Shannon, fill=factor(Volunteer)))

### # Alpha diversity: Shannon index Vs volunteer 
plot.shannon.volunteer <- ggplot(ps.sampDiv, aes(x=factor(Volunteer,level = c(1,2,3,4,6,11,12,22,23,5,7,8,9,10,13,14,15,16,17,18,19,20,21)),
                                                 y=Shannon))+
  geom_point(aes(colour=State),size=2)+
  labs(y= "Shannon index",x="Volunteer")+
  theme(axis.title = element_text(face="bold"),legend.position = "top")+
  scale_color_manual("Physiological State",values=c("#86c6be","#ff955f","#a6c64c"))+
  scale_shape_manual(values=c(15, 16, 17))
plot.shannon.volunteer

### # Alpha diversity: Observed ASVs index Vs volunteer 
ps.sampDiv <- ps.sampDiv %>% mutate(Volunteer= as.factor(Volunteer))
plot.observed.volunteer <- ggplot(ps.sampDiv, aes(x=factor(Volunteer,level = c(1,2,3,4,6,11,12,22,23,5,7,8,9,10,13,14,15,16,17,18,19,20,21)),
                                                  y=nasvs))+
  geom_point(aes(colour=State),size=2) +
  labs(y= "Observed ASVs",x="Volunteer")+
  theme(axis.title = element_text(face="bold"),legend.position = "top")+
  scale_color_manual("Physiological State",values=c("#86c6be","#ff955f","#a6c64c"))
plot.observed.volunteer

### # Alpha diversity: Faith_pd index Vs volunteer 
faith_alphaDiv <- faith_alphaDiv %>% mutate(Volunteer= as.factor(Volunteer))
plot.faithPD.volunteer <- ggplot(faith_alphaDiv, aes(x=factor(Volunteer,level = c(1,2,3,4,6,11,12,22,23,5,7,8,9,10,13,14,15,16,17,18,19,20,21)),
                                                     y=faith_pd))+
  geom_point(aes(colour=State),size=2) +
  labs(y= "Faith_pd",x="Volunteer")+
  theme(axis.title = element_text(face="bold"),legend.position = "top")+
  scale_color_manual("Physiological State",values=c("#86c6be","#ff955f","#a6c64c"))
plot.faithPD.volunteer

#join figures into the same plot_Alpha diversity Vs Volunteer
figure2A <- ggarrange(plot.shannon.volunteer, plot.observed.volunteer,plot.faithPD.volunteer,
                      labels= c("A",","),
                      ncol = 2, nrow = 2,common.legend = TRUE)
figure2A

### Alpha diversity: Shannon index Vs physiological state
plot.shannon <- alpha.shannon + geom_boxplot() +
  #bracket de significancia
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.title.x=element_blank(), legend.position="top",axis.title = element_text(face="bold"))+  #configuracion estetica de los ejes 
  labs(y="Shannon index")+ #poner titulo
  scale_fill_manual("Physiological State", values=c("#86c6be","#ff955f","#a6c64c")) #cambio de colores y leyenda 
plot.shannon

### Alpha diversity: Observed ASVs index Vs physiological state
plot.nasvs <- alpha.notus + geom_boxplot() +
  #geom_signif(annotation =c("0.021"),xmin = 1, xmax = 2,y_position = 160) + #bracket de significancia
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.title.x=element_blank(), legend.position="top",axis.title = element_text(face="bold"))+  #configuracion estetica de los ejes 
  labs(y="Observed ASVs") +
  scale_fill_manual("Physiological State", values=c("#86c6be","#ff955f","#a6c64c")) #cambio de colores y leyenda 
plot.nasvs

### Alpha diversity: Faith_pd index Vs physiological state
plot.faith <- alpha.faithPD + geom_boxplot() +
  #geom_signif(annotation =c("0.021"),xmin = 1, xmax = 2,y_position = 160) + #bracket de significancia
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.title.x=element_blank(), legend.position="top",axis.title = element_text(face="bold"))+  #configuracion estetica de los ejes 
  labs(y="Faith_pd") +
  scale_fill_manual("Physiological State", values=c("#86c6be","#ff955f","#a6c64c")) #cambio de colores y leyenda 
plot.faith

#join figures into the same plot_Alpha diversity Vs State
figure2B <- ggarrange(plot.shannon, plot.nasvs,plot.faith,
                      labels= c("B","",""),
                      ncol = 2, nrow = 2,common.legend = TRUE, tittle=FALSE)
figure2B


####  BETADIVERSITY 
out.bray <- ordinate(ps_23.rarefied, method = "MDS", distance = "bray")
beta.bray <- plot_ordination(ps_23.rarefied, out.bray, color="State")  #graficar ordination por colores

plot.beta.bray <- beta.bray + geom_point(size=4, alpha=0.75, na.rm=T) +  #sobreponer un punto mas grande 
  labs(col= "Physiological State") +  #titulo de la grafica y de la leyenda pegar , label="Volunteer" dentro del parentesis para que los numeros sean mas pequeños
  scale_colour_manual(values=c("#86c6be","#ff955f","#a6c64c")) + #cambiar colores y poner nombre a cada grupo en leyenda
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position="top")+
  geom_text(mapping = aes(label = Volunteer), size = 3, vjust = 2)
plot.beta.bray




#        Correlation matrix
#________________________________________________________
# Create a new dataframe (subsample de original dataset)
lessData <- subset(ps.sampDiv, select = c(State,Total.Colesterol, LDL,HDL,Trigliceridos,
                                          BMI,AGE,Shannon,nasvs))
# Set the variables as numeric and State as factor
lessData <- lessData %>% mutate(State= as.factor(State)) 
#%>% mutate(Trigliceridos= as.numeric(Trigliceridos)) %>% mutate(LDL= as.numeric(LDL)) %>% mutate(HDL= as.numeric(HDL))

# Createa new column indicating the state with a number
# Lactante=1 Embarazada=2 Control=3
lessData$StateN <- c(1,1,1,1,2,2,2,2,2,2,1,2,2,2,2,3,3,3,3,3,1,3,3,3,3,3,3,3,3,1,1,1,1,1,1,1,1,1,2,2,1,1,2,2,1)


# Delete first column and row 35 does not have age and BMI
lessData = lessData[,-1]
lessData = lessData[-37,]
lessData = lessData[-37,]

# Order personalizado de las columnas
lessData <- lessData[, c(9, 1,2,3,4,5,6,7,8)] 

# Matriz de correlación y su respectiva visualización
Matriz_correlación <- cor(lessData)

# Visualización
corrplot(Matriz_correlación, type="upper", number.font = 2,  order="original", diag=FALSE,
         method="color", addCoef.col="black", tl.col="black")

# Draftman´s plot (To Correlate Metadata)
library(GGally)
lessData <- lessData %>% mutate(StateN= as.factor(StateN))
ggpairs(lessData,aes(fill=StateN, color=StateN))

