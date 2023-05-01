
#Load the packages
install.packages("tidyverse")
install.packages("vegan")
install.packages("devtools")
library(devtools)
devtools::install_github("jbisanz/qiime2R")



library(tidyverse)
library(vegan)
library(qiime2R)


##############################################
#Set UP
#
#These are the things that  we need from Qiime:
#
#sample-metadata.tsv
#core-metrics-results/bray_curtis_pcoa_results.qza
#core-metrics-results/weighted_unifrac_pcoa_results.qza
#core-metrics-results/rarefied_table.qza
#rooted-tree.qza
#taxonomy.qza
#core-metrics-results/evenness_vector.qza
#core-metrics-results/faith_pd_vector.qza
#core-metrics-results/observed_otus_vector.qza
#core-metrics-results/shannon_vector.qza
#
# To get these files you need to scp them from the cluster:
#
# first on  your laptop cd to the directory where you want to save them.
# Then use this code for our example dataset today:
# mkdir core-metrics-results/
# scp john2185@bell.rcac.purdue.edu:/depot/microbiome/data/2021_ANSC595/john2185/qiime/moving_pictures_pipeline/* .
# scp john2185@bell.rcac.purdue.edu:/depot/microbiome/data/2021_ANSC595/john2185/qiime/moving_pictures_pipeline/core-metrics-results/* core-metrics-results/.
##############################################


###Set your working directory
setwd("C:/Users/alyxf/Desktop/Agolin_R/Agolin")

list.files()

if(!dir.exists("output-PiCrust2"))
  dir.create("output-PiCrust2")

#How to load a file into R
metadata16 <- read.delim("Agolin-sample-metadata.tsv", sep = "\t", header = T, quote = "", stringsAsFactors = F)
metadata16[1,]
metadata16[,1]
# When subsetting, the first number is the row and after the comma is the column
metadata16 <- metadata16[-1,]

#Now the qiime2R method
metadata16<-read_q2metadata("Agolin-sample-metadata.tsv")
str(metadata16)
levels(metadata16$`trt`)
colnames(metadata16)[5] <- "Sow_ID"
colnames(metadata16)[2] <- "Day"
colnames(metadata16)[4] <- "Treatment"
colnames(metadata16)[6] <- "Parity_Class"
str(metadata16)

row.names(metadata16) <- metadata16[,1]
row.names(metadata16) <- metadata16$SampleID
#metadata <- metadata[,-1]
row.names(metadata16)

bc_PCoAPCd16<-read_qza("q2-picrust2_output_day-16/pathabun_core_metrics_out/bray_curtis_pcoa_results.qza")
jac_PCoAPC16<-read_qza("q2-picrust2_output_day-16/pathabun_core_metrics_out/jaccard_pcoa_results.qza")

trt_colors <- c("green", "black", "red", "orange")

bc_meta16 <- bc_PCoAPCd16$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata16, by = c("SampleID" = "SampleID"))

jac_meta16 <- jac_PCoAPC16$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata16, by = c("SampleID" = "SampleID"))

# Now we are going to make an ordination plot
ggplot(bc_meta16, aes(x=PC1, y=PC2, color=Treatment)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab("PC1 (22.65%)") +
  ylab("PC2 (16.168%)") +
  scale_color_manual(values=c("black", "green"), name = "Treatment")

# Now we are going to make our code a little more re-usable
my_column <- "Treatment"
#my_column <- "DietTreatment"

ggplot(bc_meta16, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  facet_grid(~Sow_ID) +
  xlab(paste0("PC1 (", round(100*bc_PCoAPCd16$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoAPCd16$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=trt_colors, name = my_column)
ggsave(paste0("output-PiCrust2/BC2-basic_day16", my_column,".tiff"), height=4, width=8, device="tiff") # save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),bc_meta16,mean)
colnames(centroids)[1] <- "Treatment"

ggplot(bc_meta16, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoAPCd16$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoAPCd16$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=trt_colors, name = my_column)
ggsave(paste0("output-PiCrust2/BC-ellipse_day16", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches


#same thing but Jaccard
ggplot(jac_meta16, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  facet_grid(~Sow_ID) +
  xlab(paste0("PC1 (", round(100*jac_PCoAPC16$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*jac_PCoAPC16$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=trt_colors, name = my_column)
ggsave(paste0("output-PiCrust2/Jac-basic_day16", my_column,".tiff"), height=4, width=8, device="tiff") # save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),jac_meta16,mean)
colnames(centroids)[1] <- "Treatment"

ggplot(jac_meta16, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*jac_PCoAPC16$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*jac_PCoAPC16$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=trt_colors, name = my_column)
ggsave(paste0("output-PiCrust2/Jac-ellipse_day16", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches



#Bray Curtis Permonova Table
bc_dist_mat16<-read_qza("q2-picrust2_output_day-16/pathabun_core_metrics_out/bray_curtis_distance_matrix.qza")
bc_dm16 <- as.matrix(bc_dist_mat16$data) 
rownames(bc_dm16) == metadata16$SampleID ## all these values need to be "TRUE"
metadata_sub16 <- metadata16[match(rownames(bc_dm16),metadata16$SampleID),]
rownames(bc_dm16) == metadata_sub16$SampleID ## all these values need to be "TRUE"

PERMANOVA_outBC16 <- adonis2(bc_dm16 ~ Treatment, data = metadata_sub16)

write.table(PERMANOVA_outBC16,"output-PiCrust2/BCTrt_day-16_Adonis_overall.csv",sep=",", row.names = TRUE) 

#Jaccard Permonova Table
Jac_dist_mat16<-read_qza("q2-picrust2_output_day-16/pathabun_core_metrics_out/jaccard_distance_matrix.qza")
jac_dm16 <- as.matrix(Jac_dist_mat16$data) 
rownames(jac_dm16) == metadata16$SampleID ## all these values need to be "TRUE"
metadata_sub16 <- metadata16[match(rownames(jac_dm16),metadata16$SampleID),]
rownames(jac_dm16) == metadata_sub16$SampleID ## all these values need to be "TRUE"

PERMANOVA_outJAC16 <- adonis2(jac_dm16 ~ Treatment, data = metadata_sub16)

write.table(PERMANOVA_outJAC16,"output-PiCrust2/JACTrt_day-16_Adonis_overall.csv",sep=",", row.names = TRUE)



######################################################################################
##  Pairwise adonis function
##  we can also performe a pairwise comparison with the function 
##  Pairwise Adonis funtion by edro Martinez Arbizu & Sylvain Monteux
##  https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R
#######################################################################################

pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 

Treatment_PairBC16 <- pairwise.adonis2(bc_dm16 ~ Treatment, data = metadata_sub16)
write.table(Treatment_PairBC16,"output-PiCrust2/BCTrt_day-16_Adonis_pairwise.csv",sep=",", row.names = TRUE)

Treatment_PairJac16 <- pairwise.adonis2(jac_dm16 ~ Treatment, data = metadata_sub16)
write.table(Treatment_PairJac16,"output-PiCrust2/JACTrt_day-16_Adonis_pairwise.csv",sep=",", row.names = TRUE)

