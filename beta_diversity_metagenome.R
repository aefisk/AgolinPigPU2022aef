
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
setwd("C:/Users/alyxf/Desktop/Agolin_R/")

list.files()

if(!dir.exists("output"))
  dir.create("output")

#How to load a file into R
metadatagenome <- read.delim("Agolin/Agolin-sample-metadata.tsv", sep = "\t", header = T, quote = "", stringsAsFactors = F)
metadatagenome[1,]
metadatagenome[,1]
# When subsetting, the first number is the row and after the comma is the column
metadatagenome <- metadatagenome[-1,]

#Now the qiime2R method
metadatagenome<-read_q2metadata("Agolin/Agolin-sample-metadata.tsv")
str(metadatagenome)
levels(metadatagenome$`trt`)
colnames(metadatagenome)[5] <- "Sow_ID"
colnames(metadatagenome)[2] <- "Day"
colnames(metadatagenome)[4] <- "Treatment"
colnames(metadatagenome)[6] <- "Parity_Class"
str(metadatagenome)

row.names(metadatagenome) <- metadatagenome[,1]
row.names(metadatagenome) <- metadatagenome$SampleID
#metadata <- metadata[,-1]
row.names(metadatagenome)

bc_PCoAgenome<-read_qza("./q2-picrust2_output/pathabun_core_metrics_out/bray_curtis_pcoa_results.qza")
jac_PCoAgenome<-read_qza("./q2-picrust2_output/pathabun_core_metrics_out/jaccard_pcoa_results.qza")

trt_colors <- c("green", "black", "red", "orange")

bc_metagenome <- bc_PCoAgenome$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadatagenome, by = c("SampleID" = "SampleID"))

jac_metagenome <- jac_PCoAgenome$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadatagenome, by = c("SampleID" = "SampleID"))

# Now we are going to make an ordination plot
ggplot(bc_metagenome, aes(x=PC1, y=PC2, color=Treatment)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab("PC1 (22.65%)") +
  ylab("PC2 (16.148%)") +
  scale_color_manual(values=c("black", "green"), name = "Treatment")

# Now we are going to make our code a little more re-usable
my_column <- "Treatment"
#my_column <- "DietTreatment"

ggplot(bc_metagenome, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  facet_grid(~Sow_ID) +
  xlab(paste0("PC1 (", round(100*bc_PCoAgenome$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoAgenome$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=trt_colors, name = my_column)
ggsave(paste0("output/BC2-basic_daygenome", my_column,".tiff"), height=4, width=8, device="tiff") # save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),bc_metagenome,mean)
colnames(centroids)[1] <- "Treatment"

ggplot(bc_metagenome, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoAgenome$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoAgenome$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=trt_colors, name = my_column)
ggsave(paste0("output/BC-ellipse_daygenome", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(jac_metagenome, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  facet_grid(~Sow_ID) +
  xlab(paste0("PC1 (", round(100*jac_PCoAgenome$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*jac_PCoAgenome$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=trt_colors, name = my_column)
ggsave(paste0("output/Jac-basic_daygenome", my_column,".tiff"), height=4, width=8, device="tiff") # save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),jac_metagenome,mean)
colnames(centroids)[1] <- "Treatment"

ggplot(jac_metagenome, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*jac_PCoAgenome$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*jac_PCoAgenome$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=trt_colors, name = my_column)
ggsave(paste0("output/Jac-ellipse_daygenome", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches
