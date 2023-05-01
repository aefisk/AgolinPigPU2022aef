##############################################################
# title: "Alpha diversity in R - qiime2 output"
# author: "ANSC595"
# date: "March 16, 2021"
##############################################################


###Set your working directory
#setwd('C:/Users/alyxf/Desktop/Agolin_R/Agolin')

# Modified from the original online version available at 
# http://rpubs.com/dillmcfarlan/R_microbiotaSOP

# and Tutorial: Integrating QIIME2 and R for data visualization 
# and analysis using qiime2R
# https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121

##Goal
# The goal of this tutorial is to demonstrate basic analyses of microbiota data to determine if and how communities differ by variables of interest. In general, this pipeline can be used for any microbiota data set that has been clustered into operational taxonomic units (OTUs).
#
# This tutorial assumes some basic statistical knowledge. Please consider if your data fit the assumptions of each test (normality? equal sampling? Etc.). If you are not familiar with statistics at this level, we strongly recommend collaborating with someone who is. The incorrect use of statistics is a pervasive and serious problem in the sciences so don't become part of the problem! That said, this is an introductory tutorial and there are many, many further analyses that can be done with microbiota data. Hopefully, this is just the start for your data!

##Data
# The data used here are from the qiime2 moving pictures tutorial. 
# Please see their online tutorial for an explanation of the dataset.

# Data manipulation
## Load Packages

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20

library(tidyverse)
library(qiime2R)
library(ggpubr)

##Load Data
# In the code, the text before = is what the file will be called in R. 
# Make this short but unique as this is how you will tell R to use this 
# file in later commands.

# header: tells R that the first row is column names, not data
# row.names: tells R that the first column is row names, not data
# sep: tells R that the data are tab-delimited. 
# If you had a comma-delimited file, you would us sep=","


#How to load a file into R
metadata <- read.delim("Agolin-sample-metadata.tsv", sep = "\t", header = T, quote = "", stringsAsFactors = F)
metadata[1,]
metadata[,1]
# When subsetting, the first number is the row and after the comma is the column
metadata <- metadata[-1,]
# Load data

metadata<-read_q2metadata("Agolin-sample-metadata.tsv")
str(metadata)
levels(metadata$`trt`)
colnames(metadata)[5] <- "Sow_ID"
colnames(metadata)[2] <- "Day"
colnames(metadata)[4] <- "Treatment"
colnames(metadata)[6] <- "Parity_Class"
str(metadata)

evenness0 = read_qza("Results/day0-core-metrics-results/evenness_vector.qza")
evenness0<-evenness0$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

observed_features0 = read_qza("Results/day0-core-metrics-results/observed_features_vector.qza")
observed_features0<-observed_features0$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

shannon0 = read_qza("Results/day0-core-metrics-results/shannon_vector.qza")
shannon0<-shannon0$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

faith_pd = read_qza("Results/day0-core-metrics-results/faith_pd_vector.qza")
faith_pd0<-faith_pd$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged\
faith_pd0<-faith_pd0[,-1]
colnames(faith_pd0)[1]<- "SampleID"
colnames(faith_pd0)[2]<- "Faith_PD"

## Clean up the data
# You can look at your data by clicking on it in the upper-right 
# quadrant "Environment"


###Alpha Diversity tables
# These tables will be merged for convenience and added to the 
# metadata table as the original tutorial was organized.

alpha_diversity0 = merge(x=faith_pd0, y=evenness0, by.x = "SampleID", by.y = "SampleID")
alpha_diversity0 = merge(alpha_diversity0, observed_features0, by.x = "SampleID", by.y = "SampleID")
alpha_diversity0 = merge(alpha_diversity0, shannon0, by.x = "SampleID", by.y = "SampleID")
metadata0 = merge(metadata, alpha_diversity0, by.x = "SampleID", by.y = "SampleID")
row.names(metadata0) <- metadata0$SampleID
#metadata = metadata[,-1]
str(metadata0)


#Alpha-diversity
# Alpha-diversity is within sample diversity. It is how many 
# different species (OTUs) are in each sample (richness) and how 
# evenly they are distributed (evenness), which together are diversity. 
# Each sample has one value for each metric.


##Explore alpha metrics
# Now we will start to look at our data. We will first start with 
# alpha-diversity and richness. 
#
# You want the data to be roughly normal so that you can run ANOVA 
# or t-tests. If it is not normally distributed, you will need to 
# consider if you should normalize the data or usenon-parametric 
# tests such as Kruskal-Wallis.

# Here, we see that none of the data are normally distributed, 
# with the exception of "Faith" and "Observed Features".


#Plots
hist(metadata0$shannon_entropy, main="Shannon Diversity Day 0", xlab="", breaks=10)
hist(metadata0$Faith_PD, main="Faith phylogenetic Diversity Day 0", xlab="", breaks=10)
hist(metadata0$pielou_evenness, main="Evenness Day 0", xlab="", breaks=10)
hist(as.numeric(metadata0$observed_features), main="Observed Features Day 0", xlab="", breaks=10)

#Plots the qq-plot for residuals
ggqqplot(metadata0$shannon_entropy, title = "Shannon D. 0")
ggqqplot(metadata0$Faith_PD, title = "Faith PD D. 0")
ggqqplot(metadata0$pielou_evenness, title = "Evenness D. 0")
ggqqplot(metadata0$observed_features, title = "Observed Features D. 0")



# To test for normalcy statistically, we can run the Shapiro-Wilk 
# test of normality.

shapiro.test(metadata0$shannon_entropy)
shapiro.test(metadata0$Faith_PD)
shapiro.test(metadata0$pielou_evenness)
shapiro.test(metadata0$observed_features)

# The null hypothesis of these tests is that “sample distribution 
# is normal”. If the test is significant, the distribution is non-normal.

# We see that, as expected from the graphs, shannon and evenness 
# are normally distributed.


#Overall, for alpha-diversity:

# ANOVA, t-test, or general linear models with the normal distribution 
# are used when the data is roughly normal. Transforming the data to 
# achieve a normal distribution could also be completed.
#
# Kruskal-Wallis, Wilcoxon rank sum test, or general linear models 
# with another distribution are used when the data is not normal or if 
# the n is low, like less than 30.

# Our main variables of interest are

# body site: gut, tongue, right palm, left palm
# subject: 1 and 2
# month-year: 10-2008, 1-2009, 2-2009, 3-2009, 4-2009

## Categorical variables
# Now that we know which tests can be used, let's run them. 

## Normally distributed metrics

# Since it's the closest to normalcy, we will use **Evenness** as an 
#example. First, we will test body site, which is a categorical variable 
# with more than 2 levels. Thus, we run ANOVA. If age were only two 
# levels, we could run a t-test

# Does body site impact the Evenness of the microbiota?

#Run the ANOVA and save it as an object
aov.evenness.trt.d0 = aov(pielou_evenness ~ Treatment, data=metadata0)
#Call for the summary of that ANOVA, which will include P-values
summary(aov.evenness.trt.d0)

#To do all the pairwise comparisons between groups and correct for multiple comparisons, we run Tukey's honest significance test of our ANOVA.

TukeyHSD(aov.evenness.trt.d0)

# We clearly see that the evenness between hands and gut are different. 
# When we plot the data, we see that evenness decreases in the gut 
# compared to palms.

levels(metadata0$Treatment)
#Re-order the groups because the default is alphabetical order
metadata0$TreatmentD0.ord = factor(metadata0$Treatment, c("CON", "AGO"))
levels(metadata0$TreatmentD0.ord)

#Plot
boxplot(pielou_evenness ~ TreatmentD0.ord, data=metadata0, ylab="Peilou Evenness Day 0")

evenness_bp0 <- ggplot(metadata0, aes(Treatment, pielou_evenness, color=Treatment)) + 
  geom_boxplot() + 
  #facet_grid(~TreatmentD0.ord) +
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("output/evennessd0.png", evenness_bp0, height = 3, width = 3)

# Now, the above graph is kind of not correct. Our test and our graphic do not exactly match. ANOVA and Tukey are tests based on the mean, but the boxplot plots the median. Its not wrong, its just not the best method. Unfortunately plotting the average and standard deviation is a little complicated.

evenness_summary0 <- metadata0 %>% # the names of the new data frame and the data frame to be summarised
  group_by(TreatmentD0.ord) %>%   # the grouping variable
  summarise(mean_evenness0 = mean(pielou_evenness),  # calculates the mean of each group
            sd_evenness0 = sd(pielou_evenness), # calculates the standard deviation of each group
            n_evenness0 = n(),  # calculates the sample size per group
            se_evenness0 = sd(pielou_evenness)/sqrt(n())) # calculates the standard error of each group

# We can now make a bar plot of means vs body site, with standard 
# deviations or standard errors as the error bar. The following code 
# uses the standard deviations.

evenness_se0 <- ggplot(evenness_summary0, aes(TreatmentD0.ord, mean_evenness0, fill = TreatmentD0.ord)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_evenness0 - se_evenness0, ymax = mean_evenness0 + se_evenness0), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Pielou's evenness  ± s.e.", x = "") 

ggsave("output/evenness_se0.png", evenness_se0, height = 2.5, width = 3)

# or with ggplot2
boxplot(Faith_PD ~ TreatmentD0.ord, data=metadata0, ylab="Faith Phylogenetic Diversity Day 0")

faith_pd0_plot <- ggplot(metadata0, aes(TreatmentD0.ord, Faith_PD)) + 
  geom_boxplot(aes(color = Treatment)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Faith Phylogenetic Diversity D. 0", x = "") 
ggsave("output/pd0.png", faith_pd0_plot, height = 2.5, width = 3)

##Continuous variables
# For continuous variables, we use general linear models, specifying 
# the distribution that best fits our data.

# **Normally distributed metrics**

# Since days.since.experiment.start is a continuous variable, we run a 
# general linear model. We will again use evenness as our roughly normal 
# metric. The default of `glm` and `lm` is the normal distribution so we 
# don't have to specify anything.


glm.evenness.d0 = glm(pielou_evenness ~ TreatmentD0.ord, data=metadata0)
summary(glm.evenness.d0)

#The output let's us know that the intercept of our model is significantly different from 0 but our slope (*e.g.* our variable of interest) is not. This makes sense when we look at the data.

plot(pielou_evenness ~ TreatmentD0.ord, data=metadata0)
#Add the glm best fit line
plot(pielou_evenness ~ TreatmentD0.ord, data=metadata0) + abline(glm.evenness.d0)

# **Non-normally distributed metrics**

# We will again use a *general linear model* for our non-normally 
# distributed metric Faith_pd. However, this time, we change the 
# distribution from normal to something that fits the data better. 

# But which distribution should we choose? In statistics, there is no 
# one "best" model. There are only good and better models. We will use 
# the plot() function to compare two models and pick the better one.

# Quasipoisson (log) distribution
qp.faith.d0 = glm(Faith_PD ~ TreatmentD0.ord, data=metadata0, family="quasipoisson")
plot(qp.faith.d0, which=c(1,2))

# What we're looking for is no pattern in the Residuals vs. Fitted graph 
# ("stars in the sky"), which shows that we picked a good distribution 
# family to fit our data. We also want our residuals to be normally 
# distributed, which is shown by most/all of the points falling on the 
# line in the Normal Q-Q plot.

# While it's still not perfect, the quasipoisson fits much better. 
# In the residuals vs fitted graph, the y axis is from -2 to 4  whereas 
# the axis with gaussian was from -5 to 10. So, we will use quasipoisson 
# and see that ADG does not to correlate to Chao richness.
summary(qp.faith.d0)

#Plot
plot(log(Faith_PD) ~ TreatmentD0.ord, data=metadata0, ylab="ln(Faith Phylo. Diversity D. 0)")
plot(log(Faith_PD) ~ TreatmentD0.ord, data=metadata0, ylab="ln(Faith Phylo. Diversity D. 0)") + abline(qp.faith.d0)
