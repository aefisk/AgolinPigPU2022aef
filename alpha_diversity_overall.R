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

evenness = read_qza("Results/core-metrics-results/evenness_vector.qza")
evenness<-evenness$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

observed_features = read_qza("Results/core-metrics-results/observed_features_vector.qza")
observed_features<-observed_features$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

shannon = read_qza("Results/core-metrics-results/shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

faith_pd_read = read_qza("Results/core-metrics-results/faith_pd_vector.qza")
faith_pd<-faith_pd$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged\
faith_pd<-faith_pd[,-1]
colnames(faith_pd)[1]<- "SampleID"
colnames(faith_pd)[2]<- "Faith_PD"

## Clean up the data
# You can look at your data by clicking on it in the upper-right 
# quadrant "Environment"

# You always need to check the data types in your tables to make 
# sure they are what you want. We will now change some data types 
# in the meta now

str(metadata)
observed_features$observed_features_num <- lapply(observed_features$observed_features, as.numeric)
observed_features$observed_features <- as.numeric(observed_features$observed_features)
str(observed_features)


###Alpha Diversity tables
# These tables will be merged for convenience and added to the 
# metadata table as the original tutorial was organized.

alpha_diversity = merge(x=faith_pd, y=evenness, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, observed_features, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, shannon, by.x = "SampleID", by.y = "SampleID")
metadata = merge(metadata, alpha_diversity, by.x = "SampleID", by.y = "SampleID")
row.names(metadata) <- metadata$SampleID
#metadata = metadata[,-1]
str(metadata)


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
hist(metadata$shannon_entropy, main="Shannon Diversity", xlab="", breaks=10)
hist(metadata$Faith_PD, main="Faith phylogenetic Diversity", xlab="", breaks=10)
hist(metadata$pielou_evenness, main="Evenness", xlab="", breaks=10)
hist(as.numeric(metadata$observed_features), main="Observed Features", xlab="", breaks=10)

#Plots the qq-plot for residuals
ggqqplot(metadata$shannon_entropy, title = "Shannon")
ggqqplot(metadata$Faith_PD, title = "Faith PD")
ggqqplot(metadata$pielou_evenness, title = "Evenness")
ggqqplot(metadata$observed_features, title = "Observed Features")



# To test for normalcy statistically, we can run the Shapiro-Wilk 
# test of normality.

shapiro.test(metadata$shannon_entropy)
shapiro.test(metadata$Faith_PD)
shapiro.test(metadata$pielou_evenness)
shapiro.test(metadata$observed_features)

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
aov.evenness.trt = aov(pielou_evenness ~ Treatment, data=metadata)
#Call for the summary of that ANOVA, which will include P-values
summary(aov.evenness.trt)

#To do all the pairwise comparisons between groups and correct for multiple comparisons, we run Tukey's honest significance test of our ANOVA.

TukeyHSD(aov.evenness.trt)

# We clearly see that the evenness between hands and gut are different. 
# When we plot the data, we see that evenness decreases in the gut 
# compared to palms.

levels(metadata$Treatment)
#Re-order the groups because the default is alphabetical order
metadata$Treatment.ord = factor(metadata$Treatment, c("CON", "AGO"))
levels(metadata$Treatment.ord)

#Plot
boxplot(pielou_evenness ~ Treatment.ord, data=metadata, ylab="Peilou Evenness")

evenness_bp <- ggplot(metadata, aes(Treatment, pielou_evenness, color=Treatment)) + 
  geom_boxplot() + 
  #facet_grid(~Treatment.ord) +
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("output/evenness.png", evenness_bp14, height = 3, width = 3)

# Now, the above graph is kind of not correct. Our test and our graphic do not exactly match. ANOVA and Tukey are tests based on the mean, but the boxplot plots the median. Its not wrong, its just not the best method. Unfortunately plotting the average and standard deviation is a little complicated.

evenness_summary <- metadata %>% # the names of the new data frame and the data frame to be summarised
  group_by(Treatment.ord) %>%   # the grouping variable
  summarise(mean_evenness = mean(pielou_evenness),  # calculates the mean of each group
            sd_evenness = sd(pielou_evenness), # calculates the standard deviation of each group
            n_evenness = n(),  # calculates the sample size per group
            se_evenness = sd(pielou_evenness)/sqrt(n())) # calculates the standard error of each group

# We can now make a bar plot of means vs body site, with standard 
# deviations or standard errors as the error bar. The following code 
# uses the standard deviations.

evenness_se <- ggplot(evenness_summary, aes(Treatment.ord, mean_evenness, fill = Treatment.ord)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_evenness - se_evenness, ymax = mean_evenness + se_evenness), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Pielou's evenness  ± s.e.", x = "") 

ggsave("output/evenness_se.png", evenness_se, height = 2.5, width = 3)


## **Non-normally distributed metrics**

# We will use **Faith's phylogenetic diversity** here. Since body site 
# is categorical, we use Kruskal-Wallis (non-parametric equivalent of 
# ANOVA). If we have only two levels, we would run Wilcoxon rank sum 
# test (non-parametric equivalent of t-test)

kruskal.test(faith_pd ~ Treatment.ord, data= metadata)

# We can test pairwise within the age groups with Wilcoxon Rank Sum 
# Tests. This test has a slightly different syntax than our other tests

pairwise.wilcox.test(metadata$faith_pd, metadata$Treatment.ord, p.adjust.method="BH")

# Like evenness, we see that pd also increases with age.

#Plot
boxplot(faith_pd14 ~ TreatmentD14.ord, data=metadata14, ylab="Faith phylogenetic diversity")

# or with ggplot2

faith_pd14 <- ggplot(metadata14, aes(body.site.ord, faith_pd)) + 
  geom_boxplot(aes(color = body.site.ord)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Faith Phylogenetic Diversity", x = "") 
ggsave("output/pd.png", faith_pd, height = 3, width = 3)

##Continuous variables
# For continuous variables, we use general linear models, specifying 
# the distribution that best fits our data.

# **Normally distributed metrics**

# Since days.since.experiment.start is a continuous variable, we run a 
# general linear model. We will again use evenness as our roughly normal 
# metric. The default of `glm` and `lm` is the normal distribution so we 
# don't have to specify anything.

# Does days.since.experiment.start impact evenness of the microbiota?

glm.evenness.timed14 = glm(pielou_evenness ~ days.since.experiment.start, data=meta)
summary(glm.evenness.time)

#The output let's us know that the intercept of our model is significantly different from 0 but our slope (*e.g.* our variable of interest) is not. This makes sense when we look at the data.

plot(pielou_evenness ~ days.since.experiment.start, data=meta)
#Add the glm best fit line
plot(pielou_evenness ~ days.since.experiment.start, data=meta) + abline(glm.evenness.time)

# **Non-normally distributed metrics**

# We will again use a *general linear model* for our non-normally 
# distributed metric Faith_pd. However, this time, we change the 
# distribution from normal to something that fits the data better. 

# But which distribution should we choose? In statistics, there is no 
# one "best" model. There are only good and better models. We will use 
# the plot() function to compare two models and pick the better one.

# First, the Gaussian (normal) distribution, which we already know is a bad fit.

gaussian.faith.time = glm(faith_pd ~ days.since.experiment.start, data=meta, family="gaussian")
plot(gaussian.faith.time, which=c(1,2))

# Quasipoisson (log) distribution
qp.faith.time = glm(faith_pd ~ days.since.experiment.start, data=meta, family="quasipoisson")
plot(qp.faith.time, which=c(1,2))

# What we're looking for is no pattern in the Residuals vs. Fitted graph 
# ("stars in the sky"), which shows that we picked a good distribution 
# family to fit our data. We also want our residuals to be normally 
# distributed, which is shown by most/all of the points falling on the 
# line in the Normal Q-Q plot.

# While it's still not perfect, the quasipoisson fits much better. 
# In the residuals vs fitted graph, the y axis is from -2 to 4  whereas 
# the axis with gaussian was from -5 to 10. So, we will use quasipoisson 
# and see that ADG does not to correlate to Chao richness.
summary(qp.faith.time)

# Plotting this we see that, indeed, there is a trend toward correlation between Faith_pd and days.since.experiment.start.

#Plot
plot(log(faith_pd) ~ days.since.experiment.start, data=meta, ylab="ln(Faith Phylo. Diversity)")
plot(log(faith_pd) ~ days.since.experiment.start, data=meta, ylab="ln(Faith Phylo. Diversity)") + abline(qp.faith.time)

