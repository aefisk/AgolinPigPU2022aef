## This script goes through making a taxa bar plot and then running
## DESeq2 to find differentially abundant ASVs


# for help installing phyloseq, see this website
# https://bioconductor.org/packages/release/bioc/html/phyloseq.html

# to install phyloseq:
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("phyloseq")

library(qiime2R)
library(phyloseq)
library(zoo)
library(tidyverse)


##############################################
#Set UP
#
#These are the things that  we need from Qiime:
#
#sample-metadata.tsv
#core-metrics-results/rarefied_table.qza
#rooted-tree.qza
#taxonomy.qza
##############################################

setwd('C:/Users/alyxf/Desktop/Agolin_R/Agolin')
list.files()

if(!dir.exists("output"))
  dir.create("output")

metadata<-read_q2metadata("Agolin-sample-metadata.tsv")
str(metadata)
colnames(metadata)[4] = "Treatment"
levels(metadata$Treatment)
metadata$trt.ord = factor(metadata$Treatment, c("CON", "AGO"))
levels(metadata$trt.ord)
colnames(metadata)[5] <- "Sow_ID"
colnames(metadata)[2] <- "Day"
colnames(metadata)[4] <- "Treatment"
colnames(metadata)[6] <- "Parity_Class"
str(metadata)

row.names(metadata) <- metadata[ ,1]
#metadata <- metadata[,-1]
row.names(metadata)

##Qiime2r method of reading in the taxonomy files
taxonomy<-read_qza("taxonomy.qza")
head(taxonomy$data)

tax.clean<-parse_taxonomy(taxonomy$data)
head(tax.clean)

#All this is OK except that in future use of the taxonomy table, 
#these ASVs will be ignored because they are not classified. Why 
#are ASVs not classified? Its because there is not a close enough 
#match in the database. Just because there is not a good match in 
#the database does not mean they don’t exist, so I wanted to make 
#sure this data was not lost. So in my new code, from lines 200 – 224 
#I make it so that ASVs that are unclassified at any level are 
#classified as the lowest taxonomic level for which there is a 
#classification.
#Next, all these `NA` classifications with the last level that was 
#classified

tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:0] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:0] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:0] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:0] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:0] <- family
  } else if (tax.clean[i,0] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}



#################################################################
##Taxa barplot
#################################################################

physeqd8 <- qza_to_phyloseq(
  features="Results/day8-core-metrics-results/rarefied_table.qza",
  tree="rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "Agolin-sample-metadata.tsv"
)


#First get the OTU table from physeq
physeq_otu_tabled8 <- data.frame(otu_table(physeqd8), check.names = F)

tax.cleand8 = tax.clean[row.names(tax.clean) %in% rownames(physeq_otu_tabled8),]
metadata.filteredd8 = metadata[row.names(metadata) %in% colnames(physeq_otu_tabled8),]

#Assign as variables to be feed into phyloseq
OTU.physeqd8 = otu_table(as.matrix(physeq_otu_tabled8), taxa_are_rows=TRUE)

#our edited and formatted taxonomy table from the top of this script
tax.physeqd8 = tax_table(as.matrix(tax.cleand8))    
meta.physeqd8 = sample_data(metadata.filteredd8)

#We then merge these into an object of class phyloseq.

physeq_bar_plotd8 = phyloseq(OTU.physeqd8, tax.physeqd8, meta.physeqd8)


# Set colors for plotting
my_colors <- c(
  '#a6cee3','#1f88b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff0f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F0FC0", "orange","#DA5024", "#508508", "#CD9BCD",
  "#AD6F3B", "#603000","#D0285", "#652926", "#C84248", 
  "#8569D5", "#5E038F","#D1A33D", "#8A0C64", "#599861", "gray", "black"
)

#If you want different taxonomic level, find and replace the taxonomic level listed here
my_level <- c("Phylum", "Family", "Genus")
my_column <- "Treatment"  #this is the metadata column that we will use in the taxa barplot


abund_filter <- 0.05  # Our abundance threshold
#ml ="Genus"

for(ml in my_level){
  print(ml)
  
  taxa.summaryd8 <- physeq_bar_plot %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(get(my_column), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  taxa.summaryd8 <- as.data.frame(taxa.summaryd8)
  colnames(taxa.summaryd8)[1] <- my_column
  colnames(taxa.summaryd8)[2] <- ml
  
  physeq.taxa.maxd8 <- taxa.summaryd8 %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  
  physeq.taxa.maxd8 <- as.data.frame(physeq.taxa.maxd8)
  colnames(physeq.taxa.maxd8)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_metad8 <- merge(taxa.summaryd8, physeq.taxa.maxd8)
  
  
  physeq_meta_filteredd8 <- filter(physeq_metad8, overall.max>abund_filter)
  #str(physeq_meta_filtered)
  
  physeq_meta_filteredd8$Treatment= factor(physeq_meta_filteredd8$Treatment, c("CON", "AGO"))
  
  # Plot 
  ggplot(physeq_meta_filteredd8, aes(x = get(my_column), y = Abundance.average, fill = get(ml))) + 
    #facet_grid(.~LitterTreatment) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = my_colors) +
    # Remove x axis title
    #theme(axis.title.x = element_blank()) + 
    ylim(c(0,1)) +
    guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
    theme(legend.text=element_text(size=8)) +
    #theme(legend.position="bottom") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(legend.title = element_blank()) +
    ylab("Relative Abundance") +
    xlab(my_column) +
    ggtitle(paste0(ml, " (>", abund_filter * 100,"%) in at least 1 sample")) 
  ggsave(paste0("output/", ml, "d8BarPlot_", my_column, ".png"), height = 5, width = 4)
}

#################################################################
###Differential Abundance with DESeq2
#################################################################


#Adapted from https://joey011.github.io/phyloseq-extensions/DESeq2.html

#First load DESeq2.
#If you need help  with DESeq2 install, see this website
# https://bioconductor.org/packages/release/bioc/html/DESeq2.html


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library("DESeq2")


#To use DESeq, we need no zeros in our OTU table. So we will edit the table by multiplying by 2 and + 1

#First get the OTU table from physeq
physeq_otu_tabled8 <- data.frame(otu_table(physeqd8), check.names = FALSE)

OTU.cleand8 <- physeq_otu_tabled8 + 1


#Now make the phyloseq object:


OTU.physeqd8 = otu_table(as.matrix(OTU.cleand8), taxa_are_rows=TRUE)
tax.physeqd8 = tax_table(as.matrix(tax.cleand8))
meta.physeqd8 = sample_data(metadata.filteredd8)


#We then merge these into an object of class phyloseq.


physeq_deseqd8 = phyloseq(OTU.physeqd8, tax.physeqd8, meta.physeqd8)


#The following two lines actually do all the complicated DESeq2 work. The function phyloseq_to_deseq2 converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated, using the experimental design formula, also shown (the ~body.site term). The DESeq function does the rest of the testing, in this case with default testing framework, but you can actually use alternatives.


diagddsd8 = phyloseq_to_deseq2(physeq_deseqd8, ~ Treatment)
diagddsd8 = DESeq(diagddsd8, test="Wald", fitType="parametric")
#the test type of "Wald" tests for significance of coefficients in a Negative Binomial GLM. This is generally a pretty good assumption for sequencing experiments. This was designed with RNA-seq in mind, but also pretty good for 16S sequencing.


###Investigate test results table

#The following results function call creates a table of the results of the tests. Very fast. The hard work was already stored with the rest of the DESeq2-related data in our latest version of the diagdds object (see above). I then order by the adjusted p-value, removing the entries with an NA value. The rest of this example is just formatting the results table with taxonomic information for nice(ish) display in the HTML output.

#Contrast: this argument specifies what comparison to extract from the object to build a results table. There are exactly three elements:

#  1. the name of a factor in the design formula, 
#  2. the name of the numerator level for the fold change, and 
#  3. the name of the denominator level for the fold change (simplest case)

alpha = 0.05
my_contrast = c("Treatment", "CON", "AGO") 
#my_contrast = c("body.site", "gut", "right palm") 
#my_contrast = c("body.site", "gut", "tongue") 
#my_contrast = c("body.site", "tongue", "left palm") 
#my_contrast = c("body.site", "tongue", "right palm") 
#my_contrast = c("body.site", "right palm", "left palm") 

resd8 = results(diagddsd8, contrast = my_contrast, cooksCutoff = FALSE)

sigtabd8 = res[which(resd8$padj < alpha), ]
sigtab_testd8 <- as(sigtab, "data.frame")
sigtabd8 = cbind(as(sigtabd8, "data.frame"), as(tax_table(physeq_deseq)[rownames(sigtabd8), ], "matrix"))
#head(sigtab)


###Volcano Plot

with(resd8, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-15,15)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(resd8, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="purple"))
with(subset(resd8, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))


#Let's look at the OTUs that were significantly different between the two treatment groups. The following makes a nice ggplot2 summary of the results.


theme_set(theme_bw())
scale_fill_discreted8 <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
#x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
#x = sort(x, TRUE)
#sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabd8$log2FoldChange, sigtabd8$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabd8$Genus = factor(as.character(sigtabd8$Genus), levels=names(x))
DESeq_figd8 = ggplot(sigtabd8, aes(x=Genus, y = log2FoldChange, color=Phylum)) + 
  geom_point(size=3) + 
  ylab(paste0("(", my_contrast[2], "/", my_contrast[3], ")\n", "log2FoldChange")) +
  scale_color_manual(values = my_colors[c(4,6,8,10,12,14,16,18,20,22)]) +
  #ylim(0,8) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggsave(paste0("output/DESeq2-d8", my_contrast[2], "-", my_contrast[3], ".png"), DESeq_figd8, height = 5, width = 10)
